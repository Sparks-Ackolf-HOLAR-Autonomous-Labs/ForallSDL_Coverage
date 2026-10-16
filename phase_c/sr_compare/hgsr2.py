#!/usr/bin/env python3
"""HG-SR v2: the same hypergraph-provenance engine, with the four
accuracy fixes found in the v1 diagnosis (HGSR_V2.md, 2026-10-02).

Why v1 trailed PySR on the md task (see HGSR_V2.md for the numbers):

  D1 Selection starved by its own effective sample size. v1 picks the
     model by BIC with n_eff = number of regime hyperedges (5 on md),
     so every unit of complexity costs log(5) against a likelihood term
     scaled by 5: the criterion almost always stops at 3 terms.
  D2 Hyperedges too coarse. v1's regime hyperedge is an (x band, T
     band) cell and ignores box length, so a 12-row finite-size series
     counts as ONE unit of evidence; the balancing then throws away the
     size information the law needs.
  D3 Absolute loss on r = 1/kappa. Pure-element rows have r near 0.03
     and alloy rows near 1, so an absolute loss on r ignores the pure
     rows, while the benchmark scores kappa, where those rows dominate.
  D4 Hand-written basis with hand-set term costs. The PySR-discovered
     term x(1-x)/L was added at cost 7, so the selector priced it out;
     and no basis term carries the slow L^-1/2 size convergence that
     Rayleigh-scattered long-wavelength phonons give in an alloy.

v2 keeps everything that makes HG-SR different (hyperedge-balanced
loss, leave-one-hyperedge-out provenance, blast radius, exact
deterministic clean retrain) and changes only the search:

  F1 select along the complexity front by leave-one-hyperedge-out
     cross-validation (deterministic, so unlearning stays exact);
  F2 hyperedges are design cells (x, T, size), seed replicates of one
     cell stay together;
  F3 relative loss (weights 1/r^2), i.e. percentage error, which is
     what a materials user reads off a kappa law;
  F4 generated basis: products of up to two physical primitives
     (x, 1-x, x(1-x), T, 1/T, 1/L, L^-1/2), uniform term cost.
"""
from __future__ import annotations

import itertools

import numpy as np

MAX_TERMS = 4
RELTOL = 0.30
NOISE_PRIOR_DOF = 3      # I4: shrink a cell's seed scatter towards the pooled value


def primitives(features):
    p = [("x", lambda d: d["x"]),
         ("(1-x)", lambda d: 1 - d["x"]),
         ("x(1-x)", lambda d: d["x"] * (1 - d["x"])),
         ("t", lambda d: d["T"] / 1000.0),
         ("1/t", lambda d: 1000.0 / d["T"])]
    if "L" in features:
        p += [("1/l", lambda d: 100.0 / d["L"]),          # l = L / 100 A
              ("l^-1/2", lambda d: np.sqrt(100.0 / d["L"]))]
    return p


SKIP = {frozenset(("x", "x(1-x)")), frozenset(("(1-x)", "x(1-x)")),
        frozenset(("t", "1/t")), frozenset(("1/l", "l^-1/2")),
        frozenset(("x", "(1-x)"))}


def make_basis(features, generated=True, prims=None, extended=False,
               proposals=None):
    """(name, cost, fn). Generated: constant, primitives, and products of
    two distinct primitives; cost = number of primitive factors.
    extended (I3 scaling test): also squares and three-factor products.
    proposals (I2): extra (name, cost, fn) terms harvested from an open
    search engine, appended after the generated basis."""
    if not generated:
        import hgsr
        return hgsr.make_basis(features)
    prim = prims if prims is not None else primitives(features)
    b = [("1", 1, lambda d: np.ones_like(d["x"]))]
    b += [(n, 1, f) for n, f in prim]
    for (na, fa), (nb, fb) in itertools.combinations(prim, 2):
        if frozenset((na, nb)) in SKIP:
            continue
        b.append((f"{na}*{nb}", 2, (lambda fa, fb: lambda d: fa(d) * fb(d))(fa, fb)))
    if extended:
        for na, fa in prim:
            b.append((f"{na}^2", 2, (lambda fa: lambda d: fa(d) ** 2)(fa)))
        for trio in itertools.combinations(prim, 3):
            names = [t[0] for t in trio]
            if any(frozenset(p) in SKIP for p in itertools.combinations(names, 2)):
                continue
            fs = [t[1] for t in trio]
            b.append(("*".join(names), 3,
                      (lambda fs: lambda d: fs[0](d) * fs[1](d) * fs[2](d))(fs)))
    if proposals:
        b += list(proposals)
    return b


def _cols(rows, keys=("x", "T", "L", "kappa")):
    ks = set(keys) | {"x", "T", "L", "kappa"}
    for r in rows[:1]:
        ks |= {k for k, v in r.items() if isinstance(v, (int, float))}
    return {k: np.array([float(r.get(k, 0.0)) for r in rows]) for k in ks}


def cell(r, size_aware=True):
    """Design-cell hyperedge: regime plus size band (ncz ~ L / 133 A)."""
    if "cell" in r:
        return r["cell"]
    if size_aware and r.get("L", 0.0) > 0:
        return f"{r['regime']}L{int(round(r['L'] / 133.0))}"
    return r["regime"]


def _groups(rows, size_aware):
    g = {}
    for i, r in enumerate(rows):
        g.setdefault(cell(r, size_aware), []).append(i)
    return [np.array(g[k]) for k in sorted(g)]


def _solve(M, r, w, constrain=False):
    sw = np.sqrt(w)
    if constrain:
        # I1 Matthiessen: every basis term is a non-negative function on the
        # domain, so c >= 0 means every mechanism adds resistance.
        from scipy.optimize import nnls
        c, _ = nnls(M * sw[:, None], r * sw, maxiter=50 * M.shape[1])
        return c
    c, *_ = np.linalg.lstsq(M * sw[:, None], r * sw, rcond=None)
    return c


def _gloss(pred_r, r, groups, rel):
    e = (pred_r - r) / (r if rel else 1.0)
    return float(np.mean([np.mean(e[g] ** 2) for g in groups]))


def noise_weights(rows, groups, kappa):
    """I4: per-cell weight = 1 / (shrunk variance of ln kappa over the
    cell's DISTINCT values). Exact copies do not count as new seeds, so
    replicating a batch cannot raise its cell's weight (replication
    invariance is kept). Normalised to mean 1 over cells."""
    lk = np.log(kappa)
    stats = []
    for g in groups:
        u = np.unique(np.round(lk[g], 12))
        stats.append((len(u), float(np.var(u, ddof=1)) if len(u) > 1 else 0.0))
    dof = sum(m - 1 for m, _ in stats)
    pooled = (sum((m - 1) * v for m, v in stats) / dof) if dof > 0 else 1.0
    k = NOISE_PRIOR_DOF
    cw = np.array([1.0 / ((k * pooled + (m - 1) * v) / (k + m - 1))
                   for m, v in stats])
    return cw / cw.mean(), pooled


def physics_domain(x=(0.05, 0.95), T=(300.0, 1300.0), L=None, extra=None,
                   bulk=True):
    """I1 shape constraints (Matthiessen + Nordheim + finite-size), as an
    admissibility filter on candidate laws rather than a sign on every
    coefficient. A candidate is admissible when, on a grid over the
    declared domain:
      P1 r > 0 everywhere (finite positive kappa);
      P2 (if L) r does not increase with box length, from the smallest box
         through L -> infinity (a bigger box never conducts worse);
      P3 the size-free (bulk) limit is positive everywhere;
      P4 Nordheim: at every T the bulk resistance at x = 0.5 is at least
         that at x = 0.15 and x = 0.85 (where in the domain), i.e. alloy
         scattering peaks near mid-composition.
    `extra` adds fixed numeric columns (e.g. dopant flags = 0) to grid rows.
    bulk=False states P2-P4 only inside the sampled box range (the
    "limit" rows sit at the largest box), for laws declared valid only
    there; bulk=True also demands a physical L -> infinity limit."""
    xs = np.linspace(x[0], x[1], 19)
    Ts = np.linspace(T[0], T[1], 6)
    extras = extra if isinstance(extra, list) else [extra or {}]
    box, inf, size, nord = [], [], [], []
    Ls = list(np.geomspace(L[0], L[1], 6)) if L else [0.0]
    Llim = (np.inf if bulk else L[1]) if L else 0.0
    for extra in extras:
        for xv in xs:
            for Tv in Ts:
                seq = []
                for Lv in Ls:
                    box.append(dict(extra, x=xv, T=Tv, L=Lv))
                    seq.append(len(box) - 1)
                inf.append(dict(extra, x=xv, T=Tv, L=Llim))
                size.append(seq)
    lim_rows = len(inf)
    for extra in extras:
        for Tv in Ts:
            for xo in (0.15, 0.85):
                if x[0] <= xo <= x[1] and x[0] <= 0.5 <= x[1]:
                    nord.append((len(inf), len(inf) + 1))
                    inf += [dict(extra, x=0.5, T=Tv, L=Llim),
                            dict(extra, x=xo, T=Tv, L=Llim)]
    return {"box": box, "inf": inf, "size": size if L else [], "nord": nord}


def _admissible_fn(basis, phys, target="r"):
    if phys is None:
        return None
    def mat(rows):
        d = _cols([dict(r, kappa=1.0) for r in rows])
        with np.errstate(all="ignore"):
            return np.column_stack([b[2](d) for b in basis])
    Gb, Gi = mat(phys["box"]), mat(phys["inf"])
    size = np.array(phys["size"], dtype=int) if phys["size"] else None
    nord = np.array(phys["nord"], dtype=int) if phys["nord"] else None

    def ok(key, c):
        with np.errstate(all="ignore"):
            rb = Gb[:, key] @ c
            ri = Gi[:, key] @ c
            if target == "kappa":          # checks are stated on r = 1/kappa
                rb, ri = 1.0 / rb, 1.0 / ri
        if not (np.all(np.isfinite(rb)) and np.all(rb > 0)):
            return False
        if not (np.all(np.isfinite(ri)) and np.all(ri > 0)):
            return False
        tol = 1e-9 * np.max(np.abs(rb))
        if size is not None:
            # size rows i align with limit rows i (vectorised P2)
            seq = np.column_stack([rb[size], ri[:len(size)]])
            if np.any(np.diff(seq, axis=1) > tol):
                return False
        if nord is not None and np.any(ri[nord[:, 0]] <= ri[nord[:, 1]] * (1 + 1e-6)):
            return False
        return True
    return ok


def _search(F, r, w, groups, rel, max_terms, req, constrain, search, beam,
            costs, admissible=None):
    """Best subset per complexity. 'exhaustive' enumerates every subset up
    to max_terms (all subsets contain the required terms); 'beam' (I3) is a
    deterministic forward beam search: keep the `beam` lowest-loss subsets
    of each size and extend each by one term."""
    nb = F.shape[1]
    others = [i for i in range(nb) if i not in req]
    best_at, seen = {}, set()

    def evaluate(terms):
        key = tuple(sorted(terms))
        if key in seen:
            return None
        seen.add(key)
        M = F[:, key]
        c = _solve(M, r, w, constrain)
        loss = _gloss(M @ c, r, groups, rel)
        # Inadmissible subsets are not recorded, but the beam may still
        # extend them (a third term can repair a two-term law's shape).
        if constrain and len(key) > len(req) and np.any(c[[key.index(i) for i in key if i not in req]] <= 0):
            return loss          # a zero coefficient is the smaller subset
        if admissible is not None and not admissible(key, c):
            return loss
        cx = int(sum(costs[i] for i in key))
        if cx not in best_at or loss < best_at[cx]["mse"] - 1e-15:
            best_at[cx] = {"terms": list(key), "coefs": [float(x) for x in c],
                           "complexity": cx, "mse": loss}
        return loss

    kmin = max(len(req), 1)
    if search == "exhaustive":
        if req:
            evaluate(list(req))
        for k in range(1, max_terms - len(req) + 1):
            for extra in itertools.combinations(others, k):
                evaluate(list(req) + list(extra))
    else:
        level = [tuple(sorted(req))] if req else [()]
        if req:
            evaluate(list(req))
        for _ in range(max_terms - len(req)):
            cand = []
            for base in level:
                for j in others:
                    if j in base:
                        continue
                    key = tuple(sorted(base + (j,)))
                    if key in seen:
                        continue
                    loss = evaluate(list(key))
                    if loss is not None:
                        cand.append((loss, key))
            cand.sort()
            level = [k for _, k in cand[:beam]]
            if not level:
                break
    return best_at, len(seen)


def fit(rows, features, max_terms=MAX_TERMS, generated=True, rel=True,
        size_aware=True, select="logo", pooled=False, constrain=False,
        require=(), noise=False, proposals=None, prims=None, extended=False,
        search="exhaustive", beam=200, physics=None, target="r"):
    basis = make_basis(features, generated, prims, extended, proposals)
    names = [b[0] for b in basis]
    d = _cols(rows)
    # target "r" fits thermal resistivity (Matthiessen-additive scattering);
    # target "kappa" fits conductivity directly (additive accumulation,
    # e.g. a kappa_inf - A L^-1/2 finite-size tail). The variable is still
    # called r below.
    r = 1.0 / d["kappa"] if target == "r" else d["kappa"].copy()
    n = len(rows)
    groups = [np.arange(n)] if pooled else _groups(rows, size_aware)
    w = np.ones(n)
    pooled_var = None
    if not pooled:
        cw = np.ones(len(groups))
        if noise:
            cw, pooled_var = noise_weights(rows, groups, d["kappa"])
        for gi, g in enumerate(groups):
            w[g] = cw[gi] / len(g)
    if rel:
        w = w / r ** 2
    F = np.column_stack([b[2](d) for b in basis])
    req = [names.index(q) for q in require if q in names]
    best_at, n_eval = _search(F, r, w, groups, rel, max_terms, req, constrain,
                              search, beam, [b[1] for b in basis],
                              _admissible_fn(basis, physics, target))
    for m in best_at.values():
        m["names"] = [names[i] for i in m["terms"]]
    front, best = [], np.inf
    for cx in sorted(best_at):
        m = best_at[cx]
        if m["mse"] < best * (1 - 1e-9):
            front.append(m)
            best = m["mse"]

    if select == "logo" and len(groups) >= 3:
        cvg = groups if not pooled else _groups(rows, size_aware)
        for m in front:
            errs = []
            for g in cvg:
                keep = np.setdiff1d(np.arange(n), g)
                M = F[:, m["terms"]]
                if len(keep) <= len(m["terms"]):
                    errs.append(np.inf)
                    continue
                c = _solve(M[keep], r[keep], w[keep], constrain)
                p = M[g] @ c
                if target == "kappa":
                    with np.errstate(all="ignore"):
                        p = 1.0 / p
                bad = ~(p > 0)
                rg = r[g] if target == "r" else 1.0 / r[g]
                le = np.where(bad, 10.0, np.log(np.where(bad, 1.0, p) / rg))
                errs.append(np.mean(le ** 2))
            m["logo_msle"] = float(np.mean(errs))
        lo = min(m["logo_msle"] for m in front)
        # one-standard-error style tie-break: simplest within 5% of best
        sel = next(m for m in front if m["logo_msle"] <= lo * 1.05)
    else:
        n_eff = float(np.sum(w)) if not rel else float(len(groups))
        sel = min(front, key=lambda m: n_eff * np.log(m["mse"] + 1e-30)
                  + m["complexity"] * np.log(max(n_eff, 2.0)))
    return {"front": front, "selected": sel, "features": features,
            "generated": generated, "n": n, "basis": basis,
            "n_basis": len(basis), "n_evaluated": n_eval,
            "pooled_noise_var": pooled_var, "target": target}


def predict(model, rows):
    basis = model.get("basis") or make_basis(model["features"], model["generated"])
    d = _cols([dict(r, kappa=1.0) for r in rows])
    sel = model["selected"]
    M = np.column_stack([basis[i][2](d) for i in sel["terms"]])
    y = M @ np.array(sel["coefs"])
    return y if model.get("target", "r") == "kappa" else 1.0 / y


def provenance(rows, features, **kw):
    """Leave-one-hyperedge-out term provenance and blast radius, same
    semantics as hgsr.provenance, on v2's fit."""
    full = fit(rows, features, **kw)
    sel = full["selected"]
    base = dict(zip(sel["names"], sel["coefs"]))
    edges = {}
    for i, rr in enumerate(rows):
        for kind in ("regime", "party"):
            edges.setdefault(f"{kind}:{rr[kind]}", []).append(i)
    prov = {nm: [] for nm in sel["names"]}
    for eid, ix in sorted(edges.items()):
        drop = set(ix)
        keep = [rr for i, rr in enumerate(rows) if i not in drop]
        if len(keep) < MAX_TERMS + 2 or len(keep) == len(rows):
            continue
        red_sel = fit(keep, features, **kw)["selected"]
        red = dict(zip(red_sel["names"], red_sel["coefs"]))
        for nm, c0 in base.items():
            gone = nm not in red
            if gone or abs(red[nm] - c0) > RELTOL * abs(c0):
                prov[nm].append({"edge": eid,
                                 "effect": "deselected" if gone else "coef-moved"})
    blast = {}
    for nm, deps in prov.items():
        for dep in deps:
            blast.setdefault(dep["edge"], []).append(nm)
    return {"selected": sel, "term_provenance": prov, "blast_radius": blast}
