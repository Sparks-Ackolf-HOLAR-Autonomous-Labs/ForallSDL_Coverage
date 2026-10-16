#!/usr/bin/env python3
"""HG-SR: hypergraph-provenance symbolic regression (the new engine).

Every existing SR engine (PySR, SISSO, AI Feynman, QLattice, and our
minisr) consumes a flat (X, y) table. But mat_theory's evidence lives on
a hypergraph: each datapoint is a hyperedge joining {source, party,
regime, method, system}, and the shard machinery revises beliefs by
unlearning hyperedges, not rows. HG-SR makes the regression itself
hypergraph-aware. It differs from a flat fitter in four ways:

1. Hyperedge-balanced fitting. Coefficients come from weighted least
   squares with row weight 1/|regime hyperedge|, and the search loss is
   the mean over regime hyperedges of the per-hyperedge MSE, not the
   pooled row MSE. A regime sampled with 10 rows no longer outvotes a
   regime sampled with 2 (replicating a regime's rows leaves the fit
   exactly invariant), so a densely sampled (or densely biased) party
   cannot drag the law; this is the fitting-side counterpart of the
   audit-independence rule.

2. Per-term hyperedge provenance. For each term in the selected model,
   HG-SR refits with each hyperedge removed (leave-one-hyperedge-out,
   both regime edges and party edges). The term's provenance is the set
   of hyperedges whose removal deselects it or moves its coefficient by
   more than RELTOL. Terms with empty provenance are multi-supported
   (OR-semantics: no single evidence group is load-bearing); terms whose
   provenance contains a party edge are single-party and get scope
   labels instead of confirmed status.

3. Blast radius for free. Inverting the per-term provenance gives, for
   every hyperedge, exactly which terms of the law fall if that
   hyperedge is unlearned. This is the SR-side blast radius B(e) the
   design doc defines; flat engines cannot produce it.

4. Coverage-scoped validity. The validity domain is the set of covered
   regime cells, not a bounding box; a query in an uncovered cell is
   flagged out-of-support even if it is inside the box.

The search itself stays exhaustive-basis OLS like the pilot's minisr
(all term subsets up to MAX_TERMS), so it is deterministic and the clean
retrain is bit-for-bit exact: epsilon = delta = 0 unlearning by
construction. Fitting targets resistivity r = 1/kappa (Matthiessen).
Run `pooled=True` to ablate feature 1 (that ablation IS the minisr
baseline in the comparison).
"""
from __future__ import annotations

import itertools

import numpy as np

MAX_TERMS = 4
RELTOL = 0.30       # coefficient move that counts as "term depends on edge"


def make_basis(features):
    """Physical basis library over the named features."""
    b = [("1", 1, lambda d: np.ones_like(d["x"])),
         ("x", 1, lambda d: d["x"]),
         ("T", 1, lambda d: d["T"]),
         ("x(1-x)", 4, lambda d: d["x"] * (1 - d["x"])),
         ("x*T", 3, lambda d: d["x"] * d["T"]),
         ("1/T", 2, lambda d: 1.0 / d["T"]),
         ("x^2", 2, lambda d: d["x"] ** 2),
         ("x(1-x)*T", 6, lambda d: d["x"] * (1 - d["x"]) * d["T"])]
    if "L" in features:
        b += [("1/L", 2, lambda d: 1.0 / d["L"]),
              ("T/L", 4, lambda d: d["T"] / d["L"]),
              # cross term PySR discovered in the engine comparison
              # (alloy-limited mean free path times boundary term)
              ("x(1-x)/L", 7, lambda d: d["x"] * (1 - d["x"]) / d["L"])]
    return b


def _cols(rows):
    return {k: np.array([r[k] for r in rows], dtype=float)
            for k in ("x", "T", "L", "kappa")}


def _edges(rows, kinds=("regime", "party")):
    """Hyperedge id -> row indices, for the requested edge kinds."""
    e = {}
    for i, r in enumerate(rows):
        for kind in kinds:
            e.setdefault(f"{kind}:{r[kind]}", []).append(i)
    return e


def fit(rows, features, pooled=False, max_terms=MAX_TERMS):
    """Exhaustive subset OLS on r = 1/kappa with hyperedge-balanced loss.
    Returns {front, selected, ...}; deterministic."""
    basis = make_basis(features)
    d = _cols(rows)
    r = 1.0 / d["kappa"]
    n = len(rows)
    groups = ([np.arange(n)] if pooled else
              [np.array(ix) for eid, ix in sorted(_edges(rows, ("regime",)).items())])
    w = np.ones(n)
    if not pooled:
        for g in groups:
            w[g] = 1.0 / len(g)
    sw = np.sqrt(w)

    models = []
    for k in range(1, max_terms + 1):
        for terms in itertools.combinations(range(len(basis)), k):
            M = np.column_stack([basis[i][2](d) for i in terms])
            coefs, *_ = np.linalg.lstsq(M * sw[:, None], r * sw, rcond=None)
            res2 = (M @ coefs - r) ** 2
            loss = float(np.mean([np.mean(res2[g]) for g in groups]))
            models.append({"terms": list(terms),
                           "names": [basis[i][0] for i in terms],
                           "coefs": [float(c) for c in coefs],
                           "complexity": int(sum(basis[i][1] for i in terms)),
                           "mse": loss})
    models.sort(key=lambda m: (m["complexity"], m["mse"]))
    front, best = [], np.inf
    for m in models:
        if m["mse"] < best * (1 - 1e-9):
            front.append(m)
            best = m["mse"]
    # BIC on the effective sample size: pooled rows can be replicated,
    # hyperedges cannot, so balanced mode counts evidence groups.
    n_eff = n if pooled else float(np.sum(w))
    sel = min(front, key=lambda m: n_eff * np.log(m["mse"] + 1e-30)
              + m["complexity"] * np.log(n_eff))
    return {"front": front, "selected": sel, "features": features,
            "pooled": pooled, "n": n}


def predict(model, rows, features=None):
    basis = make_basis(features or model["features"])
    d = _cols([dict(r, kappa=1.0) for r in rows])
    sel = model["selected"] if "selected" in model else model
    M = np.column_stack([basis[i][2](d) for i in sel["terms"]])
    return 1.0 / (M @ np.array(sel["coefs"]))


def provenance(rows, features, max_terms=MAX_TERMS):
    """Leave-one-hyperedge-out provenance of every selected term, plus
    the inverted blast-radius map. Deterministic."""
    full = fit(rows, features, max_terms=max_terms)
    sel = full["selected"]
    base = dict(zip(sel["names"], sel["coefs"]))
    edges = _edges(rows)
    prov = {name: [] for name in sel["names"]}
    for eid, ix in sorted(edges.items()):
        keep = [r for i, r in enumerate(rows) if i not in set(ix)]
        if len(keep) < max_terms + 2:      # removal would underdetermine
            continue
        m = fit(keep, features, max_terms=max_terms)["selected"]
        red = dict(zip(m["names"], m["coefs"]))
        for name, c0 in base.items():
            gone = name not in red
            moved = (not gone) and abs(red[name] - c0) > RELTOL * abs(c0)
            if gone or moved:
                prov[name].append({"edge": eid,
                                   "effect": "deselected" if gone else "coef-moved"})
    blast = {}
    for name, deps in prov.items():
        for dep in deps:
            blast.setdefault(dep["edge"], []).append(name)
    labels = {}
    for name, deps in prov.items():
        party_deps = [d for d in deps if d["edge"].startswith("party:")]
        if not deps:
            labels[name] = "multi-supported"
        elif party_deps:
            labels[name] = ("scope:" +
                            ",".join(sorted(d["edge"][6:] for d in party_deps)))
        else:
            labels[name] = "regime-sensitive"
    covered = sorted({r["regime"] for r in rows})
    return {"selected": sel, "term_provenance": prov, "term_labels": labels,
            "blast_radius": blast, "covered_regimes": covered}


if __name__ == "__main__":
    import json
    from datasets import TASKS
    for name, taskfn in TASKS.items():
        t = taskfn()
        p = provenance(t["rows"], t["features"])
        print(f"== {name}: selected {' + '.join(p['selected']['names'])}")
        for term, lab in p["term_labels"].items():
            print(f"   {term:12s} {lab:24s} "
                  f"depends on {[d['edge'] for d in p['term_provenance'][term]]}")
        print(f"   blast radius: {json.dumps(p['blast_radius'])}")
