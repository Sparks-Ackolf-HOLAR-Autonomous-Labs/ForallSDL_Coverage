#!/usr/bin/env python3
"""How the machinery of P1, P2 and P3 scales with campaign size
(2026-09-30). CPU only, synthetic campaigns with a known law; no ledger writes.

P1  revision cost. (a) HG-SR as committed (phase_c/sr_compare/hgsr.py): fit and
    leave-one-hyperedge-out provenance time against number of rows n.
    (b) the same fit and provenance from per-cell sufficient statistics
    (cell = regime x party): a hyperedge is a set of cells, forgetting it is a
    subtraction of its cells' Gram matrices, so the cost depends on the number
    of cells, not on n. Checked against the direct fit.
P2  verification cost: Merkle root, membership proof size and check, and
    hash-chain verification against number of leaves / ledger entries.
P3  design size: strength-2 covering arrays (AETG-style greedy, scalable) and
    their 2-fold (1-loss-tolerant) extensions against number of factors k and
    levels v; pair coverage of equal-size random batches.
Writes scalability_results.json and fig_scalability.png."""
import hashlib, itertools, json, math, os, random, sys, tempfile, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(MT, "phase_c", "sr_compare"))
sys.path.insert(0, os.path.join(MT, "pilot", "src"))
import hgsr  # noqa: E402
from ledger import Ledger, leaf_hash, merkle_root, merkle_proof, verify_membership  # noqa: E402

RNG = np.random.default_rng(20260930)


# ------------------------------------------------------------------ P1 ----
def campaign(n, nx=5, nt=5, parties=5, bias_party=True):
    x = RNG.uniform(0.05, 0.9, n)
    T = RNG.uniform(300, 1100, n)
    p = RNG.integers(0, parties, n)
    r = 0.2 + 1.9 * x - 2.2 * x ** 2 + 8e-5 * T
    if bias_party:
        r = r * np.where(p == 0, 1.4, 1.0)
    r = r * np.exp(RNG.normal(0, 0.03, n))
    ix = np.minimum((x - 0.05) / 0.85 * nx, nx - 1).astype(int)
    it = np.minimum((T - 300) / 800 * nt, nt - 1).astype(int)
    return [{"x": float(a), "T": float(b), "L": 1.0, "kappa": float(1 / c), "regime": f"x{i}T{j}", "party": f"P{k}"}
            for a, b, c, i, j, k in zip(x, T, r, ix, it, p)]


BASIS = hgsr.make_basis(["x", "T"])
SUBSETS = [s for k in range(1, hgsr.MAX_TERMS + 1) for s in itertools.combinations(range(len(BASIS)), k)]


def cell_stats(rows):
    """Per (regime, party) cell: B^T B, B^T r, r^T r, count."""
    d = hgsr._cols(rows)
    B = np.column_stack([b[2](d) for b in BASIS])
    r = 1.0 / d["kappa"]
    cells = {}
    for i, row in enumerate(rows):
        cells.setdefault((row["regime"], row["party"]), []).append(i)
    return {c: (B[ix].T @ B[ix], B[ix].T @ r[ix], float(r[ix] @ r[ix]), len(ix)) for c, ix in cells.items()}


def fit_stats(S):
    """HG-SR selection from cell statistics: weights 1/|regime|, loss = mean over
    regimes of per-regime MSE, BIC on the number of regimes (as hgsr.fit)."""
    reg = {}
    for (g, _), s in S.items():
        a = reg.setdefault(g, [0, 0, 0, 0])
        for i in range(4):
            a[i] = a[i] + s[i]
    Gw = sum(a[0] / a[3] for a in reg.values())
    bw = sum(a[1] / a[3] for a in reg.values())
    models = []
    for sub in SUBSETS:
        ix = list(sub)
        c = np.linalg.lstsq(Gw[np.ix_(ix, ix)], bw[ix], rcond=None)[0]
        mse = np.mean([(c @ a[0][np.ix_(ix, ix)] @ c - 2 * c @ a[1][ix] + a[2]) / a[3] for a in reg.values()])
        models.append({"terms": ix, "coefs": c, "mse": max(float(mse), 0.0),
                       "complexity": sum(BASIS[i][1] for i in ix)})
    models.sort(key=lambda m: (m["complexity"], m["mse"]))
    front, best = [], np.inf
    for m in models:
        if m["mse"] < best * (1 - 1e-9):
            front.append(m); best = m["mse"]
    ne = float(len(reg))
    return min(front, key=lambda m: ne * np.log(m["mse"] + 1e-30) + m["complexity"] * np.log(ne))


def prov_stats(S):
    edges = {}
    for c in S:
        edges.setdefault("regime:" + c[0], []).append(c)
        edges.setdefault("party:" + c[1], []).append(c)
    base = fit_stats(S)
    out = {}
    for e, cs in edges.items():
        drop = set(cs)
        out[e] = fit_stats({c: s for c, s in S.items() if c not in drop})
    return base, out


def p1():
    res = {"direct": [], "stats": [], "agreement": []}
    for n in (250, 1000, 4000, 16000):
        rows = campaign(n)
        t = time.time(); hgsr.fit(rows, ["x", "T"]); tf = time.time() - t
        t = time.time(); hgsr.provenance(rows, ["x", "T"]); tp = time.time() - t
        res["direct"].append({"n": n, "fit_s": tf, "provenance_s": tp})
        S = cell_stats(rows)
        a = hgsr.fit(rows, ["x", "T"])["selected"]
        b = fit_stats(S)
        same = sorted(a["terms"]) == sorted(b["terms"])
        rel = float(np.max(np.abs(np.array(a["coefs"]) - b["coefs"]) / np.abs(a["coefs"]))) if same else None
        res["agreement"].append({"n": n, "same_terms": same, "max_rel_coef_diff": rel})
    for n in (1000, 10000, 100000, 1000000):
        rows = campaign(n)
        t = time.time(); S = cell_stats(rows); ts = time.time() - t
        t = time.time(); prov_stats(S); tp = time.time() - t
        res["stats"].append({"n": n, "cells": len(S), "build_stats_s": ts, "provenance_s": tp})
    for nx in (5, 10, 20):   # more regimes: cost follows the cell count
        rows = campaign(20000, nx=nx, nt=nx)
        S = cell_stats(rows)
        t = time.time(); prov_stats(S); tp = time.time() - t
        res["stats"].append({"n": 20000, "cells": len(S), "regimes": nx * nx, "provenance_s": tp})
    return res


# ------------------------------------------------------------------ P2 ----
def p2():
    res = {"merkle": [], "chain": []}
    for n in (10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6):
        leaves = [hashlib.sha256(str(i).encode()).hexdigest() for i in range(n)]
        t = time.time(); root = merkle_root(leaves); tr = time.time() - t
        pr = merkle_proof(leaves, n // 3)
        t = time.time(); ok = verify_membership(leaves[n // 3], pr, root); tv = time.time() - t
        res["merkle"].append({"leaves": n, "root_s": tr, "proof_len": len(pr), "verify_s": tv, "ok": ok})
    for n in (100, 1000, 10000):
        with tempfile.TemporaryDirectory() as d:
            led = Ledger(os.path.join(d, "l.jsonl"))
            t = time.time()
            for i in range(n):
                led.append("synthetic", {"i": i})
            ta = time.time() - t
            led = Ledger(os.path.join(d, "l.jsonl"))
            t = time.time(); ok = led.verify_chain(); tv = time.time() - t
        res["chain"].append({"entries": n, "append_total_s": ta, "verify_s": tv, "ok": ok})
    return res


# ------------------------------------------------------------------ P3 ----
def aetg(k, v, mult=1, base=None, ncand=30, seed=0):
    rng = random.Random(seed)
    rows = [list(r) for r in (base or [])]
    need = {}
    for a, b in itertools.combinations(range(k), 2):
        for u in range(v):
            for w in range(v):
                need[(a, b, u, w)] = mult
    for r in rows:
        for a, b in itertools.combinations(range(k), 2):
            need[(a, b, r[a], r[b])] -= 1
    left = {p for p, c in need.items() if c > 0}
    while left:
        best, bs = None, -1
        seedpair = rng.choice(sorted(left))
        for _ in range(ncand):
            r = [None] * k
            a, b, u, w = seedpair
            r[a], r[b] = u, w
            order = [j for j in range(k) if r[j] is None]
            rng.shuffle(order)
            for j in order:
                sc = []
                for l in range(v):
                    s = 0
                    for i in range(k):
                        if r[i] is not None and i != j:
                            p = (min(i, j), max(i, j), r[i] if i < j else l, l if i < j else r[i])
                            s += p in left
                    sc.append(s)
                m = max(sc)
                r[j] = rng.choice([l for l in range(v) if sc[l] == m])
            s = sum((a_, b_, r[a_], r[b_]) in left for a_, b_ in itertools.combinations(range(k), 2))
            if s > bs:
                best, bs = r, s
        rows.append(best)
        for a_, b_ in itertools.combinations(range(k), 2):
            p = (a_, b_, best[a_], best[b_])
            if p in need:
                need[p] -= 1
                if need[p] <= 0:
                    left.discard(p)
    return rows


def pair_cov(rows, k, v):
    have = {(a, b, r[a], r[b]) for r in rows for a, b in itertools.combinations(range(k), 2)}
    return len(have) / (math.comb(k, 2) * v * v)


def p3():
    res = []
    for v in (3, 5):
        for k in (4, 6, 8, 12, 16, 24, 32):
            t = time.time()
            ca1 = aetg(k, v)
            t1 = time.time() - t
            ca2 = aetg(k, v, mult=2, base=ca1)
            rnd = [[random.Random(i * 7 + j).randrange(v) for j in range(k)] for i in range(len(ca1))]
            res.append({"k": k, "v": v, "full_factorial": v ** k, "CA1": len(ca1), "CA2": len(ca2),
                        "ratio": len(ca2) / len(ca1), "lower_bound": v * v, "build_s": t1,
                        "random_same_n_pair_coverage": pair_cov(rnd, k, v)})
            print("P3", res[-1], flush=True)
    return res


def figure(R):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.4))
    d = R["P1"]["direct"]; s = [r for r in R["P1"]["stats"] if "regimes" not in r]
    ax[0].loglog([r["n"] for r in d], [r["provenance_s"] for r in d], "o-", label="refit per hyperedge")
    ax[0].loglog([r["n"] for r in s], [r["provenance_s"] for r in s], "s-", label="cell statistics")
    ax[0].set_xlabel("evidence rows"); ax[0].set_ylabel("provenance time (s)"); ax[0].legend(fontsize=6)
    ax[0].set_title("P1 revision", fontsize=8)
    m = R["P2"]["merkle"]
    ax[1].loglog([r["leaves"] for r in m], [r["root_s"] for r in m], "o-", label="commit (root)")
    ax[1].loglog([r["leaves"] for r in m], [max(r["verify_s"], 1e-6) for r in m], "s-", label="check one record")
    ax[1].set_xlabel("records"); ax[1].set_ylabel("time (s)"); ax[1].legend(fontsize=6); ax[1].set_title("P2 verification", fontsize=8)
    for v, mk in ((3, "o-"), (5, "s-")):
        q = [r for r in R["P3"] if r["v"] == v]
        ax[2].semilogy([r["k"] for r in q], [r["CA1"] for r in q], mk, label=f"covering, {v} levels")
        ax[2].semilogy([r["k"] for r in q], [r["full_factorial"] for r in q], mk.replace("-", ":"), alpha=0.5, label=f"full factorial, {v} levels")
    ax[2].set_xlabel("factors"); ax[2].set_ylabel("runs"); ax[2].legend(fontsize=5); ax[2].set_title("P3 design size", fontsize=8)
    from matplotlib.ticker import NullFormatter
    for a in ax:
        a.tick_params(labelsize=6)
        a.yaxis.set_minor_formatter(NullFormatter())
    ax[2].legend(fontsize=5, loc="center left")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig_scalability.png"), dpi=200)


if __name__ == "__main__" and "--fig" in sys.argv:
    figure(json.load(open(os.path.join(HERE, "scalability_results.json"))))
elif __name__ == "__main__":
    R = {}
    t0 = time.time()
    R["P3"] = p3()
    R["P2"] = p2()
    R["P1"] = p1()
    R["wall_s"] = time.time() - t0
    json.dump(R, open(os.path.join(HERE, "scalability_results.json"), "w"), indent=1, default=float)
    figure(R)
    print(json.dumps({k: R[k] for k in ("P1", "P2")}, indent=1, default=float))
