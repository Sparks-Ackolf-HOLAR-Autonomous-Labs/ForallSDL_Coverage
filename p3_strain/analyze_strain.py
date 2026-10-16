#!/usr/bin/env python3
"""Pre-registered analysis of the P3 strain seed batch
(strain_batch_design.json, 'analysis_preregistered'). Reads only runs
whose hashes were committed to the ledger (commit_strain.py).

Model (weighted equally, OLS on ln kappa):
  ln k = b0 + b1 x + b2 x^2 + b3 t + (a + b x + c t) eps + d eps*dir,
  t = T/1000, dir = 0 cross-plane / 1 in-plane.
Bootstrap (2000 resamples, seeded) for CIs. Erasure analysis: pair
coverage of C1-only and C2 after removing never-completed runs.
Materials call: kappa change vs relaxed bulk for pseudomorphic films on
Si(001) (eps = -0.042 x) and on relaxed Ge (eps = +0.04 (1-x)), both
geometries, 300 and 1100 K; |eps| > 0.03 is flagged as extrapolation.
Usage: python3 analyze_strain.py
"""
import itertools
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, "strain_batch_design.json")))
LEDGER = os.path.join(HERE, "..", "phase_c", "ledger_campaign.jsonl")
NAMES = ["b0", "b1_x", "b2_x2", "b3_t", "a_eps", "b_eps_x", "c_eps_t", "d_eps_dir"]


def committed():
    for l in open(LEDGER):
        e = json.loads(l)
        if e["payload"].get("op_id") == "commit-" + D["campaign"]:
            return set(e["payload"].get("never_completed", []))
    raise SystemExit("batch not committed to the ledger yet")


def X(rows):
    return np.array([[1, r["x"], r["x"] ** 2, r["T"] / 1000, r["eps"], r["eps"] * r["x"],
                      r["eps"] * r["T"] / 1000, r["eps"] * r["dir"]] for r in rows])


def fit(rows):
    y = np.log([r["kappa"] for r in rows])
    c, *_ = np.linalg.lstsq(X(rows), y, rcond=None)
    res = y - X(rows) @ c
    return c, float(np.sqrt(np.mean(res ** 2)))


def dlnk(c, x, T, eps, d):
    t = T / 1000
    return (c[4] + c[5] * x + c[6] * t) * eps + c[7] * eps * d


def coverage(rows):
    lv = [sorted({p[k] for p in D["points"]}) for k in ("x", "eps", "T", "dir")]
    keys = ("x", "eps", "T", "dir")
    need = {(a, b, u, v) for a, b in itertools.combinations(range(4), 2)
            for u in lv[a] for v in lv[b]}
    have = {(a, b, r[keys[a]], r[keys[b]]) for r in rows
            for a, b in itertools.combinations(range(4), 2)}
    return len(need & have), len(need)


def main():
    missing = committed()
    rows = []
    for p in D["points"]:
        if p["run_id"] in missing:
            continue
        j = json.load(open(os.path.join(HERE, "results", p["run_id"] + ".json")))
        rows.append(dict(p, kappa=j["kappa_W_mK"]))
    c, rms = fit(rows)
    rng = np.random.default_rng(20260929)
    boot = np.array([fit([rows[i] for i in rng.integers(0, len(rows), len(rows))])[0]
                     for _ in range(2000)])
    ci = {n: [round(float(np.percentile(boot[:, i], 2.5)), 4),
              round(float(np.percentile(boot[:, i], 97.5)), 4)] for i, n in enumerate(NAMES)}
    c1_rows = [r for r in rows if r["in_C1"]]
    calls = {}
    for sub, epsf in (("on_Si001", lambda x: -0.042 * x), ("on_relaxed_Ge", lambda x: 0.04 * (1 - x))):
        for T in (300, 1100):
            for d, dn in ((0, "cross-plane"), (1, "in-plane")):
                for x in (0.1, 0.3, 0.5, 0.7, 0.9):
                    e = epsf(x)
                    k = f"{sub}|T{T}|{dn}|x{x}"
                    samp = np.array([dlnk(b, x, T, e, d) for b in boot])
                    calls[k] = {"eps": round(e, 4), "pct_change": round(100 * (np.exp(dlnk(c, x, T, e, d)) - 1), 1),
                                "ci95_pct": [round(100 * (np.exp(np.percentile(samp, q)) - 1), 1) for q in (2.5, 97.5)],
                                "extrapolated": abs(e) > 0.03 + 1e-9}
    out = {"n_completed": len(rows), "never_completed": sorted(missing),
           "coef": {n: round(float(v), 4) for n, v in zip(NAMES, c)}, "ci95": ci,
           "rms_ln": round(rms, 4),
           "coverage_C2_completed": coverage(rows), "coverage_C1_completed": coverage(c1_rows),
           "predictions": {
               "compressive_raises_kappa_cross_plane_300K_x0.5": bool(dlnk(c, 0.5, 300, -0.02, 0) > 0),
               "sensitivity_smaller_in_alloy_than_dilute": bool(abs(c[4] + c[5] * 0.5) < abs(c[4] + c[5] * 0.1))},
           "materials_calls": calls}
    json.dump(out, open(os.path.join(HERE, "analysis_strain.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "materials_calls"}, indent=1))
    for k, v in calls.items():
        if "x0.5" in k or "x0.1" in k:
            print(k, v)


if __name__ == "__main__":
    main()
