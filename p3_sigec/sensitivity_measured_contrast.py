#!/usr/bin/env python3
"""P3 design comparisons under the measured, radius-dependent SiC contrast.

Not pre-registered: a reporting check added after ledger seq 105. The
registered P3 study (seq 56-57) scored designs with one SiC contrast for
every radius (stiffness edge primary, density edge as the check). The
small-SiC MD checks then measured the contrast falling with radius
(about R^-4 between 2 and 3 nm, nano_sigec/small_sic/r30). This script
re-scores the same designs (same seeds, same random batches, same
erasure trials) with lam2(R) taken from the measured curve by the rule
in r30/p1_numbers.py (log interpolation up to 3 nm, local slope beyond),
and reports whether the P3 orderings change. No new simulation.
Output: results/p3_sigec_measured_contrast.json.
"""
import json, os, random, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import run_p3_sigec as P  # noqa: E402  (checks the P3 design hash on import)

r = json.load(open(os.path.join(MT, "nano_sigec", "small_sic", "r30", "r30_results.json")))
xs = np.log([c[0] for c in r["measured_curve"]])
ys = np.log([c[1] for c in r["measured_curve"]])


def lamfun(Rn):  # same rule as r30/analyze_r30.py and r30/p1_numbers.py
    if Rn <= 3.0:
        return float(np.exp(np.interp(np.log(Rn), xs, ys)))
    s = min(0.0, (ys[-1] - ys[-2]) / (xs[-1] - xs[-2]))
    return float(np.exp(ys[-1] + s * (np.log(Rn) - xs[-1])))


_orig = P.A.contrast
P.A.contrast = lambda x, kind: lamfun(kind[1]) if isinstance(kind, tuple) else _orig(x, kind)
_ft = P.form_tuple


def form_tuple(name, contrast):
    f = _ft(name, contrast)
    return f if f[0] == "sub" else (f[0], f[1], f[2], ("measured", f[2] * 1e9))


P.form_tuple = form_tuple


def main():
    t0 = time.time()
    y = P.response("measured")
    ys_ = np.load(os.path.join(HERE, "results", "response_stiffness.npy"))
    res = {"not_preregistered": True, "after_ledger_seq": 105,
           "lam2_nm": {f"R{R:g}": lamfun(R) for R in (2.5, 5, 10, 25, 50, 100)},
           "true_min_candidate": {P.NAMES[j]: P.LEV[j][P.FULL[P.CAND[int(np.argmin(y[P.CAND]))]][j]] for j in range(6)},
           "true_min_kappa": float(y[P.CAND].min()),
           "rank_corr_vs_stiffness_all9000": float(np.corrcoef(np.argsort(np.argsort(y)), np.argsort(np.argsort(ys_)))[0, 1])}
    ca1 = [P.covering(1, s) for s in range(20)]
    ref = min(range(20), key=lambda s: len(ca1[s]))
    ca2 = [P.covering(2, 1000 + s, ca1[s]) for s in range(20)]
    n1, n2 = len(ca1[ref]), len(ca2[ref])
    rng = random.Random(7)
    rnd1 = [rng.sample(P.FULL, n1) for _ in range(200)]
    rnd2 = [rng.sample(P.FULL, n2) for _ in range(200)]
    ev = {}
    for name, sets in (("CA1", ca1), ("CA2", ca2), ("RND_nCA1", rnd1), ("RND_nCA2", rnd2)):
        rows = [P.evaluate(s, y)[0] for s in sets]
        reg = [q["regret"] for q in rows]
        ev[name] = {"n_ref": len(sets[ref]) if name.startswith("CA") else len(sets[0]), "regret_median": P.med(reg),
                    "regret_q25_q75": [float(np.quantile(reg, q)) for q in (0.25, 0.75)],
                    "log_rmse_median": P.med([q["log_rmse"] for q in rows])}
    o, _ = P.evaluate(P.ofat(), y)
    ev["OFAT"] = {"n": len(P.ofat()), "regret": o["regret"], "pick": o["pick"]}
    ref2, _ = P.evaluate(ca2[ref], y, importance=True)
    ev["CA2_ref"] = {"regret": ref2["regret"], "pick": ref2["pick"], "importance": ref2["importance"]}
    rng = random.Random(11)
    er = {}
    for name, rows in (("CA1", ca1[ref]), ("CA2", ca2[ref])):
        er[name] = {}
        for k in (1, 2, 3):
            regs = []
            for t in range(500):
                keep = rows[:]
                for i in sorted(rng.sample(range(len(keep)), k), reverse=True):
                    keep.pop(i)
                if t < 60:
                    regs.append(P.evaluate(keep, y)[0]["regret"])
            er[name][f"k{k}"] = {"regret_median_60trials": P.med(regs)}
    res["designs"], res["erasure"] = ev, er
    res["wall_seconds"] = time.time() - t0
    json.dump(res, open(os.path.join(HERE, "results", "p3_sigec_measured_contrast.json"), "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
