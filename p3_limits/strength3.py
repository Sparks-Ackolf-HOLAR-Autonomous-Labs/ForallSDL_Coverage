#!/usr/bin/env python3
"""Limitation 2: two-factor coverage says nothing about three-factor interactions.

Builds strength-3 covering batches for the registered P3 SiGeC space with the same
greedy rule as the registered strength-2 batches (most deficient three-factor cells
touched; seeds 0-19), then removes redundant runs (latest first) while every
three-factor cell stays sampled. Scores them with the registered surrogate and
screening score and compares with the registered strength-2 batches and with random
batches of equal size. Writes results/strength3.json.
"""
import time
import numpy as np
import common as C

SEEDS = range(20)


def prune(rows):
    rows = list(rows)
    cnt = C.cell_counts(rows, C.KI3, C.NK3)
    for j in range(len(rows) - 1, -1, -1):
        k = C.KI3[C.R.IDX[tuple(rows[j])]]
        if (cnt[k] > 1).all():
            np.add.at(cnt, k, -1)
            rows.pop(j)
    return rows


def summ(v):
    v = np.asarray(v, float)
    return {"median": float(np.median(v)), "q25": float(np.quantile(v, .25)), "q75": float(np.quantile(v, .75)),
            "min": float(v.min()), "max": float(v.max())}


def score(sets, y, yd):
    ev = [C.R.evaluate(s, y)[0] for s in sets]
    return {"runs": summ([len(s) for s in sets]),
            "coverage2": summ([C.coverage_t(s, 2) for s in sets]),
            "coverage3": summ([C.coverage_t(s, 3) for s in sets]),
            "regret": summ([e["regret"] for e in ev]), "log_rmse": summ([e["log_rmse"] for e in ev]),
            "found_best_frac": float(np.mean([e["regret"] == 0 for e in ev])),
            "regret_density_contrast": summ([C.R.evaluate(s, yd)[0]["regret"] for s in sets])}


def main():
    t0 = time.time()
    y = C.load_response()
    yd = C.R.response(C.R.D["model_settings"]["SiC_contrast_check"])
    nl = sorted(C.NL)
    res = {"n_three_factor_cells": int(C.NK3), "n_two_factor_cells": int(C.NK2),
           "lower_bound_strength3": nl[-1] * nl[-2] * nl[-3], "full_factorial": C.N}
    raw = [C.greedy_cover(3, 1, s) for s in SEEDS]
    S3 = [prune(r) for r in raw]
    assert all(C.coverage_t(s, 3) == 1 for s in S3)
    res["sizes_greedy"] = [len(r) for r in raw]
    res["sizes_pruned"] = [len(r) for r in S3]
    ref = int(np.argmin([len(r) for r in S3]))
    U1 = [C.R.covering(1, s) for s in SEEDS]
    U2 = [C.R.covering(2, 1000 + s, U1[s]) for s in SEEDS]
    rng = np.random.default_rng(23)
    n3 = int(np.median([len(r) for r in S3]))
    RND3 = [[C.FULL[i] for i in rng.choice(C.N, n3, replace=False)] for _ in range(200)]
    res["designs"] = {"strength2_single": score(U1, y, yd), "strength2_doubled": score(U2, y, yd),
                      "strength3": score(S3, y, yd), f"random_{n3}": score(RND3, y, yd)}
    res["random_size_for_full_strength3"] = {}
    for n in (300, 600, 1200):
        cov = [C.coverage_t([C.FULL[i] for i in rng.choice(C.N, n, replace=False)], 3) for _ in range(50)]
        res["random_size_for_full_strength3"][n] = {"coverage3_median": float(np.median(cov)),
                                                     "complete_frac": float(np.mean(np.array(cov) == 1))}
    e, _ = C.R.evaluate(S3[ref], y, importance=True)
    imp = e["importance"]
    tot = sum(imp.values())
    res["reference_strength3"] = {"n": len(S3[ref]), "regret": e["regret"], "pick": e["pick"],
                                  "log_rmse": e["log_rmse"],
                                  "importance_share": {k: v / tot for k, v in imp.items()},
                                  "rows": [list(map(int, r)) for r in S3[ref]]}
    res["wall_seconds"] = time.time() - t0
    C.jdump(res, "strength3.json")
    print({k: (v["runs"]["median"], round(v["regret"]["median"], 3), round(v["log_rmse"]["median"], 3),
               v["found_best_frac"]) for k, v in res["designs"].items()}, res["sizes_pruned"])


if __name__ == "__main__":
    main()
