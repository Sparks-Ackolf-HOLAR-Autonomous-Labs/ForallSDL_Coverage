#!/usr/bin/env python3
"""Post-hoc (not registered) significance test for the P3 SiGeC design study.

Rebuilds the registered covering batches (20 constructions each, seeds as in
run_p3_sigec.py) and the 200 random batches of equal size (same RNG seed), on
the committed stiffness-contrast response (results/response_stiffness.npy),
and tests whether covering batches pick a better low-Ge recipe (lower
screening regret) or fit better (lower log-RMSE) than random batches of the
same size: two-sided permutation test on the difference in medians (20000
label permutations) and the Mann-Whitney U test. Writes only
results/posthoc_permutation.json; the registered results are not touched.
"""
import json, os, random
import numpy as np
from scipy.stats import mannwhitneyu
import run_p3_sigec as R  # checks the design hash on import

y = np.load(os.path.join(R.OUT, "response_stiffness.npy"))
ca1 = [R.covering(1, s) for s in range(20)]
ref = min(range(20), key=lambda s: len(ca1[s]))
ca2 = [R.covering(2, 1000 + s, ca1[s]) for s in range(20)]
n1, n2 = len(ca1[ref]), len(ca2[ref])
rng = random.Random(7)
rnd1 = [rng.sample(R.FULL, n1) for _ in range(200)]
rnd2 = [rng.sample(R.FULL, n2) for _ in range(200)]


def scores(sets):
    rows = [R.evaluate(s, y)[0] for s in sets]
    return np.array([r["regret"] for r in rows]), np.array([r["log_rmse"] for r in rows])


def perm(a, b, n=20000, seed=0):
    g = np.random.default_rng(seed)
    obs = np.median(a) - np.median(b)
    pool = np.concatenate([a, b])
    k = len(a)
    hits = 0
    for _ in range(n):
        g.shuffle(pool)
        if abs(np.median(pool[:k]) - np.median(pool[k:])) >= abs(obs) - 1e-12:
            hits += 1
    return float(obs), (hits + 1) / (n + 1)


out = {"note": "post hoc, not registered", "n_runs": {"single": n1, "doubled": n2}}
for name, cov, rnd in (("single_vs_random", ca1, rnd1), ("doubled_vs_random", ca2, rnd2)):
    rc, lc = scores(cov)
    rr, lr = scores(rnd)
    d = {}
    for metric, a, b in (("regret", rc, rr), ("log_rmse", lc, lr)):
        diff, p = perm(a, b)
        d[metric] = {"median_covering": float(np.median(a)), "median_random": float(np.median(b)),
                     "diff_median": diff, "p_permutation": p,
                     "p_mannwhitney": float(mannwhitneyu(a, b, alternative="two-sided").pvalue)}
    out[name] = d
json.dump(out, open(os.path.join(R.OUT, "posthoc_permutation.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
