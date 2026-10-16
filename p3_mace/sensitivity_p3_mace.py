#!/usr/bin/env python3
"""POST-HOC (not registered): context for the registered Q1-Q5 results. Run after
analyze_p3_mace.py. (a) percentile of CA1 and OFAT within the 200 registered random draws,
(b) spread of the hidden log kappa, (c) an x_Ge-only surrogate on each design, to test whether
the MD-box response is effectively one-dimensional."""
import json, os, random, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import analyze_p3_mace as A  # noqa: E402
from sklearn.ensemble import ExtraTreesRegressor  # noqa: E402

k, _ = A.load()
des = {n: [tuple(r) for r in v] for n, v in A.DES.items()}
hidden = [tuple(r) for r in A.HID["hidden"] if "H_" + A.D.pid(tuple(r)) in k]
yh = np.log([k["H_" + A.D.pid(r)] for r in hidden])
rng = random.Random(64)
pool = [r for r in des["RND_pool"] if A.D.pid(r) in k]
out = {"hidden_log_kappa_sd": float(np.std(yh)), "hidden_kappa_range_W_mK": [float(np.exp(yh.min())), float(np.exp(yh.max()))]}
for n, size in (("RND15", 15), ("RND12", 12)):
    r = [A.score(rng.sample(pool, size), k, hidden)["log_rmse"] for _ in range(200)]
    ref = {"RND15": "CA1", "RND12": "OFAT"}[n]
    v = A.score(des[ref], k, hidden)["log_rmse"]
    out[f"{ref}_percentile_in_{n}"] = float(np.mean(np.array(r) < v))


def xonly(rows):
    rows = [r for r in rows if A.D.pid(r) in k]
    m = ExtraTreesRegressor(n_estimators=300, random_state=0).fit(
        np.array([[A.feat(r)[0], A.feat(r)[5]] for r in rows]), np.log([k[A.D.pid(r)] for r in rows]))
    p = m.predict(np.array([[A.feat(r)[0], A.feat(r)[5]] for r in hidden]))
    return float(np.sqrt(np.mean((p - yh) ** 2)))


out["x_and_T_only_log_rmse"] = {n: xonly(des[n]) for n in ("CA1", "CA2", "OFAT", "RND_pool")}
json.dump(out, open(os.path.join(HERE, "results", "p3_mace_posthoc.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
