#!/usr/bin/env python3
"""P3 hidden-subset analysis (p3_mace_design.json, ledger seq 64-65). Written
and pushed before any MACE-MD output existed. Reads raw/out/<id>.json, applies
the registered QC (profile R2 >= 0.9 on both halves), and scores every design
on the hidden points with the registered surrogate and metrics Q1-Q5.
Writes results/p3_mace_results.json."""
import hashlib, itertools, json, os, random, sys
import numpy as np
from scipy.stats import spearmanr
from sklearn.ensemble import ExtraTreesRegressor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
reg = open(os.path.join(HERE, "PREREGISTERED_SHA256")).read().split()[0]
if hashlib.sha256(open(os.path.join(HERE, "p3_mace_design.json"), "rb").read()).hexdigest() != reg:
    sys.exit("p3_mace_design.json changed after registration")
import design_p3_mace as D  # noqa: E402

DES = json.load(open(os.path.join(HERE, "planner_designs.json")))["designs"]
HID = json.load(open(os.path.join(HERE, "hidden_set.json")))
OUT = os.environ.get("P3M_OUT", os.path.join(HERE, "raw", "out"))  # override only for dry runs


def load():
    k, lost = {}, []
    P = json.load(open(os.path.join(HERE, "chpc", "points.json")))["points"]
    for pid in P:
        f = os.path.join(OUT, pid + ".json")
        if not os.path.exists(f):
            lost.append((pid, "no output"))
            continue
        d = json.load(open(f))
        if min(d["profile_R2"]) < 0.9 or not np.isfinite(d["kappa_W_mK"]) or d["kappa_W_mK"] <= 0:
            lost.append((pid, f"QC R2={d['profile_R2']}"))
            continue
        k[pid] = d["kappa_W_mK"]
    return k, lost


def feat(r):
    x, y, f, d, T = (D.F[D.NAMES[j]][r[j]] for j in range(5))
    return [x, y, float(y > 0 and f == "clustered"), float(d == "B"), float(d == "P"), T / 1000]


def fit(rows, k):
    rows = [r for r in rows if D.pid(r) in k]
    m = ExtraTreesRegressor(n_estimators=300, random_state=0, n_jobs=8)
    m.fit(np.array([feat(r) for r in rows]), np.log([k[D.pid(r)] for r in rows]))
    return m, rows


def score(rows, k, hidden):
    m, used = fit(rows, k)
    Xh = np.array([feat(r) for r in hidden])
    yh = np.log([k["H_" + D.pid(r)] for r in hidden])
    p = m.predict(Xh)
    pick = int(np.argmin(p))
    return {"n_used": len(used), "log_rmse": float(np.sqrt(np.mean((p - yh) ** 2))),
            "spearman": float(spearmanr(p, yh).correlation),
            "pick_regret": float(np.exp(yh[pick] - yh.min()) - 1),
            "pair_coverage_after_losses": coverage(used)}


def coverage(rows):
    c = D.counts(rows)
    return sum(1 for key in D.KEYS if c.get(key, 0) >= 1) / len(D.KEYS)


def main():
    k, lost = load()
    hidden = [tuple(r) for r in HID["hidden"] if "H_" + D.pid(tuple(r)) in k]
    res = {"design_sha256": reg, "n_measured": len(k), "lost": lost, "n_hidden_ok": len(hidden)}
    reps = []
    for r in HID["replicates"]:
        a, b = k.get("H_" + D.pid(tuple(r))), k.get("R_" + D.pid(tuple(r)))
        if a and b:
            reps.append(abs(np.log(a / b)))
    res["noise_floor_abs_log_diff"] = reps
    des = {n: [tuple(r) for r in v] for n, v in DES.items()}
    S = {n: score(des[n], k, hidden) for n in ("CA1", "CA2", "OFAT", "RND_pool")}
    rng = random.Random(64)
    pool = [r for r in des["RND_pool"] if D.pid(r) in k]
    for n, size in (("RND15", 15), ("RND12", 12)):
        sc = [score(rng.sample(pool, size), k, hidden) for _ in range(200)]
        S[n] = {m: float(np.median([s[m] for s in sc])) for m in ("log_rmse", "spearman", "pick_regret")}
        S[n]["log_rmse_q25_q75"] = [float(np.quantile([s["log_rmse"] for s in sc], q)) for q in (0.25, 0.75)]
    boot = [score([rng.choice(pool) for _ in pool], k, hidden)["log_rmse"] for _ in range(200)]
    S["RND_pool_bootstrap_log_rmse_median"] = float(np.median(boot))
    res["scores"] = S
    allrows = list(dict.fromkeys(des["CA2"] + des["OFAT"] + des["RND_pool"]))
    m, used = fit(allrows, k)
    X = np.array([feat(r) for r in used]); y = np.log([k[D.pid(r)] for r in used])
    base = np.mean((m.predict(X) - y) ** 2)
    groups = {"x_Ge": [0], "y_C": [1], "carbon_form": [2], "dopant": [3, 4], "T_K": [5]}
    rs = np.random.default_rng(0)
    imp = {}
    for g, cols in groups.items():
        Xp = X.copy(); Xp[:, cols] = X[rs.permutation(len(X))][:, cols]
        imp[g] = float(np.mean((m.predict(Xp) - y) ** 2) - base)
    allk = {tuple(r): k[D.pid(r)] for r in used}
    allk.update({r: k["H_" + D.pid(r)] for r in hidden})
    pairs = []
    for r, v in allk.items():
        if r[1] > 0 and r[2] == 0:
            q = (r[0], r[1], 1, r[3], r[4])
            if q in allk:
                pairs.append(v / allk[q])
    grid = np.array([[x, 0, 0, 0, 0, 0.3] for x in np.linspace(0.1, 0.7, 61)])
    xmin = float(grid[int(np.argmin(m.predict(grid))), 0])
    Xd = np.array([feat(r) for r in D.FULL if r[1] > 0])
    Xc = Xd.copy(); Xd[:, 2] = 0; Xc[:, 2] = 1
    contrast = float(np.mean(np.exp(m.predict(Xd) - m.predict(Xc))))
    res["materials"] = {"importance": imp, "ranking": sorted(imp, key=imp.get, reverse=True),
                        "dispersed_over_clustered_matched_pairs": pairs,
                        "dispersed_over_clustered_surrogate": contrast, "argmin_x_yC0_300K": xmin}
    P = {}
    P["Q1"] = {"CA1": S["CA1"]["log_rmse"], "RND15_median": S["RND15"]["log_rmse"],
               "pass": S["CA1"]["log_rmse"] <= S["RND15"]["log_rmse"]}
    P["Q2"] = {"OFAT": S["OFAT"]["log_rmse"], "RND12_median": S["RND12"]["log_rmse"],
               "pass": S["OFAT"]["log_rmse"] > S["RND12"]["log_rmse"]}
    P["Q3"] = {"ranking": res["materials"]["ranking"], "pass": res["materials"]["ranking"][-1] == "dopant"}
    q4 = float(np.mean(pairs)) if len(pairs) >= 4 else contrast
    P["Q4"] = {"dispersed_over_clustered": q4, "basis": "matched pairs" if len(pairs) >= 4 else "surrogate", "pass": q4 < 1}
    P["Q5"] = {"argmin_x": xmin, "pass": 0.3 <= xmin <= 0.7}
    res["predictions"] = P
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    json.dump(res, open(os.path.join(HERE, "results", "p3_mace_results.json"), "w"), indent=1)
    print(json.dumps({q: v["pass"] for q, v in P.items()}), {n: round(v["log_rmse"], 3) for n, v in S.items() if isinstance(v, dict)})


if __name__ == "__main__":
    main()
