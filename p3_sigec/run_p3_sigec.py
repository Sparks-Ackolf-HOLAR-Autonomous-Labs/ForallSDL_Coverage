#!/usr/bin/env python3
"""P3 SiGeC + dopant campaign-design study (p3_sigec_design.json, ledger seq 56).

Refuses to run if the design changed after registration. Evaluates the
committed SiGeC response model on the full factorial, builds the designs,
and measures coverage, erasure tolerance, screening regret and transversal
follow-up exactly as registered. Writes results/p3_sigec_results.json.
"""
import hashlib, itertools, json, os, random, sys, time
import numpy as np
from sklearn.ensemble import ExtraTreesRegressor

HERE = os.path.dirname(os.path.abspath(__file__))
MT = os.path.dirname(HERE)
DESIGN = os.path.join(HERE, "p3_sigec_design.json")
reg = open(os.path.join(HERE, "PREREGISTERED_SHA256")).read().split()[0]
if hashlib.sha256(open(DESIGN, "rb").read()).hexdigest() != reg:
    sys.exit("p3_sigec_design.json changed after registration")
D = json.load(open(DESIGN))
sys.path.insert(0, os.path.join(MT, "nano_sigec"))
os.environ.setdefault("SIGEC_M1_SUMMARY", os.path.join(MT, "nano_sigec", "m1c", "m1c_summary.json"))
import analyze_sigec as A  # noqa: E402  (checks the SiGeC design hash on import)

F = D["factors"]
NAMES = list(F)
LEV = [F[k] for k in NAMES]
NL = [len(l) for l in LEV]
PAIRS = list(itertools.combinations(range(len(NAMES)), 2))
FULL = list(itertools.product(*[range(n) for n in NL]))  # level indices
OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)


def form_tuple(name, contrast):
    if name.startswith("substitutional"):
        return ("sub", 0.0, None, None)
    R = float(name.split("_R")[1].replace("nm", "")) * 1e-9
    return ("SiC", None, R, contrast)


def response(contrast):
    modes = {x: A.SModes(x) for x in F["x_Ge"]}
    y = np.empty(len(FULL))
    for i, r in enumerate(FULL):
        x, yc, form, d, dop, T = (LEV[j][r[j]] for j in range(6))
        g = 0.0 if dop == "none" else A.AN.gamma_dop(x, dop.split("_")[0], float(dop.split("_")[1][:-4]), 0.0)
        y[i] = A.kap(modes[x], x, T, d * 1e-6, yc, form_tuple(form, contrast), g)
    return y


def features(rows):
    X = []
    for r in rows:
        x, yc, form, d, dop, T = (LEV[j][r[j]] for j in range(6))
        if yc == 0:
            sub, lr = 0.0, -1.0
        elif form.startswith("substitutional"):
            sub, lr = 1.0, 0.0
        else:
            sub, lr = 0.0, np.log10(float(form.split("_R")[1].replace("nm", "")))
        X.append([x, yc, sub, lr, np.log10(d), float(dop.startswith("B")), float(dop.startswith("P")), T / 1000])
    return np.array(X)


XF = features(FULL)
IDX = {r: i for i, r in enumerate(FULL)}
GROUPS = {"x_Ge": [0], "y_C_at_frac": [1], "carbon_form": [2, 3], "d_eff_um": [4], "dopant": [5, 6], "T_K": [7]}
CAND = [i for i, r in enumerate(FULL) if LEV[0][r[0]] <= 0.2 and LEV[5][r[5]] == 1100]


def pair_keys():
    return {(a, b, u, v) for a, b in PAIRS for u in range(NL[a]) for v in range(NL[b])}


KEYS = pair_keys()


def counts(rows):
    c = {}
    for r in rows:
        for a, b in PAIRS:
            k = (a, b, r[a], r[b])
            c[k] = c.get(k, 0) + 1
    return c


def coverage(rows):
    return len(set(counts(rows)) & KEYS) / len(KEYS)


def complete(rows):
    return coverage(rows) >= 1.0


def covering(mult, seed, base=None):
    rng = random.Random(seed)
    rows = list(base or [])
    full = list(FULL)
    while True:
        c = counts(rows)
        deficit = {k for k in KEYS if c.get(k, 0) < mult}
        if not deficit:
            return rows
        rng.shuffle(full)
        best, bs = None, -1
        for cand in full:
            s = sum(1 for a, b in PAIRS if (a, b, cand[a], cand[b]) in deficit)
            if s > bs:
                best, bs = cand, s
        rows.append(best)


def fit(rows, y):
    ii = [IDX[r] for r in rows]
    m = ExtraTreesRegressor(n_estimators=300, min_samples_leaf=1, random_state=0, n_jobs=8)
    m.fit(XF[ii], np.log(y[ii]))
    return m


def evaluate(rows, y, importance=False):
    m = fit(rows, y)
    p = m.predict(XF)
    out = {"n": len(rows), "coverage": coverage(rows),
           "log_rmse": float(np.sqrt(np.mean((p - np.log(y)) ** 2)))}
    c = CAND
    pick = c[int(np.argmin(p[c]))]
    out["regret"] = float(y[pick] / y[c].min() - 1)
    out["pick"] = {NAMES[j]: LEV[j][FULL[pick][j]] for j in range(6)}
    if importance:
        rng = np.random.default_rng(0)
        base = np.mean((p - np.log(y)) ** 2)
        imp = {}
        for g, cols in GROUPS.items():
            Xp = XF.copy()
            perm = rng.permutation(len(Xp))
            Xp[:, cols] = XF[perm][:, cols]
            imp[g] = float(np.mean((m.predict(Xp) - np.log(y)) ** 2) - base)
        out["importance"] = imp
    return out, m


def ofat():
    base = [F["x_Ge"].index(0.2), F["y_C_at_frac"].index(0.01), F["carbon_form"].index("SiC_R25nm"),
            F["d_eff_um"].index(0.3), 0, F["T_K"].index(700)]
    rows = [tuple(base)]
    for j in range(6):
        for l in range(NL[j]):
            if l != base[j]:
                r = list(base)
                r[j] = l
                rows.append(tuple(r))
    return rows


def med(v):
    return float(np.median(v))


def main():
    t0 = time.time()
    res = {"design_sha256": reg, "ledger_prereg_seq": 56}
    y = response(D["model_settings"]["SiC_contrast_primary"])
    yd = response(D["model_settings"]["SiC_contrast_check"])
    res["response"] = {"min": float(y.min()), "max": float(y.max()), "eval_seconds": time.time() - t0,
                       "true_min_candidate": {NAMES[j]: LEV[j][FULL[CAND[int(np.argmin(y[CAND]))]][j]] for j in range(6)},
                       "true_min_kappa": float(y[CAND].min())}
    np.save(os.path.join(OUT, "response_stiffness.npy"), y)

    # ---- designs
    ca1 = [covering(1, s) for s in range(20)]
    ref = min(range(20), key=lambda s: len(ca1[s]))
    ca2 = [covering(2, 1000 + s, ca1[s]) for s in range(20)]
    n1, n2 = len(ca1[ref]), len(ca2[ref])
    res["sizes"] = {"CA1_all": [len(r) for r in ca1], "CA2_all": [len(r) for r in ca2], "CA1_ref": n1, "CA2_ref": n2,
                    "lower_bound": sorted(NL)[-1] * sorted(NL)[-2], "OFAT": len(ofat()), "full": len(FULL)}
    json.dump({"CA1_ref": ca1[ref], "CA2_ref": ca2[ref], "OFAT": ofat()}, open(os.path.join(OUT, "designs.json"), "w"))
    rng = random.Random(7)
    rnd1 = [rng.sample(FULL, n1) for _ in range(200)]
    rnd2 = [rng.sample(FULL, n2) for _ in range(200)]

    ev = {}
    for name, sets in (("CA1", ca1), ("CA2", ca2), ("RND_nCA1", rnd1), ("RND_nCA2", rnd2)):
        rows = [evaluate(s, y)[0] for s in sets]
        ev[name] = {"coverage_median": med([r["coverage"] for r in rows]),
                    "log_rmse_median": med([r["log_rmse"] for r in rows]),
                    "regret_median": med([r["regret"] for r in rows]),
                    "regret_q25_q75": [float(np.quantile([r["regret"] for r in rows], q)) for q in (0.25, 0.75)],
                    "log_rmse_q25_q75": [float(np.quantile([r["log_rmse"] for r in rows], q)) for q in (0.25, 0.75)]}
        if name.startswith("CA"):
            dens = [evaluate(s, yd)[0] for s in sets]
            ev[name]["density_contrast_regret_median"] = med([r["regret"] for r in dens])
    o, _ = evaluate(ofat(), y)
    ev["OFAT"] = {"coverage": o["coverage"], "log_rmse": o["log_rmse"], "regret": o["regret"], "pick": o["pick"]}
    ref2, mref = evaluate(ca2[ref], y, importance=True)
    ev["CA2_ref"] = ref2
    res["designs"] = ev

    # ---- erasure
    er = {}
    rng = random.Random(11)
    for name, rows in (("CA1", ca1[ref]), ("CA2", ca2[ref]), ("RND_nCA2", rnd2[0])):
        er[name] = {}
        for k in (1, 2, 3):
            ok, reg_k = 0, []
            for t in range(500):
                keep = rows[:]
                for i in sorted(rng.sample(range(len(keep)), k), reverse=True):
                    keep.pop(i)
                ok += complete(keep)
                if t < 60:
                    reg_k.append(evaluate(keep, y)[0]["regret"])
            er[name][f"k{k}"] = {"complete_frac": ok / 500, "regret_median_60trials": med(reg_k)}
    res["erasure"] = er

    # ---- transversal follow-up after CA2_ref
    base = ca2[ref]
    trees = np.stack([t.predict(XF) for t in mref.estimators_])
    sd = trees.std(0)
    mat = [0, 1, 2, 3]
    cells = {}
    for i, r in enumerate(FULL):
        for a, b in itertools.combinations(mat, 2):
            if 2 in (a, b) and LEV[1][r[1]] == 0:
                continue
            cells.setdefault((a, b, r[a], r[b]), []).append(i)
    cu = {k: float(sd[v].mean()) for k, v in cells.items()}
    cut = np.quantile(list(cu.values()), 0.8)
    unc = {k for k, v in cu.items() if v >= cut}
    hit, fu = set(), []
    have = set(base)
    while hit != unc:
        best, bs = None, (-1, -1.0)
        for i, r in enumerate(FULL):
            if r in have:
                continue
            h = {k for k in unc - hit if k in cells and k[0] in mat and r[k[0]] == k[2] and r[k[1]] == k[3]
                 and not (2 in (k[0], k[1]) and LEV[1][r[1]] == 0)}
            s = (len(h), sd[i])
            if s > bs:
                best, bs = (i, h), s
        fu.append(FULL[best[0]])
        have.add(FULL[best[0]])
        hit |= best[1]
    m_ = len(fu)
    tv, _ = evaluate(base + fu, y)
    rng = random.Random(13)
    pool = [r for r in FULL if r not in set(base)]
    rr = [evaluate(base + rng.sample(pool, m_), y)[0] for _ in range(200)]
    order = [FULL[i] for i in np.argsort(-sd) if FULL[i] not in set(base)][:m_]
    tv_top, _ = evaluate(base + order, y)
    res["transversal"] = {"n_uncertain_cells": len(unc), "n_cells": len(cu), "m": m_,
                          "uncertain_cells_by_pair": {f"{NAMES[a]}|{NAMES[b]}": sum(1 for k in unc if k[:2] == (a, b))
                                                      for a, b in itertools.combinations(mat, 2)},
                          "base": {k: ref2[k] for k in ("log_rmse", "regret")},
                          "transversal": {k: tv[k] for k in ("log_rmse", "regret")},
                          "random_median": {"log_rmse": med([r["log_rmse"] for r in rr]), "regret": med([r["regret"] for r in rr])},
                          "top_sd": {k: tv_top[k] for k in ("log_rmse", "regret")},
                          "followup_points": [{NAMES[j]: LEV[j][r[j]] for j in range(6)} for r in fu]}

    # ---- predictions
    P = {}
    P["P3-1"] = {"CA1": n1, "rnd_cov_median": ev["RND_nCA1"]["coverage_median"],
                 "pass": n1 <= 40 and ev["RND_nCA1"]["coverage_median"] < 0.9}
    P["P3-2"] = {"CA1": ev["CA1"]["regret_median"], "RND": ev["RND_nCA1"]["regret_median"],
                 "pass": ev["CA1"]["regret_median"] <= ev["RND_nCA1"]["regret_median"]}
    P["P3-3"] = {"CA2_k1": er["CA2"]["k1"]["complete_frac"], "CA1_k1": er["CA1"]["k1"]["complete_frac"],
                 "pass": er["CA2"]["k1"]["complete_frac"] == 1.0 and er["CA1"]["k1"]["complete_frac"] < 0.5}
    P["P3-4"] = {"ratio": n2 / n1, "pass": 1.6 <= n2 / n1 <= 2.2}
    P["P3-5"] = {"transversal": tv["log_rmse"], "random_median": res["transversal"]["random_median"]["log_rmse"],
                 "pass": tv["log_rmse"] < res["transversal"]["random_median"]["log_rmse"]}
    imp = ref2["importance"]
    P["P3-6"] = {"importance": imp, "ranking": sorted(imp, key=imp.get, reverse=True),
                 "pass": min(imp, key=imp.get) == "dopant"}
    P["P3-7"] = {"OFAT": o["regret"], "RND_nCA1_median": ev["RND_nCA1"]["regret_median"],
                 "pass": o["regret"] > ev["RND_nCA1"]["regret_median"]}
    res["predictions"] = P
    res["wall_seconds"] = time.time() - t0
    json.dump(res, open(os.path.join(OUT, "p3_sigec_results.json"), "w"), indent=1, default=str)
    print(json.dumps({k: v["pass"] for k, v in P.items()}), res["sizes"], res["wall_seconds"])


if __name__ == "__main__":
    main()
