#!/usr/bin/env python3
"""P3 replication-reserve study (reserve_design.json; ledger entry in PREREGISTERED_SHA256).

Refuses to run if the design or the response file changed after registration.
Writes results/reserve_results.json.
"""
import hashlib, json, os, sys, time
import numpy as np
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.model_selection import KFold

HERE = os.path.dirname(os.path.abspath(__file__))
MT = os.path.dirname(HERE)
DESIGN = os.path.join(HERE, "reserve_design.json")
reg = open(os.path.join(HERE, "PREREGISTERED_SHA256")).read().split()[0]
if hashlib.sha256(open(DESIGN, "rb").read()).hexdigest() != reg:
    sys.exit("reserve_design.json changed after registration")
D = json.load(open(DESIGN))
RESP = os.path.join(MT, "p3_sigec", "results", "response_stiffness.npy")
if hashlib.sha256(open(RESP, "rb").read()).hexdigest() != D["response_sha256"]:
    sys.exit("response file changed")
sys.path.insert(0, os.path.join(MT, "p3_sigec"))
import run_p3_sigec as P  # noqa: E402  (checks the P3 design hash on import)

Y = np.load(RESP)
LNY = np.log(Y)
XF, FULL, IDX, CAND = P.XF, P.FULL, P.IDX, np.array(P.CAND)
CANDSET = set(CAND.tolist())
BASE = [IDX[tuple(r)] for r in json.load(open(os.path.join(MT, "p3_sigec", "results", "designs.json")))["CA2_ref"]]
TRUE_MIN = Y[CAND].min()
SIG, TAU, CONF, MAXREP, EXTRA = 0.05, 0.15, 0.141, 6, D["budget_extra"]
NT = 60


def model(idx, lny):
    m = ExtraTreesRegressor(n_estimators=NT, random_state=0, n_jobs=4)
    m.fit(XF[idx], lny)
    return m


def observe(i, rng, f):
    v = LNY[i] + rng.normal(0, SIG)
    fault = rng.random() < f
    if fault:
        v += np.log(rng.uniform(0.4, 0.7))
    return v, fault


def neighbours(i):
    r = FULL[i]
    out = []
    for j in range(len(r)):
        for l in range(P.NL[j]):
            if l != r[j]:
                q = list(r); q[j] = l
                out.append(IDX[tuple(q)])
    return out


def score(acc, faults):
    """acc: dict idx -> accepted ln value; faults: dict idx -> whether the accepted value is a fault."""
    idx = np.array(sorted(acc)); v = np.array([acc[i] for i in idx])
    m = model(idx, v)
    p = m.predict(XF)
    c = [i for i in idx if i in CANDSET]
    best = min(c, key=lambda i: acc[i]) if c else None
    out = {"log_rmse": float(np.sqrt(np.mean((p - LNY) ** 2))),
           "surrogate_regret": float(Y[CAND[int(np.argmin(p[CAND]))]] / TRUE_MIN - 1)}
    if best is None:
        out.update(false_claim=None, report_regret=None)
    else:
        out["false_claim"] = bool(np.exp(acc[best]) / Y[best] < 0.8)
        out["report_regret"] = float(Y[best] / TRUE_MIN - 1)
        out["reported_is_fault"] = bool(faults[best])
    return out


def trial(f, t):
    rng = np.random.default_rng(1000 * int(round(f * 100)) + t)
    base = {}; bfault = {}
    for i in BASE:
        base[i], bfault[i] = observe(i, rng, f)
    res = {}
    # policy streams use their own generators so extra runs do not shift each other
    m0 = model(np.array(BASE), np.array([base[i] for i in BASE]))
    p0 = m0.predict(XF)
    unmeasured = [i for i in range(len(FULL)) if i not in base]
    # CA2+random
    r1 = np.random.default_rng(50000 + t + int(f * 1e4))
    acc, fl = dict(base), dict(bfault)
    for i in r1.choice(unmeasured, EXTRA, replace=False):
        acc[i], fl[i] = observe(i, r1, f)
    res["CA2+random"] = score(acc, fl)
    # CA2+exploit
    r2 = np.random.default_rng(60000 + t + int(f * 1e4))
    acc, fl = dict(base), dict(bfault)
    cand_un = [i for i in CAND if i not in base]
    for i in sorted(cand_un, key=lambda i: p0[i])[:EXTRA]:
        acc[i], fl[i] = observe(i, r2, f)
    res["CA2+exploit"] = score(acc, fl)
    # CA2+reserve
    r3 = np.random.default_rng(70000 + t + int(f * 1e4))
    acc, fl = dict(base), dict(bfault)
    idx = np.array(BASE); v = np.array([base[i] for i in BASE])
    cv = np.empty(len(idx))
    for tr, te in KFold(5, shuffle=True, random_state=0).split(idx):
        cv[te] = model(idx[tr], v[tr]).predict(XF[idx[te]])
    resid = v - cv
    flagged = [int(idx[k]) for k in np.argsort(resid) if resid[k] < -TAU][:MAXREP]
    left = EXTRA
    n_rep = n_conf = n_aside = conf_true = aside_fault = 0
    probes = []
    for i in flagged:
        if left <= 0:
            break
        rep, rf = observe(i, r3, f); left -= 1; n_rep += 1
        if abs(rep - base[i]) <= CONF:
            n_conf += 1; conf_true += (not bfault[i])
            acc[i] = 0.5 * (rep + base[i]); fl[i] = bfault[i] and rf
            if left > 0:
                pm = model(np.array(sorted(acc)), np.array([acc[k] for k in sorted(acc)])).predict(XF)
                nb = [k for k in neighbours(i) if k not in acc]
                if nb:
                    k = min(nb, key=lambda k: pm[k])
                    acc[k], fl[k] = observe(k, r3, f); left -= 1; probes.append(k)
        else:
            n_aside += 1; aside_fault += bfault[i]
            acc[i], fl[i] = rep, rf
    if left > 0:
        pm = model(np.array(sorted(acc)), np.array([acc[k] for k in sorted(acc)])).predict(XF)
        for k in sorted([k for k in CAND if k not in acc], key=lambda k: pm[k])[:left]:
            acc[k], fl[k] = observe(k, r3, f)
    s = score(acc, fl)
    s.update(n_flagged=len(flagged), n_replicated=n_rep, n_confirmed=n_conf, n_set_aside=n_aside,
             n_confirmed_not_fault=int(conf_true), n_set_aside_fault=int(aside_fault),
             n_base_faults=int(sum(bfault.values())), n_base_faults_in_region=int(sum(bfault[i] for i in BASE if i in CANDSET)))
    res["CA2+reserve"] = s
    return res


def summarise(rows, pol):
    g = [r[pol] for r in rows]
    fc = [x["false_claim"] for x in g if x["false_claim"] is not None]
    out = {"false_claim_rate": float(np.mean(fc)) if fc else None,
           "report_regret_median": float(np.median([x["report_regret"] for x in g if x["report_regret"] is not None])),
           "surrogate_regret_median": float(np.median([x["surrogate_regret"] for x in g])),
           "log_rmse_median": float(np.median([x["log_rmse"] for x in g]))}
    if pol == "CA2+reserve":
        rep = sum(x["n_replicated"] for x in g)
        conf = sum(x["n_confirmed"] for x in g)
        out.update(flagged_mean=float(np.mean([x["n_flagged"] for x in g])),
                   replicated_total=rep, confirmed_total=conf,
                   set_aside_frac=(sum(x["n_set_aside"] for x in g) / rep) if rep else 0.0,
                   confirmed_true_frac=(sum(x["n_confirmed_not_fault"] for x in g) / conf) if conf else None,
                   set_aside_that_were_faults=sum(x["n_set_aside_fault"] for x in g),
                   trials_with_confirmed=float(np.mean([x["n_confirmed"] > 0 for x in g])),
                   base_faults_mean=float(np.mean([x["n_base_faults"] for x in g])))
    return out


def main():
    t0 = time.time()
    out = {"design_sha256": reg, "base_n": len(BASE), "base_in_region": sum(i in CANDSET for i in BASE),
           "true_min": float(TRUE_MIN), "by_f": {}, "trials": {}}
    for f in D["measurement_model"]["f_values"]:
        rows = [trial(f, t) for t in range(D["measurement_model"]["trials_per_f"])]
        out["trials"][str(f)] = rows
        out["by_f"][str(f)] = {p: summarise(rows, p) for p in ("CA2+random", "CA2+exploit", "CA2+reserve")}
        print(f, json.dumps(out["by_f"][str(f)]), time.time() - t0, flush=True)
    B = out["by_f"]
    e5, r5 = B["0.05"]["CA2+exploit"], B["0.05"]["CA2+reserve"]
    pr = {
        "PR1": {"exploit": e5["false_claim_rate"], "reserve": r5["false_claim_rate"],
                "pass": e5["false_claim_rate"] >= 0.15 and r5["false_claim_rate"] <= 0.5 * e5["false_claim_rate"]},
        "PR2": {"exploit": e5["report_regret_median"], "reserve": r5["report_regret_median"],
                "pass": r5["report_regret_median"] <= e5["report_regret_median"]},
        "PR3": {"random": B["0.0"]["CA2+random"]["log_rmse_median"], "reserve": B["0.0"]["CA2+reserve"]["log_rmse_median"],
                "pass": B["0.0"]["CA2+reserve"]["log_rmse_median"] <= 1.10 * B["0.0"]["CA2+random"]["log_rmse_median"]},
        "PR4": {"set_aside_frac": B["0.0"]["CA2+reserve"]["set_aside_frac"],
                "pass": B["0.0"]["CA2+reserve"]["set_aside_frac"] <= 0.10},
        "PR5": {"exploit": B["0.1"]["CA2+exploit"]["false_claim_rate"], "reserve": B["0.1"]["CA2+reserve"]["false_claim_rate"],
                "pass": B["0.1"]["CA2+reserve"]["false_claim_rate"] <= 0.7 * B["0.1"]["CA2+exploit"]["false_claim_rate"]},
    }
    out["predictions"] = pr
    out["wall_seconds"] = time.time() - t0
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    json.dump(out, open(os.path.join(HERE, "results", "reserve_results.json"), "w"), indent=1, default=str)
    print(json.dumps({k: v["pass"] for k, v in pr.items()}), out["wall_seconds"])


if __name__ == "__main__":
    main()
