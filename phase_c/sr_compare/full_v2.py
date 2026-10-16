#!/usr/bin/env python3
"""HG-SR full v2: the five improvement items of HGSR_V2.md applied on top
of the v2 of HGSR_V2.md, scored on the same yardsticks plus a second system
(2026-10-02). Writes results/full_v2_<task>.json.

Ladder (cumulative):
  v2              HGSR_V2.md (F1-F4)
  +I1 shape       physics admissibility filter inside the declared domain
                  (r > 0; a bigger box never conducts worse; Nordheim peak
                  near x = 0.5); see hgsr2.physics_domain
  +I4 noise       per-cell weights from seed scatter (distinct values only)
  +I2 proposals   PySR hall-of-fame terms harvested per training fold
  full v2         + I3 deterministic beam search (beam 200) instead of
                  exhaustive enumeration; equality with exhaustive recorded
Rejected I1 variants (scored on the external hold-out only, for the
record): term-wise Matthiessen signs (NNLS) and a demanded bulk
(L -> infinity) limit.

Tasks:
  md, md-alloy   as in improve_hgsr.py: cv5, LOCO, external hold-out of
                 111 unseen production runs.
  mace           I5, a second system and method family: SiGeC + dopant,
                 MACE-MP MD (P3's hidden-set study). Train on the 71
                 planner-design points, score on the 16 hidden points
                 drawn from the ledger hash (seq 65), which no engine and
                 no design choice here has seen.

Usage: full_v2.py md|md-alloy|mace
"""
import glob
import json
import os
import sys
import time

import numpy as np

import hgsr2
import propose
from datasets import TASKS, regime
from improve_hgsr import alloy_task, holdout_rows, loco_folds
from run_compare import errors, folds

HERE = os.path.dirname(os.path.abspath(__file__))
MACE = os.path.join(HERE, "..", "..", "p3_mace")

# ---------------------------------------------------------------- domains
MD_PHYS = hgsr2.physics_domain(x=(0.05, 0.95), T=(300, 1300), L=(130, 800),
                               bulk=False)
MD_PHYS_BULK = hgsr2.physics_domain(x=(0.05, 0.95), T=(300, 1300),
                                    L=(130, 800), bulk=True)
MD_VARS = {"x": lambda r: r["x"], "t": lambda r: r["T"] / 1000.0,
           "il": lambda r: (100.0 / r["L"] if 0 < r["L"] < np.inf else 0.0)}

MACE_EXTRA = [{"c": 0.0, "cc": 0.0, "B": 0.0, "P": 0.0},
              {"c": 3.0, "cc": 3.0, "B": 1.0, "P": 0.0},
              {"c": 3.0, "cc": 0.0, "B": 0.0, "P": 1.0}]
MACE_PHYS = hgsr2.physics_domain(x=(0.1, 0.7), T=(300, 1100), L=None,
                                 extra=MACE_EXTRA)
MACE_VARS = {"x": lambda r: r["x"], "t": lambda r: r["T"] / 1000.0,
             "c": lambda r: r["c"], "cc": lambda r: r["cc"],
             "B": lambda r: r["B"], "P": lambda r: r["P"]}


def mace_prims():
    return [("x", lambda d: d["x"]),
            ("(1-x)", lambda d: 1 - d["x"]),
            ("x(1-x)", lambda d: d["x"] * (1 - d["x"])),
            ("t", lambda d: d["T"] / 1000.0),
            ("1/t", lambda d: 1000.0 / d["T"]),
            ("c", lambda d: d["c"]),          # carbon, at%
            ("cc", lambda d: d["cc"]),        # carbon in a SiC cluster, at%
            ("B", lambda d: d["B"]),          # 0.2 at% boron present
            ("P", lambda d: d["P"])]          # 0.2 at% phosphorus present


hgsr2.SKIP |= {frozenset(("c", "cc")), frozenset(("B", "P"))}


def mace_rows():
    out = []
    for f in sorted(glob.glob(os.path.join(MACE, "raw", "out", "*.json"))):
        j = json.load(open(f))
        pid = j["pid"]
        x, T = float(j["x_Ge"]), float(j["T_K"])
        c = 100.0 * float(j["y_C"])
        out.append({"x": x, "T": T, "L": 0.0, "c": c,
                    "cc": c if j["carbon_form"] == "clustered" else 0.0,
                    "B": 1.0 if j["dopant"] == "B" else 0.0,
                    "P": 1.0 if j["dopant"] == "P" else 0.0,
                    "kappa": float(j["kappa_W_mK"]), "party": "MACE-MD",
                    "regime": regime(x, T), "source": pid,
                    "cell": pid.split("_", 1)[1] if pid[:2] in ("H_", "R_") else pid,
                    "split": pid[0]})
    return out


# ---------------------------------------------------------------- engines
class HG:
    def __init__(self, name, prop=False, vars_=None, phys_dom=None, **kw):
        self.name, self.kw, self.prop = name, kw, prop
        self.vars, self.phys_dom = vars_, phys_dom
        self.cache = {}

    def fit(self, rows, feats):
        kw = dict(self.kw)
        if self.prop:
            dom = self.phys_dom["box"] + self.phys_dom["inf"]
            kw["proposals"] = propose.proposals_for(
                rows, feats, self.vars, dom,
                basis_kw={"prims": kw.get("prims")}, cache=self.cache)
        return hgsr2.fit(rows, feats, **kw)

    def predict(self, m, rows):
        return hgsr2.predict(m, rows)

    def describe(self, m):
        s = m["selected"]
        return "r = " + " + ".join(f"{c:.4g}*{n}" for n, c in zip(s["names"], s["coefs"]))


class ExtraTrees:
    """P3's registered surrogate (ExtraTrees, 300 trees, seed 0) on
    [x, y_C, is_clustered, is_B, is_P, T/1000], target ln kappa."""
    name = "extratrees"

    def _X(self, rows):
        return np.array([[r["x"], r["c"] / 100.0, 1.0 if r["cc"] > 0 else 0.0,
                          r["B"], r["P"], r["T"] / 1000.0] for r in rows])

    def fit(self, rows, feats):
        from sklearn.ensemble import ExtraTreesRegressor
        m = ExtraTreesRegressor(n_estimators=300, random_state=0)
        m.fit(self._X(rows), np.log([r["kappa"] for r in rows]))
        return m

    def predict(self, m, rows):
        return np.exp(m.predict(self._X(rows)))

    def describe(self, m):
        return "ExtraTrees(300)"


def ladder(base, vars_, dom):
    v2 = dict(base)
    s1 = dict(v2, physics=dom)
    s4 = dict(s1, noise=True)
    return [HG("v2", **v2),
            HG("+I1 shape filter", **s1),
            HG("+I4 noise weights", **s4),
            HG("+I2 PySR proposals", prop=True, vars_=vars_, phys_dom=dom, **s4),
            HG("full v2 (+I3 beam)", prop=True, vars_=vars_, phys_dom=dom,
               search="beam", beam=200, **s4)]


def cv(eng, rows, feats, fl):
    yt, yp = [], []
    for hold in fl:
        hs = set(hold)
        tr = [r for i, r in enumerate(rows) if i not in hs]
        m = eng.fit(tr, feats)
        yt += [rows[i]["kappa"] for i in hold]
        yp += list(eng.predict(m, [rows[i] for i in hold]))
    return {k: round(v, 4) for k, v in errors(yt, yp).items()}


def hidden_scores(yt, yp):
    from scipy.stats import spearmanr
    yt, yp = np.asarray(yt), np.asarray(yp)
    e = errors(yt, yp)
    ok = np.isfinite(yp) & (yp > 0)
    lr = float(np.sqrt(np.mean((np.log(np.where(ok, yp, 1e-9)) - np.log(yt)) ** 2)))
    pick = int(np.argmin(np.where(ok, yp, np.inf)))
    return {"log_rmse": round(lr, 4), "medape_pct": round(e["medape_pct"], 3),
            "spearman": round(float(spearmanr(yt, yp).correlation), 4),
            "pick_regret": round(float(yt[pick] / yt.min() - 1), 4)}


def exhaustive_match(eng, rows, feats):
    """I3 check: does beam search select exactly what exhaustive does?"""
    kw = dict(eng.kw)
    kw.pop("search", None)
    kw.pop("beam", None)
    if eng.prop:
        dom = eng.phys_dom["box"] + eng.phys_dom["inf"]
        kw["proposals"] = propose.proposals_for(
            rows, feats, eng.vars, dom, basis_kw={"prims": kw.get("prims")},
            cache=eng.cache)
    t0 = time.time()
    a = hgsr2.fit(rows, feats, **kw)
    ta = time.time() - t0
    t0 = time.time()
    b = hgsr2.fit(rows, feats, search="beam", beam=200, **kw)
    tb = time.time() - t0
    return {"same_selected": a["selected"]["names"] == b["selected"]["names"]
            and np.allclose(a["selected"]["coefs"], b["selected"]["coefs"], rtol=1e-9),
            "same_front": [m["names"] for m in a["front"]] == [m["names"] for m in b["front"]],
            "n_basis": a["n_basis"], "evaluated_exhaustive": a["n_evaluated"],
            "evaluated_beam": b["n_evaluated"], "s_exhaustive": round(ta, 2),
            "s_beam": round(tb, 2)}


def run_md(task_name):
    t = TASKS["md"]() if task_name == "md" else alloy_task()
    rows, feats = t["rows"], t["features"]
    hold = holdout_rows()
    engs = ladder({}, MD_VARS, MD_PHYS)
    out = {"task": task_name, "n": len(rows), "n_holdout": len(hold), "engines": []}
    for e in engs:
        t0 = time.time()
        m = e.fit(rows, feats)
        rec = {"engine": e.name, "expression": e.describe(m),
               "fit_s": round(time.time() - t0, 2), "n_basis": m["n_basis"],
               "proposals": [b[0] for b in m["basis"] if b[0].startswith("pysr:")],
               "cv5_rows": cv(e, rows, feats, folds(len(rows))),
               "loco": cv(e, rows, feats, loco_folds(rows)),
               "holdout": {k: round(v, 4) for k, v in errors(
                   [h["kappa"] for h in hold], e.predict(m, hold)).items()}}
        inf = e.predict(m, [{"x": x, "T": 300.0, "L": np.inf} for x in (0.15, 0.5, 0.85)])
        rec["kappa_Linf_300K_x0.15_0.5_0.85"] = [round(float(v), 3) for v in inf]
        out["engines"].append(rec)
        print(f"{task_name:9s} {e.name:24s} cv5 {rec['cv5_rows']['medape_pct']:6.1f}% "
              f"| loco {rec['loco']['medape_pct']:6.1f}% | hold "
              f"{rec['holdout']['medape_pct']:6.1f}% {rec['holdout']['rmse_kappa']:.3f} "
              f"| {rec['expression'][:90]}", flush=True)
    out["I3_beam_vs_exhaustive"] = exhaustive_match(engs[-1], rows, feats)
    # rejected I1 variants, external hold-out only
    rej = []
    for name, kw in [("I1a term-wise signs (NNLS)", dict(constrain=True, require=("x(1-x)",))),
                     ("I1b demanded bulk limit", dict(physics=MD_PHYS_BULK))]:
        m = hgsr2.fit(rows, feats, **kw)
        inf = hgsr2.predict(m, [{"x": x, "T": 300.0, "L": np.inf} for x in (0.15, 0.5, 0.85)])
        rej.append({"variant": name, "expression": HG("").describe(m),
                    "holdout": {k: round(v, 4) for k, v in errors(
                        [h["kappa"] for h in hold], hgsr2.predict(m, hold)).items()},
                    "kappa_Linf_300K_x0.15_0.5_0.85": [round(float(v), 3) for v in inf]})
        print("rejected", name, rej[-1]["holdout"]["medape_pct"], rej[-1]["kappa_Linf_300K_x0.15_0.5_0.85"])
    out["rejected_I1_variants"] = rej
    return out


def run_mace():
    allr = mace_rows()
    train = [r for r in allr if r["split"] == "x"]
    hidden = [r for r in allr if r["split"] == "H"]
    feats = ["x", "T", "c", "cc", "B", "P"]
    base = {"prims": mace_prims()}
    engs = ladder(base, MACE_VARS, MACE_PHYS)
    import engines as E
    others = [("pysr", E.REGISTRY["pysr"]), ("sisso", E.REGISTRY["sisso"])]
    out = {"task": "mace", "n_train": len(train), "n_hidden": len(hidden),
           "engines": []}
    yt = [r["kappa"] for r in hidden]
    allengs = [(e.name, e) for e in engs] + [("extratrees", ExtraTrees())]
    allengs += [(n, ent["cls"]()) for n, ent in others if ent["available"]]
    for name, e in allengs:
        t0 = time.time()
        try:
            m = e.fit(train, feats)
            yp = e.predict(m, hidden)
            rec = {"engine": name, "expression": e.describe(m),
                   "fit_s": round(time.time() - t0, 2),
                   "hidden": hidden_scores(yt, yp)}
            if isinstance(m, dict) and "n_basis" in m:
                rec["n_basis"] = m["n_basis"]
                rec["proposals"] = [b[0] for b in m["basis"] if b[0].startswith("pysr:")]
        except Exception as ex:                       # report, do not hide
            rec = {"engine": name, "error": repr(ex)[:300]}
        out["engines"].append(rec)
        h = rec.get("hidden", {})
        print(f"mace {name:24s} hidden log-RMSE {h.get('log_rmse')} MedAPE {h.get('medape_pct')} "
              f"rho {h.get('spearman')} regret {h.get('pick_regret')} | {rec.get('expression', rec.get('error'))[:90]}",
              flush=True)
    out["I3_beam_vs_exhaustive"] = exhaustive_match(engs[-1], train, feats)
    return out


if __name__ == "__main__":
    task = sys.argv[1]
    res = run_mace() if task == "mace" else run_md(task)
    with open(os.path.join(HERE, "results", f"full_v2_{task}.json"), "w") as f:
        json.dump(res, f, indent=1)
