#!/usr/bin/env python3
"""Ablation ladder: why HG-SR v1 trailed PySR/SISSO and what each fix
recovers (2026-10-02). Writes results/improve_hgsr.json.

Variants (cumulative, one fix at a time on top of v1):
  v1                 hgsr.py as benchmarked in COMPARISON.md (2026-09-27)
  +F1 logo-select    select on the front by leave-one-hyperedge-out CV
  +F2 cells          hyperedges are (x, T, size) design cells
  +F3 relative loss  percentage-error fitting
  +F4 generated      products of physical primitives incl. L^-1/2 (= v2)
plus PySR and SISSO with the original settings (2026-09-27).

Tasks: the original md task (24 rows, includes the ballistic pure
endpoints) and md-alloy (x in [0.05, 0.95], the ledgered scope-limit
of 2026-09-27). Protocols: the original 5-fold round-robin row CV,
and leave-one-design-cell-out (LOCO; seed replicates held out together,
the stricter test of extrapolation to an unseen composition-size cell).
"""
import json
import os
import sys
import time

import numpy as np

import engines
import hgsr
import hgsr2
from datasets import TASKS
from run_compare import errors, folds

HERE = os.path.dirname(os.path.abspath(__file__))


class V:
    def __init__(self, name, **kw):
        self.name, self.kw = name, kw

    def fit(self, rows, feats):
        if self.kw.get("v1"):
            return ("v1", hgsr.fit(rows, feats))
        return ("v2", hgsr2.fit(rows, feats, **self.kw))

    def predict(self, m, rows):
        return hgsr.predict(m[1], rows) if m[0] == "v1" else hgsr2.predict(m[1], rows)

    def describe(self, m):
        s = m[1]["selected"]
        return "r = " + " + ".join(f"{c:.4g}*{n}" for n, c in zip(s["names"], s["coefs"]))


VARIANTS = [
    V("hgsr-v1", v1=True),
    V("+F1 logo-select", generated=False, rel=False, size_aware=False),
    V("+F2 size cells", generated=False, rel=False, size_aware=True),
    V("+F3 relative loss", generated=False, rel=True, size_aware=True),
    V("+F4 generated basis (v2)", generated=True, rel=True, size_aware=True),
]


def alloy_task():
    t = TASKS["md"]()
    return {"name": "md-alloy", "features": t["features"],
            "rows": [r for r in t["rows"] if 0.05 <= r["x"] <= 0.95]}


def loco_folds(rows):
    g = {}
    for i, r in enumerate(rows):
        g.setdefault(hgsr2.cell(r), []).append(i)
    return [g[k] for k in sorted(g)]


def score(eng, rows, feats, fl):
    yt, yp = [], []
    for hold in fl:
        hs = set(hold)
        tr = [r for i, r in enumerate(rows) if i not in hs]
        m = eng.fit(tr, feats)
        yt += [rows[i]["kappa"] for i in hold]
        yp += list(eng.predict(m, [rows[i] for i in hold]))
    e = errors(yt, yp)
    return {k: round(v, 4) for k, v in e.items()}


def holdout_rows():
    """External hold-out: the 2026-09-28 production campaign (126
    standard-deck runs at 40 compositions/temperatures the benchmark
    never saw, same 4x4 / swap-50 protocol). Restricted to the law's
    declared domain: x in [0.05, 0.90], T below the near-melting probes."""
    pd = os.path.join(HERE, "..", "production")
    des = json.load(open(os.path.join(pd, "production_design.json")))
    out = []
    for d in des["runs"]:
        f = os.path.join(pd, "results", d["run_id"] + ".json")
        if d["deck"] != "std" or not os.path.exists(f):
            continue
        j = json.load(open(f))
        if not (0.05 <= d["x"] <= 0.90) or d["T"] >= 1400:
            continue
        x, T = float(d["x"]), float(d["T"])
        from datasets import regime
        out.append({"x": x, "T": T, "L": float(j["box_A"][2]),
                    "kappa": float(j["kappa_W_mK"]), "party": "CHPC-MD-C",
                    "regime": regime(x, T), "source": d["run_id"]})
    return out


def main():
    engs = [(v.name, v) for v in VARIANTS]
    want_ext = "--no-ext" not in sys.argv
    if want_ext:
        for en in ("pysr", "sisso"):
            ent = engines.REGISTRY[en]
            if ent["available"]:
                engs.append((en, ent["cls"]()))
    tasks = [TASKS["md"](), alloy_task()]
    hold = holdout_rows()
    print(f"external hold-out: {len(hold)} runs")
    out = []
    for t in tasks:
        rows, feats = t["rows"], t["features"]
        for name, e in engs:
            t0 = time.time()
            m = e.fit(rows, feats)
            fs = time.time() - t0
            rec = {"task": t["name"], "engine": name, "n": len(rows),
                   "expression": e.describe(m), "fit_s": round(fs, 2),
                   "cv5_rows": score(e, rows, feats, folds(len(rows))),
                   "loco": score(e, rows, feats, loco_folds(rows)),
                   "holdout": {k: round(v, 4) for k, v in errors(
                       [h["kappa"] for h in hold], e.predict(m, hold)).items()}}
            out.append(rec)
            print(f"{t['name']:9s} {name:26s} cv5 {rec['cv5_rows']['medape_pct']:6.1f}% "
                  f"{rec['cv5_rows']['rmse_kappa']:7.3f} | loco "
                  f"{rec['loco']['medape_pct']:6.1f}% {rec['loco']['rmse_kappa']:7.3f} "
                  f"| hold {rec['holdout']['medape_pct']:6.1f}% {rec['holdout']['rmse_kappa']:6.3f} "
                  f"| {rec['expression'][:60]}", flush=True)
    with open(os.path.join(HERE, "results", "improve_hgsr.json"), "w") as f:
        json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()
