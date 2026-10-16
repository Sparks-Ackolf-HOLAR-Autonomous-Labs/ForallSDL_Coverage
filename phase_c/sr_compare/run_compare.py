#!/usr/bin/env python3
"""Run every available SR engine on every benchmark task with one
protocol and write results/results_<engine>_<task>.json.

Protocol per (engine, task):
  - full fit: expression, engine complexity, wall time, in-sample errors;
  - determinism: fit twice with identical settings, compare expressions;
  - 5-fold CV (deterministic round-robin folds): out-of-sample errors.
Errors are on kappa (not r): RMSE (dominated by the pure endpoints in
the md task, where kappa spans 0.7 to 35 W/mK) and median absolute
percentage error (robust, scale-free).

Usage: run_compare.py [engine ...]   (default: all available)
"""
import json
import os
import sys
import time

import numpy as np

import engines
import hgsr
from datasets import TASKS

HERE = os.path.dirname(os.path.abspath(__file__))
K = 5


def errors(kappa_true, kappa_pred):
    t, p = np.asarray(kappa_true), np.asarray(kappa_pred)
    ok = np.isfinite(p)
    if not ok.all():
        p = np.where(ok, p, 0.0)      # a non-finite prediction is just wrong
    ape = np.abs(p - t) / t
    return {"rmse_kappa": float(np.sqrt(np.mean((p - t) ** 2))),
            "medape_pct": float(100 * np.median(ape)),
            "n_nonfinite_pred": int((~ok).sum())}


def folds(n, k=K):
    return [list(range(i, n, k)) for i in range(k)]


def run(engine_name, task_name):
    ent = engines.REGISTRY[engine_name]
    task = TASKS[task_name]()
    rows, feats = task["rows"], task["features"]
    if not ent["available"]:
        return {"engine": engine_name, "task": task_name,
                "available": False, "reason": ent["reason"]}
    eng = ent["cls"]()

    t0 = time.time()
    m1 = eng.fit(rows, feats)
    fit_s = time.time() - t0
    m2 = eng.fit(rows, feats)
    expr1, expr2 = eng.describe(m1), eng.describe(m2)

    ins = errors([r["kappa"] for r in rows], eng.predict(m1, rows))
    yt, yp = [], []
    for hold in folds(len(rows)):
        train = [r for i, r in enumerate(rows) if i not in set(hold)]
        test = [rows[i] for i in hold]
        mf = eng.fit(train, feats)
        yt += [r["kappa"] for r in test]
        yp += list(eng.predict(mf, test))
    cv = errors(yt, yp)

    res = {"engine": engine_name, "task": task_name, "available": True,
           "n_rows": len(rows), "expression": expr1,
           "complexity": eng.complexity(m1),
           "fit_seconds": round(fit_s, 2),
           "deterministic": expr1 == expr2,
           "expression_second_run": None if expr1 == expr2 else expr2,
           "in_sample": ins, "cv_5fold": cv}
    if engine_name == "hgsr":
        res["provenance"] = hgsr.provenance(rows, feats)
    return res


def main():
    want = sys.argv[1:] or [n for n, e in engines.REGISTRY.items()]
    for en in want:
        for tn in TASKS:
            res = run(en, tn)
            out = os.path.join(HERE, "results", f"results_{en}_{tn}.json")
            with open(out, "w") as f:
                json.dump(res, f, indent=1)
            if res.get("available"):
                print(f"{en:14s} {tn:11s} cv-MedAPE "
                      f"{res['cv_5fold']['medape_pct']:6.1f}%  rmse "
                      f"{res['cv_5fold']['rmse_kappa']:8.3f}  "
                      f"det={res['deterministic']}  {res['expression'][:70]}")
            else:
                print(f"{en:14s} {tn:11s} UNAVAILABLE: {res['reason'][:80]}")


if __name__ == "__main__":
    main()
