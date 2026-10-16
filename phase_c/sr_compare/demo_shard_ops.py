#!/usr/bin/env python3
"""Demo: the non-forget shard operations on HG-SR, on the real md task
evidence plus the demo_bias biased party. Writes results/demo_shard_ops.json.

Answers the question (2026-09-26): how does the
hypergraph machinery do the SISA-like process, and what can it do
besides unlearning? Each block below is one operator from the U family
executed through shard_ops.apply, with the belief change reported by
the same clean-retrain comparison that verifies forget.
"""
import copy
import json
import os

import numpy as np

import hgsr
import shard_ops
from datasets import TASKS

HERE = os.path.dirname(os.path.abspath(__file__))


def taint(rows):
    rng = np.random.default_rng(7)
    biased = []
    for r in rows:
        if r["regime"] == "x3T3":
            for _ in range(2):
                b = copy.deepcopy(r)
                b["kappa"] = round(b["kappa"] * 1.4 *
                                   (1 + 0.02 * rng.standard_normal()), 4)
                b["party"] = "BIASED-MD"
                b["source"] = "biased-" + b["source"]
                biased.append(b)
    return rows + biased


def sel(model):
    s = model["selected"]
    return dict(zip(s["names"], [round(c, 6) for c in s["coefs"]]))


def main():
    task = TASKS["md"]()
    rows, feats = task["rows"], task["features"]
    clean = hgsr.fit(rows, feats)
    tainted = taint(rows)
    taint_fit = hgsr.fit(tainted, feats)
    out = {"clean_model": sel(clean), "tainted_model": sel(taint_fit)}

    # 1. quarantine: same model as forget, but reversible
    q = shard_ops.apply(tainted, feats, "quarantine", edge="party:BIASED-MD")
    forgot = shard_ops.apply(tainted, feats, "forget", edge="party:BIASED-MD")
    out["quarantine"] = {
        "model": sel(q["model"]),
        "equals_forget_model": sel(q["model"]) == sel(forgot["model"]),
        "equals_clean_bitwise": (
            q["model"]["selected"]["names"] == clean["selected"]["names"]
            and np.allclose(q["model"]["selected"]["coefs"],
                            clean["selected"]["coefs"], rtol=0, atol=0)),
    }
    # restore lifts it: back to the tainted fit bitwise (nothing was lost)
    r = shard_ops.apply(tainted, feats, "restore", edge="party:BIASED-MD")
    out["restore"] = {
        "equals_tainted_bitwise": (
            r["model"]["selected"]["names"] == taint_fit["selected"]["names"]
            and np.allclose(r["model"]["selected"]["coefs"],
                            taint_fit["selected"]["coefs"], rtol=0, atol=0))}

    # 2. downweight: continuous trust dial between retain and quarantine
    sweep = {}
    for g in (1.0, 0.5, 0.2, 0.05, 0.0):
        m = shard_ops.apply(tainted, feats, "downweight",
                            edge="party:BIASED-MD", gamma=g)["model"]
        sweep[str(g)] = {"terms": m["selected"]["names"],
                         "spurious_T_term": "T" in m["selected"]["names"]}
    out["downweight_sweep"] = sweep

    # 3. correct: the biased party's calibration is fixed (kappa / 1.4)
    # instead of thrown away; evidence is kept, bias removed
    fixed_rows = copy.deepcopy(tainted)
    for row in fixed_rows:
        if row["party"] == "BIASED-MD":
            row["kappa"] = round(row["kappa"] / 1.4, 4)
    corr = hgsr.fit(fixed_rows, feats)
    out["correct"] = {
        "model": sel(corr),
        "terms_match_clean": corr["selected"]["names"] == clean["selected"]["names"]}

    # 4. split: instead of unlearning, recognize the biased rows as their
    # own sub-regime; the balanced loss then caps their vote (a
    # world-is-richer response, evidence kept, no term deleted)
    def retag(row):
        return row["regime"] + ("b" if row["party"] == "BIASED-MD" else "")
    sp = shard_ops.apply(tainted, feats, "split", edge="regime:x3T3",
                         retag=retag)
    out["split"] = {
        "model": sel(sp["model"]),
        "delta_vs_tainted": sp["delta"],
        "spurious_T_term_after_split": "T" in sp["model"]["selected"]["names"]}

    # 5. replace: swap the biased rows for their corrected re-measurement
    repl = [dict(r, kappa=round(r["kappa"] / 1.4, 4), party="FIXED-MD")
            for r in tainted if r["party"] == "BIASED-MD"]
    rp = shard_ops.apply(tainted, feats, "replace", edge="party:BIASED-MD",
                         replacement=repl)
    out["replace"] = {
        "model": sel(rp["model"]),
        "terms_match_clean": rp["model"]["selected"]["names"] == clean["selected"]["names"]}

    with open(os.path.join(HERE, "results", "demo_shard_ops.json"), "w") as f:
        json.dump(out, f, indent=1)
    for k, v in out.items():
        print(f"== {k}: {json.dumps(v)[:160]}")


if __name__ == "__main__":
    main()
