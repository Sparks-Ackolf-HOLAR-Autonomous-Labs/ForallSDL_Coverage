#!/usr/bin/env python3
"""Extended basis (squares + three-factor products, 43 terms on the MD
features) through the same comparison as full_v2.py (2026-10-02).

Engines: the adopted HG-SR v2 (shape filter + noise weights + beam 200)
with the standard generated basis vs the same with extended=True.
Yardsticks, unchanged from full_v2.py: md and md-alloy -> 5-fold CV over
rows, leave-one-composition/size-out (LOCO), external hold-out of 111
unseen production runs, plus the L -> infinity readout at 300 K; mace ->
the 16 hidden points (ledger seq 65). Then the four HG-SR properties on
md-alloy for the extended engine.

Disclosure: the md-alloy hold-out number for the extended basis (3.5%)
was seen earlier, in the full v2 study (adopted_v2.json, check B). CV, LOCO, the md task,
the mace hidden set and the properties are new evidence; the hold-out on
md-alloy is not blind for this choice.

Adoption rule (stated before this script ran): extended replaces standard
only if it is no worse on LOCO and on the hold-out for BOTH md tasks, no
worse on the mace hidden log-RMSE by more than the 0.08 seed scatter, and
keeps all four properties.
Usage: ~/mat_theory_pilot/venv/bin/python extended_compare.py
"""
import json
import os
import time

import numpy as np

import hgsr2
from full_v2 import (HG, MACE_PHYS, MD_PHYS, cv, hidden_scores, mace_prims,
                     mace_rows)
from datasets import TASKS
from improve_hgsr import alloy_task, holdout_rows, loco_folds
from run_compare import errors, folds

HERE = os.path.dirname(os.path.abspath(__file__))
ADOPT = dict(noise=True, search="beam", beam=200)


def law(m):
    s = m["selected"]
    return dict(zip(s["names"], s["coefs"]))


def md_block(name, t, hold):
    rows, feats = t["rows"], t["features"]
    out = []
    for ext in (False, True):
        e = HG(f"adopted v2{' + extended' if ext else ''}", physics=MD_PHYS,
               extended=ext, **ADOPT)
        t0 = time.time()
        m = e.fit(rows, feats)
        inf = e.predict(m, [{"x": x, "T": 300.0, "L": np.inf} for x in (0.15, 0.5, 0.85)])
        rec = {"engine": e.name, "n_basis": m["n_basis"], "law": law(m),
               "fit_s": round(time.time() - t0, 2),
               "cv5_rows": cv(e, rows, feats, folds(len(rows))),
               "loco": cv(e, rows, feats, loco_folds(rows)),
               "holdout": {k: round(v, 4) for k, v in errors(
                   [h["kappa"] for h in hold], e.predict(m, hold)).items()},
               "kappa_Linf_300K_x0.15_0.5_0.85": [round(float(v), 3) for v in inf]}
        out.append(rec)
        print(f"{name:9s} {e.name:24s} cv5 {rec['cv5_rows']['medape_pct']:5.1f}% | loco "
              f"{rec['loco']['medape_pct']:5.1f}% | hold {rec['holdout']['medape_pct']:5.1f}% "
              f"| Linf {rec['kappa_Linf_300K_x0.15_0.5_0.85']} | {list(rec['law'])}", flush=True)
    return out


def mace_block():
    allr = mace_rows()
    tr = [r for r in allr if r["split"] == "x"]
    hid = [r for r in allr if r["split"] == "H"]
    feats = ["x", "T", "c", "cc", "B", "P"]
    out = []
    for ext in (False, True):
        e = HG(f"adopted v2{' + extended' if ext else ''}", physics=MACE_PHYS,
               prims=mace_prims(), extended=ext, **ADOPT)
        m = e.fit(tr, feats)
        rec = {"engine": e.name, "n_basis": m["n_basis"], "law": law(m),
               "hidden": hidden_scores([r["kappa"] for r in hid], e.predict(m, hid))}
        out.append(rec)
        print("mace", e.name, rec["hidden"], list(rec["law"]), flush=True)
    return out


def properties(rows, F):
    KW = dict(physics=MD_PHYS, extended=True, **ADOPT)
    a, b = hgsr2.fit(rows, F, **KW), hgsr2.fit(rows, F, **KW)
    c0 = hgsr2.cell(rows[0])
    s = hgsr2.fit(rows + [dict(r) for r in rows if hgsr2.cell(r) == c0] * 4, F, **KW)
    bias = [dict(r, kappa=r["kappa"] * 1.4, party="BIASED-MD", source="b" + r["source"])
            for r in rows if r["T"] > 1000]
    pv = hgsr2.provenance(rows + bias, F, **KW)
    u = hgsr2.fit([r for r in rows + bias if r["party"] != "BIASED-MD"], F, **KW)
    la = law(a)
    same = list(law(s)) == list(la)
    return {"clean_law": la, "deterministic": la == law(b),
            "replication_same_terms": same,
            "replication_max_rel_coef_drift": max(abs(law(s)[k] - v) / abs(v)
                                                  for k, v in la.items()) if same else None,
            "biased_law_terms": pv["selected"]["names"],
            "blast_radius_biased_party": pv["blast_radius"].get("party:BIASED-MD"),
            "unlearn_bitwise_exact": law(u) == la}


def main():
    hold = holdout_rows()
    al = alloy_task()
    out = {"config": "adopted v2 (shape + noise + beam 200), standard vs extended basis",
           "md": md_block("md", TASKS["md"](), hold),
           "md-alloy": md_block("md-alloy", al, hold),
           "mace": mace_block(),
           "properties_extended_md-alloy": properties(al["rows"], al["features"])}
    print(json.dumps(out["properties_extended_md-alloy"], indent=1))
    s, x = {}, {}
    for t in ("md", "md-alloy"):
        s[t], x[t] = out[t]
    ok = all(x[t]["loco"]["medape_pct"] <= s[t]["loco"]["medape_pct"] and
             x[t]["holdout"]["medape_pct"] <= s[t]["holdout"]["medape_pct"]
             for t in ("md", "md-alloy"))
    ok &= out["mace"][1]["hidden"]["log_rmse"] <= out["mace"][0]["hidden"]["log_rmse"] + 0.08
    p = out["properties_extended_md-alloy"]
    ok &= p["deterministic"] and p["replication_same_terms"] and p["unlearn_bitwise_exact"] \
        and (p["replication_max_rel_coef_drift"] or 0) < 1e-9
    # fourth property, localized blast radius (standard engine: exactly one
    # term, adopted_v2.json C). Missing from the first run of this script,
    # which printed adopt_extended: True; fixed before any write-up.
    ok &= len(p["blast_radius_biased_party"] or []) <= 1
    out["adopt_extended"] = bool(ok)
    print("adopt_extended:", out["adopt_extended"])
    json.dump(out, open(os.path.join(HERE, "results", "extended_compare.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
