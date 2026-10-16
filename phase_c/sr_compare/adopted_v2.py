#!/usr/bin/env python3
"""The adopted HG-SR v2 configuration after the improvement items
(2026-10-02): shape filter (I1) + noise weights (I4) + beam
search (I3); PySR proposals (I2) stay an optional discovery channel, off
by default, because on the external hold-out they raised the md-alloy
error from 4.3% to 11.6% (results/full_v2_md-alloy.json).

Checks, written to results/adopted_v2.json:
  A  beam (200) selects exactly what exhaustive search selects, per task;
  B  scaling: extended basis (squares and three-factor products) on the
     md-alloy task, exhaustive vs beam, at 4 and 5 terms;
  C  the four HG-SR properties on md-alloy (as demo_v2_props.py).
"""
import json
import os
import time

import hgsr2
from full_v2 import MACE_PHYS, MD_PHYS, mace_prims, mace_rows
from datasets import TASKS
from improve_hgsr import alloy_task, holdout_rows
from run_compare import errors

HERE = os.path.dirname(os.path.abspath(__file__))
ADOPT = dict(noise=True)


def law(m):
    s = m["selected"]
    return dict(zip(s["names"], s["coefs"]))


def match(rows, feats, **kw):
    t0 = time.time()
    a = hgsr2.fit(rows, feats, **kw)
    ta = time.time() - t0
    t0 = time.time()
    b = hgsr2.fit(rows, feats, search="beam", beam=200, **kw)
    tb = time.time() - t0
    return {"same_selected": law(a) == law(b), "n_basis": a["n_basis"],
            "evaluated_exhaustive": a["n_evaluated"], "evaluated_beam": b["n_evaluated"],
            "s_exhaustive": round(ta, 2), "s_beam": round(tb, 2), "law": law(b)}


out = {"config": "shape filter + noise weights + beam 200 (no PySR proposals)"}
md, al = TASKS["md"](), alloy_task()
mace = [r for r in mace_rows() if r["split"] == "x"]
out["A"] = {
    "md": match(md["rows"], md["features"], physics=MD_PHYS, **ADOPT),
    "md-alloy": match(al["rows"], al["features"], physics=MD_PHYS, **ADOPT),
    "mace": match(mace, ["x", "T", "c", "cc", "B", "P"], physics=MACE_PHYS,
                  prims=mace_prims(), **ADOPT),
}
print(json.dumps(out["A"], indent=1), flush=True)

hold = holdout_rows()
B = []
for k in (4, 5):
    rec = {"max_terms": k}
    for mode in ("beam", "exhaustive"):
        if mode == "exhaustive" and k == 5:
            rec["exhaustive"] = "skipped: est. > 1 h (subsets grow as n^5)"
            continue
        t0 = time.time()
        m = hgsr2.fit(al["rows"], al["features"], physics=MD_PHYS, extended=True,
                      max_terms=k, search=mode, beam=200, **ADOPT)
        rec[mode] = {"s": round(time.time() - t0, 1), "n_basis": m["n_basis"],
                     "evaluated": m["n_evaluated"], "law": law(m),
                     "holdout_medape": round(errors([h["kappa"] for h in hold],
                                                    hgsr2.predict(m, hold))["medape_pct"], 2)}
        print(k, mode, rec[mode], flush=True)
    if isinstance(rec.get("exhaustive"), dict):
        rec["same_selected"] = rec["exhaustive"]["law"] == rec["beam"]["law"]
    B.append(rec)
out["B"] = B

rows, F = al["rows"], al["features"]
KW = dict(physics=MD_PHYS, search="beam", beam=200, **ADOPT)
a, b = hgsr2.fit(rows, F, **KW), hgsr2.fit(rows, F, **KW)
c0 = hgsr2.cell(rows[0])
s = hgsr2.fit(rows + [dict(r) for r in rows if hgsr2.cell(r) == c0] * 4, F, **KW)
bias = [dict(r, kappa=r["kappa"] * 1.4, party="BIASED-MD", source="b" + r["source"])
        for r in rows if r["T"] > 1000]
pv = hgsr2.provenance(rows + bias, F, **KW)
u = hgsr2.fit([r for r in rows + bias if r["party"] != "BIASED-MD"], F, **KW)
la = law(a)
out["C"] = {"clean_law": la, "deterministic": la == law(b),
            "replication_same_terms": list(law(s)) == list(la),
            "replication_max_rel_coef_drift": max(abs(law(s)[k] - v) / abs(v) for k, v in la.items())
            if list(law(s)) == list(la) else None,
            "biased_law_terms": pv["selected"]["names"],
            "blast_radius_biased_party": pv["blast_radius"].get("party:BIASED-MD"),
            "unlearn_bitwise_exact": law(u) == la}
print(json.dumps(out["C"], indent=1))
json.dump(out, open(os.path.join(HERE, "results", "adopted_v2.json"), "w"), indent=1)
