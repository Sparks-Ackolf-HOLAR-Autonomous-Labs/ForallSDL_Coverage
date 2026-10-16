#!/usr/bin/env python3
"""Check that the full v2 configuration (shape filter, noise weights,
PySR proposals, beam search) keeps the four HG-SR properties on the
md-alloy task (2026-10-02). PySR is re-run from scratch for
every fit here (no proposal cache), so determinism and exact unlearning
are tested end to end, open search included. Writes
results/demo_full_props.json."""
import json
import os

import hgsr2
import propose
from full_v2 import MD_PHYS, MD_VARS
from improve_hgsr import alloy_task

HERE = os.path.dirname(os.path.abspath(__file__))
KW = dict(physics=MD_PHYS, noise=True, search="beam", beam=200)
DOM = MD_PHYS["box"] + MD_PHYS["inf"]


def fit(rows, feats):
    P = propose.proposals_for(rows, feats, MD_VARS, DOM)
    return hgsr2.fit(rows, feats, proposals=P, **KW), [p[0] for p in P]


def sel(m):
    s = m["selected"]
    return s["names"], s["coefs"]


t = alloy_task()
rows, F = t["rows"], t["features"]
a, pa = fit(rows, F)
b, pb = fit(rows, F)
c0 = hgsr2.cell(rows[0])
rep = rows + [dict(r) for r in rows if hgsr2.cell(r) == c0] * 4
s, ps = fit(rep, F)
bias = [dict(r, kappa=r["kappa"] * 1.4, party="BIASED-MD", source="b" + r["source"])
        for r in rows if r["T"] > 1000]
Pb = propose.proposals_for(rows + bias, F, MD_VARS, DOM)
pv = hgsr2.provenance(rows + bias, F, proposals=Pb, **KW)
u, pu = fit([r for r in rows + bias if r["party"] != "BIASED-MD"], F)
na, ca = sel(a)
out = {
    "config": "full v2: shape filter + noise weights + PySR proposals + beam 200",
    "clean_law": dict(zip(na, ca)),
    "proposals_clean": pa,
    "deterministic_including_pysr": sel(a) == sel(b) and pa == pb,
    "replicated_cell": c0,
    "replication_same_terms": sel(s)[0] == na,
    "replication_max_rel_coef_drift": max(abs(x - y) / abs(y) for x, y in zip(sel(s)[1], ca))
    if sel(s)[0] == na else None,
    "replication_proposals_changed": ps != pa,
    "biased_law_terms": pv["selected"]["names"],
    "blast_radius_biased_party": pv["blast_radius"].get("party:BIASED-MD"),
    "unlearn_bitwise_exact": sel(u) == sel(a) and pu == pa,
}
print(json.dumps(out, indent=1))
json.dump(out, open(os.path.join(HERE, "results", "demo_full_props.json"), "w"), indent=1)
