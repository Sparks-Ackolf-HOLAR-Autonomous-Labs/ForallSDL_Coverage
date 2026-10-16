#!/usr/bin/env python3
"""Check that HG-SR v2 keeps the properties that make HG-SR different
(determinism, replication invariance, blast-radius localisation, exact
unlearning) on the md-alloy task. Writes results/demo_v2_props.json."""
import json
import os

import hgsr2
from improve_hgsr import alloy_task

HERE = os.path.dirname(os.path.abspath(__file__))
t = alloy_task()
rows, F = t["rows"], t["features"]
a = hgsr2.fit(rows, F)["selected"]
b = hgsr2.fit(rows, F)["selected"]
c0 = hgsr2.cell(rows[0])
rep = rows + [dict(r) for r in rows if hgsr2.cell(r) == c0] * 4
s = hgsr2.fit(rep, F)["selected"]
bias = [dict(r, kappa=r["kappa"] * 1.4, party="BIASED-MD", source="b" + r["source"])
        for r in rows if r["T"] > 1000]
pv = hgsr2.provenance(rows + bias, F)
u = hgsr2.fit([r for r in rows + bias if r["party"] != "BIASED-MD"], F)["selected"]
out = {
    "clean_law": dict(zip(a["names"], a["coefs"])),
    "deterministic": a["names"] == b["names"] and a["coefs"] == b["coefs"],
    "replicated_cell": c0,
    "replication_same_terms": s["names"] == a["names"],
    "replication_max_rel_coef_drift": max(abs(x - y) / abs(y)
                                          for x, y in zip(s["coefs"], a["coefs"])),
    "biased_law_terms": pv["selected"]["names"],
    "blast_radius_biased_party": pv["blast_radius"].get("party:BIASED-MD"),
    "unlearn_bitwise_exact": u["names"] == a["names"] and u["coefs"] == a["coefs"],
}
print(json.dumps(out, indent=1))
json.dump(out, open(os.path.join(HERE, "results", "demo_v2_props.json"), "w"), indent=1)
