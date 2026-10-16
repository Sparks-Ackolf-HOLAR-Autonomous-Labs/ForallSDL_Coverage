#!/usr/bin/env python3
"""Two demos of what hypergraph-awareness buys that flat (X, y) engines
structurally cannot provide. Writes results/demo_bias.json.

Demo 1 (oversampling): replicate one regime's rows 5x, as happens when
one party samples its favorite corner densely. A flat pooled loss lets
the oversampled regime drag the law; the hyperedge-balanced loss keeps
coefficients close to the un-replicated fit.

Demo 2 (biased party + exact unlearning): add a biased evidence party
(+40% kappa, mid-composition regimes). HG-SR's leave-one-hyperedge-out
provenance names which terms of the law now depend on that party
(localization), and unlearning the party is a clean retrain that
restores the original model bit-for-bit (epsilon = delta = 0). A flat
engine has no party concept: it can neither localize the damage nor
certify its removal.
"""
import copy
import json
import os

import numpy as np

import hgsr
from datasets import TASKS

HERE = os.path.dirname(os.path.abspath(__file__))


def coefdrift(a, b):
    """Relative drift between two selected models' shared terms."""
    ca, cb = dict(zip(a["names"], a["coefs"])), dict(zip(b["names"], b["coefs"]))
    shared = set(ca) & set(cb)
    if not shared:
        return float("inf")
    return max(abs(cb[t] - ca[t]) / (abs(ca[t]) + 1e-12) for t in shared)


def main():
    task = TASKS["md"]()
    rows, feats = task["rows"], task["features"]
    out = {}

    # Demo 1: oversample regime x3T0 (x=0.5, 300 K) 5x
    dense = rows + [copy.deepcopy(r) for r in rows if r["regime"] == "x3T0"
                    for _ in range(4)]
    base_pool = hgsr.fit(rows, feats, pooled=True)["selected"]
    base_bal = hgsr.fit(rows, feats, pooled=False)["selected"]
    dens_pool = hgsr.fit(dense, feats, pooled=True)["selected"]
    dens_bal = hgsr.fit(dense, feats, pooled=False)["selected"]
    out["demo1_oversampling"] = {
        "replicated_regime": "x3T0 (x=0.5, 300 K) 5x",
        "pooled_drift": round(coefdrift(base_pool, dens_pool), 4),
        "balanced_drift": round(coefdrift(base_bal, dens_bal), 4),
        "pooled_terms_before": base_pool["names"],
        "pooled_terms_after": dens_pool["names"],
        "balanced_terms_before": base_bal["names"],
        "balanced_terms_after": dens_bal["names"],
    }

    # Demo 2: biased party, provenance localization, exact unlearning
    # the dangerous case: the biased party dominates a thinly covered
    # regime (high-T has 3 clean rows; the biased party adds 6 there)
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
    tainted = rows + biased

    clean = hgsr.fit(rows, feats)
    taint = hgsr.fit(tainted, feats)
    prov = hgsr.provenance(tainted, feats)
    dep_terms = prov["blast_radius"].get("party:BIASED-MD", [])
    unlearned = hgsr.fit([r for r in tainted if r["party"] != "BIASED-MD"], feats)
    exact = (unlearned["selected"]["names"] == clean["selected"]["names"] and
             np.allclose(unlearned["selected"]["coefs"], clean["selected"]["coefs"],
                         rtol=0, atol=0))
    out["demo2_biased_party"] = {
        "injected": f"{len(biased)} rows, +40% kappa, party BIASED-MD",
        "clean_model": dict(zip(clean["selected"]["names"], clean["selected"]["coefs"])),
        "tainted_model": dict(zip(taint["selected"]["names"], taint["selected"]["coefs"])),
        "coef_drift_from_bias": round(coefdrift(clean["selected"], taint["selected"]), 4),
        "provenance_blast_radius_of_biased_party": dep_terms,
        "term_labels_under_taint": prov["term_labels"],
        "unlearn_exact_bitwise": bool(exact),
    }

    with open(os.path.join(HERE, "results", "demo_bias.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
