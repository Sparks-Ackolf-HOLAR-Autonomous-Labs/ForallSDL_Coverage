#!/usr/bin/env python3
"""Shard operations on HG-SR: the SISA-like mechanics, beyond forget.

SISA (Bourtoule et al., IEEE S&P 2021) gets cheap unlearning from
structure: partition the training set into shards, train an isolated
sub-model per shard, aggregate; deleting a point retrains only its
shard. HG-SR keeps the two useful SISA properties but relocates them:

  - Isolation: SISA isolates *compute* (one sub-model per shard).
    HG-SR isolates *influence*: the hyperedge-balanced loss caps every
    evidence hyperedge's vote at 1/|edge| per row (a shard cannot
    outvote the rest by being big), and leave-one-hyperedge-out
    provenance records exactly which terms each hyperedge supports.
  - Bounded recompute: SISA bounds it at one sub-model. HG-SR's full
    clean retrain is an exhaustive deterministic OLS costing seconds,
    so recompute needs no bounding; what stays bounded is the *belief
    change*, which the blast radius B(e) states in advance.

One more structural difference: SISA shards are a disjoint random
partition of rows. HG-SR shards are semantic and overlapping (every row
sits in a regime edge AND a party edge AND a source edge), so an
operation can target "the literature party" or "the x2T0 regime"
without renumbering anything, and no aggregate ensemble is needed: the
law stays one global fit, and every operation is exact
(epsilon = delta = 0) because the refit is the mechanism.

The operator family U (00_terminology_shard_operations.md) maps onto
four row/weight primitives, all going through the same deterministic
refit; `apply()` below implements them:

  forget      remove the edge's rows, refit (bitwise-exact; in hgsr.py)
  quarantine  weight 0 on the edge's rows, rows kept: reversible forget
  downweight  weight gamma in (0,1) on the edge's rows (continuous
              trust; quarantine is gamma=0, retain is gamma=1)
  restore     lift a quarantine/downweight (gamma back to 1)
  correct     edit fields of the edge's rows (e.g. recalibrated kappa),
              refit; ledger keeps old rows
  split       re-tag the edge's rows into finer edges (e.g. literature
              party split by microstructure: bulk/film/nano) so the
              balanced loss and provenance see them separately; the
              fitted-law change quantifies what pooling had hidden
  replace     forget + add replacement rows in one transaction
  scope-limit no refit: restrict a claim's validity to covered regimes
              or attach a party scope label (hgsr.provenance already
              emits the labels; the ledger records the restriction)
  retain /    no refit: status bookkeeping on the ledger, backed by the
  strengthen  provenance labels (a term going multi-supported after new
              evidence lands IS the strengthen witness)

Verification for every operation is the clean-retrain comparison the
project already uses: rerun `hgsr.fit` on the operated evidence and
compare selected terms/coefficients; determinism makes it exact.
"""
from __future__ import annotations

import copy

import numpy as np

import hgsr


def _edge_rows(rows, edge):
    kind, _, val = edge.partition(":")
    return [i for i, r in enumerate(rows) if str(r.get(kind)) == val]


def apply(rows, features, op, edge=None, gamma=None, patch=None,
          retag=None, replacement=None, max_terms=hgsr.MAX_TERMS):
    """Apply one shard operation and return {rows, weights_scale, model,
    delta} where delta reports the belief change vs the pre-op fit.

    op: forget | quarantine | downweight | restore | correct | split |
        replace. (retain/strengthen/scope-limit are ledger operations:
        they change no rows and need no refit, so they are not here.)
    edge: hyperedge id, e.g. "party:LIT" or "regime:x2T0".
    gamma: downweight factor (downweight only).
    patch: {field: new_value} applied to the edge's rows (correct only).
    retag: callable row -> new edge value for the edge's kind (split).
    replacement: list of new rows (replace only).
    """
    before = hgsr.fit(rows, features, max_terms=max_terms)["selected"]
    new_rows = copy.deepcopy(rows)
    scale = {}                       # edge id -> extra weight multiplier

    if op == "forget":
        new_rows = [r for i, r in enumerate(new_rows)
                    if i not in set(_edge_rows(rows, edge))]
    elif op in ("quarantine", "downweight", "restore"):
        g = {"quarantine": 0.0, "restore": 1.0}.get(op, gamma)
        if g is None:
            raise ValueError("downweight needs gamma")
        scale[edge] = float(g)
    elif op == "correct":
        for i in _edge_rows(rows, edge):
            new_rows[i].update(patch)
    elif op == "split":
        kind = edge.split(":")[0]
        for i in _edge_rows(rows, edge):
            new_rows[i][kind] = retag(new_rows[i])
    elif op == "replace":
        new_rows = [r for i, r in enumerate(new_rows)
                    if i not in set(_edge_rows(rows, edge))] + list(replacement)
    else:
        raise ValueError(f"unknown shard operation: {op}")

    model = fit_weighted(new_rows, features, scale, max_terms=max_terms)
    after = model["selected"]
    delta = {
        "terms_dropped": [t for t in before["names"] if t not in after["names"]],
        "terms_added": [t for t in after["names"] if t not in before["names"]],
        "coef_moves": {
            t: {"before": b, "after": a}
            for t, b, a in (
                (t, before["coefs"][before["names"].index(t)],
                 after["coefs"][after["names"].index(t)])
                for t in before["names"] if t in after["names"])
            if abs(a - b) > 1e-12},
    }
    return {"rows": new_rows, "weights_scale": scale,
            "model": model, "delta": delta}


def fit_weighted(rows, features, scale=None, max_terms=hgsr.MAX_TERMS):
    """hgsr.fit with per-hyperedge trust multipliers on top of the
    balanced weights (quarantine/downweight/restore). scale={} or None
    reproduces hgsr.fit exactly."""
    if not scale:
        return hgsr.fit(rows, features, max_terms=max_terms)
    mult = np.ones(len(rows))
    for edge, g in scale.items():
        ix = _edge_rows(rows, edge)
        mult[ix] = mult[ix] * g
    kept = [i for i in range(len(rows)) if mult[i] > 0]
    if len(kept) < len(rows):
        # gamma=0 rows carry zero weight; drop them from the refit but
        # NOT from the ledger (that asymmetry is what makes quarantine
        # reversible and forget not).
        sub = [rows[i] for i in kept]
        m = np.array([mult[i] for i in kept])
    else:
        sub, m = rows, mult
    if np.allclose(m, m[0] if len(m) else 1.0):
        # uniform multiplier: identical to the plain balanced fit
        return hgsr.fit(sub, features, max_terms=max_terms)
    return _fit_scaled(sub, features, m, max_terms)


def _fit_scaled(rows, features, mult, max_terms):
    """Balanced fit with non-uniform trust multipliers (downweight)."""
    import itertools
    basis = hgsr.make_basis(features)
    d = hgsr._cols(rows)
    r = 1.0 / d["kappa"]
    n = len(rows)
    groups = [np.array(ix) for _, ix in
              sorted(hgsr._edges(rows, ("regime",)).items())]
    w = np.ones(n)
    for g in groups:
        w[g] = 1.0 / len(g)
    w = w * mult
    sw = np.sqrt(w)
    models = []
    for k in range(1, max_terms + 1):
        for terms in itertools.combinations(range(len(basis)), k):
            M = np.column_stack([basis[i][2](d) for i in terms])
            coefs, *_ = np.linalg.lstsq(M * sw[:, None], r * sw, rcond=None)
            res2 = (M @ coefs - r) ** 2
            # weighted per-hyperedge MSE keeps downweighted rows from
            # steering the selection, mirroring hgsr.fit's loss
            loss = float(np.mean(
                [np.average(res2[g], weights=mult[g]) if mult[g].sum() > 0
                 else 0.0 for g in groups]))
            models.append({"terms": list(terms),
                           "names": [basis[i][0] for i in terms],
                           "coefs": [float(c) for c in coefs],
                           "complexity": int(sum(basis[i][1] for i in terms)),
                           "mse": loss})
    models.sort(key=lambda m: (m["complexity"], m["mse"]))
    front, best = [], np.inf
    for m in models:
        if m["mse"] < best * (1 - 1e-9):
            front.append(m)
            best = m["mse"]
    n_eff = float(np.sum(w))
    sel = min(front, key=lambda m: n_eff * np.log(m["mse"] + 1e-30)
              + m["complexity"] * np.log(n_eff))
    return {"front": front, "selected": sel, "features": features,
            "pooled": False, "n": n}
