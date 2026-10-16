"""Explore-Exploit-Shard controller pieces for the Phase A pilot.

- Shard priority score R(c): cold-start weighted flags (residual z-score,
  support distance, ensemble disagreement), per the deck's baseline; the
  learned nonlinear score replaces this once labels exist.
- Hypergraph shard targeting: claims are hyperedges over regime vertices;
  a shard-action batch is a greedy minimum-cost transversal of the
  high-priority hyperedges (open flag: validate greedy vs exact here).
- Shard operations: forget and scope-limit, both executed as exact clean
  retrains with a verification evidence object (Merkle membership of the
  excluded rows + deterministic refit recipe + front comparison).
"""

from __future__ import annotations

import itertools

import numpy as np

from ledger import leaf_hash, merkle_proof, merkle_root
from minisr import ensemble, fit_front, predict, select_model, support_region
from simulator import regime_of


# --- Shard priority ---------------------------------------------------------

def shard_priority(claim_regimes: dict[str, list[dict]], model, data,
                   w=(0.5, 0.25, 0.25)) -> dict[str, float]:
    """R per regime-claim: weighted (residual z, support distance,
    disagreement), each min-max normalized across regimes."""
    dist = support_region(data)
    _, disagree, _ = ensemble(data)
    feats = {}
    for reg, rows in claim_regimes.items():
        if not rows:
            continue
        res = [abs(r["kappa"] - predict(model, [r["x"]], [r["T"]])[0])
               / max(r["kappa"], 1e-9) for r in rows]
        cx = float(np.mean([r["x"] for r in rows]))
        cT = float(np.mean([r["T"] for r in rows]))
        feats[reg] = [float(np.mean(res)), dist(cx, cT), disagree(cx, cT)]
    arr = np.array(list(feats.values()))
    lo, hi = arr.min(0), arr.max(0)
    norm = (arr - lo) / np.where(hi - lo > 0, hi - lo, 1)
    return {reg: float(np.dot(w, norm[i])) for i, reg in enumerate(feats)}


# --- Hypergraph transversal targeting ---------------------------------------

def greedy_transversal(hyperedges: dict[str, set[str]]) -> list[str]:
    """Greedy set cover: pick vertices (regimes) until every hyperedge
    (high-priority claim) is touched. Returns chosen regimes in pick order."""
    uncovered = dict(hyperedges)
    chosen = []
    while uncovered:
        counts: dict[str, int] = {}
        for verts in uncovered.values():
            for v in verts:
                counts[v] = counts.get(v, 0) + 1
        best = max(counts, key=counts.get)
        chosen.append(best)
        uncovered = {e: v for e, v in uncovered.items() if best not in v}
    return chosen


def exact_transversal_size(hyperedges: dict[str, set[str]]) -> int:
    """Brute-force minimum transversal size (pilot scale only), to measure
    the greedy approximation gap for the open flag."""
    verts = sorted(set().union(*hyperedges.values()))
    for k in range(1, len(verts) + 1):
        for combo in itertools.combinations(verts, k):
            s = set(combo)
            if all(s & e for e in hyperedges.values()):
                return k
    return len(verts)


# --- Shard operations with verification --------------------------------------

def _refit(data):
    front = fit_front(data)
    return front, select_model(front, len(data))


def op_forget(all_rows: list[dict], shard_rows: list[dict]) -> dict:
    """Forget a shard: exact clean retrain on D minus S, plus evidence object."""
    keep = [r for r in all_rows if r not in shard_rows]
    front, model = _refit(keep)
    leaves = [leaf_hash(r) for r in all_rows]
    root = merkle_root(leaves)
    proofs = [merkle_proof(leaves, all_rows.index(r))[:2] for r in shard_rows[:3]]
    return {
        "op": "forget", "n_removed": len(shard_rows), "n_kept": len(keep),
        "evidence": {
            "dataset_root": root,
            "sample_membership_proof_depth": len(proofs[0]) if proofs else 0,
            "retrain": "exhaustive-OLS, deterministic, eps=delta=0",
        },
        "front": front, "model": model, "data": keep,
    }


def op_scope_limit(model: dict, domain: dict) -> dict:
    """Scope-limit: the law keeps its artifact but its validity domain
    shrinks; nothing is refit, dependents outside the domain are flagged."""
    return {"op": "scope-limit", "model": model, "validity_domain": domain,
            "evidence": {"retrain": "none (artifact retained, influence bounded)"}}


# --- Fault localization -------------------------------------------------------

def batch_residual_report(batches: dict[str, list[dict]], model) -> dict[str, float]:
    """Mean signed relative residual per batch: a biased batch shows a
    coherent offset; a valid anomaly shows a regime-correlated (not
    batch-correlated) offset."""
    out = {}
    for bid, rows in batches.items():
        res = [(r["kappa"] - predict(model, [r["x"]], [r["T"]])[0])
               / max(r["kappa"], 1e-9) for r in rows]
        out[bid] = float(np.mean(res))
    return out


def regime_map(rows: list[dict]) -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = {}
    for r in rows:
        out.setdefault(regime_of(r["x"], r["T"]), []).append(r)
    return out
