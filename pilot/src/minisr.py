"""Deterministic mini symbolic regression for the Phase A pilot.

We regress thermal RESISTIVITY r = 1/kappa on a fixed physical basis
library and enumerate all term subsets up to MAX_TERMS, exactly. The
Pareto front over (complexity, validation MSE) plays the role of the PySR
front; because the search is exhaustive and deterministic, the clean
retrain is bit-for-bit reproducible, which makes shard-operation
verification exact (epsilon = delta = 0 by construction).

PySR replaces this on CHPC in Phase B; passports (support region,
ensemble disagreement) and the claim format are the stable interface.
"""

from __future__ import annotations

import itertools

import numpy as np

# (name, complexity, callable) — complexity ~ node count of the term
BASIS = [
    ("1",       1, lambda x, T: np.ones_like(x)),
    ("x",       1, lambda x, T: x),
    ("T",       1, lambda x, T: T),
    ("x(1-x)",  4, lambda x, T: x * (1 - x)),
    ("x*T",     3, lambda x, T: x * T),
    ("1/T",     2, lambda x, T: 1.0 / T),
    ("x^2",     2, lambda x, T: x**2),
    ("x(1-x)*T", 6, lambda x, T: x * (1 - x) * T),
]
MAX_TERMS = 4


def _design_matrix(x, T, terms):
    return np.column_stack([BASIS[i][2](x, T) for i in terms])


def fit_front(data: list[dict]) -> list[dict]:
    """Exhaustive subset OLS on r = 1/kappa; returns the Pareto front,
    each entry: {terms, names, coefs, complexity, mse}."""
    x = np.array([d["x"] for d in data])
    T = np.array([d["T"] for d in data])
    r = 1.0 / np.array([d["kappa"] for d in data])
    models = []
    for k in range(1, MAX_TERMS + 1):
        for terms in itertools.combinations(range(len(BASIS)), k):
            M = _design_matrix(x, T, terms)
            coefs, *_ = np.linalg.lstsq(M, r, rcond=None)
            mse = float(np.mean((M @ coefs - r) ** 2))
            models.append({
                "terms": list(terms),
                "names": [BASIS[i][0] for i in terms],
                "coefs": [float(c) for c in coefs],
                "complexity": int(sum(BASIS[i][1] for i in terms)),
                "mse": mse,
            })
    models.sort(key=lambda m: (m["complexity"], m["mse"]))
    front, best = [], np.inf
    for m in models:
        if m["mse"] < best * (1 - 1e-9):
            front.append(m)
            best = m["mse"]
    return front


def predict(model: dict, x, T):
    M = _design_matrix(np.asarray(x, float), np.asarray(T, float), model["terms"])
    r = M @ np.array(model["coefs"])
    return 1.0 / r


def select_model(front: list[dict], n: int) -> dict:
    """BIC-style selection over the front: n*ln(mse) + complexity*ln(n)."""
    return min(front, key=lambda m: n * np.log(m["mse"] + 1e-30)
               + m["complexity"] * np.log(n))


# Passports -----------------------------------------------------------------

def support_region(data: list[dict]):
    """Returns a callable dist(x, T): normalized distance to nearest
    training point. Passport 1: using the law far from support is flagged."""
    pts = np.array([[d["x"], d["T"] / 1000.0] for d in data])

    def dist(x: float, T: float) -> float:
        q = np.array([x, T / 1000.0])
        return float(np.min(np.linalg.norm(pts - q, axis=1)))
    return dist


def ensemble(data: list[dict], m: int = 12, seed: int = 0):
    """Passport 2: m bootstrap refits. Returns per-term front-membership
    frequency and a prediction-disagreement callable."""
    rng = np.random.default_rng(seed)
    n = len(data)
    fronts, picks = [], []
    for _ in range(m):
        idx = rng.integers(0, n, n)
        f = fit_front([data[i] for i in idx])
        fronts.append(f)
        picks.append(select_model(f, n))
    term_freq = {}
    for p in picks:
        for name in p["names"]:
            term_freq[name] = term_freq.get(name, 0) + 1
    term_freq = {k: v / m for k, v in term_freq.items()}

    def disagreement(x: float, T: float) -> float:
        preds = [predict(p, [x], [T])[0] for p in picks]
        return float(np.std(preds) / (np.mean(preds) + 1e-12))
    return term_freq, disagreement, picks
