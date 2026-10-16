"""Shared pieces for the P3 limitation studies (p3_limits/).

Reuses the registered P3 SiGeC design study (../p3_sigec/run_p3_sigec.py):
same factors, full factorial, response model, surrogate and scores. The only
change is the thread count of the tree ensemble (4 instead of 8, a shared
machine); the ensemble's random_state is unchanged, so scores are identical.
"""
import itertools, json, os, sys
import numpy as np
from sklearn.ensemble import ExtraTreesRegressor

HERE = os.path.dirname(os.path.abspath(__file__))
MT = os.path.dirname(HERE)
OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, os.path.join(MT, "p3_sigec"))
import run_p3_sigec as R  # noqa: E402  (checks both registered design hashes on import)

NJOBS = 4


def _fit(rows, y):
    ii = [R.IDX[r] for r in rows]
    m = ExtraTreesRegressor(n_estimators=300, min_samples_leaf=1, random_state=0, n_jobs=NJOBS)
    m.fit(R.XF[ii], np.log(y[ii]))
    return m


R.fit = _fit  # evaluate() looks fit up as a module global

FULL = R.FULL
LEV, NL, NAMES = R.LEV, R.NL, R.NAMES
NF = len(NL)
N = len(FULL)
FULLA = np.array(FULL)


def load_response():
    """Response of the registered model (stiffness contrast), checked against the saved copy."""
    y = np.load(os.path.join(MT, "p3_sigec", "results", "response_stiffness.npy"))
    return y


def key_index(t):
    """Index array (N, C(NF,t)) of the t-factor level cells each full-factorial row touches."""
    combos = list(itertools.combinations(range(NF), t))
    cols, off = [], 0
    for c in combos:
        sizes = [NL[j] for j in c]
        idx = np.zeros(N, dtype=np.int64)
        mul = 1
        for j, s in zip(reversed(c), reversed(sizes)):
            idx += FULLA[:, j] * mul
            mul *= s
        cols.append(idx + off)
        off += int(np.prod(sizes))
    return np.stack(cols, 1), off


KI2, NK2 = key_index(2)
KI3, NK3 = key_index(3)


def cell_counts(rows, ki, nk):
    c = np.zeros(nk, dtype=np.int64)
    for r in rows:
        np.add.at(c, ki[R.IDX[tuple(r)]], 1)
    return c


def coverage_t(rows, t):
    ki, nk = (KI2, NK2) if t == 2 else (KI3, NK3)
    return float((cell_counts(rows, ki, nk) > 0).mean())


def greedy_cover(t, mult, seed, cost=None, base=None):
    """Greedy covering: add the row with the largest (deficient cells touched) / (marginal cost)
    until every t-factor cell is sampled `mult` times. cost=None gives the uniform-cost rule.
    cost is an object with .marginal(state) -> array(N) and .add(state, i)."""
    ki, nk = (KI2, NK2) if t == 2 else (KI3, NK3)
    rng = np.random.default_rng(seed)
    rows = [tuple(r) for r in (base or [])]
    cnt = cell_counts(rows, ki, nk)
    state = cost.init(rows) if cost else None
    while True:
        need = cnt < mult
        if not need.any():
            return rows
        gain = need[ki].sum(1).astype(float)
        score = gain / cost.marginal(state) if cost else gain
        perm = rng.permutation(N)
        i = perm[int(np.argmax(score[perm]))]
        rows.append(FULL[i])
        np.add.at(cnt, ki[i], 1)
        if cost:
            cost.add(state, i)


def jdump(obj, name):
    json.dump(obj, open(os.path.join(OUT, name), "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
