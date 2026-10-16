"""Run-cost models for the P3 SiGeC factor space (x_Ge, y_C, carbon form, d_eff, dopant, T).

LAB_SHARED (primary): a run is one thermal-diffusivity measurement of a sintered
pellet at one temperature. A pellet recipe (x_Ge, y_C, carbon form, d_eff,
dopant) is paid once per batch, however many temperatures it is measured at.
Hours of instrument and operator time, assumed (not measured in this project):
  milling to reach the target grain size  d_eff 0.1/0.3/0.5/0.9 um -> 20/10/5/2 h
  hot pressing (sinter, cool)             3 h per pellet
  carbon addition, only when y_C > 0      SiC 2.5 or 5 nm: 4 h (nanoparticle dispersion)
                                          SiC 25 or 100 nm: 1 h (powder blending)
                                          substitutional: 8 h (rapid-solidification route)
  dopant master alloy (B or P)            1 h
  germanium powder (material cost)        4 x x_Ge h-equivalent
  measurement at T (laser flash, ramp)    1 + (T - 300 K)/400 K h  (1 h at 300 K, 3 h at 1100 K)
LAB_RUN: as LAB_SHARED, but every run pays its own pellet (no reuse across temperatures).
SIM: evidence cost when each run is scored by the model and checked atomistically
where the project's own validity limits require it; measured per-point costs:
  almaBTE host phonons                    75 min per distinct x_Ge (party4_bte/RESULTS.md)
  T >= 1100 K (BTE not accepted there)    + 69 min: one usable Tersoff rNEMD estimate,
                                            3 box lengths x 23 min (chpc/RESULTS.md,
                                            p3_prework/covering scenario B)
  SiC radius <= 5 nm with y_C > 0         + 246 min: MACE-MP-0 rNEMD check point
                                            (party4_mace/RESULTS.md, 4.1 h per point)
  floor                                   1 min per run (job overhead)
"""
import numpy as np
import common as C

MILL = {0.1: 20.0, 0.3: 10.0, 0.5: 5.0, 0.9: 2.0}
CARB = {"substitutional_eps0": 8.0, "SiC_R2.5nm": 4.0, "SiC_R5nm": 4.0, "SiC_R25nm": 1.0, "SiC_R100nm": 1.0}


def _vals(i):
    return [C.LEV[j][C.FULL[i][j]] for j in range(C.NF)]


PELLET = np.empty(C.N)
MEAS = np.empty(C.N)
SIMRUN = np.empty(C.N)
RECIPE = np.empty(C.N, dtype=np.int64)
XLEV = C.FULLA[:, 0].copy()
_rec = {}
for i in range(C.N):
    x, yc, form, d, dop, T = _vals(i)
    PELLET[i] = MILL[d] + 3.0 + (CARB[form] if yc > 0 else 0.0) + (0.0 if dop == "none" else 1.0) + 4.0 * x
    MEAS[i] = 1.0 + (T - 300) / 400.0
    small_sic = yc > 0 and form in ("SiC_R2.5nm", "SiC_R5nm")
    SIMRUN[i] = 1.0 + (69.0 if T >= 1100 else 0.0) + (246.0 if small_sic else 0.0)
    RECIPE[i] = _rec.setdefault(tuple(C.FULL[i][:5]), len(_rec))
NREC = len(_rec)


class LabShared:
    name, unit = "lab_shared_pellet", "h"

    def init(self, rows):
        s = np.zeros(NREC, bool)
        for r in rows:
            s[RECIPE[C.R.IDX[tuple(r)]]] = True
        return s

    def marginal(self, s):
        return MEAS + np.where(s[RECIPE], 0.0, PELLET)

    def add(self, s, i):
        s[RECIPE[i]] = True

    def marg1(self, s, i):
        return MEAS[i] + (0.0 if s[RECIPE[i]] else PELLET[i])

    def total(self, rows):
        ii = [C.R.IDX[tuple(r)] for r in rows]
        recs = {RECIPE[i]: PELLET[i] for i in ii}
        return float(MEAS[ii].sum() + sum(recs.values()))


class LabRun:
    name, unit = "lab_pellet_per_run", "h"
    COST = PELLET + MEAS

    def init(self, rows):
        return None

    def marginal(self, s):
        return self.COST

    def add(self, s, i):
        pass

    def marg1(self, s, i):
        return self.COST[i]

    def total(self, rows):
        return float(sum(self.COST[C.R.IDX[tuple(r)]] for r in rows))


class Sim:
    name, unit = "simulation_evidence", "min"

    def init(self, rows):
        s = np.zeros(C.NL[0], bool)
        for r in rows:
            s[r[0]] = True
        return s

    def marginal(self, s):
        return SIMRUN + np.where(s[XLEV], 0.0, 75.0)

    def add(self, s, i):
        s[XLEV[i]] = True

    def marg1(self, s, i):
        return SIMRUN[i] + (0.0 if s[XLEV[i]] else 75.0)

    def total(self, rows):
        ii = [C.R.IDX[tuple(r)] for r in rows]
        return float(SIMRUN[ii].sum() + 75.0 * len({C.FULL[i][0] for i in ii}))


MODELS = [LabShared(), LabRun(), Sim()]


# ---------------------------------------------------------------- pellet-aware greedy
# Rows grouped by pellet recipe: ROWS_RT[r, t] is the full-factorial row of recipe r at level t of T.
ROWS_RT = np.empty((NREC, C.NL[5]), dtype=np.int64)
for i in range(C.N):
    ROWS_RT[RECIPE[i], C.FULL[i][5]] = i
_pairs = list(__import__("itertools").combinations(range(C.NF), 2))
NONT_COLS = [k for k, (a, b) in enumerate(_pairs) if 5 not in (a, b)]
T_COLS = [k for k, (a, b) in enumerate(_pairs) if 5 in (a, b)]


def greedy_pellet(mult, seed, base=None):
    """Cost-weighted greedy whose unit is one pellet measured at a set of temperatures.
    Each step takes the (recipe, temperature set) with the largest
    (deficient pair-cell samples added) / (pellet cost if new + measurement costs)."""
    rng = np.random.default_rng(seed)
    rows = [tuple(r) for r in (base or [])]
    cnt = C.cell_counts(rows, C.KI2, C.NK2)
    made = np.zeros(NREC, bool)
    for r in rows:
        made[RECIPE[C.R.IDX[r]]] = True
    taken = np.zeros(C.N, bool)
    for r in rows:
        taken[C.R.IDX[r]] = True
    KN = C.KI2[ROWS_RT[:, 0]][:, NONT_COLS]          # (NREC, 10), same for every T
    KT = C.KI2[ROWS_RT][:, :, T_COLS]                 # (NREC, 5, 5)
    meas_rt = MEAS[ROWS_RT]                           # (NREC, 5)
    pel = PELLET[ROWS_RT[:, 0]]
    while True:
        deficit = np.clip(mult - cnt, 0, None)
        if not deficit.any():
            return rows
        dn = deficit[KN]                               # (NREC, 10)
        gT = (deficit[KT] > 0).sum(2).astype(float)    # (NREC, 5)
        gT[taken[ROWS_RT]] = -1e9                      # no duplicate rows inside a step
        order = np.argsort(-(gT / meas_rt), axis=1, kind="stable")
        gs = np.take_along_axis(gT, order, 1)
        ms = np.take_along_axis(meas_rt, order, 1)
        best_ratio = np.full(NREC, -np.inf)
        best_s = np.zeros(NREC, int)
        csum_g = np.cumsum(gs, 1)
        csum_m = np.cumsum(ms, 1)
        pc = np.where(made, 0.0, pel)
        for s in range(1, C.NL[5] + 1):
            gain = np.minimum(dn, s).sum(1) + csum_g[:, s - 1]
            ratio = gain / (pc + csum_m[:, s - 1])
            ratio[gs[:, s - 1] < 0] = -np.inf
            better = ratio > best_ratio + 1e-12
            best_ratio[better] = ratio[better]
            best_s[better] = s
        perm = rng.permutation(NREC)
        r = perm[int(np.argmax(best_ratio[perm]))]
        for t in order[r, :best_s[r]]:
            i = ROWS_RT[r, t]
            rows.append(C.FULL[i])
            taken[i] = True
            np.add.at(cnt, C.KI2[i], 1)
        made[r] = True


# ---------------------------------------------------------------- LP / ILP bounds
def lp_bound(model, mult, integer=False, time_limit=120):
    """Lowest cost of any batch sampling every pair cell >= mult times (LP relaxation, or ILP
    with a time limit; for the ILP the best bound and the incumbent are both returned)."""
    from scipy.optimize import milp, LinearConstraint, Bounds
    from scipy.sparse import csr_matrix, hstack, vstack, identity
    rows_i = np.repeat(np.arange(C.N), C.KI2.shape[1])
    A = csr_matrix((np.ones(rows_i.size), (C.KI2.ravel(), rows_i)), shape=(C.NK2, C.N))
    if model.name == "lab_shared_pellet":
        P = csr_matrix((np.ones(C.N), (np.arange(C.N), RECIPE)), shape=(C.N, NREC))
        c = np.concatenate([MEAS, PELLET[ROWS_RT[:, 0]]])
        Acov = hstack([A, csr_matrix((C.NK2, NREC))])
        Alink = hstack([identity(C.N), -P])            # z_i <= p_recipe
        cons = [LinearConstraint(Acov, mult, np.inf), LinearConstraint(Alink, -np.inf, 0)]
        n = C.N + NREC
        ub = np.concatenate([np.full(C.N, 1.0), np.ones(NREC)])  # distinct rows (the greedy may repeat none here)
        const = 0.0
    else:
        c = (PELLET + MEAS) if model.name == "lab_pellet_per_run" else SIMRUN.copy()
        cons = [LinearConstraint(A, mult, np.inf)]
        n = C.N
        ub = np.full(C.N, float(mult))
        const = 75.0 * C.NL[0] if model.name == "simulation_evidence" else 0.0  # every x level must appear
    res = milp(c, constraints=cons, bounds=Bounds(0, ub), integrality=np.ones(n) if integer else None,
               options={"time_limit": time_limit, "disp": False})
    out = {"status": int(res.status), "message": str(res.message)}
    if res.x is not None:
        out["objective"] = float(res.fun) + const
    if integer and getattr(res, "mip_dual_bound", None) is not None:
        out["dual_bound"] = float(res.mip_dual_bound) + const
    if integer and res.x is not None:
        z = res.x[:C.N]
        rows = []
        for i in np.nonzero(z > 0.5)[0]:
            rows += [C.FULL[i]] * int(round(z[i]))
        out["rows"] = rows
    return out
