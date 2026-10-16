"""Mock evidence party for the Phase A pilot: kappa_L(x, T) of Si1-xGex.

Ground truth follows the standard alloy picture (series resistance of the
pure-crystal 1/T law plus an x(1-x) alloy-scattering term), with two
deliberate wrinkles that the physics of the real testbed provides for free:

- FAULTY BATCH: one batch is generated with a calibration bias
  (kappa scaled by 0.72), standing in for a biased MD batch. The right
  shard operation is forget (or correct if the bias were identified).
- VALID ANOMALY: above T_ANOM the world is richer than the model
  (an extra defect/near-melting term suppresses kappa). Data are CORRECT;
  the right shard operation is scope-limit, not forget.

Units: kappa in W/(m K), T in K, x = Ge fraction in [0, 1].
Coefficients are order-of-magnitude physical, tuned for pilot legibility,
not literature-fitted. The real party (LAMMPS/Tersoff) replaces this file
in Phase B; the interface (run_batch) is what carries over.
"""

from __future__ import annotations

import numpy as np

# 1/kappa = (A + B/T contributions) + alloy scattering
A_PURE = 1.0e-3      # residual resistivity-like term, (m K)/W
B_PURE = 2.0        # Umklapp: contributes B/T... modeled as C*T form below
C_UMKLAPP = 2.2e-6   # * T   -> pure Si ~ 140 W/mK at 300 K
D_ALLOY = 6.0e-2     # * x(1-x): alloy scattering dominates mid-composition
E_GE = 1.5e-3        # * x: Ge is softer, higher baseline resistivity

T_ANOM = 1000.0      # K; above this the richer-world term kicks in
F_ANOM = 2.0e-5      # * (T - T_ANOM): extra suppression (defects/anharmonicity)

NOISE_REL = 0.04     # 4 percent relative noise, Green-Kubo-like scatter


def kappa_true(x: float, T: float) -> float:
    r = A_PURE + C_UMKLAPP * T + D_ALLOY * x * (1 - x) + E_GE * x
    if T > T_ANOM:
        r += F_ANOM * (T - T_ANOM)          # world is richer up here
    return 1.0 / r


def kappa_model_ideal(x: float, T: float) -> float:
    """What the in-domain law (no anomaly term) predicts. Reference only."""
    return 1.0 / (A_PURE + C_UMKLAPP * T + D_ALLOY * x * (1 - x) + E_GE * x)


def run_batch(rows: list[dict], seed: int, fault_scale: float = 1.0) -> list[dict]:
    """Execute a committed batch design. Each row: {x, T}. Returns results.

    fault_scale != 1.0 simulates a biased evidence party for that batch.
    """
    rng = np.random.default_rng(seed)
    out = []
    for r in rows:
        k = kappa_true(r["x"], r["T"])
        k_obs = k * fault_scale * (1.0 + NOISE_REL * rng.standard_normal())
        out.append({**r, "kappa": float(k_obs), "seed": seed})
    return out


# Coverage-designed seeding (design element 2 of 04_experimental_design_v2) --

X_BANDS = [(0.0, 0.2), (0.2, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 1.0)]
T_BANDS = [(300, 500), (500, 700), (700, 900), (900, 1100), (1100, 1300)]


def regime_of(x: float, T: float) -> str:
    xi = min(int(x * 5), 4)
    ti = min(int((T - 300) / 200), 4)
    return f"x{xi}T{ti}"


def seed_design(reps: int = 2, rng_seed: int = 7) -> list[dict]:
    """Full-coverage design: every (x-band, T-band) regime sampled `reps`
    times (strength-2 coverage of the two factors, trivially, since we take
    the full grid of bands). reps >= 2 gives k-loss tolerance k = 1 per
    regime: any single failed run leaves every regime still covered."""
    rng = np.random.default_rng(rng_seed)
    rows = []
    for xlo, xhi in X_BANDS:
        for tlo, thi in T_BANDS:
            for _ in range(reps):
                rows.append({
                    "x": float(round(rng.uniform(xlo, xhi), 4)),
                    "T": float(round(rng.uniform(tlo, thi), 1)),
                })
    return rows
