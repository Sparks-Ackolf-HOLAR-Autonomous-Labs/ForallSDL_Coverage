#!/usr/bin/env python3
"""Does the factor-of-2 uncertainty on small-SiC scattering change any decision?

Sweeps the particle contrast parameter lam2 over every value the model form
admits: from the density-only contrast (smallest bracket) through the
stiffness-only contrast (largest bracket) to lam2 -> infinity, where the
cross-section saturates at the geometric limit 2*pi*R^2. For each value it
recomputes the three decision-bearing outputs of the SiGeC campaign:
  (a) kappa_L reduction per 1 at% C (x=0.2, 0.3 um effective grains),
  (b) the Ge fraction that matches carbon-free Si0.8Ge0.2 (the Ge saving),
  (c) the screening pick (argmin kappa_L over x, C content, SiC radius).
Reads only committed evidence (almaBTE modes, converged MACE summary).
No new simulation; runtime is seconds."""
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault("SIGEC_M1_SUMMARY", os.path.join(HERE, "..", "m1c", "m1c_summary.json"))
sys.path.insert(0, os.path.dirname(HERE))
import analyze_sigec as S  # noqa: E402  (checks the design hash on import)

D_EFF = 0.3e-6
Y = 0.01
LAM = {"density": S.contrast(0.2, "density"), "stiffness": S.contrast(0.2, "stiffness")}
SWEEP = np.logspace(np.log10(LAM["density"]), 4, 41)
modes = {x: S.SModes(x) for x in S.XS}


def k(x, T, y, R, lam2):
    return modes[x].kappa2(T, D_EFF, 0.0, np_R=R, np_phi=S.phi(x, y), lam2=lam2) if y else modes[x].kappa2(T, D_EFF)


def ge_match(T, R, lam2):
    tgt = k(0.2, T, 0, None, 0)
    return S.xmin_for([k(x, T, Y, R, lam2) for x in S.XS], tgt)


def pick(T, lam2):
    best = None
    for x in S.XS:
        for y in S.YS:
            for R in S.RS:
                v = k(x, T, y, R, lam2)
                if best is None or v < best[0]:
                    best = (v, x, y, R * 1e9)
    return best


out = {"lam2_brackets_x0.2": LAM, "sweep": []}
for lam2 in SWEEP:
    row = {"lam2": lam2}
    for T in (300, 1100):
        k0 = k(0.2, T, 0, None, 0)
        for Rn in (2.5, 5.0):
            row[f"red_pct_R{Rn:g}_T{T}"] = 100 * (1 - k(0.2, T, Y, Rn * 1e-9, lam2) / k0)
            xm = ge_match(T, Rn * 1e-9, lam2)
            row[f"ge_saving_pct_R{Rn:g}_T{T}"] = None if xm is None else 100 * (1 - xm / 0.2)
    out["sweep"].append(row)
out["picks"] = {}
for name, lam2 in [("density", LAM["density"]), ("stiffness", LAM["stiffness"]), ("geometric_cap", 1e6)]:
    out["picks"][name] = {f"T{T}": dict(zip(("kappa", "x", "y_C", "R_nm"), pick(T, lam2))) for T in (300, 1100)}
json.dump(out, open(os.path.join(HERE, "decision_robustness.json"), "w"), indent=1, default=float)

f = lambda r, key: r[key]
lo, st, cap = out["sweep"][0], min(out["sweep"], key=lambda r: abs(r["lam2"] - LAM["stiffness"])), out["sweep"][-1]
for key in [kk for kk in lo if kk != "lam2"]:
    print(f"{key:28s} density {lo[key]!s:>8.6} | stiffness {st[key]!s:>8.6} | cap {cap[key]!s:>8.6}")
print(json.dumps(out["picks"], indent=1, default=float))
