#!/usr/bin/env python3
"""Figure data for P1 Fig. SiGeC at the measured radius-dependent SiC
contrast adopted at ledger entry 110 (r50 decision). Same lam2(R) rule as
analyze_r50.py (log-log interpolation of the measured curve up to 5 nm,
non-increasing extrapolation of the last segment beyond), same scattering
model and grain size as ../../radius_contrast/fit_radius_contrast.py.
Derived reporting only; no new decision. Output: measured_curves.json."""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(HERE)), "radius_contrast"))
import fit_radius_contrast as FR  # noqa: E402
S = FR.S
curve = json.load(open(os.path.join(HERE, "r50_results.json")))["measured_curve"]
lx = np.log([c[0] for c in curve]); ly = np.log([c[1] for c in curve])


def lamfun(Rn):  # same rule as analyze_r50.py
    if Rn <= 5.0:
        return float(np.exp(np.interp(np.log(Rn), lx, ly)))
    s = min(0.0, (ly[-1] - ly[-2]) / (lx[-1] - lx[-2]))
    return float(np.exp(ly[-1] + s * (np.log(Rn) - lx[-1])))


RADII = (2.5, 5.0, 10.0, 100.0)
out = {"source": "r50_results.json measured_curve (ledger 110)", "T": 300, "d_eff_m": FR.D_EFF,
       "lam2": {f"R{r:g}": lamfun(r) for r in RADII}, "panel_a_x0.2": {}, "kappa_ref_x0.2_noC": FR.K0[300], "panel_b_x_equal_kappa": {}}
for Rn in RADII:
    out["panel_a_x0.2"][f"R{Rn:g}"] = [[y, FR.kC(0.2, 300, y, Rn * 1e-9, lamfun(Rn))] for y in S.YS]
    out["panel_b_x_equal_kappa"][f"R{Rn:g}"] = [
        [y, S.xmin_for([FR.kC(x, 300, y, Rn * 1e-9, lamfun(Rn)) for x in S.XS], FR.K0[300])] for y in S.YS if y > 0]
json.dump(out, open(os.path.join(HERE, "measured_curves.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
