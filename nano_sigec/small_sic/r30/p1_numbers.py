#!/usr/bin/env python3
"""Numbers quoted in P1 after the r30 decision (ledger 105), from the
registered measured curve in r30_results.json: recipe kappa_L at the two
radii and both compositions, and the carbon-grain-boundary comparison at
20 nm and 0.5 um grains. Derived reporting only; no new decision.
Output: p1_numbers.json."""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(HERE)), "radius_contrast"))
import fit_radius_contrast as FR  # noqa: E402
r = json.load(open(os.path.join(HERE, "r30_results.json")))
xs = np.log([c[0] for c in r["measured_curve"]]); ys = np.log([c[1] for c in r["measured_curve"]])


def lamfun(Rn):  # same rule as analyze_r30.py
    if Rn <= 3.0:
        return float(np.exp(np.interp(np.log(Rn), xs, ys)))
    s = min(0.0, (ys[-1] - ys[-2]) / (xs[-1] - xs[-2]))
    return float(np.exp(ys[-1] + s * (np.log(Rn) - xs[-1])))


m, S = FR.MODES[0.2], FR.S
gb = {}
for d in (20e-9, 0.5e-6):
    gb[f"d{d:g}"] = {"no_C": m.kappa2(300, d),
                     **{f"C1_R{Rn:g}": m.kappa2(300, d, 0.0, np_R=Rn * 1e-9, np_phi=S.phi(0.2, 0.01), lam2=lamfun(Rn))
                        for Rn in (2.5, 5.0)}}
out = {"slope_2_to_3nm": float((ys[-1] - ys[-2]) / (xs[-1] - xs[-2])),
       "lam2": {"R2.5": lamfun(2.5), "R5": lamfun(5.0)},
       "recipe_kappa": FR.recipe_kappa(lamfun), "carbon_vs_grain_x0.2_300K": gb}
json.dump(out, open(os.path.join(HERE, "p1_numbers.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
