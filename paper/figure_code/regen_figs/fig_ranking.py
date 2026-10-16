#!/usr/bin/env python3
"""Fig. ranking (P3): factor shares of ln kappa_L under the registered and the
measured SiC contrast, for the full factorial and for the 62-run doubled
covering batch, at one-column print width.

Inputs (committed files only): p3_sigec/results/p3_sigec_results.json
(registered contrast, batch importance P3-6) and
p3_sigec/results/p3_sigec_measured_contrast.json (main-effect shares of the
full factorial under both contrasts; batch importance under the measured one).
"""
import json

import numpy as np

import paper_fig_style as S

S.apply_style()
import matplotlib.pyplot as plt  # noqa: E402

R = S.PROJECT / "p3_sigec" / "results"
reg = json.load(open(R / "p3_sigec_results.json"))
mea = json.load(open(R / "p3_sigec_measured_contrast.json"))
fac = ["d_eff_um", "x_Ge", "carbon_form", "y_C_at_frac", "T_K", "dopant"]
lab = ["grain size", "Ge fraction", "carbon form", "carbon content", "temperature", "dopant"]
imp = reg["predictions"]["P3-6"]["importance"]
tot = sum(imp.values())
series = [
    ("full factorial, registered contrast", [mea["main_effect_pct_of_var_lnkappa_full_factorial"]["stiffness_registered"][f] for f in fac], "#9ECAE1", None),
    ("62-run batch, registered contrast", [100 * imp[f] / tot for f in fac], "#0072B2", None),
    ("full factorial, measured contrast", [mea["main_effect_pct_of_var_lnkappa_full_factorial"]["measured_curve"][f] for f in fac], "#D8C3E3", None),
    ("62-run batch, measured contrast", [mea["importance_pct_CA2_ref"][f] for f in fac], "#6A3D9A", None),
]
fig, ax = plt.subplots(figsize=(S.COL, 2.9))
y = np.arange(len(fac))
h = 0.2
for i, (name, v, c, _) in enumerate(series):
    ax.barh(y + (i - 1.5) * h, v, height=h, color=c, label=name)
ax.set_yticks(y); ax.set_yticklabels(lab)
ax.invert_yaxis()
ax.set_xlabel("share of $\\ln\\kappa_L$ (%)")
ax.set_xlim(0, 85)
ax.grid(axis="y", visible=False)
ax.legend(loc="lower right", fontsize=7, handlelength=1.0, borderaxespad=0.1)
S.save(fig, "fig_ranking")
