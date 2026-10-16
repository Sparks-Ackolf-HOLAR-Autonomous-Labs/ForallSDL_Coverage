#!/usr/bin/env python3
"""Figure for radius_contrast_results.json: 3.3 in wide, fonts >= 8 pt.
(a) effective contrast lam2 per MD seed against particle radius, with the
fitted radius forms (x = 0.2 calibration, registered frame);
(b) Ge saving at equal kappa_L (x = 0.2, 0.3 um grains, 1 at% C) at 300 K
per fitted form, R = 2.5 and 5 nm, bars = seed recombinations."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import fit_radius_contrast as F

HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, "radius_contrast_results.json")))
plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
                     "legend.fontsize": 8, "font.family": "DejaVu Sans", "axes.linewidth": 0.6})
INK, MUTED = "#222222", "#888888"
COL = {"const": "#555555", "a_core_fraction": "#0072B2", "b_core_shell": "#009E73", "c_power_law": "#D55E00"}
LAB = {"const": "constant", "a_core_fraction": "(a) core fraction", "b_core_shell": "(b) core-shell",
       "c_power_law": "(c) power law"}
LS = {"const": ":", "a_core_fraction": "--", "b_core_shell": "-.", "c_power_law": "-"}

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(3.3, 4.6), gridspec_kw={"height_ratios": [1.15, 1]})

# (a)
Rg = np.linspace(0.8, 6.0, 300)
for m in COL:
    p = D["fits"][f"x0.2_registered_{m}"]["params"]
    lam = np.clip(F.MODELS[m]["f"](Rg, p, F.SHELL_NM), 1e-7, None)
    ax1.plot(Rg, lam, LS[m], color=COL[m], lw=1.2, label=LAB[m])
for key, mk, fc, lab in (("x0.2_registered", "o", INK, "x = 0.2"), ("x0.2_host_matched", "o", "white", "x = 0.2, host"),
                         ("x0.3_registered", "^", INK, "x = 0.3")):
    Rs, ls = [], []
    for R in (10, 15, 20):
        e = D["effective_lam2"].get(f"{key}_R{R}A")
        if e:
            Rs += [R / 10] * 2
            ls += e["lam2_seeds"]
    ax1.plot(Rs, ls, mk, ms=4.5, mfc=fc, mec=INK, mew=0.8, ls="none", label=lab)
ax1.axhline(D["lam2_density_x0.2"], color="#AAAAAA", lw=2.0, zorder=0)
for R in (2.5, 5.0):
    ax1.axvline(R, color="#CCCCCC", lw=0.6, zorder=0)
ax1.set_yscale("log")
ax1.set_xlim(0.8, 6.0)
ax1.set_ylim(3e-5, 0.4)
ax1.set_xlabel("SiC particle radius R (nm)")
ax1.set_ylabel(r"effective contrast $\lambda^2$")
h, l = ax1.get_legend_handles_labels()
from matplotlib.lines import Line2D
h1 = h[:4] + [Line2D([], [], color="#AAAAAA", lw=2.0)]
l1 = l[:4] + ["P1 (density edge)"]
leg = ax1.legend(h1, l1, loc="lower left", frameon=False, handlelength=1.8, borderaxespad=0.1,
                 labelspacing=0.15, prop={"size": 8})
ax1.add_artist(leg)
ax1.legend(h[4:], l[4:], loc="upper right", title="MD seeds", title_fontsize=8, frameon=False,
           handletextpad=0.2, borderaxespad=0.1, labelspacing=0.15, prop={"size": 8})
ax1.text(0.02, 0.97, "(a)", transform=ax1.transAxes, va="top", fontweight="bold")

# (b)
models = list(COL)
y = np.arange(len(models))[::-1]
ax2.axvspan(27, 33, color="#E8E8E8", zorder=0)
ax2.text(30, len(models) - 0.45, "P1: 27-33%", ha="center", va="bottom", fontsize=8, color="#444444")
for i, m in zip(y, models):
    f = D["fits"][f"x0.2_registered_{m}"]
    for Rn, dy, mk in (("2.5", 0.14, "s"), ("5", -0.14, "o")):
        k = f"ge_saving_pct_R{Rn}_T300"
        v, (lo, hi) = f["claims"][k], f["seed_recombinations"][k]
        ax2.errorbar(v, i + dy, xerr=[[v - lo], [hi - v]], fmt=mk, ms=4.5, color=COL[m], mfc=COL[m] if Rn == "5" else "white",
                     mew=0.9, elinewidth=0.8, capsize=0)
ax2.plot([], [], "s", mfc="white", mec=INK, ms=4.5, label="R = 2.5 nm")
ax2.plot([], [], "o", color=INK, ms=4.5, label="R = 5 nm")
ax2.set_yticks(y)
ax2.set_yticklabels([LAB[m] for m in models])
ax2.set_ylim(-0.6, len(models) - 0.1)
ax2.set_xlim(0, 45)
ax2.set_xlabel("Ge saving at equal $\\kappa_L$, 300 K (%)")
ax2.legend(loc="lower right", frameon=False, prop={"size": 8}, handletextpad=0.2, borderaxespad=0.1)
ax2.text(-0.02, 1.02, "(b)", transform=ax2.transAxes, va="bottom", ha="right", fontweight="bold")
for ax in (ax1, ax2):
    ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout(pad=0.3, h_pad=0.8)
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(HERE, f"fig_radius_contrast.{ext}"), dpi=300)
print("saved")
