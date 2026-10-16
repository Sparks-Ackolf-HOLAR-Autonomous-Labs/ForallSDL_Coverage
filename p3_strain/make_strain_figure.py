#!/usr/bin/env python3
"""Figure for the P3 strain seed batch: strain part of ln kappa (fitted
composition/temperature baseline removed) vs in-plane strain, by
transport geometry. Usage: python3 make_strain_figure.py"""
import json
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import analyze_strain as a

HERE = os.path.dirname(os.path.abspath(__file__))
res = json.load(open(os.path.join(HERE, "analysis_strain.json")))
c = [res["coef"][n] for n in a.NAMES]
fig, ax = plt.subplots(figsize=(6.4, 4.4))
cols = {0: "#1b6ca8", 1: "#e07b39"}
lab = {0: "cross-plane (heat along the film normal)", 1: "in-plane (heat along the film)"}
for p in a.D["points"]:
    k = json.load(open(os.path.join(HERE, "results", p["run_id"] + ".json")))["kappa_W_mK"]
    base = c[0] + c[1] * p["x"] + c[2] * p["x"] ** 2 + c[3] * p["T"] / 1000
    y = 100 * (np.exp(np.log(k) - base) - 1)
    ax.scatter(p["eps"] * 100 + (0.08 if p["dir"] else -0.08), y, color=cols[p["dir"]],
               s=18 + 0.02 * p["T"], alpha=0.8, edgecolor="white", linewidth=0.5)
e = np.linspace(-0.03, 0.01, 50)
for d in (0, 1):
    ax.plot(100 * e, 100 * (np.exp(a.dlnk(c, 0.5, 700, e, d)) - 1), color=cols[d], lw=2, label=lab[d] + ", fit at x=0.5, 700 K")
ax.axhline(0, color="0.6", lw=0.8)
ax.axvline(0, color="0.6", lw=0.8)
ax.set_xlabel("biaxial in-plane strain (%)   [SiGe on Si(001): about -4.2 x %]")
ax.set_ylabel("kappa change vs unstrained (%)")
ax.set_title("Strained SiGe films (Tersoff-1989 MD, 50 runs): the response is anisotropic", fontsize=9.5)
ax.legend(fontsize=7.5, frameon=False, loc="lower left")
ax.text(0.99, 0.97, "marker size grows with T (300 to 1100 K)", transform=ax.transAxes,
        ha="right", va="top", fontsize=7, color="0.4")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "fig_strain.png"), dpi=200)
