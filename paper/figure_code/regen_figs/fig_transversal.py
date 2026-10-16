#!/usr/bin/env python3
"""ESI Fig. transversal (P3): where the nine transversal follow-up runs went
and what they changed, at full two-column print width.

Inputs (committed files only): p3_sigec/results/p3_sigec_results.json
(transversal block) and p3_sigec/results/designs.json (62-run batch).
"""
import json
import sys

import numpy as np

import paper_fig_style as S

S.apply_style()
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(S.PROJECT / "p3_sigec"))
import run_p3_sigec as R  # noqa: E402

res = json.load(open(S.PROJECT / "p3_sigec" / "results" / "p3_sigec_results.json"))
des = json.load(open(S.PROJECT / "p3_sigec" / "results" / "designs.json"))
tv = res["transversal"]
F = R.F
nice = {"x_Ge": "Ge", "y_C_at_frac": "C content", "carbon_form": "C form", "d_eff_um": "grain"}

fig, axs = plt.subplots(1, 3, figsize=(S.FULL, 2.4), gridspec_kw={"width_ratios": [1.0, 1.15, 0.95]})
# (a) uncertain cells by factor pair
a = axs[0]
items = sorted(tv["uncertain_cells_by_pair"].items(), key=lambda kv: kv[1])
a.barh(range(len(items)), [v for _, v in items], color="#0072B2")
a.set_yticks(range(len(items)))
a.set_yticklabels([" × ".join(nice[p] for p in k.split("|")) for k, _ in items])
a.set_xlabel(f"uncertain cells ({tv['n_uncertain_cells']} of {tv['n_cells']})")
a.grid(axis="y", visible=False)

# (b) follow-up runs on the carbon form x grain size plane
b = axs[1]
fi, gi = R.NAMES.index("carbon_form"), R.NAMES.index("d_eff_um")
g = np.zeros((len(F["carbon_form"]), len(F["d_eff_um"])), int)
for r in des["CA2_ref"]:
    g[r[fi], r[gi]] += 1
b.imshow(g.T, cmap="Greys", vmin=0, vmax=12, origin="lower", aspect="auto")
for i in range(g.shape[0]):
    for j in range(g.shape[1]):
        b.text(i, j + 0.32, str(g[i, j]), ha="center", va="center", fontsize=6, color=S.GREY)
tcol = {300: "#56B4E9", 500: "#009E73", 700: "#E69F00", 900: "#CC79A7", 1100: "#D55E00"}
seen = {}
slots = [(-0.22, -0.05), (0.0, -0.05), (0.22, -0.05), (-0.22, -0.3), (0.0, -0.3), (0.22, -0.3)]
used = {}
for p in tv["followup_points"]:
    i, j = F["carbon_form"].index(p["carbon_form"]), F["d_eff_um"].index(p["d_eff_um"])
    dx, dy = slots[used.get((i, j), 0)]
    used[(i, j)] = used.get((i, j), 0) + 1
    lab = f"{p['T_K']} K" if p["T_K"] not in seen else None
    seen[p["T_K"]] = 1
    b.plot(i + dx, j + dy, "o", ms=4.5, color=tcol[p["T_K"]],
           mec="white", mew=0.5, label=lab)
b.set_xticks(range(5)); b.set_xticklabels(["on site", "2.5", "5", "25", "100"])
b.set_yticks(range(4)); b.set_yticklabels([str(v) for v in F["d_eff_um"]])
b.set_xlabel("carbon form (SiC radius, nm)")
b.set_ylabel("effective grain size ($\\mu$m)")
b.grid(False)
b.legend(fontsize=6, loc="lower left", bbox_to_anchor=(0.45, 1.0), ncol=2, handletextpad=0.1, borderpad=0.1, columnspacing=0.6)

# (c) model error and decision after nine follow-up runs
c = axs[2]
rows = [("none", tv["base"]), ("transversal", tv["transversal"]), ("random (median)", tv["random_median"]),
        ("most uncertain", tv["top_sd"])]
x = np.arange(len(rows))
c.barh(x, [v["log_rmse"] for _, v in rows], color=["#BBBBBB", "#0072B2", S.GREY, "#D55E00"], height=0.6)
for i, (_, v) in enumerate(rows):
    c.text(v["log_rmse"] + 0.0004, i, f"pick +{100 * v['regret']:.0f}%", va="center", fontsize=6)
c.set_yticks(x); c.set_yticklabels([n for n, _ in rows])
c.invert_yaxis()
c.set_xlim(0.118, 0.141)
c.set_xlabel("model error ($\\ln\\kappa_L$)")
c.set_title("nine follow-up runs", fontsize=7.5)
c.grid(axis="y", visible=False)
for ax, l in zip(axs, "abc"):
    ax.text(-0.12 if l == "b" else (-0.55 if l == "a" else -0.6), 1.03, l, transform=ax.transAxes, fontweight="bold", fontsize=9)
S.save(fig, "fig_transversal")
