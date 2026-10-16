#!/usr/bin/env python3
"""Fig. campaign loop (P3, Methods): the batch-design loop of the SiGeC
campaigns and what three batches of comparable size sample in one two-factor
slice, at full two-column print width.

a  Schematic of the loop as described in Methods: six processing factors,
   a covering batch (every pair of levels at least once, or twice for
   tolerance to one lost run), runs of the response model or molecular
   dynamics, lost runs, the follow-up and replication-reserve rules, and the
   reported recipe.
b  Runs per cell in the carbon content x carbon form slice for the
   one-factor-at-a-time batch, one random batch and the single covering
   batch (replaces the former separate slice figure).

Inputs (committed files only): the reference batches in
p3_sigec/results/designs.json and the registered design (run_p3_sigec.py
refuses to import if the design changed after registration). The random
batch is one seeded draw of 33 recipes from the full factorial.
"""
import json
import random
import sys

import numpy as np

import paper_fig_style as S

S.apply_style()
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

sys.path.insert(0, str(S.PROJECT / "p3_sigec"))
import run_p3_sigec as R  # noqa: E402

C_CA, C_OF = "#0072B2", "#D55E00"
DASH = "#444444"

HA = 1.45  # height of the schematic row, inches
fig = plt.figure(figsize=(S.FULL, HA + 3.0))
sf_a, sf_b = fig.subfigures(2, 1, height_ratios=[HA, 3.0], hspace=0.0)

# ---------------------------------------------------------------- (a) schematic
ax = sf_a.add_axes([0, 0, 1, 1])
ax.set_xlim(0, S.FULL)
ax.set_ylim(0, HA)
ax.axis("off")
boxes = [
    ("Six factors", "Ge fraction (6)\nC content (5)\ncarbon form (5)\ngrain size (4)\nB or P dopant (3)\ntemperature (5)", "#F2F2F2", S.INK),
    ("Covering batch", "each level pair\nonce (1-fold) or\ntwice (2-fold),\nfixed before\nthe first run", "#DCEAF5", C_CA),
    ("Runs", "phonon scattering\nmodel, or\nmolecular\ndynamics", "#F2F2F2", S.INK),
    ("Lost runs", "cracked pellet,\ntripped furnace,\ncancelled job;\n2-fold keeps all\npairs after one", "#F2F2F2", S.INK),
    ("Follow-up", "tree model; runs\ntouch each\nuncertain region;\nreserve repeats\nlow readings", "#F2F2F2", S.INK),
    ("Report", "lowest predicted\n$\\kappa_L$ and the\nlevel pairs\nbehind it", "#F2F2F2", S.INK),
]
n = len(boxes)
gap = 0.13
w = (S.FULL - 0.06 - gap * (n - 1)) / n
y0, h = 0.33, 1.08
xs = [0.03 + i * (w + gap) for i in range(n)]
for x, (title, body, fc, ec) in zip(xs, boxes):
    ax.add_patch(FancyBboxPatch((x, y0), w, h, boxstyle="round,pad=0,rounding_size=0.05",
                                fc=fc, ec=ec, lw=0.9))
    ax.text(x + w / 2, y0 + h - 0.1, title, ha="center", va="top", fontsize=7.5, fontweight="bold",
            color=ec if ec != S.INK else S.INK)
    ax.text(x + w / 2, y0 + h - 0.3, body, ha="center", va="top", fontsize=7, linespacing=1.15)
for i in range(n - 1):
    ax.add_patch(FancyArrowPatch((xs[i] + w + 0.015, y0 + h / 2), (xs[i + 1] - 0.015, y0 + h / 2),
                                 arrowstyle="-|>", mutation_scale=8, color=S.INK, lw=0.9))
# feedback: model -> next batch of runs
xa, xb = xs[4] + w / 2, xs[2] + w / 2
ax.add_patch(FancyArrowPatch((xa, y0 - 0.01), (xb, y0 - 0.01), connectionstyle="bar,fraction=-0.12",
                             arrowstyle="-|>", mutation_scale=8, color=S.INK, lw=0.9))
ax.text((xa + xb) / 2, 0.2, "next batch of runs", ha="center", va="center", fontsize=7)
ax.text(0.0, HA, "a", ha="left", va="top", fontweight="bold", fontsize=9)

# ---------------------------------------------------------------- (b) one slice
des = json.load(open(S.PROJECT / "p3_sigec" / "results" / "designs.json"))
rnd = random.Random(3).sample(R.FULL, len(des["CA1_ref"]))
batches = [("one factor at a time", des["OFAT"]), ("random", rnd), ("covering, 1-fold", des["CA1_ref"])]
A, B = R.NAMES.index("y_C_at_frac"), R.NAMES.index("carbon_form")
c_lab = ["0", "0.5", "1", "2", "3"]
f_lab = ["on\nsite", "2.5", "5", "25", "100"]
best = (4, 1)  # 3 at% C as 2.5 nm SiC: the carbon levels of the best recipe (response minimum)

cmap = ListedColormap(["#FFFFFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"])
axs = sf_b.subplots(1, 3)
for k, (axb, (name, rows)) in enumerate(zip(axs, batches)):
    g = np.zeros((5, 5), int)
    for r in rows:
        g[r[A], r[B]] += 1
    axb.imshow(np.minimum(g, 4), cmap=cmap, vmin=0, vmax=4, origin="lower")
    for i in range(5):
        for j in range(5):
            if g[i, j] == 0:
                axb.text(j, i, "–", ha="center", va="center", color=DASH, fontsize=9, fontweight="bold")
            else:
                axb.text(j, i, str(g[i, j]), ha="center", va="center", fontsize=7,
                         color="white" if g[i, j] >= 3 else S.INK)
    # separator between carbon on lattice sites and the SiC particle radii
    axb.axvline(0.5, color=S.INK, lw=1.6)
    axb.add_patch(plt.Rectangle((best[1] - 0.5, best[0] - 0.5), 1, 1, fill=False, ec="black", lw=1.6,
                                zorder=4))
    axb.set_xticks(range(5)); axb.set_xticklabels(f_lab)
    axb.set_yticks(range(5)); axb.set_yticklabels(c_lab)
    axb.grid(False)
    for s in ("top", "right"):
        axb.spines[s].set_visible(True)
    cov = 100 * R.coverage(rows)
    axb.set_title(f"{name}, {len(rows)} runs\nall pairs sampled: {cov:.0f}%", fontsize=8)
axs[0].set_ylabel("carbon content (at%)")
sf_b.supxlabel("carbon form: on lattice sites, or SiC particles of the given radius (nm)", fontsize=8)
sf_b.text(0.0, 1.0, "b", ha="left", va="top", fontweight="bold", fontsize=9)
S.save(fig, "fig_campaign_loop")
