#!/usr/bin/env python3
"""Fig. slice (P3): what three batches of the same size sample in one
two-factor slice of the SiGeC design space (carbon content x carbon form),
at full two-column print width.

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

sys.path.insert(0, str(S.PROJECT / "p3_sigec"))
import run_p3_sigec as R  # noqa: E402

des = json.load(open(S.PROJECT / "p3_sigec" / "results" / "designs.json"))
rnd = random.Random(3).sample(R.FULL, len(des["CA1_ref"]))
batches = [("one factor at a time", des["OFAT"]), ("random", rnd), ("covering", des["CA1_ref"])]
A, B = R.NAMES.index("y_C_at_frac"), R.NAMES.index("carbon_form")
c_lab = ["0", "0.5", "1", "2", "3"]
f_lab = ["on site", "2.5", "5", "25", "100"]
best = (4, 1)  # 3 at% C as 2.5 nm SiC: the carbon levels of the best recipe (response min)

cmap = ListedColormap(["#FFFFFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"])
fig, axs = plt.subplots(1, 3, figsize=(S.FULL, 2.45))
for ax, (name, rows) in zip(axs, batches):
    g = np.zeros((5, 5), int)
    for r in rows:
        g[r[A], r[B]] += 1
    ax.imshow(np.minimum(g, 4), cmap=cmap, vmin=0, vmax=4, origin="lower")
    for i in range(5):
        for j in range(5):
            if g[i, j] == 0:
                ax.text(j, i, "–", ha="center", va="center", color="#B22222", fontsize=8)
            else:
                ax.text(j, i, str(g[i, j]), ha="center", va="center", fontsize=7,
                        color="white" if g[i, j] >= 3 else S.INK)
    ax.add_patch(plt.Rectangle((best[1] - 0.5, best[0] - 0.5), 1, 1, fill=False, ec="#D55E00", lw=1.6))
    ax.set_xticks(range(5)); ax.set_xticklabels(f_lab)
    ax.set_yticks(range(5)); ax.set_yticklabels(c_lab)
    ax.set_xlabel("carbon form (SiC radius, nm)")
    ax.grid(False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(True)
    cov = 100 * R.coverage(rows)
    ax.set_title(f"{name}, {len(rows)} runs\nall pairs sampled: {cov:.0f}%", fontsize=8)
axs[0].set_ylabel("carbon content (at%)")
for ax, l in zip(axs, "abc"):
    ax.text(-0.22, 1.13, l, transform=ax.transAxes, fontweight="bold", fontsize=9)
S.save(fig, "fig_slice")
