#!/usr/bin/env python3
"""ESI Fig. assignment (P3): the lowest-cost assignment of evidence sources to
composition-temperature cells in scenario A, at full two-column print width.
(The costs of the five scenarios are given once, in ESI Table S4.)

Input (committed file only): p3_prework/covering/results.json (integer
program and greedy rule, scenarios A-E).
"""
import json

import numpy as np

import paper_fig_style as S

S.apply_style()
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

c = json.load(open(S.PROJECT / "p3_prework" / "covering" / "results.json"))
XC, TB = c["grid"]["x_columns"], c["grid"]["T_bands"]
xl = ["Si", "0.05–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–0.95", "Ge"]
tl = ["300–500", "500–700", "700–900", "900–1100", "≥1100"]
A = c["scenarios"]["A_base"]
src = [("DGX-MD", "Tersoff MD", S.SOURCE_COLOURS["tersoff_md"]),
       ("GPUMD-HNEMD", "Tersoff HNEMD (GPU)", "#F0E442"),
       ("MACE-MD", "MACE-MP MD", S.SOURCE_COLOURS["mace_md"]),
       ("BTE-DFT", "DFT phonon BTE", S.SOURCE_COLOURS["bte_dft"]),
       ("CHPC-MD", "Tersoff MD (CHPC)", "#A6611A")]
fixed = c["fixed_coverage"]

fig, a = plt.subplots(figsize=(S.FULL, 2.6))
for i, xc in enumerate(XC):
    for j, tb in enumerate(TB):
        cell = f"{xc}/{tb}"
        a.add_patch(Rectangle((i, j), 1, 1, fc="white", ec="#CCCCCC", lw=0.6))
        if cell in fixed:
            a.add_patch(Rectangle((i, j), 1, 1, fc=S.SOURCE_COLOURS["database"], alpha=0.22, lw=0))
        if cell in A["infeasible_cells"]:
            a.text(i + 0.8, j + 0.5, "×", ha="center", va="center", fontsize=12, color="black")
        got = [k for k, _, _ in src if cell in A["optimal"]["assignment"].get(k, [])]
        for n, k in enumerate(got):
            col = dict((s[0], s[2]) for s in src)[k]
            a.add_patch(Rectangle((i + 0.12 + 0.38 * n, j + 0.3), 0.34, 0.4, fc=col, ec="white", lw=0.4))
a.set_xlim(0, len(XC)); a.set_ylim(0, len(TB))
a.set_xticks(np.arange(len(XC)) + 0.5); a.set_xticklabels(xl)
a.set_yticks(np.arange(len(TB)) + 0.5); a.set_yticklabels(tl)
a.set_xlabel("Ge fraction $x$"); a.set_ylabel("temperature (K)")
a.grid(False)
h = [Rectangle((0, 0), 1, 1, fc=col) for k, _, col in src if A["optimal"]["assignment"].get(k)]
l = [lab for k, lab, _ in src if A["optimal"]["assignment"].get(k)]
h.append(Rectangle((0, 0), 1, 1, fc=S.SOURCE_COLOURS["database"], alpha=0.22)); l.append("held experimental values")
h.append(plt.Line2D([], [], ls="", marker="x", color="black", ms=6)); l.append("no source pair can cover")
a.legend(h, l, fontsize=7, loc="center left", bbox_to_anchor=(1.02, 0.5), ncol=1, handlelength=0.9)
S.save(fig, "fig_assignment")
