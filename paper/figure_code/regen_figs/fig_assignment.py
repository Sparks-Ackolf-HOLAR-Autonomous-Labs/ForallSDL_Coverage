#!/usr/bin/env python3
"""ESI Fig. assignment (P3): the lowest-cost assignment of evidence sources to
composition-temperature cells, and its cost under five scenarios, at full
two-column print width.

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

fig, axs = plt.subplots(1, 2, figsize=(S.FULL, 2.75), gridspec_kw={"width_ratios": [1.55, 1.0]})
a = axs[0]
for i, xc in enumerate(XC):
    for j, tb in enumerate(TB):
        cell = f"{xc}/{tb}"
        a.add_patch(Rectangle((i, j), 1, 1, fc="white", ec="#CCCCCC", lw=0.6))
        if cell in fixed:
            a.add_patch(Rectangle((i, j), 1, 1, fc=S.SOURCE_COLOURS["database"], alpha=0.22, lw=0))
        if cell in A["infeasible_cells"]:
            a.text(i + 0.5, j + 0.5, "×", ha="center", va="center", fontsize=12, color="#B22222")
        got = [k for k, _, _ in src if cell in A["optimal"]["assignment"].get(k, [])]
        for n, k in enumerate(got):
            col = dict((s[0], s[2]) for s in src)[k]
            a.add_patch(Rectangle((i + 0.12 + 0.38 * n, j + 0.3), 0.34, 0.4, fc=col, ec="white", lw=0.4))
a.set_xlim(0, len(XC)); a.set_ylim(0, len(TB))
a.set_xticks(np.arange(len(XC)) + 0.5); a.set_xticklabels(xl, fontsize=6.3)
a.set_yticks(np.arange(len(TB)) + 0.5); a.set_yticklabels(tl)
a.set_xlabel("Ge fraction $x$"); a.set_ylabel("temperature (K)")
a.grid(False)
h = [Rectangle((0, 0), 1, 1, fc=col) for k, _, col in src if A["optimal"]["assignment"].get(k)]
l = [lab for k, lab, _ in src if A["optimal"]["assignment"].get(k)]
h.append(Rectangle((0, 0), 1, 1, fc=S.SOURCE_COLOURS["database"], alpha=0.22)); l.append("held experimental values")
a.legend(h, l, fontsize=6, loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=3, handlelength=0.9, columnspacing=0.8)

b = axs[1]
names = {"A_base": "base", "B_fidelity": "converged\nestimates", "C_measured_gpumd": "measured HNEMD\nvalidity",
         "D_bte_T4_ok": "BTE accepted\n≥1100 K", "E_campaign_T_bands": "campaign bands\nonly (21 cells)"}
keys = list(names)
y = np.arange(len(keys))
b.barh(y - 0.18, [c["scenarios"][k]["optimal"]["cost_h"] for k in keys], height=0.34, color="#0072B2", label="exact")
b.barh(y + 0.18, [c["scenarios"][k]["greedy"]["cost_h"] for k in keys], height=0.34, color=S.GREY, label="greedy")
for i, k in enumerate(keys):
    g = c["scenarios"][k]["greedy"]
    b.text(g["cost_h"] + 0.8, i + 0.18, f"+{100 * g['gap_rel']:.0f}%", va="center", fontsize=6)
b.set_yticks(y); b.set_yticklabels([names[k] for k in keys], fontsize=6.3)
b.invert_yaxis()
b.set_xlabel("computation (h)")
b.set_xlim(0, 75)
b.grid(axis="y", visible=False)
b.legend(fontsize=6, loc="lower left", bbox_to_anchor=(0.1, 1.0), ncol=2, handlelength=0.9)
for ax, lab, xo in zip(axs, "ab", (-0.1, -0.5)):
    ax.text(xo, 1.0, lab, transform=ax.transAxes, fontweight="bold", fontsize=9)
S.save(fig, "fig_assignment")
