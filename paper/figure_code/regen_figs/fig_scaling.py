#!/usr/bin/env python3
"""ESI Fig. scaling (P3): covering batch size and the coverage of an
equal-size random batch against the number of factors, at full two-column
print width.

Input (committed file only): scalability/scalability_results.json (greedy
constructor that never lists the full factorial; 3 and 5 levels per factor).
"""
import json

import numpy as np

import paper_fig_style as S

S.apply_style()
import matplotlib.pyplot as plt  # noqa: E402

d = json.load(open(S.PROJECT / "scalability" / "scalability_results.json"))["P3"]
col = {3: "#009E73", 5: "#0072B2"}
fig, axs = plt.subplots(1, 2, figsize=(S.FULL, 2.3))
a, b = axs
for v in (3, 5):
    r = [x for x in d if x["v"] == v]
    k = [x["k"] for x in r]
    a.plot(k, [x["CA1"] for x in r], "o-", color=col[v], ms=3.5, label=f"{v} levels, every pair once")
    a.plot(k, [x["CA2"] for x in r], "s--", color=col[v], ms=3.5, mfc="white", label=f"{v} levels, every pair twice")
    a.axhline(r[0]["lower_bound"], color=col[v], lw=0.8, ls=":")
    b.plot(k, [100 * x["random_same_n_pair_coverage"] for x in r], "o-", color=col[v], ms=3.5, label=f"{v} levels")
a.set_xscale("log", base=2)
a.set_xticks([4, 8, 16, 32]); a.set_xticklabels(["4", "8", "16", "32"])
a.set_xlabel("number of factors")
a.set_ylabel("runs in the covering batch")
a.set_ylim(0, 132)
a.legend(fontsize=6, loc="upper left", ncol=2, columnspacing=0.8)
b.set_xscale("log", base=2)
b.set_xticks([4, 8, 16, 32]); b.set_xticklabels(["4", "8", "16", "32"])
b.set_xlabel("number of factors")
b.set_ylabel("pairs sampled by a random\nbatch of the same size (%)")
b.set_ylim(70, 100)
b.axhline(100, color=S.INK, lw=0.8)
b.text(4, 98.3, "covering batch", fontsize=6.5, color=S.INK)
b.legend(fontsize=6.5, loc="lower right")
for ax, l in zip(axs, "ab"):
    ax.text(-0.16, 1.0, l, transform=ax.transAxes, fontweight="bold", fontsize=9)
S.save(fig, "fig_scaling")
