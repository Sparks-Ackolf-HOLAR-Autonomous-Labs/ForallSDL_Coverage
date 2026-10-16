#!/usr/bin/env python3
"""Fig. reserve (P3): how three rules for spending nine spare runs protect the
reported best recipe when some runs read low, at one-column print width.

Input (committed file only): p3_reserve/results/reserve_results.json
(60 trials per fault rate, registered design reserve_design.json).
Bars are Wilson 95% intervals over the 60 trials.
"""
import json

import numpy as np

import paper_fig_style as S

S.apply_style()
import matplotlib.pyplot as plt  # noqa: E402

r = json.load(open(S.PROJECT / "p3_reserve" / "results" / "reserve_results.json"))
rules = [("CA2+random", "random recipes", S.GREY, "o"),
         ("CA2+exploit", "lowest-predicted recipes", "#6A3D9A", "^"),
         ("CA2+reserve", "replication reserve", "#0072B2", "s")]
fs = sorted(r["trials"], key=float)


def wilson(k, n, z=1.96):
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


fig, axs = plt.subplots(2, 1, figsize=(S.COL, 3.6), sharex=True)
off = {"CA2+random": -0.006, "CA2+exploit": 0.0, "CA2+reserve": 0.006}
for key, lab, col, mk in rules:
    x, y, lo, hi, reg = [], [], [], [], []
    for f in fs:
        t = r["trials"][f]
        k = sum(1 for tr in t if tr[key]["false_claim"])
        a, b = wilson(k, len(t))
        x.append(100 * float(f) + 100 * off[key]); y.append(100 * k / len(t))
        lo.append(100 * a); hi.append(100 * b)
        reg.append(100 * np.quantile([tr[key]["surrogate_regret"] for tr in t], [0.25, 0.5, 0.75]))
    y = np.array(y)
    axs[0].errorbar(x, y, yerr=[y - np.array(lo), np.array(hi) - y], fmt=mk + "-", color=col, ms=4,
                    capsize=2, lw=1.2, label=lab)
    reg = np.array(reg)
    axs[1].errorbar(x, reg[:, 1], yerr=[reg[:, 1] - reg[:, 0], reg[:, 2] - reg[:, 1]], fmt=mk + "-",
                    color=col, ms=4, capsize=2, lw=1.2)
    print(lab, "median and IQR of excess (%):", np.round(reg, 1).tolist())
axs[0].set_ylabel("trials reporting a\nwrong best value (%)")
axs[0].legend(fontsize=7, loc="upper left")
axs[0].set_ylim(0, 100)
axs[0].text(0.99, 0.99, "n = 60 trials per point", transform=axs[0].transAxes, ha="right", va="top",
            fontsize=7)
axs[1].set_ylabel("excess $\\kappa_L$ of the\nrecommended recipe (%)")
axs[1].set_xlabel("share of runs that read 30–60% low (%)")
axs[1].set_xticks([0, 5, 10])
axs[1].set_ylim(0, max(40, float(np.ceil(axs[1].get_ylim()[1] / 10) * 10)))
for ax, l in zip(axs, "ab"):
    ax.text(-0.25, 1.0, l, transform=ax.transAxes, fontweight="bold", fontsize=9)
S.save(fig, "fig_reserve")
