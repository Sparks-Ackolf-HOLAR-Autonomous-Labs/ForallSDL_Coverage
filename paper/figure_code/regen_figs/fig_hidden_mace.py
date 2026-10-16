#!/usr/bin/env python3
"""Fig. hidden MACE-MD test (P3): the 89 MACE-MP-0 molecular dynamics runs of
the SiGe(C) box and how each batch scores on the 16 hidden recipes, at full
two-column print width.

Inputs (committed files only): p3_mace/raw/out/*.json (89 runs), the planner
designs and hidden set, and p3_mace/results/p3_mace_results.json. Panel c
re-draws the 200 random 15-run batches with the seeded generator of
analyze_p3_mace.py and checks that its median matches the committed result.
"""
import itertools
import json
import random
import sys

import numpy as np

import paper_fig_style as S

S.apply_style()
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(S.PROJECT / "p3_mace"))
import analyze_p3_mace as A  # noqa: E402  (checks the registered design hash)
D = A.D

res = json.load(open(S.PROJECT / "p3_mace" / "results" / "p3_mace_results.json"))
post = json.load(open(S.PROJECT / "p3_mace" / "results" / "p3_mace_posthoc.json"))
k, lost = A.load()
runs = []
for pid in k:
    d = json.load(open(A.OUT + "/" + pid + ".json"))
    runs.append(d)

des = {n: [tuple(r) for r in v] for n, v in A.DES.items()}
hidden = [tuple(r) for r in A.HID["hidden"] if "H_" + D.pid(tuple(r)) in k]
rng = random.Random(64)
pool = [r for r in des["RND_pool"] if D.pid(r) in k]
rnd15 = [A.score(rng.sample(pool, 15), k, hidden)["log_rmse"] for _ in range(200)]
assert abs(np.median(rnd15) - res["scores"]["RND15"]["log_rmse"]) < 1e-9, "random batches differ from the committed result"

T_COL = {300: "#56B4E9", 700: "#E69F00", 1100: "#D55E00"}
fig, axs = plt.subplots(1, 3, figsize=(S.FULL, 2.45), gridspec_kw={"width_ratios": [1.1, 1.0, 1.0]})

# (a) all runs: kappa against Ge content
a = axs[0]
jit = {"none": -0.012, "dispersed": 0.0, "clustered": 0.012}
for d in runs:
    form = d["carbon_form"] if d["y_C"] > 0 else "none"
    mk = "o" if d["y_C"] == 0 else ("^" if form == "dispersed" else "s")
    a.plot(d["x_Ge"] + jit[form], d["kappa_W_mK"], mk, ms=3.2, mfc=T_COL[d["T_K"]] if d["y_C"] == 0 else "white",
           mec=T_COL[d["T_K"]], mew=0.8, alpha=0.9)
for T, c in T_COL.items():
    a.plot([], [], "o", color=c, ms=3.5, label=f"{T} K")
a.plot([], [], "o", color=S.GREY, ms=3.5, label="no carbon")
a.plot([], [], "^", mfc="white", mec=S.GREY, ms=3.5, label="C on sites")
a.plot([], [], "s", mfc="white", mec=S.GREY, ms=3.5, label="C cluster")
a.set_xlabel("Ge fraction $x$")
a.set_ylabel("$\\kappa$ (W m$^{-1}$ K$^{-1}$)")
a.legend(fontsize=6, ncol=2, loc="upper right", handletextpad=0.2, columnspacing=0.6)
a.set_ylim(0.25, 1.05)

# (b) matched-pair effect of each factor against the seed-to-seed scatter
b = axs[1]
allk = {}
for d in runs:
    r = (D.F["x_Ge"].index(d["x_Ge"]), D.F["y_C"].index(d["y_C"]),
         (D.F["carbon_form"].index(d["carbon_form"]) if d["y_C"] > 0 else 0),
         D.F["dopant"].index(d["dopant"]), D.F["T_K"].index(d["T_K"]))
    if d["pid"].startswith("R_"):
        continue
    allk[r] = np.log(d["kappa_W_mK"])
names = ["Ge fraction", "temperature", "carbon content", "carbon form", "dopant"]
idx = [0, 4, 1, 2, 3]
effects = []
for j in idx:
    v = []
    for r1, r2 in itertools.combinations(allk, 2):
        diff = [i for i in range(5) if r1[i] != r2[i]]
        if diff == [j] and not (j == 2 and r1[1] == 0):
            v.append(abs(allk[r1] - allk[r2]))
    effects.append(v)
noise = res["noise_floor_abs_log_diff"]
b.axvspan(0, max(noise), color=S.LIGHT, zorder=0)
b.text(max(noise) / 2, -0.75, "repeat\nseeds", fontsize=6, color=S.GREY, ha="center", va="center")
rs = np.random.default_rng(1)
for i, v in enumerate(effects):
    yy = i + rs.uniform(-0.18, 0.18, len(v))
    b.plot(v, yy, "o", ms=2.6, color="#0072B2", alpha=0.6, mew=0)
    b.plot([np.median(v)] * 2, [i - 0.3, i + 0.3], color=S.INK, lw=1.4)
b.set_yticks(range(5)); b.set_yticklabels([f"{n} ({len(v)})" for n, v in zip(names, effects)])
b.set_ylim(4.5, -1.1)
b.set_xlabel("$|\\Delta \\ln\\kappa|$, one factor changed")
b.grid(axis="y", visible=False)

# (c) hidden-recipe error of each batch
c = axs[2]
hn, _, _ = c.hist(rnd15, bins=18, color=S.GREY, alpha=0.55, label="random, 15 runs (200)")
top = hn.max() * 1.1
sc = res["scores"]
marks = [("covering, 15", sc["CA1"]["log_rmse"], "#0072B2", "-"),
         ("covering, 30", sc["CA2"]["log_rmse"], "#0072B2", "--"),
         ("one factor at a time, 12", sc["OFAT"]["log_rmse"], "#D55E00", "-"),
         ("covering 15, Ge and T only", post["x_and_T_only_log_rmse"]["CA1"], "#009E73", ":")]
for lab, v, col, ls in marks:
    c.vlines(v, 0, top, color=col, ls=ls, lw=1.4, label=lab)
c.set_xlabel("hidden-recipe RMSE in $\\ln\\kappa$")
c.set_ylabel("random batches")
c.legend(fontsize=5.8, loc="upper right", handlelength=1.6)
c.set_ylim(0, top * 1.75)
for ax, l in zip(axs, "abc"):
    ax.text(-0.18 if l != "b" else -0.42, 1.04, l, transform=ax.transAxes, fontweight="bold", fontsize=9)
S.save(fig, "fig_hidden_mace")
print("effects n:", [len(v) for v in effects], "medians:", [round(float(np.median(v)), 3) for v in effects])
