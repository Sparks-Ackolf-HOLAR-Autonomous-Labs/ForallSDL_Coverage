#!/usr/bin/env python3
"""Fig. SiGeC design (P3): batch designs on the SiGeC + dopant design space,
at full two-column print width.

a  pair coverage against batch size; b  coverage after lost runs;
c  screening excess by design. (Factor importance of the 62-run batch is
   drawn once, in fig_ranking.)

Inputs (committed files only): p3_sigec/results/p3_sigec_results.json and
the registered design (p3_sigec/run_p3_sigec.py refuses to import if the
design changed after registration). Panel a draws 200 random batches per
size with the same seeded generator as the original figure script.
"""
import json
import random
import sys

import numpy as np

import paper_fig_style as S

S.apply_style()
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(S.PROJECT / "p3_sigec"))
import run_p3_sigec as R  # noqa: E402

res = json.load(open(S.PROJECT / "p3_sigec" / "results" / "p3_sigec_results.json"))
dd = res["designs"]
C_CA, C_RND, C_OF = "#0072B2", S.GREY, "#D55E00"

fig, ax = plt.subplots(1, 3, figsize=(S.FULL, 2.55), gridspec_kw={"width_ratios": [1.0, 1.0, 1.2]})

# (a) pair coverage against batch size
rng = random.Random(3)
ns = [10, 15, 20, 23, 30, 33, 40, 50, 62, 80, 100, 130]
q = np.array([np.quantile([R.coverage(rng.sample(R.FULL, n)) for _ in range(200)],
                          [0.1, 0.5, 0.9]) for n in ns]) * 100
a = ax[0]
a.fill_between(ns, q[:, 0], q[:, 2], color=C_RND, alpha=0.3, lw=0)
a.plot(ns, q[:, 1], color=C_RND, lw=1.2, label="random batches")
a.plot([res["sizes"]["CA1_ref"]], [100], "o", color=C_CA, ms=5, label="covering, 1-fold")
a.plot([res["sizes"]["CA2_ref"]], [100], "s", color=C_CA, ms=5, mfc="white", mew=1.1,
       label="covering, 2-fold")
a.plot([res["sizes"]["OFAT"]], [100 * dd["OFAT"]["coverage"]], "^", color=C_OF, ms=5.5,
       label="one factor at a time")
a.set_xlabel("Runs in the batch")
a.set_ylabel("Two-factor combinations\nsampled (%)")
a.set_xlim(0, 135)
a.set_ylim(0, 105)
a.legend(loc="lower right", bbox_to_anchor=(1.04, 0.0), handlelength=1.0, handletextpad=0.4,
         borderaxespad=0.1)
a.set_title("a", loc="left", fontweight="bold")

# (b) erasure tolerance
c = ax[1]
er = res["erasure"]
ks = np.array([0, 1, 2, 3])
for k, col, mk, ls, lb in (("CA1", C_CA, "o", "-", "covering, 1-fold"),
                           ("CA2", C_CA, "s", "--", "covering, 2-fold")):
    v = [100.0] + [100 * er[k][f"k{j}"]["complete_frac"] for j in (1, 2, 3)]
    c.plot(ks, v, marker=mk, ls=ls, color=col, ms=4, lw=1.1,
           mfc="white" if k == "CA2" else col, label=lb, zorder=3)
# random batch of the 2-fold size: never complete, before or after losses.
# Drawn just right of the covering points so it is not hidden under them.
vr = [0.0] + [100 * er["RND_nCA2"][f"k{j}"]["complete_frac"] for j in (1, 2, 3)]
c.plot(ks + 0.12, vr, marker="d", ls=":", color=C_RND, mec=S.INK, mew=0.5, ms=4.5, lw=1.0,
       zorder=4, clip_on=False)
c.annotate("random, size of 2-fold:\n0% in every trial", xy=(2.12, 3.0), xytext=(1.9, 24),
           ha="center", fontsize=7, color=S.INK,
           arrowprops=dict(arrowstyle="-", color=S.INK, lw=0.6))
c.set_xticks(ks)
c.set_xlabel("Runs lost")
c.set_ylabel("Batches still sampling every\ntwo-factor combination (%)")
c.set_ylim(-5, 105)
c.set_xlim(-0.2, 3.3)
c.legend(loc="center right", handlelength=1.4, handletextpad=0.4, bbox_to_anchor=(1.04, 0.6),
         borderaxespad=0.0)
c.set_title("b", loc="left", fontweight="bold")

# (c) screening excess
b = ax[2]
n1, n2, no = res["sizes"]["CA1_ref"], res["sizes"]["CA2_ref"], res["sizes"]["OFAT"]
lab = [f"one factor\nat a time, {no}", f"random, {n1}", f"covering\n1-fold, {n1}",
       f"random, {n2}", f"covering\n2-fold, {n2}"]
keys = [None, "RND_nCA1", "CA1", "RND_nCA2", "CA2"]
for i, k in enumerate(keys):
    if k is None:
        b.plot([100 * dd["OFAT"]["regret"]], [i], "^", color=C_OF, ms=5.5)
        continue
    lo, hi = dd[k]["regret_q25_q75"]
    col = C_CA if k.startswith("CA") else C_RND
    b.hlines(i, 100 * lo, 100 * hi, color=col, lw=5, alpha=0.5)
    b.plot([100 * dd[k]["regret_median"]], [i], "o", color=col, ms=4.5)
b.set_yticks(range(5))
b.set_yticklabels(lab)
b.invert_yaxis()
b.grid(axis="y", visible=False)
b.tick_params(axis="y", length=0)
b.set_xlabel("Excess $\\kappa_L$ of pick (%)")
for yv in (0.5, 2.5):
    b.axhline(yv, color="#BBBBBB", lw=0.6)
b.set_xlim(-1, 37)
b.set_title("c", loc="left", fontweight="bold")

S.save(fig, "fig_sigec_design")
