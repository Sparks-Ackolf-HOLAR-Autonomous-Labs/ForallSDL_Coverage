#!/usr/bin/env python3
"""Fig. SiGeC design (P3): batch designs on the SiGeC + dopant design space,
at full two-column print width.

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

fig, axs = plt.subplots(2, 2, figsize=(S.FULL, 4.3))
ax = axs.ravel()

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
a.legend(loc="lower right", handlelength=1.4)
a.set_title("a", loc="left", fontweight="bold")

# (b) screening regret
b = ax[2]
lab = ["one factor\nat a time", "random", "covering\n1-fold", "random", "covering\n2-fold"]
keys = [None, "RND_nCA1", "CA1", "RND_nCA2", "CA2"]
for i, k in enumerate(keys):
    if k is None:
        b.plot([i], [100 * dd["OFAT"]["regret"]], "^", color=C_OF, ms=5.5)
        continue
    lo, hi = dd[k]["regret_q25_q75"]
    col = C_CA if k.startswith("CA") else C_RND
    b.vlines(i, 100 * lo, 100 * hi, color=col, lw=5, alpha=0.5)
    b.plot([i], [100 * dd[k]["regret_median"]], "o", color=col, ms=4.5)
b.set_xticks(range(5))
b.set_xticklabels(lab)
b.grid(axis="x", visible=False)
b.set_ylabel("Excess $\\kappa_L$ of the chosen\nlow-Ge composition (%)")
b.text(1.5, 37, f"{res['sizes']['CA1_ref']} runs", ha="center", fontsize=7)
b.text(3.5, 37, f"{res['sizes']['CA2_ref']} runs", ha="center", fontsize=7)
for xv in (0.5, 2.5):
    b.axvline(xv, color="#BBBBBB", lw=0.6)
b.set_ylim(-2, 40)
b.set_title("c", loc="left", fontweight="bold")

# (c) erasure tolerance
c = ax[1]
er = res["erasure"]
ks = [0, 1, 2, 3]
for k, col, mk, ls, lb in (("CA1", C_CA, "o", "-", "covering, 1-fold"),
                           ("CA2", C_CA, "s", "--", "covering, 2-fold"),
                           ("RND_nCA2", C_RND, "d", "-", "random, size of 2-fold")):
    v = [100.0 * (1.0 if k != "RND_nCA2" else 0.0)] + \
        [100 * er[k][f"k{j}"]["complete_frac"] for j in (1, 2, 3)]
    c.plot(ks, v, marker=mk, ls=ls, color=col, ms=4, lw=1.1,
           mfc="white" if k == "CA2" else col, label=lb,
           zorder=2 if k == "RND_nCA2" else 3)
c.set_xticks(ks)
c.set_xlabel("Runs lost")
c.set_ylabel("Batches still sampling every\ntwo-factor combination (%)")
c.set_ylim(-5, 105)
c.set_xlim(-0.2, 3.2)
c.legend(loc="center right", handlelength=2.0)
c.set_title("b", loc="left", fontweight="bold")

# (d) factor importance
d = ax[3]
imp = dd["CA2_ref"]["importance"]
names = {"d_eff_um": "grain size", "carbon_form": "carbon form", "x_Ge": "Ge content",
         "y_C_at_frac": "C content", "T_K": "temperature", "dopant": "B/P dopant"}
order = sorted(imp, key=imp.get)
tot = sum(imp.values())
d.barh(range(len(order)), [100 * imp[k] / tot for k in order], color=C_CA, height=0.6)
d.set_yticks(range(len(order)))
d.set_yticklabels([names[k] for k in order])
d.grid(axis="y", visible=False)
d.tick_params(axis="y", length=0)
d.set_xlabel("Relative importance for $\\ln\\kappa_L$ (%)")
d.set_title("d", loc="left", fontweight="bold")

S.save(fig, "fig_sigec_design")
