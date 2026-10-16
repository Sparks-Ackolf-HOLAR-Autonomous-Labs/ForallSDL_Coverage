#!/usr/bin/env python3
"""One figure at print size (17.1 cm wide, two-column RSC width) from the results JSON files."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import common as C

plt.rcParams.update({"font.size": 7, "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
                     "legend.fontsize": 6, "axes.spines.top": False, "axes.spines.right": False})
CC = json.load(open(os.path.join(C.OUT, "cost_covering.json")))
S3 = json.load(open(os.path.join(C.OUT, "strength3.json")))
LR = json.load(open(os.path.join(C.OUT, "loss_replay.json")))
fig, ax = plt.subplots(1, 3, figsize=(17.1 / 2.54, 5.6 / 2.54), constrained_layout=True)

# a: cost against screening score, single batches
cols = {"lab_shared_pellet": "#4C72B0", "lab_pellet_per_run": "#DD8452", "simulation_evidence": "#55A868"}
lab = {"lab_shared_pellet": "lab, pellet shared", "lab_pellet_per_run": "lab, pellet per run",
       "simulation_evidence": "simulation"}
mk = {"uniform": "o", "cost_weighted": "s", "ilp": "D", "random_at_cost_weighted_cost": "x"}
a = ax[0]
for m, o in CC["models"].items():
    d = o["single"]
    base = d["uniform"]["cost"]["median"]
    for k, s in mk.items():
        if k not in d:
            continue
        c = d[k]["cost"]["median"] if "cost" in d[k] and isinstance(d[k]["cost"], dict) else d[k].get("cost", d[k].get("budget"))
        a.scatter(c / base, d[k]["regret"]["median"] * 100, marker=s, color=cols[m], s=16,
                  facecolor="none" if k == "uniform" else cols[m], lw=0.9)
for m in cols:
    a.plot([], [], color=cols[m], lw=3, label=lab[m])
for k, s in (("uniform", "uniform covering"), ("cost_weighted", "cost-weighted"), ("ilp", "minimum cost"),
             ("random_at_cost_weighted_cost", "random, same cost")):
    a.scatter([], [], marker=mk[k], color="0.3", s=14, facecolor="none" if k == "uniform" else "0.3", label=s)
a.set_xlabel("batch cost relative to uniform covering")
a.set_ylabel("median excess $\\kappa_L$ of pick (%)")
a.legend(frameon=False, loc="upper left", ncol=2, handlelength=1.0, columnspacing=0.6, borderaxespad=0.1)
a.set_xlim(0.68, 1.05)
a.set_ylim(8, 38)
a.set_title("a", loc="left", fontweight="bold")

# b: three-factor coverage against batch size
b = ax[1]
rng = np.random.default_rng(3)
ns = [20, 35, 63, 100, 150, 191, 300, 450, 600]
q = []
for n in ns:
    v = [C.coverage_t([C.FULL[i] for i in rng.choice(C.N, n, replace=False)], 3) for _ in range(30)]
    q.append(np.quantile(v, [0.1, 0.5, 0.9]))
q = np.array(q)
b.fill_between(ns, q[:, 0] * 100, q[:, 2] * 100, color="0.85", lw=0)
b.plot(ns, q[:, 1] * 100, color="0.5", lw=1, label="random batches")
D = S3["designs"]
for k, lbl, s in (("strength2_single", "pair covering", "o"), ("strength2_doubled", "doubled pair covering", "s"),
                  ("strength3", "three-factor covering", "D")):
    b.scatter(D[k]["runs"]["median"], D[k]["coverage3"]["median"] * 100, marker=s, color="#C44E52", s=16, zorder=3, label=lbl)
b.axvline(S3["lower_bound_strength3"], color="0.3", lw=0.6, ls=":")
b.set_xlim(0, 620)
b.set_xlabel("runs in batch")
b.set_ylabel("three-factor cells sampled (%)")
b.legend(frameon=False, loc="center right", handlelength=1.0)
b.set_title("b", loc="left", fontweight="bold")

# c: replayed losses, no requeue
c = ax[2]
pats = [("midrun_preemption (per-run records + zero-loss aggregates)", "mid-run preemption\n(0.5%)", "bernoulli"),
        ("p3_mace:node_fault (aggregate only)", "MACE, one faulty node\n(21%)", "bernoulli"),
        ("gpumd_endpoint_v2:node_fault_only", "endpoint, node faults\n(51%, replayed)", "replay"),
        ("gpumd_sizefree:first_attempt", "size-free, node faults\n(53%, replayed)", "replay"),
        ("gpumd_endpoint_v2:first_attempt", "endpoint, all first-\nattempt losses (87%)", "replay")]
yy = np.arange(len(pats))
for j, (dn, colr) in enumerate((("single", "#8172B3"), ("doubled", "#64B5CD"))):
    v = []
    for k, _, kind in pats:
        e = LR[kind][k][dn]
        v.append((e["construction_order"] if kind == "replay" else e)["coverage_median"] * 100)
    c.barh(yy + (0.2 if j == 0 else -0.2), v, height=0.38, color=colr,
           label={"single": "covering, 33 to 36 runs", "doubled": "doubled, 61 to 66 runs"}[dn])
c.set_yticks(yy)
c.set_yticklabels([p[1] for p in pats])
c.invert_yaxis()
c.set_xlabel("pair cells kept, median (%)")
c.set_xlim(0, 100)
c.legend(frameon=False, loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=1, borderaxespad=0.1)
c.set_title("c", loc="left", fontweight="bold", x=-0.75)
fig.savefig(os.path.join(C.OUT, "fig_p3_limits.png"), dpi=300)
fig.savefig(os.path.join(C.OUT, "fig_p3_limits.pdf"))
