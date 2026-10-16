#!/usr/bin/env python3
"""Fig. strain (P3): what the covering-designed strained-film campaign
measured, and how far a half-size batch drifts from it, at full two-column
print width.

a  The anisotropic strain response. Each run's ln kappa is moved to the
   reference condition (x = 0.5, 700 K) with the pre-registered fit: its
   composition and temperature baseline is removed, and so is the difference
   between its own strain term and the strain term at the reference
   condition. What is left is the run's strain effect at the reference
   condition plus its residual, so the runs can be read against one line.
   Filled markers: the 25 runs of the single covering subset; open markers:
   the 25 runs that complete the doubled batch.
b  The six substrate calls of the main text (Si0.5Ge0.5 on Si(001) and on
   relaxed Ge, 300 and 1100 K, cross-plane and in-plane) re-estimated from
   25 runs. Histogram: largest deviation from the 50-run values over 1000
   random 25-run subsets of the same 50 runs; line: the 25-run covering
   subset; diamond and bar: median and range of the covering subset after any one
   of its runs is lost.
   This comparison was made after the fact (2026-10-11); its numbers are
   written to p3_strain/halfbatch_posthoc.json.

Inputs (committed files only): p3_strain/analysis_strain.json,
p3_strain/strain_batch_design.json and p3_strain/results/ (per-run kappa).
"""
import json
import sys

import numpy as np

import paper_fig_style as S

S.apply_style()
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

HERE = S.PROJECT / "p3_strain"
sys.path.insert(0, str(HERE))
import analyze_strain as a  # noqa: E402

res = json.load(open(HERE / "analysis_strain.json"))
c = np.array([res["coef"][n] for n in a.NAMES])

missing = a.committed()
rows = []
for p in a.D["points"]:
    if p["run_id"] in missing:
        continue
    j = json.load(open(HERE / "results" / (p["run_id"] + ".json")))
    rows.append(dict(p, kappa=j["kappa_W_mK"]))
assert len(rows) == res["n_completed"]
c_check, _ = a.fit(rows)
assert np.allclose(np.round(c_check, 4), c), "fit does not reproduce analysis_strain.json"

X0, T0 = 0.5, 700
COLS = {0: "#0072B2", 1: "#D55E00"}
LAB = {0: "cross-plane", 1: "in-plane"}

# ---- panel b numbers --------------------------------------------------------
CALLS = [(-0.021, 300, 0), (-0.021, 1100, 0), (0.02, 300, 0), (0.02, 1100, 0),
         (-0.021, 300, 1), (-0.021, 1100, 1), (0.02, 300, 1), (0.02, 1100, 1)]


def calls(rs):
    Xm = a.X(rs)
    if np.linalg.matrix_rank(Xm) < Xm.shape[1]:
        return None
    cc, *_ = np.linalg.lstsq(Xm, np.log([r["kappa"] for r in rs]), rcond=None)
    return np.array([100 * (np.exp(a.dlnk(cc, X0, T, e, d)) - 1) for e, T, d in CALLS])


ref = calls(rows)
c1 = [r for r in rows if r["in_C1"]]
dev_c1 = float(np.abs(calls(c1) - ref).max())
dev_c1_loss = np.array([np.abs(calls(c1[:i] + c1[i + 1:]) - ref).max() for i in range(len(c1))])
rng = np.random.default_rng(20261011)
dev_rand, cov_rand = [], []
for _ in range(1000):
    idx = rng.choice(len(rows), len(c1), replace=False)
    sub = [rows[i] for i in idx]
    cv = calls(sub)
    dev_rand.append(np.nan if cv is None else float(np.abs(cv - ref).max()))
    cov_rand.append(a.coverage(sub)[0])
dev_rand = np.array(dev_rand)
post = {
    "note": "after-the-fact comparison (2026-10-11), not pre-registered",
    "calls_eps_T_dir": CALLS, "calls_50run_pct": np.round(ref, 2).tolist(),
    "covering25_max_abs_dev_pp": round(dev_c1, 2),
    "covering25_one_loss_max_abs_dev_pp": [round(float(dev_c1_loss.min()), 2),
                                           round(float(np.median(dev_c1_loss)), 2),
                                           round(float(dev_c1_loss.max()), 2)],
    "random25_n": 1000, "random25_unidentifiable": int(np.isnan(dev_rand).sum()),
    "random25_max_abs_dev_pp_median": round(float(np.nanmedian(dev_rand)), 2),
    "random25_max_abs_dev_pp_p90": round(float(np.nanpercentile(dev_rand, 90)), 2),
    "share_random25_worse_than_covering25": round(float(np.nanmean(dev_rand > dev_c1)), 3),
    "share_random25_worse_than_covering25_one_loss_median":
        round(float(np.nanmean(dev_rand > np.median(dev_c1_loss))), 3),
    "random25_pair_coverage_median": int(np.median(cov_rand)),
    "pairs_total": a.coverage(rows)[1],
}
json.dump(post, open(HERE / "halfbatch_posthoc.json", "w"), indent=1)
print(json.dumps(post, indent=1))

# ---- figure ------------------------------------------------------------------
fig, (ax, bx) = plt.subplots(1, 2, figsize=(S.FULL, 2.6),
                             gridspec_kw={"width_ratios": [1.25, 1]})

e = np.linspace(-0.03, 0.01, 81)
boot_rng = np.random.default_rng(20260929)
boot = np.array([a.fit([rows[i] for i in boot_rng.integers(0, len(rows), len(rows))])[0]
                 for _ in range(2000)])
ax.axhline(0, color=S.GREY, lw=0.7)
ax.axvline(0, color=S.GREY, lw=0.7)
for d in (0, 1):
    curves = np.array([100 * (np.exp(a.dlnk(cb, X0, T0, e, d)) - 1) for cb in boot])
    lo, hi = np.percentile(curves, [2.5, 97.5], axis=0)
    ax.fill_between(100 * e, lo, hi, color=COLS[d], alpha=0.18, lw=0)
    ax.plot(100 * e, 100 * (np.exp(a.dlnk(c, X0, T0, e, d)) - 1), color=COLS[d], lw=1.4)
for r in rows:
    base = c[0] + c[1] * r["x"] + c[2] * r["x"] ** 2 + c[3] * r["T"] / 1000
    shift = a.dlnk(c, r["x"], r["T"], r["eps"], r["dir"]) - a.dlnk(c, X0, T0, r["eps"], r["dir"])
    y = 100 * (np.exp(np.log(r["kappa"]) - base - shift) - 1)
    col = COLS[r["dir"]]
    ax.scatter(r["eps"] * 100 + (0.07 if r["dir"] else -0.07), y, s=16, zorder=3,
               facecolor=col if r["in_C1"] else "white", edgecolor=col, linewidth=0.9)
ax.set_xlabel("Biaxial in-plane strain (%)")
ax.set_ylabel(r"$\kappa$ change from unstrained (%)")
ax.set_title("a", loc="left", fontweight="bold")
hd = [Line2D([], [], color=COLS[d], lw=1.4, label=LAB[d]) for d in (0, 1)]
hd += [Line2D([], [], ls="", marker="o", ms=4, mfc=S.GREY, mec=S.GREY, label="covering subset (25)"),
       Line2D([], [], ls="", marker="o", ms=4, mfc="white", mec=S.GREY, label="added for doubling (25)")]
ax.set_ylim(-13, 41)
ax.legend(handles=hd, loc="upper right", frameon=False, ncol=2, columnspacing=1.0)

ok = dev_rand[~np.isnan(dev_rand)]
bx.hist(ok, bins=np.arange(0, np.ceil(ok.max()) + 0.5, 0.5), color=S.LIGHT, edgecolor=S.GREY,
        linewidth=0.5, label="random 25 of the 50 runs")
ytop = np.histogram(ok, bins=np.arange(0, np.ceil(ok.max()) + 0.5, 0.5))[0].max() * 1.08
bx.errorbar([np.median(dev_c1_loss)], [ytop], xerr=[[np.median(dev_c1_loss) - dev_c1_loss.min()],
            [dev_c1_loss.max() - np.median(dev_c1_loss)]], fmt="D", ms=4, color="#009E73",
            mfc="white", capsize=2.5, lw=1.0, label="covering subset, one run lost")
bx.axvline(dev_c1, color="#009E73", lw=1.6, label="covering subset (25 runs)")
bx.set_xlabel("Largest error of the eight calls (% points)")
bx.set_ylim(0, ytop * 1.35)
bx.set_ylabel("Random batches")
bx.set_title("b", loc="left", fontweight="bold")
bx.legend(loc="upper right", frameon=False, bbox_to_anchor=(1.0, 0.66))
S.save(fig, "fig_strain")
