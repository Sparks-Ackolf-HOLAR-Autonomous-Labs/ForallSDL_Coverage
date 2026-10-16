#!/usr/bin/env python3
"""Registered analysis of the R = 5 nm fixed-volume-fraction SiC check
(r50_design.json). Refuses to run if the design changed after registration.

Same inversion as ../r30/analyze_r30.py: each seed-paired reduction
r = 1 - kappa(particle)/kappa(pocket-only) -> lam2 with the model side of
../../radius_contrast/fit_radius_contrast.py (registered frame). Then

  ratio_R5 = lam2(R = 5.0 nm, phi 10%) / lam2(R = 3.0 nm, phi 10%, ledger 104)

over the 4 seed pairings, the branch of the decision rule, and the claims
recomputed with the measured curve extended to 5 nm. Output: r50_results.json."""
import hashlib, itertools, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
reg = open(os.path.join(HERE, "PREREGISTERED_SHA256")).read().split()[0]
if hashlib.sha256(open(os.path.join(HERE, "r50_design.json"), "rb").read()).hexdigest() != reg:
    sys.exit("r50_design.json changed after registration")
NS = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(NS, "radius_contrast"))
import fit_radius_contrast as FR  # noqa: E402

T, X = 300, 0.2
SEEDS_NEW, SEEDS_R30 = (20261007, 20261008), (20261006, 20261007)
THRESH = float(np.sqrt((5 / 3) ** -3.9))
RES = os.path.join(HERE, "results")
RES30 = os.path.join(os.path.dirname(HERE), "r30", "results")


def read_run(res, tag, seed):
    base = os.path.join(res, f"{tag}_s{seed}")
    kv = {}
    for line in open(base + ".log"):
        p = line.split()
        if p and p[0] in ("SIC_G1_START", "SIC_G1_END", "SIC_G3_CCOORD", "MPT_BOX"):
            kv[p[0]] = [float(v) for v in p[1:]]
    return {"kappa": json.load(open(base + ".json"))["kappa_W_mK"], **kv}


runs = {(t, s): read_run(RES, t, s) for t in ("sh50", "p50") for s in SEEDS_NEW}
runs30 = {(t, s): read_run(RES30, t, s) for t in ("sh30", "p30") for s in SEEDS_R30}
out = {"design_sha256": reg, "runs": {f"{t}_s{s}": v for (t, s), v in runs.items()}, "gates": {}}
parts = [("p50", s) for s in SEEDS_NEW]
g1 = {f"{t}_s{s}": [runs[(t, s)]["SIC_G1_START"][0] - runs[(t, s)]["SIC_G1_START"][1],
                    runs[(t, s)]["SIC_G1_END"][0] - runs[(t, s)]["SIC_G1_END"][1]] for t, s in parts}
out["gates"]["G1"] = {"gaps_A": g1, "pass": all(min(v) >= 5.7 for v in g1.values())}
out["gates"]["G2"] = {"pass": all(np.isfinite(v["kappa"]) and v["kappa"] > 0 for v in runs.values())}
cc = {f"{t}_s{s}": runs[(t, s)]["SIC_G3_CCOORD"][0] for t, s in parts}
out["gates"]["G3"] = {"mean_C_coordination": cc, "pass": all(v >= 3.3 for v in cc.values())}
g4 = [runs[("sh50", s)]["kappa"] for s in SEEDS_NEW]
out["gates"]["G4"] = {"pocket_kappa": g4, "pass": abs(g4[0] - g4[1]) / np.mean(g4) <= 0.10}

base, _ = FR.scaled_modes(X)


def lam2_points(r_seeds, phi, k_shell, R_A):
    d = FR.bisect(lambda dd: base.kappa2(T, dd) - k_shell, 1e-11, 1e-2)
    k0 = base.kappa2(T, d)
    tab = np.array([1 - base.kappa2(T, d, np_R=R_A * 1e-10, np_phi=phi, lam2=l) / k0 for l in FR.LAM_GRID])
    return [FR.lam_of({"k": tab}, "k", r) for r in r_seeds], d


def point(rr, R, seeds):
    rs = [1 - rr[(f"p{R}", s)]["kappa"] / rr[(f"sh{R}", s)]["kappa"] for s in seeds]
    box = rr[(f"p{R}", seeds[0])]["MPT_BOX"]
    phi = 4 / 3 * np.pi * R ** 3 / (box[0] * box[1] * 7 * box[2] / 20)  # analyze_md.py convention
    ksh = float(np.mean([rr[(f"sh{R}", s)]["kappa"] for s in seeds]))
    return {"r_seeds": rs, "phi_eff": float(phi), "k_shell": ksh}


old, _ = FR.md_points()
pts = {"R50_phi10_new": (point(runs, 50, SEEDS_NEW), 50), "R30_phi10_ledger104": (point(runs30, 30, SEEDS_R30), 30),
       "R20_phi10_seq59": (old[X][20], 20), "R15_phi4_seq59": (old[X][15], 15), "R10_phi1_seq59": (old[X][10], 10)}
lam = {}
for k, (p, R) in pts.items():
    ls, d = lam2_points(p["r_seeds"], p["phi_eff"], p["k_shell"], R)
    lam[k] = ls
    out.setdefault("points", {})[k] = {**p, "R_A": R, "lam2_seeds": ls, "d_registered_m": d}

rR5 = [None if (x is None or y is None or y == 0 or not np.isfinite(x)) else x / y
       for x, y in itertools.product(lam["R50_phi10_new"], lam["R30_phi10_ledger104"])]
out["ratios"] = {"ratio_R5": rR5, "threshold": THRESH}


def branch_of(rr):
    v = [x for x in rr if x is not None]
    if len(v) < len(rr) or not v:
        return "undecided"
    if all(x < THRESH for x in v):
        return "keeps_falling"
    if all(x >= THRESH for x in v):
        return "stops_falling"
    return "undecided"


b = branch_of(rR5)
out["branch"] = b


def gmean(v):
    v = [x for x in v if x is not None and np.isfinite(x) and x > 0]
    return float(np.exp(np.mean(np.log(v)))) if v else None


curve = sorted([(1.0, gmean(lam["R10_phi1_seq59"])), (1.5, gmean(lam["R15_phi4_seq59"])),
                (2.0, gmean(lam["R20_phi10_seq59"])), (3.0, gmean(lam["R30_phi10_ledger104"])),
                (5.0, gmean(lam["R50_phi10_new"]))])
out["measured_curve"] = curve


def lamfun(Rn):
    xs = np.log([c[0] for c in curve]); ys = np.log([c[1] for c in curve])
    if Rn <= 5.0:
        return float(np.exp(np.interp(np.log(Rn), xs, ys)))
    slope = min(0.0, (ys[-1] - ys[-2]) / (xs[-1] - xs[-2]))
    return float(np.exp(ys[-1] + slope * (np.log(Rn) - xs[-1])))


if all(c[1] is not None for c in curve):
    out["claims_measured_curve"] = FR.claims(lamfun)
    out["pick_measured_curve"] = FR.pick(lamfun)
s5 = out.get("claims_measured_curve", {}).get("ge_saving_pct_R5_T300")
out["pred"] = {
    "PrF1": {"ratio_R5": rR5, "pass": b == "keeps_falling"},
    "PrF2": {"r_R50": out["points"]["R50_phi10_new"]["r_seeds"],
             "pass": bool(np.mean(out["points"]["R50_phi10_new"]["r_seeds"]) >= 0.05)},
    "PrF3": {"ge_saving_pct_R5_T300": s5, "pass": s5 is not None and abs(s5 - 8.695288209342621) <= 3.0},
}
fails = [k for k, g in out["gates"].items() if not g["pass"]]
if fails:  # registered gates_fail rule
    bad = set()
    if "G1" in fails:
        bad |= {k for k, v in g1.items() if min(v) < 5.7}
    if "G3" in fails:
        bad |= {k for k, v in cc.items() if v < 3.3}
    keep = [i for i, (s_new, _) in enumerate(itertools.product(SEEDS_NEW, SEEDS_R30)) if f"p50_s{s_new}" not in bad]
    b2 = branch_of([rR5[i] for i in keep]) if bad and keep else ("undecided" if bad else None)
    out["gates_fail"] = {"failed": fails, "excluded_runs": sorted(bad), "branch_without": b2}
    out["decision_taken"] = (b2 == b) if bad else ("G2" not in fails and "G4" not in fails)
else:
    out["decision_taken"] = True
json.dump(out, open(os.path.join(HERE, "r50_results.json"), "w"), indent=1, default=float)
print(json.dumps({"branch": b, "decision_taken": out["decision_taken"], "ratio_R5": rR5}, default=float))
print({k: g["pass"] for k, g in out["gates"].items()}, {k: v["pass"] for k, v in out["pred"].items()})
