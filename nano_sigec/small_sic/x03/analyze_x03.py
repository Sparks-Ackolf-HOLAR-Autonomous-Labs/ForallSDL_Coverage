#!/usr/bin/env python3
"""Registered analysis of the x = 0.3 embedded-SiC MD check (x03_design.json).
Refuses to run if the design changed after registration. Output: x03_results.json."""
import hashlib, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
reg = open(os.path.join(HERE, "PREREGISTERED_SHA256")).read().split()[0]
if hashlib.sha256(open(os.path.join(HERE, "x03_design.json"), "rb").read()).hexdigest() != reg:
    sys.exit("x03_design.json changed after registration")
RES = {0.3: os.path.join(HERE, "results"), 0.2: os.path.join(HERE, "..", "md", "results")}
SEEDS, RADII = (20260930, 20260931), (15, 20)


def read_run(x, tag, seed):
    base = os.path.join(RES[x], f"{tag}_s{seed}")
    kv = {}
    for line in open(base + ".log"):
        p = line.split()
        if p and p[0] in ("SIC_G1_START", "SIC_G1_END", "SIC_G3_CCOORD"):
            kv[p[0]] = [float(v) for v in p[1:]]
    return {"kappa": json.load(open(base + ".json"))["kappa_W_mK"], **kv}


runs = {x: {(t, s): read_run(x, t, s) for t in ["host"] + [f"{k}{R}" for R in RADII for k in ("sh", "p")]
            for s in SEEDS} for x in (0.2, 0.3)}
out = {"design_sha256": reg, "runs_x03": {f"{t}_s{s}": v for (t, s), v in runs[0.3].items()}, "gates": {}}
r3 = runs[0.3]
g1 = {f"p{R}_s{s}": [r3[(f"p{R}", s)]["SIC_G1_START"][0] - r3[(f"p{R}", s)]["SIC_G1_START"][1],
                     r3[(f"p{R}", s)]["SIC_G1_END"][0] - r3[(f"p{R}", s)]["SIC_G1_END"][1]]
      for R in RADII for s in SEEDS}
out["gates"]["G1"] = {"gaps_A": g1, "pass": all(min(v) >= 5.7 for v in g1.values())}
out["gates"]["G2"] = {"pass": all(np.isfinite(v["kappa"]) and v["kappa"] > 0 for v in r3.values())}
cc = {f"p{R}_s{s}": r3[(f"p{R}", s)]["SIC_G3_CCOORD"][0] for R in RADII for s in SEEDS}
out["gates"]["G3"] = {"mean_C_coordination": cc, "pass": all(v >= 3.3 for v in cc.values())}
kh = [r3[("host", s)]["kappa"] for s in SEEDS]
out["gates"]["G4"] = {"host_kappa": kh, "pass": bool(abs(kh[0] - kh[1]) / np.mean(kh) <= 0.08)}


def red(x, R, skip=()):
    rs = [1 - runs[x][(f"p{R}", s)]["kappa"] / runs[x][(f"sh{R}", s)]["kappa"] for s in SEEDS
          if (x, f"p{R}", s) not in skip]
    return {"r_seeds": rs, "r_mean": float(np.mean(rs)), "half_range": float(np.ptp(rs) / 2)}


def branch(skip=()):
    rows = {R: {"x0.2": red(0.2, R), "x0.3": red(0.3, R, skip)} for R in RADII}
    rho = float(np.mean([rows[R]["x0.3"]["r_mean"] / rows[R]["x0.2"]["r_mean"] for R in RADII]))
    b = "density_tracks" if rho < 0.25 else ("insensitive" if rho >= 0.75 else "intermediate")
    return rows, rho, b


rows, rho, b = branch()
out["reductions"] = {f"R{R}A": v for R, v in rows.items()}
out["rho"], out["branch"] = rho, b
kh2 = np.mean([runs[0.2][("host", s)]["kappa"] for s in SEEDS])
out["pred"] = {"PrX1": {"r_R20": rows[20]["x0.3"]["r_mean"], "pass": rows[20]["x0.3"]["r_mean"] >= 0.10},
               "PrX2": {"host_x03": float(np.mean(kh)), "host_x02": float(kh2), "pass": bool(np.mean(kh) < kh2)},
               "PrX3": {"rho": rho, "pass": rho >= 0.5}}
fails = [k for k, g in out["gates"].items() if not g["pass"]]
if fails:  # registered gates_fail rule: decision only if the branch is the same without the failing runs
    bad = set()
    if "G1" in fails:
        bad |= {(0.3, k.split("_s")[0], int(k.split("_s")[1])) for k, v in g1.items() if min(v) < 5.7}
    if "G3" in fails:
        bad |= {(0.3, k.split("_s")[0], int(k.split("_s")[1])) for k, v in cc.items() if v < 3.3}
    b2 = branch(tuple(bad))[2] if bad else None
    out["gates_fail"] = {"failed": fails, "excluded_runs": sorted(map(str, bad)), "branch_without": b2}
    out["decision_taken"] = b2 == b
else:
    out["decision_taken"] = True
json.dump(out, open(os.path.join(HERE, "x03_results.json"), "w"), indent=1, default=float)
print(json.dumps({k: out[k] for k in ("rho", "branch", "decision_taken")}, default=float))
print({k: g["pass"] for k, g in out["gates"].items()}, {k: v["pass"] for k, v in out["pred"].items()})
