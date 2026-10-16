#!/usr/bin/env python3
"""Analysis of the embedded-SiC MD check (md_design.json, ledger seq 58).

Run only after commit_md.py has ledgered md/results/. Refuses to run if
md_design.json changed after registration. Applies the gates, the three
predictions and the decision rule exactly as registered; every choice the
design leaves open is listed under "analysis_time_choices" in the output.
Output: md_results.json."""
import hashlib, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DESIGN = os.path.join(HERE, "md_design.json")
reg = open(os.path.join(HERE, "PREREGISTERED_SHA256")).read().split()[0]
if hashlib.sha256(open(DESIGN, "rb").read()).hexdigest() != reg:
    sys.exit("md_design.json changed after registration")
sys.path.insert(0, os.path.dirname(HERE))
os.environ.setdefault("SIGEC_M1_SUMMARY", os.path.join(HERE, "..", "m1c", "m1c_summary.json"))
import analyze_sigec as S  # noqa: E402

RES = os.environ.get("SMALL_SIC_RES", os.path.join(HERE, "md", "results"))
X, T = 0.2, 300
LAM_DEN = S.contrast(X, "density")
LAM_STF = S.contrast(X, "stiffness")
LAM_CAP = 1e6
K_TERSOFF_BULK = 1 / (0.209 + 1.919 * X - 2.244 * X ** 2 + 8.35e-5 * T)  # 1.894 W/mK, P1 final law
CHOICES = ["first analysis pass (md_results_firstpass_scalebug.json) did not apply the trial factor inside the intrinsic-rate bisection, so f hit its bound (1e3); fixed after that pass, before any interpretation; kept for the record"]


def read_run(tag, seed):
    base = os.path.join(RES, f"{tag}_s{seed}")
    j = json.load(open(base + ".json"))
    kv = {}
    for line in open(base + ".log"):
        p = line.split()
        if p and p[0] in ("SIC_G1_START", "SIC_G1_END", "SIC_G3_CCOORD", "MPT_BOX"):
            kv[p[0]] = [float(v) for v in p[1:]]
    return {"kappa": j["kappa_W_mK"], "box": kv["MPT_BOX"], "g1s": kv.get("SIC_G1_START"),
            "g1e": kv.get("SIC_G1_END"), "ccoord": kv.get("SIC_G3_CCOORD", [None])[0],
            "dT_half": j["deltaT_across_half_box_K"]}


pts = [l.split() for l in open(os.path.join(HERE, "md", "points.txt")) if l.strip()]
runs = {}
for tag, R, _, sic, seed in pts:
    runs[(tag, int(seed))] = dict(read_run(tag, seed), R_A=float(R), sic=int(sic))
SEEDS = sorted({k[1] for k in runs})
RADII = [10, 15, 20]
out = {"design_sha256": reg, "runs": {f"{k[0]}_s{k[1]}": v for k, v in runs.items()}, "gates": {}, "pred": {}}

# ---- gates
g1 = {}
for (tag, s), r in runs.items():
    if r["sic"]:
        gaps = [r["g1s"][0] - r["g1s"][1], r["g1e"][0] - r["g1e"][1]]
        g1[f"{tag}_s{s}"] = {"gap_start_A": gaps[0], "gap_end_A": gaps[1]}
out["gates"]["G1"] = {"rows": g1, "pass": all(min(v.values()) >= 5.7 for v in g1.values())}
out["gates"]["G2"] = {"pass": all(np.isfinite(r["kappa"]) and r["kappa"] > 0 for r in runs.values())}
cc = {f"{k[0]}_s{k[1]}": r["ccoord"] for k, r in runs.items() if r["sic"]}
out["gates"]["G3"] = {"mean_C_coordination": cc, "pass": all(v is not None and v >= 3.3 for v in cc.values())}
kh = [runs[("host", s)]["kappa"] for s in SEEDS]
out["gates"]["G4"] = {"host_kappa": kh, "rel_diff": abs(kh[0] - kh[1]) / np.mean(kh),
                      "pass": bool(abs(kh[0] - kh[1]) / np.mean(kh) <= 0.08)}
CHOICES.append("G4 'agree within 8%' read as |k1 - k2| / mean(k1, k2) <= 0.08")
out["gates_pass"] = all(g["pass"] for g in out["gates"].values())

# ---- MD reductions, seed-paired
md = {}
for R in RADII:
    rs = [1 - runs[(f"p{R}", s)]["kappa"] / runs[(f"sh{R}", s)]["kappa"] for s in SEEDS]
    ksh = [runs[(f"sh{R}", s)]["kappa"] for s in SEEDS]
    box = runs[(f"p{R}", SEEDS[0])]["box"]
    v_p = 4 / 3 * np.pi * R ** 3
    phi_eff = v_p / (box[0] * box[1] * 7 * box[2] / 20)
    md[R] = {"r_seeds": rs, "r_mean": float(np.mean(rs)), "half_range": abs(rs[0] - rs[1]) / 2,
             "kappa_shell_mean": float(np.mean(ksh)), "phi_eff": phi_eff}

# ---- model side: Tersoff-scaled modes, box-matched boundary length
base = S.SModes(X)
tau0 = base.tau[T].copy()


def set_f(f):
    base.tau[T] = tau0 / f


def bisect(fun, lo, hi, n=200, log=True):
    for _ in range(n):
        mid = np.sqrt(lo * hi) if log else 0.5 * (lo + hi)
        if fun(mid) > 0:
            hi = mid
        else:
            lo = mid
    return np.sqrt(lo * hi) if log else 0.5 * (lo + hi)


def _fk(f):
    set_f(f)
    return K_TERSOFF_BULK - base.kappa2(T)


f_scale = bisect(_fk, 1e-3, 1e3)
set_f(f_scale)
out["model"] = {"intrinsic_rate_factor": f_scale, "bulk_kappa_scaled": base.kappa2(T), "target": K_TERSOFF_BULK}


def r_model(R_A, d, phi, lam2, reff=1.0):
    k0 = base.kappa2(T, d)
    # effective-radius form: cross-section of radius reff*R at the same particle count
    k1 = base.kappa2(T, d, np_R=reff * R_A * 1e-10, np_phi=phi * reff ** 3, lam2=lam2)
    return 1 - k1 / k0


band = {}
for R in RADII:
    ksh = md[R]["kappa_shell_mean"]
    d = bisect(lambda dd: base.kappa2(T, dd) - ksh, 1e-11, 1e-2)
    lo, cap = r_model(R, d, md[R]["phi_eff"], LAM_DEN), r_model(R, d, md[R]["phi_eff"], LAM_CAP)
    stf = r_model(R, d, md[R]["phi_eff"], LAM_STF)
    lam_eff = []
    for r in md[R]["r_seeds"]:
        if r <= lo:
            lam_eff.append(LAM_DEN if np.isclose(r, lo) else None)
        elif r >= cap:
            lam_eff.append(None)
        else:
            lam_eff.append(bisect(lambda l: r_model(R, d, md[R]["phi_eff"], l) - r, 1e-8, LAM_CAP))
    m, h = md[R]["r_mean"], md[R]["half_range"]
    where = "above_cap" if m - h > cap else ("below_density" if m + h < lo else "inside_band")
    band[R] = {"boundary_length_m": d, "r_model_density": lo, "r_model_stiffness": stf, "r_model_cap": cap,
               "lam2_eff_seeds": lam_eff, "classification": where}
    md[R].update(band[R])
out["md_vs_model"] = {f"R{R}A": md[R] for R in RADII}
CHOICES.append("lam2_eff is solved per seed; a seed outside [density, cap] gets no lam2_eff (None)")
CHOICES.append("classification uses the seed mean with the two-seed half range as the tolerance")

# ---- predictions
out["pred"]["PrS1"] = {"r_R20": md[20]["r_mean"], "pass": md[20]["r_mean"] >= 0.05}
out["pred"]["PrS2"] = {"pass": all(md[R]["classification"] == "inside_band" for R in RADII)}
ksh10 = md[10]["kappa_shell_mean"]
out["pred"]["PrS3"] = {"shell_over_host": ksh10 / np.mean(kh), "pass": bool(abs(1 - ksh10 / np.mean(kh)) < 0.05)}

# ---- decision rule, on the unscaled almaBTE modes of the SiGeC campaign
modes = {x: S.SModes(x) for x in S.XS}
D_EFF, Y = 0.3e-6, 0.01


def claims(R_nm, lam2, reff=1.0):
    row = {}
    for TT in (300, 1100):
        k0 = modes[0.2].kappa2(TT, D_EFF)
        kC = lambda x: modes[x].kappa2(TT, D_EFF, 0.0, np_R=reff * R_nm * 1e-9,
                                       np_phi=S.phi(x, Y) * reff ** 3, lam2=lam2)
        row[f"red_pct_T{TT}"] = 100 * (1 - kC(0.2) / k0)
        xm = S.xmin_for([kC(x) for x in S.XS], k0)
        row[f"ge_saving_pct_T{TT}"] = None if xm is None else 100 * (1 - xm / 0.2)
    return row


cls = [md[R]["classification"] for R in RADII]
dec = {"classifications": dict(zip([f"R{R}A" for R in RADII], cls))}
if "above_cap" in cls:
    dec["branch"] = "above_cap"
    err = lambda c: sum((r_model(R, md[R]["boundary_length_m"], md[R]["phi_eff"], LAM_STF, c) - md[R]["r_mean"]) ** 2
                        for R in RADII)
    cs = np.linspace(1.0, 4.0, 601)
    c = float(cs[np.argmin([err(c) for c in cs])])
    CHOICES.append("above_cap: R_eff/R fitted with the stiffness contrast (the physical upper bracket), "
                   "same particle count, least squares over the three radii, grid 1.0 to 4.0")
    dec["R_eff_over_R"] = c
    dec["claims"] = {f"R{Rn:g}nm": claims(Rn, LAM_STF, c) for Rn in (2.5, 5.0)}
elif "below_density" in cls:
    dec["branch"] = "below_density"
    dec["claims"] = {f"R{Rn:g}nm": claims(Rn, LAM_DEN) for Rn in (2.5, 5.0)}
    CHOICES.append("below_density: MD-bounded saving reported at the density contrast (the model's lower edge)")
else:
    dec["branch"] = "inside_band"
    ls = [l for R in RADII for l in md[R]["lam2_eff_seeds"] if l is not None]
    dec["lam2_range"] = [min(ls), max(ls)]
    dec["claims"] = {f"R{Rn:g}nm": {"lam2_min": claims(Rn, min(ls)), "lam2_max": claims(Rn, max(ls))}
                     for Rn in (2.5, 5.0)}
dec["prior_band"] = {f"R{Rn:g}nm": {"density": claims(Rn, LAM_DEN), "stiffness": claims(Rn, LAM_STF)}
                     for Rn in (2.5, 5.0)}
out["decision"] = dec
out["analysis_time_choices"] = CHOICES
json.dump(out, open(os.path.join(HERE, "md_results.json"), "w"), indent=1, default=float)
print(json.dumps({"gates": {k: v["pass"] for k, v in out["gates"].items()},
                  "pred": {k: v["pass"] for k, v in out["pred"].items()},
                  "md": {R: (md[R]["r_mean"], md[R]["half_range"], md[R]["r_model_density"], md[R]["r_model_cap"],
                             md[R]["classification"]) for R in RADII},
                  "decision": {k: v for k, v in dec.items() if k != "prior_band"}}, indent=1, default=float))
