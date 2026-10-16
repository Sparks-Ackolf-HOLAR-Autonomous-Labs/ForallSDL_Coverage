#!/usr/bin/env python3
"""Radius-dependent effective SiC contrast (POST HOC, not pre-registered).

The P1 decision analysis uses one effective contrast lam2 for every SiC
radius (the density edge, lam2 = 9.3e-3, after the registered MD decision,
ledger seq 61). The bond-angle analysis (../../bond_angles) found that a
nominal R = 2 nm particle is crystalline only to 1.45 nm, so lam2 may depend
on R. This script:

 1. inverts each committed MD reduction (x = 0.2: R = 1.0, 1.5, 2.0 nm;
    x = 0.3: R = 1.5, 2.0 nm; two seeds each) into an effective lam2, with
    the model side built exactly as in ../small_sic/analyze_md.py
    (Tersoff-scaled almaBTE modes, box-matched boundary length d), in two
    frames: d matched to the shell-only control per radius (registered
    frame) and d matched once to the bare host (post-hoc frame of
    ../small_sic/sensitivity_md.py);
 2. fits radius-dependent lam2(R) forms to the reductions themselves
    (least squares in r, all seeds), and to every seed recombination;
 3. propagates lam2(2.5 nm), lam2(5 nm), and lam2(R) on the full grid
    through the decision analysis of ../small_sic/decision_robustness.py:
    Ge saving at equal kappa_L (x = 0.2, 0.3 um grains, 1 at% C, 300 and
    1100 K) and the screening pick (argmin kappa_L over x, C content, R).

Reads committed files only; writes only into this folder. No MD.
Runtime: about a minute on 2 cores."""
import itertools, json, os, sys
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "2")
import numpy as np
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
NS = os.path.dirname(HERE)
SS = os.path.join(NS, "small_sic")
MT = os.path.dirname(NS)
os.environ.setdefault("SIGEC_M1_SUMMARY", os.path.join(NS, "m1c", "m1c_summary.json"))
sys.path.insert(0, NS)
import analyze_sigec as S  # noqa: E402  (checks the design hash on import; writes nothing)

T = 300
SEEDS = (20260930, 20260931)
LAM_DEN = S.contrast(0.2, "density")
LAM_STF = S.contrast(0.2, "stiffness")
LAM_GRID = np.logspace(-8, 4, 481)
SHELL_NM = 0.55          # bond_angles: crystalline core to 1.45 nm at nominal R = 2.0 nm
SHELL_ALT = (0.45, 0.65)  # sensitivity on the shell thickness

MDR = json.load(open(os.path.join(SS, "md_results.json")))
BAD = json.load(open(os.path.join(MT, "bond_angles", "bad_results.json")))


def k_tersoff(x):  # P1 final law, as in analyze_md.py
    return 1 / (0.209 + 1.919 * x - 2.244 * x ** 2 + 8.35e-5 * T)


def bisect(fun, lo, hi, n=200):
    for _ in range(n):
        mid = np.sqrt(lo * hi)
        if fun(mid) > 0:
            hi = mid
        else:
            lo = mid
    return np.sqrt(lo * hi)


def scaled_modes(x):
    m = S.SModes(x)
    tau0 = m.tau[T].copy()

    def fk(f):
        m.tau[T] = tau0 / f
        return k_tersoff(x) - m.kappa2(T)
    f = bisect(fk, 1e-3, 1e3)
    m.tau[T] = tau0 / f
    return m, f


def read_x03():
    rd = os.path.join(SS, "x03", "results")
    runs = {}
    for tag in ("host", "sh15", "p15", "sh20", "p20"):
        for s in SEEDS:
            j = json.load(open(os.path.join(rd, f"{tag}_s{s}.json")))
            box = None
            for line in open(os.path.join(rd, f"{tag}_s{s}.log")):
                p = line.split()
                if p and p[0] == "MPT_BOX":
                    box = [float(v) for v in p[1:]]
            runs[(tag, s)] = {"kappa": j["kappa_W_mK"], "box": box}
    return runs


def md_points():
    """Per-seed MD reductions with the geometry needed by the model side."""
    pts = {0.2: {}, 0.3: {}}
    for R in (10, 15, 20):
        m = MDR["md_vs_model"][f"R{R}A"]
        pts[0.2][R] = {"r_seeds": m["r_seeds"], "phi_eff": m["phi_eff"], "k_shell": m["kappa_shell_mean"]}
    kh02 = float(np.mean([MDR["runs"][f"host_s{s}"]["kappa"] for s in SEEDS]))
    x3 = read_x03()
    for R in (15, 20):
        rs = [1 - x3[(f"p{R}", s)]["kappa"] / x3[(f"sh{R}", s)]["kappa"] for s in SEEDS]
        box = x3[(f"p{R}", SEEDS[0])]["box"]
        phi = 4 / 3 * np.pi * R ** 3 / (box[0] * box[1] * 7 * box[2] / 20)  # analyze_md.py convention
        pts[0.3][R] = {"r_seeds": rs, "phi_eff": phi,
                       "k_shell": float(np.mean([x3[(f"sh{R}", s)]["kappa"] for s in SEEDS]))}
    kh03 = float(np.mean([x3[("host", s)]["kappa"] for s in SEEDS]))
    return pts, {0.2: kh02, 0.3: kh03}


def response_tables(pts, khost):
    """r_model(R, lam2) on LAM_GRID for each composition, radius and frame."""
    tab, meta = {}, {}
    for x in (0.2, 0.3):
        base, f = scaled_modes(x)
        d_host = bisect(lambda dd: base.kappa2(T, dd) - khost[x], 1e-11, 1e-2)
        meta[x] = {"intrinsic_rate_factor": f, "d_host_m": d_host}
        for R, p in pts[x].items():
            d_sh = bisect(lambda dd: base.kappa2(T, dd) - p["k_shell"], 1e-11, 1e-2)
            for frame, d in (("registered", d_sh), ("host_matched", d_host)):
                k0 = base.kappa2(T, d)
                r = np.array([1 - base.kappa2(T, d, np_R=R * 1e-10, np_phi=p["phi_eff"], lam2=l) / k0
                              for l in LAM_GRID])
                tab[(x, R, frame)] = r
                meta[x][f"d_{frame}_R{R}A_m"] = d
    return tab, meta


def r_of(tab, key, lam2):
    lam2 = np.clip(lam2, LAM_GRID[0], LAM_GRID[-1])
    return np.interp(np.log(lam2), np.log(LAM_GRID), tab[key])


def lam_of(tab, key, r):
    rr = tab[key]
    if r <= rr[0]:
        return None
    if r >= rr[-1]:
        return float("inf")
    return float(np.exp(np.interp(r, rr, np.log(LAM_GRID))))


# ---- radius-dependent forms; R in nm; s = interface-shell thickness
def fcore(R, s):
    return np.clip((R - s) / R, 0, 1) ** 3


MODELS = {
    "const": {"n": 1, "desc": "lam2 = L (one contrast for all radii, the current P1 assumption, refitted)",
              "f": lambda R, p, s: np.exp(p[0]) + 0 * R, "p0": [np.log(1e-2)]},
    "a_core_fraction": {"n": 1, "desc": "lam2 = L_inf * f_core, f_core = ((R - s)/R)^3 (contrast proportional to crystalline-core volume fraction)",
                        "f": lambda R, p, s: np.exp(p[0]) * fcore(R, s), "p0": [np.log(1e-2)]},
    "a2_core_amplitude": {"n": 1, "desc": "lam2 = L_inf * f_core^2 (core-only scattering amplitude, shell inert)",
                          "f": lambda R, p, s: np.exp(p[0]) * fcore(R, s) ** 2, "p0": [np.log(1e-2)]},
    "b_core_shell": {"n": 2, "desc": "lam2 = (a_core f_core + a_shell (1 - f_core))^2 (volume-averaged Rayleigh amplitude; shell is a disordered layer with its own signed contrast)",
                     "f": lambda R, p, s: (p[0] * fcore(R, s) + p[1] * (1 - fcore(R, s))) ** 2, "p0": [0.1, 0.1]},
    "c_power_law": {"n": 2, "desc": "lam2 = A (R / 2 nm)^p",
                    "f": lambda R, p, s: np.exp(p[0]) * (R / 2.0) ** p[1], "p0": [np.log(1e-2), 0.0]},
}


def fit(tab, x, frame, model, data, s=SHELL_NM):
    """data: list of (R_A, r_obs). Least squares in the observed reduction r."""
    M = MODELS[model]

    def res(p):
        return [r_of(tab, (x, R, frame), max(M["f"](R / 10, p, s), 1e-12)) - r for R, r in data]
    best = None
    starts = [M["p0"]]
    if model == "b_core_shell":
        starts = [[a, b] for a in (-0.4, -0.1, 0.1, 0.4) for b in (-0.4, -0.1, 0.1, 0.4)]
    if model == "c_power_law":
        starts = [[np.log(l), q] for l in (1e-3, 1e-2, 1e-1) for q in (-4, 0, 2)]
    for p0 in starts:
        sol = least_squares(res, p0, method="lm" if len(data) > M["n"] else "trf")
        if best is None or sol.cost < best.cost:
            best = sol
    rms = float(np.sqrt(np.mean(np.square(best.fun))))
    return list(map(float, best.x)), rms


# ---- decision side: unscaled almaBTE modes, as in decision_robustness.py / analyze_md.claims
D_EFF, Y = 0.3e-6, 0.01
MODES = {x: S.SModes(x) for x in S.XS}
K0 = {TT: MODES[0.2].kappa2(TT, D_EFF) for TT in (300, 1100)}


def kC(x, TT, y, Rm, lam2):
    if y == 0:
        return MODES[x].kappa2(TT, D_EFF)
    return MODES[x].kappa2(TT, D_EFF, 0.0, np_R=Rm, np_phi=S.phi(x, y), lam2=lam2)


def claims(lamfun):
    row = {}
    for Rn in (2.5, 5.0):
        lam2 = float(lamfun(Rn))
        row[f"lam2_R{Rn:g}"] = lam2
        for TT in (300, 1100):
            row[f"red_pct_R{Rn:g}_T{TT}"] = 100 * (1 - kC(0.2, TT, Y, Rn * 1e-9, lam2) / K0[TT])
            xm = S.xmin_for([kC(x, TT, Y, Rn * 1e-9, lam2) for x in S.XS], K0[TT])
            row[f"ge_saving_pct_R{Rn:g}_T{TT}"] = None if xm is None else 100 * (1 - xm / 0.2)
    return row


def pick(lamfun):
    out = {}
    for TT in (300, 1100):
        best = None
        for x in S.XS:
            for y in S.YS:
                for Rm in S.RS:
                    v = kC(x, TT, y, Rm, float(lamfun(Rm * 1e9)))
                    if best is None or v < best[0]:
                        best = (v, x, y, Rm * 1e9)
        out[f"T{TT}"] = dict(zip(("kappa", "x", "y_C", "R_nm"), best))
    return out


def recipe_kappa(lamfun):
    """kappa_L of the P1 recipe and its nearest alternatives, to show pick margins."""
    return {f"x{x:g}_C3_R{Rn:g}_T{TT}": kC(x, TT, 0.03, Rn * 1e-9, float(lamfun(Rn)))
            for x in (0.3, 0.5) for Rn in (2.5, 5.0) for TT in (300, 1100)}


def lamfun_for(model, p, s=SHELL_NM):
    f = MODELS[model]["f"]
    # physical guard for the far extrapolation in the pick: 0 < lam2 <= geometric-cap proxy 1e4
    return lambda Rn: float(np.clip(f(Rn, p, s), 1e-12, 1e4))


def main():
    pts, khost = md_points()
    tab, meta = response_tables(pts, khost)
    out = {"status": "POST HOC: not pre-registered, no ledger entry; analysis of committed MD and bond-angle results",
           "inputs": ["../small_sic/md_results.json", "../small_sic/x03/results/*", "../../bond_angles/bad_results.json",
                      "../analyze_sigec.py (modes, contrast, phi, xmin_for)"],
           "shell_nm": SHELL_NM, "shell_source": "bond_angles sic_radial: angle s.d. < 0.1 deg to r = 1.45 nm, jumps to 4.8 deg at 1.5-1.6 nm (nominal R = 2.0 nm, x = 0.3)",
           "lam2_density_x0.2": LAM_DEN, "lam2_stiffness_x0.2": LAM_STF, "model_side": meta,
           "models": {k: v["desc"] for k, v in MODELS.items()},
           "current_P1": {"lam2": LAM_DEN, "claims": claims(lambda R: LAM_DEN), "pick": pick(lambda R: LAM_DEN),
                                  "recipe_kappa": recipe_kappa(lambda R: LAM_DEN)}}
    # bond-angle profile used for the shell
    prof = BAD["sic_radial"]["profile"]
    out["bond_angle_profile"] = [{"r_mid_nm": (q["r_lo_A"] + q["r_hi_A"]) / 20, "std_deg": q["std_deg"]} for q in prof]

    # 1. per-seed effective lam2
    inv = {}
    for x in (0.2, 0.3):
        for frame in ("registered", "host_matched"):
            for R, p in pts[x].items():
                inv[f"x{x}_{frame}_R{R}A"] = {"r_seeds": p["r_seeds"],
                                              "lam2_seeds": [lam_of(tab, (x, R, frame), r) for r in p["r_seeds"]],
                                              "lam2_mean_r": lam_of(tab, (x, R, frame), float(np.mean(p["r_seeds"]))),
                                              "r_at_density_edge": float(r_of(tab, (x, R, frame), LAM_DEN))}
    out["effective_lam2"] = inv

    # 2.-3. fits and propagation
    fits = {}
    for x in (0.2, 0.3):
        Rs = sorted(pts[x])
        for frame in ("registered", "host_matched"):
            data = [(R, r) for R in Rs for r in pts[x][R]["r_seeds"]]
            for model in MODELS:
                if MODELS[model]["n"] > len(Rs) - (0 if x == 0.2 else 0):
                    continue
                key = f"x{x}_{frame}_{model}"
                p, rms = fit(tab, x, frame, model, data)
                lf = lamfun_for(model, p)
                row = {"params": p, "rms_resid_r": rms, "n_radii": len(Rs), "n_params": MODELS[model]["n"],
                       "lam2_at": {f"{R/10:g}nm": lf(R / 10) for R in Rs} | {"2.5nm": lf(2.5), "5nm": lf(5.0)},
                       "claims": claims(lf), "pick": pick(lf), "recipe_kappa": recipe_kappa(lf)}
                # seed recombinations (2^n_radii)
                combo = []
                for c in itertools.product(range(2), repeat=len(Rs)):
                    dat = [(R, pts[x][R]["r_seeds"][i]) for R, i in zip(Rs, c)]
                    pc, _ = fit(tab, x, frame, model, dat)
                    combo.append(claims(lamfun_for(model, pc)) | {"params": pc})
                row["seed_recombinations"] = {
                    k: [min(cc[k] for cc in combo if cc[k] is not None), max(cc[k] for cc in combo if cc[k] is not None)]
                    for k in combo[0] if k != "params"}
                row["seed_recombination_picks"] = sorted({json.dumps(pick(lamfun_for(model, cc["params"]))["T300"]
                                                                     | {"kappa": None}) for cc in combo})
                if model in ("a_core_fraction", "a2_core_amplitude", "b_core_shell"):
                    row["shell_sensitivity"] = {}
                    for s in SHELL_ALT:
                        ps, rs = fit(tab, x, frame, model, data, s)
                        row["shell_sensitivity"][f"s{s:g}nm"] = {"params": ps, "rms_resid_r": rs,
                                                                "claims": claims(lamfun_for(model, ps, s))}
                fits[key] = row
                print(key, "rms", round(rms, 4), "lam2(2.5,5)", f"{lf(2.5):.3g} {lf(5.0):.3g}",
                      "Ge300", [round(row["claims"][f"ge_saving_pct_R{r}_T300"] or -1, 1) for r in ("2.5", "5")],
                      "Ge1100", [round(row["claims"][f"ge_saving_pct_R{r}_T1100"] or -1, 1) for r in ("2.5", "5")],
                      "pick", [(v["x"], v["y_C"], v["R_nm"]) for v in row["pick"].values()], flush=True)
    out["fits"] = fits
    json.dump(out, open(os.path.join(HERE, "radius_contrast_results.json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
