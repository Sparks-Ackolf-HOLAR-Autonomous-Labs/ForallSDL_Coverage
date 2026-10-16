#!/usr/bin/env python3
"""Stage A of the pre-registered SiGeC campaign (sigec_design.json, ledger seq 46).

Reads only committed evidence: ../nano/modes/x*.npz (almaBTE per-mode data,
ledger seq 45) and m1/m1_summary.json (MACE gate, ledgered before this runs).
Refuses to run if sigec_design.json changed after registration."""
import hashlib, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(MT, "nano"))
import analyze_nano as AN  # noqa: E402  (checks the nano design hash on import)

DESIGN = os.path.join(HERE, "sigec_design.json")
reg = open(os.path.join(HERE, "PREREGISTERED_SHA256")).read().split()[0]
if hashlib.sha256(open(DESIGN, "rb").read()).hexdigest() != reg:
    sys.exit("sigec_design.json changed after registration")
D = json.load(open(DESIGN))
G = D["stage_A_analysis"]["grid"]
XS, TS, YS, RS = G["x_Ge"], G["T_K"], G["y_C"], [r * 1e-9 for r in G["SiC_radius_nm"]]
DEFF = [0.1e-6, 0.3e-6, 0.5e-6]
EPS = (0.0, 39.0)
# M1c decision rule (m1c/m1c_design.json): re-run on the converged summary via env overrides
M1 = json.load(open(os.environ.get("SIGEC_M1_SUMMARY", os.path.join(HERE, "m1", "m1_summary.json"))))
OUT_NAME = os.environ.get("SIGEC_OUT_NAME", "sigec_results.json")
OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)

# ---- constants from the design
M_C = 12.011
A_SIC = 4.360
RHO = {"Si": 2.33, "Ge": 5.32, "SiC": 3.21}
KB = {"Si": 98.0, "Ge": 75.0, "SiC": 225.0}
XM1 = [0.0, 0.2, 0.5]
DC = [M1["per_x"][f"{x:g}"]["delta_C"] for x in XM1] if M1["gates_pass"] else [M1["modified_vegard_delta_C"]] * 3
DPREC = [M1["per_x"][f"{x:g}"]["dE_prec_mean_eV"] for x in XM1]


def delta_C(x):
    return float(np.interp(min(x, 0.5), XM1, DC))


def dE_prec(x):
    return float(np.interp(min(x, 0.5), XM1, DPREC))


def contrast(x, kind):
    if kind == "density":
        h = (1 - x) * RHO["Si"] + x * RHO["Ge"]
        return ((RHO["SiC"] - h) / h) ** 2
    h = (1 - x) * KB["Si"] + x * KB["Ge"]
    return ((KB["SiC"] - h) / h) ** 2


def phi(x, y):
    a = (1 - x) * 5.431 + x * 5.658
    return y * (A_SIC ** 3 / 4) / (a ** 3 / 8)


class SModes(AN.Modes):
    def kappa2(self, T, d=np.inf, gamma=0.0, np_R=None, np_phi=0.0, lam2=0.0):
        tau = self.tau[T]
        ok = tau > 0
        rate = np.zeros_like(tau)
        rate[ok] = 1.0 / tau[ok]
        if np.isfinite(d):
            rate = rate + self.v / d
        if gamma:
            rate = rate + (np.pi / 6) * gamma * self.w ** 2 * self.Dt
        if np_R and np_phi > 0:
            R = np_R
            q = self.w / np.maximum(self.v, 1000.0)
            sg = 2 * np.pi * R ** 2
            sr = (4.0 / 9.0) * np.pi * R ** 2 * lam2 * (q * R) ** 4
            sig = np.where(sr > 0, 1.0 / (1.0 / np.maximum(sr, 1e-300) + 1.0 / sg), 0.0)
            n = np_phi / (4.0 / 3.0 * np.pi * R ** 3)
            rate = rate + self.v * n * sig
        m = ok & (rate > 0)
        return float((self.C[T][m] * self.v2[m] / rate[m]).sum() / 3.0)


def gamma_sub(x, y, eps):
    mbar, _, _ = AN.host(x)
    return y * ((1 - M_C / mbar) ** 2 + eps * delta_C(x) ** 2)


FORMS = [("sub", e, None, None) for e in EPS] + [("SiC", None, R, c) for R in RS for c in ("density", "stiffness")]


def form_label(f):
    return f"sub_eps{f[1]:g}" if f[0] == "sub" else f"SiC_R{f[2]*1e9:g}nm_{f[3]}"


def kap(m, x, T, d, y, f, gdop=0.0):
    if y == 0 or f is None:
        return m.kappa2(T, d, gdop)
    if f[0] == "sub":
        return m.kappa2(T, d, gamma_sub(x, y, f[1]) + gdop)
    return m.kappa2(T, d, gdop, np_R=f[2], np_phi=phi(x, y), lam2=contrast(x, f[3]))


def lnT_interp(k900, k1100, T=1023.0):
    return float(np.exp(np.interp(T, [900, 1100], [np.log(k900), np.log(k1100)])))


def xmin_for(ks, tgt):
    if ks[0] <= tgt:
        return 0.0
    for i in range(1, len(XS)):
        if ks[i] <= tgt < ks[i - 1]:
            f = (np.log(ks[i - 1]) - np.log(tgt)) / (np.log(ks[i - 1]) - np.log(ks[i]))
            return XS[i - 1] + f * (XS[i] - XS[i - 1])
    return None


def main():
    modes = {x: SModes(x) for x in XS}
    res = {"design_sha256": reg, "M1_gates_pass": M1["gates_pass"], "delta_C_used": dict(zip(map(str, XM1), DC)),
           "dE_prec_used_eV": dict(zip(map(str, XM1), DPREC)), "checks": {}, "pred": {}}
    # implementation checks (not pre-registered): y=0 reproduces the nano campaign
    chk = []
    for x in (0.0, 0.2, 0.5):
        for T in (300, 1100):
            chk.append(abs(modes[x].kappa2(T, 20e-9) / modes[x].kappa(T, 20e-9) - 1))
    res["checks"]["y0_reproduces_nano"] = max(chk)

    # ---- D1 table (compact: d in {20 nm, 100 nm, d_eff band, bulk})
    DS1 = [20e-9, 100e-9] + DEFF + [np.inf]
    d1 = []
    for x in XS:
        m = modes[x]
        for T in TS:
            for d in DS1:
                k0 = kap(m, x, T, d, 0, None)
                for y in YS[1:]:
                    for f in FORMS:
                        k = kap(m, x, T, d, y, f)
                        d1.append({"x": x, "T": T, "d_m": None if not np.isfinite(d) else d, "y_C": y,
                                   "form": form_label(f), "kappa": k, "kappa_noC": k0})
    res["D1"] = d1
    idx = {(r["x"], r["T"], r["d_m"], r["y_C"], r["form"]): r for r in d1}

    def K(x, T, d, y, f):
        if y == 0:
            return kap(modes[x], x, T, d, 0, None)
        return idx[(x, T, d, y, form_label(f))]["kappa"]

    # ---- D3 per-carbon lever at y=0.01
    d3 = []
    for x in XS:
        for T in (300, 1100):
            for d in DEFF:
                k0 = K(x, T, d, 0, None)
                row = {"x": x, "T": T, "d_m": d, "kappa_noC": k0}
                for f in FORMS:
                    row[form_label(f)] = 100 * (K(x, T, d, 0.01, f) / k0 - 1)
                sic = {k: v for k, v in row.items() if k.startswith("SiC")}
                row["best_SiC"] = min(sic, key=sic.get)
                d3.append(row)
    res["D3_lever_pct_at_yC_0.01"] = d3

    # ---- D2 carbon-for-germanium frontier
    d2 = []
    for T in (300, 1100):
        for d in DEFF:
            for tgt in (2.5, 4.0):
                for y in YS:
                    for f in ([None] if y == 0 else FORMS):
                        ks = [K(x, T, d, y, f) for x in XS]
                        d2.append({"T": T, "d_m": d, "target": tgt, "y_C": y,
                                   "form": "none" if f is None else form_label(f), "x_min": xmin_for(ks, tgt)})
    res["D2_frontier"] = d2

    # ---- D4 screening at d_eff = 0.3 um
    d4 = {}
    for T in (300, 1100):
        cand = []
        for x in XS:
            for y in YS:
                for f in ([None] if y == 0 else FORMS):
                    cand.append((K(x, T, 0.3e-6, y, f), x, y, "none" if f is None else form_label(f), f))
        best = min(cand, key=lambda c: c[0])
        eq = [c for c in cand if c[4] is None or c[4][0] == "SiC" or dE_prec(c[1]) >= 0]
        beq = min(eq, key=lambda c: c[0])
        nox = min((c for c in cand if c[2] == 0), key=lambda c: c[0])
        d4[f"T{T}"] = {"unconstrained": best[:4], "equilibrium_constrained": beq[:4], "carbon_free": nox[:4]}
    res["D4_screening"] = d4

    # ---- D5 where carbon goes
    res["D5_M1"] = {x: {k: M1["per_x"][x][k] for k in ("delta_C", "R2", "dE_prec_mean_eV", "dE_site_mean_eV")}
                    for x in M1["per_x"]}

    # ---- predictions
    P = res["pred"]
    dm = M1["modified_vegard_delta_C"]
    d0 = M1["per_x"]["0"]["delta_C"]
    P["PrC1"] = {"delta_C_x0": d0, "modified_vegard": dm, "rel": d0 / dm - 1, "pass": abs(d0 / dm - 1) <= 0.3}
    dp = {x: M1["per_x"][x]["dE_prec_mean_eV"] for x in M1["per_x"]}
    P["PrC2"] = {"dE_prec": dp, "pass": all(v < -0.5 for v in dp.values()) and dp["0.5"] < dp["0"]}
    ds = {x: M1["per_x"][x]["dE_site_mean_eV"] for x in ("0.2", "0.5")}
    P["PrC2b"] = {"dE_site": ds, "pass": all(v is not None and v > 0 for v in ds.values())}
    k_np = K(0.2, 300, 0.3e-6, 0.01, ("SiC", None, 2.5e-9, "stiffness"))
    k_s0 = K(0.2, 300, 0.3e-6, 0.01, ("sub", 0.0, None, None))
    k_s39 = K(0.2, 300, 0.3e-6, 0.01, ("sub", 39.0, None, None))
    P["PrC3"] = {"k_SiC_R2.5_stiff": k_np, "k_sub_eps0": k_s0, "k_sub_eps39": k_s39,
                 "kappa_noC": K(0.2, 300, 0.3e-6, 0, None), "pass": k_s39 < k_np < k_s0}
    rows = []
    for R in (50e-9, 100e-9, 150e-9):
        for c in ("density", "stiffness"):
            for d in DEFF:
                m = modes[0.2]
                f = ("SiC", None, R, c)
                k0 = lnT_interp(kap(m, 0.2, 900, d, 0, None), kap(m, 0.2, 1100, d, 0, None))
                k1 = lnT_interp(kap(m, 0.2, 900, d, 0.03, f), kap(m, 0.2, 1100, d, 0.03, f))
                rows.append({"R_nm": R * 1e9, "contrast": c, "d_eff_m": d, "pct": 100 * (k1 / k0 - 1)})
    worst = min(r["pct"] for r in rows)
    P["PrC4"] = {"rows": rows, "largest_reduction_pct": worst, "measured_total_kappa_drop_pct": -20.9,
                 "pass": worst > -10}
    tgt = K(0.2, 300, 0.3e-6, 0, None)
    best_x = None
    for f in FORMS:
        if f[0] == "sub" and f[1] != 39.0:
            continue
        if f[0] == "SiC" and f[3] != "stiffness":
            continue
        xm = xmin_for([K(x, 300, 0.3e-6, 0.01, f) for x in XS], tgt)
        if xm is not None and (best_x is None or xm < best_x[0]):
            best_x = (xm, form_label(f))
    P["PrC5"] = {"target_kappa": tgt, "x_min": best_x[0], "form": best_x[1], "pass": best_x[0] <= 0.15}
    fC = ("SiC", None, 2.5e-9, "stiffness")
    ch = []
    for x in (0.15, 0.2, 0.3, 0.5):
        m = modes[x]
        for T in TS:
            k0 = kap(m, x, T, 0.3e-6, 0.01, fC)
            for dop in ("B", "P"):
                for n in (1e19, 3e19, 1e20):
                    k = kap(m, x, T, 0.3e-6, 0.01, fC, AN.gamma_dop(x, dop, n, 0.0))
                    ch.append(100 * (k / k0 - 1))
    P["PrC6"] = {"largest_abs_change_pct": max(abs(c) for c in ch), "pass": max(abs(c) for c in ch) < 5}
    json.dump(res, open(os.path.join(OUT, OUT_NAME), "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    print(json.dumps({"checks": res["checks"], "pred": {k: {kk: vv for kk, vv in v.items() if kk != "rows"} for k, v in P.items()},
                      "D4": d4}, indent=1, default=str))


if __name__ == "__main__":
    main()
