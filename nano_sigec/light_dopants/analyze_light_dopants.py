#!/usr/bin/env python3
"""Analysis of the pre-registered light-dopant extension (light_dopants_design.json,
ledger seq 74). Reads only committed evidence: the almaBTE host modes (ledger 45)
via analyze_sigec, and the converged MACE summary (ledger 54-55). Refuses to run if
the design changed after registration. Runtime: seconds."""
import hashlib, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
reg = open(os.path.join(HERE, "PREREGISTERED_SHA256")).read().split()[0]
DES = os.path.join(HERE, "light_dopants_design.json")
if hashlib.sha256(open(DES, "rb").read()).hexdigest() != reg:
    sys.exit("light_dopants_design.json changed after registration")
D = json.load(open(DES))
os.environ.setdefault("SIGEC_M1_SUMMARY", os.path.join(HERE, "..", "m1c", "m1c_summary.json"))
sys.path.insert(0, os.path.dirname(HERE))
import analyze_sigec as S  # noqa: E402
AN = S.AN

C = D["constants"]
MX = C["mass_amu"]
EPS_U, DEL_U = C["eps_upper"], C["delta_upper"]
LAM_V = C["void_lam2"]
XS, TS, DEFF = S.XS, sorted(S.TS), D["grid"]["d_eff_m"]
D0, X0 = 0.3e-6, 0.2
RV = [2.5e-9, 5e-9, 10e-9, 25e-9]
PHIV = [0.001, 0.003, 0.01, 0.03]
NI = [1e18, 1e19, 1e20]
modes = {x: S.SModes(x) for x in XS}


def incl_rate(m, R, phi, lam2):
    q = m.w / np.maximum(m.v, 1000.0)
    sg = 2 * np.pi * R ** 2
    sr = (4.0 / 9.0) * np.pi * R ** 2 * lam2 * (q * R) ** 4
    sig = np.where(sr > 0, 1.0 / (1.0 / np.maximum(sr, 1e-300) + 1.0 / sg), 0.0)
    return m.v * (phi / (4.0 / 3.0 * np.pi * R ** 3)) * sig


def kap(x, T, d=D0, gamma=0.0, incl=(), porosity=0.0):
    """kappa_L with boundary d, point-defect Gamma and any list of (R, phi, lam2) inclusions."""
    m = modes[x]
    tau = m.tau[T]
    ok = tau > 0
    rate = np.zeros_like(tau)
    rate[ok] = 1.0 / tau[ok]
    rate = rate + m.v / d
    if gamma:
        rate = rate + (np.pi / 6) * gamma * m.w ** 2 * m.Dt
    for R, phi, lam2 in incl:
        rate = rate + incl_rate(m, R, phi, lam2)
    k = (m.C[T][ok & (rate > 0)] * m.v2[ok & (rate > 0)] / rate[ok & (rate > 0)]).sum() / 3.0
    return float(k * (1 - porosity) / (1 + porosity / 2))


def gamma_int(x, X, n, bound):
    mbar, _, nat = AN.host(x)
    c = n / nat
    return c * ((MX[X] / mbar) ** 2 + (EPS_U * DEL_U ** 2 if bound == "upper" else 0.0))


def sic(x, y, R=5e-9):
    return (R, S.phi(x, y), S.contrast(x, "density"))


def void(R, phi, lam2):
    return (R, phi, lam2)


def pct(a, b):
    return 100 * (a / b - 1)


def xmatch(T, **kw):
    tgt = kap(X0, T)
    return S.xmin_for([kap(x, T, **kw) for x in XS], tgt)


def main():
    res = {"design_sha256": reg, "checks": {}, "pred": {}}
    res["checks"]["reproduces_sigec_noC"] = max(abs(kap(x, T) / modes[x].kappa2(T, D0) - 1) for x in XS for T in TS)
    res["checks"]["reproduces_sigec_SiC5nm_1pct"] = max(
        abs(kap(x, T, incl=[sic(x, 0.01)]) / S.kap(modes[x], x, T, D0, 0.01, ("SiC", None, 5e-9, "density")) - 1)
        for x in XS for T in TS)

    # interstitial H, D, He (isolated)
    inter = []
    for X in ("H", "D", "He"):
        for x in XS:
            for T in TS:
                for d in DEFF:
                    for n in NI:
                        for b in ("mass", "upper"):
                            inter.append({"X": X, "x": x, "T": T, "d": d, "n": n, "bound": b,
                                          "pct": pct(kap(x, T, d, gamma_int(x, X, n, b)), kap(x, T, d))})
    res["interstitial"] = inter
    # substitutional B, P on top of the SiC recipe
    subs = []
    for x in (0.2, 0.3, 0.5):
        for T in TS:
            k0 = kap(x, T, incl=[sic(x, 0.01)])
            for X in ("B", "P"):
                for n in (1e19, 3e19, 1e20):
                    subs.append({"X": X, "x": x, "T": T, "n": n,
                                 "pct": pct(kap(x, T, gamma=AN.gamma_dop(x, X, n, 0.0), incl=[sic(x, 0.01)]), k0)})
    res["substitutional_on_SiC"] = subs
    # cavities
    cav = []
    for x in XS:
        for T in TS:
            for d in DEFF:
                k0 = kap(x, T, d)
                for R in RV:
                    for phi in PHIV:
                        for l in LAM_V:
                            cav.append({"x": x, "T": T, "d": d, "R_nm": R * 1e9, "phi": phi, "lam2": l,
                                        "pct": pct(kap(x, T, d, incl=[void(R, phi, l)], porosity=phi), k0),
                                        "pct_scatter_only": pct(kap(x, T, d, incl=[void(R, phi, l)]), k0)})
    res["cavities"] = cav

    P = res["pred"]
    sel = lambda X, b: [r for r in inter if r["X"] == X and r["bound"] == b and r["x"] == X0 and r["d"] == D0]
    m1 = max(abs(r["pct"]) for X in ("H", "D") for r in sel(X, "mass"))
    P["PrL1"] = {"max_abs_pct": m1, "pass": m1 < 1}
    m2 = max(abs(r["pct"]) for X in ("H", "D", "He") for r in sel(X, "upper"))
    P["PrL2"] = {"max_abs_pct": m2, "pass": m2 < 5}
    hd = {}
    for r in inter:
        if r["X"] in ("H", "D"):
            hd.setdefault((r["x"], r["T"], r["d"], r["n"], r["bound"]), {})[r["X"]] = r["pct"]
    m3 = max(abs(v["H"] - v["D"]) for v in hd.values())
    P["PrL3"] = {"max_abs_H_minus_D_pct_points": m3, "pass": m3 < 0.5}
    rows4 = {}
    for T in (300, 1100):
        k0 = kap(X0, T)
        rows4[T] = {"SiC_density": pct(kap(X0, T, incl=[sic(X0, 0.01)]), k0),
                    **{f"void_lam{l:g}": pct(kap(X0, T, incl=[void(5e-9, S.phi(X0, 0.01), l)],
                                                 porosity=S.phi(X0, 0.01)), k0) for l in LAM_V}}
    P["PrL4"] = {"rows": rows4, "pass": all(min(r["void_lam1"], r["void_lam4"]) < r["SiC_density"] for r in rows4.values())}
    sav = {}
    for T in (300, 1100):
        for l in LAM_V:
            xm = xmatch(T, incl=[void(5e-9, 0.01, l)], porosity=0.01)
            sav[f"T{T}_lam{l:g}"] = None if xm is None else 100 * (1 - xm / X0)
        xm = xmatch(T, incl=[sic(0.2, 0.01)])  # reference: SiC at 1 at% C, but phi follows x
        sav[f"T{T}_SiC_ref"] = None if xm is None else 100 * (1 - xm / X0)
    P["PrL5"] = {"ge_saving_pct": sav, "pass": all(v is not None and v >= 25 for k, v in sav.items() if k.startswith("T300_lam"))}
    rows6 = {}
    for T in (300, 1100):
        k0 = kap(X0, T)
        a = kap(X0, T, incl=[sic(X0, 0.01)])
        b = kap(X0, T, incl=[void(5e-9, 0.01, 1.0)], porosity=0.01)
        ab = kap(X0, T, incl=[sic(X0, 0.01), void(5e-9, 0.01, 1.0)], porosity=0.01)
        rows6[T] = {"SiC": pct(a, k0), "void": pct(b, k0), "both": pct(ab, k0), "sum": pct(a, k0) + pct(b, k0)}
    P["PrL6"] = {"rows": rows6, "pass": all(abs(r["both"]) < abs(r["sum"]) for r in rows6.values())}
    m7 = max(abs(r["pct"]) for r in subs)
    P["PrL7"] = {"max_abs_pct": m7, "pass": m7 < 1}

    # decision rule: lever / inert at the reference point
    lev = {}
    for X in ("H", "D", "He"):
        v = max(abs(r["pct"]) for r in inter if r["X"] == X and r["x"] == X0 and r["d"] == D0 and r["n"] == 1e20)
        lev[f"{X} isolated"] = v
    lev["B, P substitutional"] = max(abs(r["pct"]) for r in subs if r["x"] == X0)
    lev["cavities 1 vol% R 5 nm"] = max(abs(r["pct"]) for r in cav if r["x"] == X0 and r["d"] == D0
                                        and r["phi"] == 0.01 and r["R_nm"] == 5)
    res["decision"] = {k: {"max_abs_pct": v, "class": "lever" if v >= 5 else ("inert" if v < 1 else "minor")}
                       for k, v in lev.items()}
    # screening: does adding cavities move x or the carbon part of the pick?
    picks = {}
    for T in (300, 1100):
        for label, extra in (("no_cavities", None), ("cavities_1pct_5nm", 0.01), ("cavities_3pct_5nm", 0.03)):
            best = None
            for x in XS:
                for y in S.YS:
                    for R in S.RS:
                        inc = [sic(x, y, R)] if y else []
                        if extra:
                            inc.append(void(5e-9, extra, 1.0))
                        v = kap(x, T, incl=inc, porosity=extra or 0.0)
                        if best is None or v < best[0]:
                            best = (v, x, y, R * 1e9 if y else None)
            picks[f"T{T}_{label}"] = {"kappa": best[0], "x": best[1], "y_C": best[2], "SiC_R_nm": best[3]}
    res["screening"] = picks
    ch = any((picks[f"T{T}_{l}"]["x"], picks[f"T{T}_{l}"]["y_C"], picks[f"T{T}_{l}"]["SiC_R_nm"]) !=
             (picks[f"T{T}_no_cavities"]["x"], picks[f"T{T}_no_cavities"]["y_C"], picks[f"T{T}_no_cavities"]["SiC_R_nm"])
             for T in (300, 1100) for l in ("cavities_1pct_5nm", "cavities_3pct_5nm"))
    res["recipe_changed"] = ch
    json.dump(res, open(os.path.join(HERE, "light_dopants_results.json"), "w"), indent=1, default=float)
    print(json.dumps({"checks": res["checks"], "pred": P, "decision": res["decision"], "screening": picks,
                      "recipe_changed": ch}, indent=1, default=float))


if __name__ == "__main__":
    main()
