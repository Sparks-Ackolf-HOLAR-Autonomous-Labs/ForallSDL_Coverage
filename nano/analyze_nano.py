#!/usr/bin/env python3
"""Analyse the pre-registered nano campaign (nano_design.json, ledger seq 40).

Reads only committed evidence: modes/x*.npz, modes/x*_tsweep.json and the
open-data party. Writes results/nano_results.json, results/RESULTS_TABLES.md
and figures. Refuses to run if nano_design.json no longer matches the hash
recorded in PREREGISTERED_SHA256.
"""
import hashlib, json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MT = os.path.dirname(HERE)
DESIGN = os.path.join(HERE, "nano_design.json")
reg = open(os.path.join(HERE, "PREREGISTERED_SHA256")).read().split()[0]
if hashlib.sha256(open(DESIGN, "rb").read()).hexdigest() != reg:
    sys.exit("nano_design.json changed after registration")
D = json.load(open(DESIGN))
K = D["physics_constants"]
M, R = K["mass_amu"], K["tetrahedral_covalent_radius_A"]
XS, TS = D["grid"]["x_Ge"], D["grid"]["T_K"]
DS = np.logspace(np.log10(5e-9), np.log10(1e-5), 25)
EPS = [K["strain_parameter_eps"]["lower_bound"], K["strain_parameter_eps"]["upper_bound"]]
OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)


def load(x):
    z = np.load(os.path.join(HERE, "modes", f"x{x:g}.npz"))
    ts = json.load(open(os.path.join(HERE, "modes", f"x{x:g}_tsweep.json")))
    return z, ts


def host(x):
    mbar = (1 - x) * M["Si"] + x * M["Ge"]
    rbar = (1 - x) * R["Si"] + x * R["Ge"]
    a_cm = ((1 - x) * K["lattice_A"]["Si"] + x * K["lattice_A"]["Ge"]) * 1e-8
    return mbar, rbar, 8.0 / a_cm ** 3


def dos_tilde(omega):
    """Per-atom DOS normalised to integral 3, Gaussian-broadened, evaluated at each mode."""
    w = omega[omega > 0].astype(np.float64)
    sig = 0.005 * w.max()  # dos_broadening in the design
    grid = np.linspace(0, w.max() * 1.05, 1200)
    hist = np.zeros_like(grid)
    for c in np.array_split(w, 40):
        hist += np.exp(-0.5 * ((grid[:, None] - c[None, :]) / sig) ** 2).sum(1)
    hist *= 3.0 / getattr(np, "trapezoid", getattr(np, "trapz", None))(hist, grid)
    return np.interp(omega, grid, hist)


def gamma_dop(x, dop, amount, eps):
    if dop == "none" or amount == 0:
        return 0.0
    mbar, rbar, nat = host(x)
    c = amount if dop == "C" else amount / nat
    return c * ((1 - M[dop] / mbar) ** 2 + eps * (1 - R[dop] / rbar) ** 2)


class Modes:
    def __init__(self, x):
        z, self.ts = load(x)
        self.x = x
        self.w = z["omega"].astype(np.float64)
        self.v2 = z["v2"].astype(np.float64)
        self.v = np.sqrt(self.v2)
        self.C = {T: z[f"C_{T}"].astype(np.float64) for T in TS}
        self.tau = {T: z[f"tau_{T}"].astype(np.float64) for T in TS}
        self.Dt = dos_tilde(self.w)

    def kappa(self, T, d=np.inf, gamma=0.0):
        tau = self.tau[T]
        ok = tau > 0
        rate = np.zeros_like(tau)
        rate[ok] = 1.0 / tau[ok]
        if np.isfinite(d):
            rate = rate + self.v / d
        if gamma:
            rate = rate + (np.pi / 6) * gamma * self.w ** 2 * self.Dt
        m = ok & (rate > 0)
        return float((self.C[T][m] * self.v2[m] / rate[m]).sum() / 3.0)


def main():
    modes = {x: Modes(x) for x in XS}
    res = {"design_sha256": reg, "grain_sizes_m": DS.tolist(), "gates": {}, "pred": {}}

    # ---- N0 gates
    g1 = []
    for x, m in modes.items():
        for T in TS:
            if str(T) in m.ts:
                k, ref = m.kappa(T), m.ts[str(T)]["rta_xx"]
                g1.append({"x": x, "T": T, "per_mode": k, "tsweep_rta": ref, "rel": k / ref - 1})
    res["gates"]["G1"] = {"rows": g1, "max_abs_rel": max(abs(r["rel"]) for r in g1),
                          "pass": all(abs(r["rel"]) <= 0.005 for r in g1)}
    mbar = 0.8 * M["Si"] + 0.2 * M["Ge"]
    g_all = 0.8 * (1 - M["Si"] / mbar) ** 2 + 0.2 * (1 - M["Ge"] / mbar) ** 2
    k_t = modes[0.0].kappa(300, gamma=g_all)
    ref = modes[0.2].ts["300"]["rta_xx"]
    res["gates"]["G2"] = {"g_alloy_x0.2": g_all, "kappa_tamura_on_Si": k_t, "kappa_vca_x0.2": ref,
                          "ratio": k_t / ref, "pass": 1 / 1.5 <= k_t / ref <= 1.5}
    g3 = True
    for x, m in modes.items():
        for T in TS:
            ks = [m.kappa(T, d) for d in DS]
            kb = m.kappa(T)
            g3 &= all(np.diff(ks) > 0)
            if 0 < x < 1:
                g3 &= ks[-1] >= 0.7 * kb
    ksi = modes[0.0].kappa(300, DS[0]) / modes[0.0].kappa(300)
    res["gates"]["G3"] = {"pure_Si_5nm_over_bulk": ksi, "pass": bool(g3 and ksi < 0.25)}
    res["gates_pass"] = all(g["pass"] for g in res["gates"].values())

    # ---- D1 grid (undoped)
    grid = {}
    for x, m in modes.items():
        for T in TS:
            grid[f"{x:g}_{T}"] = {"bulk": m.kappa(T), "d": [m.kappa(T, d) for d in DS]}
    res["D1_undoped"] = grid

    def kd(x, T, d):  # kappa at an arbitrary grain size by log-log interpolation
        if not np.isfinite(d):
            return grid[f"{x:g}_{T}"]["bulk"]
        return float(np.exp(np.interp(np.log(d), np.log(DS), np.log(grid[f"{x:g}_{T}"]["d"]))))

    # ---- D2 Ge-savings frontier (ln kappa linear in x between grid points)
    front = {}
    for T in (300, 1100):
        for tgt in (2.5, 4.0):
            row = []
            for d in list(DS) + [np.inf]:
                ks = [kd(x, T, d) for x in XS]
                xmin = None
                for i in range(1, len(XS)):
                    if ks[i] <= tgt < ks[i - 1] or ks[0] <= tgt:
                        if ks[0] <= tgt:
                            xmin = 0.0
                        else:
                            f = (np.log(ks[i - 1]) - np.log(tgt)) / (np.log(ks[i - 1]) - np.log(ks[i]))
                            xmin = XS[i - 1] + f * (XS[i] - XS[i - 1])
                        break
                row.append({"d_m": None if not np.isfinite(d) else float(d), "x_min": xmin})
            front[f"T{T}_k{tgt}"] = row
    res["D2_frontier"] = front

    # ---- D3 dopant lever
    lever = []
    for dop, spec in D["grid"]["dopants"].items():
        if dop == "none":
            continue
        amounts = spec["at_frac"] if dop == "C" else spec["n_cm3"]
        for x in XS:
            m = modes[x]
            for T in (300, 1100):
                for d in (np.inf, 100e-9, 20e-9):
                    k0 = m.kappa(T, d)
                    for a in amounts:
                        for eps in EPS:
                            g = gamma_dop(x, dop, a, eps)
                            k = m.kappa(T, d, g)
                            lever.append({"dop": dop, "amount": a, "x": x, "T": T,
                                          "d_m": None if not np.isfinite(d) else d, "eps": eps,
                                          "Gamma": g, "kappa0": k0, "kappa": k,
                                          "pct_change": 100 * (k / k0 - 1)})
    res["D3_lever"] = lever

    # ---- D4 screening
    scr = {}
    for T in (300, 1100):
        best = min(((kd(x, T, 10e-9), x) for x in XS))
        bulkbest = min(((grid[f"{x:g}_{T}"]["bulk"], x) for x in XS))
        # grain size delivering half of the bulk -> 10 nm reduction at the best x
        xb = best[1]
        kb, k10 = grid[f"{xb:g}_{T}"]["bulk"], kd(xb, T, 10e-9)
        half = np.exp(0.5 * (np.log(kb) + np.log(k10)))
        dd = np.logspace(np.log10(10e-9), -5, 400)
        d_half = float(dd[np.argmin([abs(kd(xb, T, q) - half) for q in dd])])
        scr[f"T{T}"] = {"best_x_at_10nm": xb, "kappa_at_10nm": best[0],
                        "best_x_bulk": bulkbest[1], "kappa_bulk_best": bulkbest[0],
                        "d_half_reduction_m": d_half}
    res["D4_screening"] = scr

    # ---- predictions
    k20 = kd(0.2, 300, 20e-9)
    res["pred"]["Pr1"] = {"kappa": k20, "pass": 2 <= k20 <= 5}
    s02 = k20 / grid["0.2_300"]["bulk"]
    s05 = kd(0.5, 300, 20e-9) / grid["0.5_300"]["bulk"]
    res["pred"]["Pr2"] = {"supp_x0.2": s02, "supp_x0.5": s05, "pass": s02 < s05}
    worst = min(r["pct_change"] for r in lever if r["dop"] in ("B", "P") and 0.15 <= r["x"] <= 0.85
                and r["amount"] <= 1e20 and r["eps"] == 0.0)
    res["pred"]["Pr3"] = {"largest_reduction_pct": worst, "pass": worst > -10}
    m = modes[0.2]
    g = gamma_dop(0.2, "B", 1e20, EPS[1])
    kb_, kgb, kdp, kboth = m.kappa(300), m.kappa(300, 20e-9), m.kappa(300, gamma=g), m.kappa(300, 20e-9, g)
    S = kboth * kb_ / (kdp * kgb)
    res["pred"]["Pr4"] = {"S": S, "eps": EPS[1], "pass": S >= 1}
    pts = json.load(open(os.path.join(MT, "party5_databases", "open_db", "open_db_points.json")))
    sel = [p for p in pts if p["system"] == "nanostructured" and 0.15 <= p["x_ge"] <= 0.25
           and abs(p["T_K"] - 300) <= 25]
    inf = []
    for p in sel:
        x = p["x_ge"]
        # ln kappa linear in x between 0.15/0.2/0.3 grid rows, then solve for d
        dd = np.logspace(np.log10(5e-9), -5, 400)
        xs = [q for q in XS if 0.1 <= q <= 0.3]

        def kx(d):
            return float(np.exp(np.interp(x, xs, [np.log(kd(q, 300, d)) for q in xs])))
        ks = np.array([kx(q) for q in dd])
        if p["kappa_W_mK"] <= ks[0]:
            d_eff = dd[0]
        elif p["kappa_W_mK"] >= ks[-1]:
            d_eff = None
        else:
            d_eff = float(np.exp(np.interp(np.log(p["kappa_W_mK"]), np.log(ks), np.log(dd))))
        inf.append({"x": x, "kappa_total": p["kappa_W_mK"], "d_eff_m": d_eff,
                    "doi": p["upstream"].get("doi"), "sample": p["upstream"].get("sample_id")})
    dvals = [r["d_eff_m"] for r in inf if r["d_eff_m"] is not None]
    med = float(np.median(dvals)) if dvals else None
    res["pred"]["Pr5"] = {"n_points": len(inf), "rows": inf, "median_d_eff_m": med,
                          "pass": bool(med and 10e-9 <= med <= 300e-9)}
    json.dump(res, open(os.path.join(OUT, "nano_results.json"), "w"), indent=1)
    print(json.dumps({"gates": {k: v["pass"] for k, v in res["gates"].items()},
                      "G1max": res["gates"]["G1"]["max_abs_rel"], "G2": res["gates"]["G2"]["ratio"],
                      "G3": res["gates"]["G3"]["pure_Si_5nm_over_bulk"],
                      "pred": {k: {kk: vv for kk, vv in v.items() if kk != "rows"} for k, v in res["pred"].items()},
                      "D4": scr}, indent=1))


if __name__ == "__main__":
    main()
