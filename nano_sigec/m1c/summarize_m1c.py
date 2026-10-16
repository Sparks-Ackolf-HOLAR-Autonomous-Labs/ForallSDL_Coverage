#!/usr/bin/env python3
"""Reduce the ledgered M1c outputs (out/*.json, ledger seq 54) to
m1c_summary.json in the same schema as ../m1/m1_summary.json (so
analyze_sigec.py can consume it), then evaluate gates M1c-a/b/c and
predictions PM1-PM4 against the M1 values exactly as pre-registered
(m1c_design.json, ledger seq 53)."""
import glob, json, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
R = {json.load(open(p))["task"]: json.load(open(p)) for p in glob.glob(os.path.join(HERE, "out", "*.json"))}
M1 = json.load(open(os.path.join(HERE, "..", "m1", "m1_summary.json")))
expt = {"Si": 5.431, "Ge": 5.658, "SiC": 4.360}
b = {k: R[f"bulk_{k}"] for k in expt}
m1a = {k: {"a_mace": b[k]["a_A"], "a_expt": v, "rel": b[k]["a_A"] / v - 1} for k, v in expt.items()}
S = {"M1a": {"rows": m1a, "pass": all(abs(r["rel"]) <= 0.025 for r in m1a.values())}}
mu_si = b["Si"]["E_eV"] / b["Si"]["n_atoms"]
e_sic = b["SiC"]["E_eV"] / (b["SiC"]["n_atoms"] / 2)
gate_a = all(r["gate_pass"] for r in R.values())
per_x, r2ok = {}, True
for xs in ("0", "0.2", "0.5"):
    ys, lna, dprec, dsite = [], [], [], []
    for s in (1, 2, 3):
        g = lambda lab: R[f"{xs}_{s}_{lab}"]
        for nc in (0, 2, 4):
            q = g(f"nC{nc}"); ys.append(nc / q["n_sites"]); lna.append(np.log(q["a_A"]))
        dprec.append(g("nC0")["E_eV"] + e_sic - g("nC1")["E_eV"] - 2 * mu_si)
        if xs != "0":
            dsite.append({"dE_site": g("maxGe")["E_eV"] - g("minGe")["E_eV"],
                          "nGe_nn_max": g("maxGe")["n_Ge_nn"], "nGe_nn_min": g("minGe")["n_Ge_nn"]})
    ys, lna = np.array(ys), np.array(lna)
    A = np.vstack([ys, np.ones_like(ys)]).T
    coef, *_ = np.linalg.lstsq(A, lna, rcond=None)
    r2 = 1 - ((lna - A @ coef) ** 2).sum() / ((lna - lna.mean()) ** 2).sum()
    r2ok &= r2 >= 0.9
    per_x[xs] = {"n_configs": 3, "delta_C": float(coef[0]), "R2": float(r2), "a0_A": float(np.exp(coef[1])),
                 "dE_prec_eV": dprec, "dE_prec_mean_eV": float(np.mean(dprec)), "site_pref": dsite,
                 "dE_site_mean_eV": float(np.mean([d["dE_site"] for d in dsite])) if dsite else None}
S["M1b"] = {"all_converged": bool(gate_a), "all_R2_ge_0.9": bool(r2ok), "pass": bool(gate_a and r2ok),
            "max_fmax_eV_A": max(r["fmax_norm"] for r in R.values()),
            "max_abs_pressure_GPa": max(abs(r["pressure_GPa"]) for r in R.values()),
            "max_steps": max(r["steps"] for r in R.values()), "n_relaxations": len(R)}
S["per_x"], S["mu_Si_eV"], S["E_SiC_fu_eV"] = per_x, mu_si, e_sic
S["modified_vegard_delta_C"] = (4.360 - 5.431) / (0.5 * 5.431)
S["gates_pass"] = S["M1a"]["pass"] and S["M1b"]["pass"]
S["source"] = "M1c converged rerun (ledger seq 53 pre-registration, seq 54 raw)"

P, cmp = {}, {}
for xs in per_x:
    cmp[xs] = {k: {"M1": M1["per_x"][xs][k], "M1c": per_x[xs][k]} for k in ("delta_C", "dE_prec_mean_eV", "dE_site_mean_eV")}
P["PM1"] = {"pass": S["gates_pass"], "M1c-a": gate_a, "M1c-b": bool(r2ok), "M1c-c": S["M1a"]["pass"]}
dp = {x: per_x[x]["dE_prec_mean_eV"] for x in per_x}
P["PM2"] = {"shift_eV": {x: dp[x] - M1["per_x"][x]["dE_prec_mean_eV"] for x in dp},
            "pass": all(abs(dp[x] - M1["per_x"][x]["dE_prec_mean_eV"]) <= 0.10 and dp[x] < -0.5 for x in dp) and dp["0.5"] < dp["0"]}
ds = {x: per_x[x]["dE_site_mean_eV"] for x in ("0.2", "0.5")}
P["PM3"] = {"shift_eV": {x: ds[x] - M1["per_x"][x]["dE_site_mean_eV"] for x in ds},
            "pass": all(ds[x] > 0 and abs(ds[x] - M1["per_x"][x]["dE_site_mean_eV"]) <= 0.05 for x in ds)}
dc = {x: per_x[x]["delta_C"] for x in per_x}
P["PM4"] = {"rel_shift": {x: dc[x] / M1["per_x"][x]["delta_C"] - 1 for x in dc},
            "rel_vs_vegard_x0": dc["0"] / S["modified_vegard_delta_C"] - 1,
            "pass": all(abs(dc[x] / M1["per_x"][x]["delta_C"] - 1) <= 0.05 for x in dc)
                    and abs(dc["0"] / S["modified_vegard_delta_C"] - 1) <= 0.30}
S["comparison_vs_M1"], S["predictions"] = cmp, P
json.dump(S, open(os.path.join(HERE, "m1c_summary.json"), "w"), indent=1)
print(json.dumps({"gates": S["M1b"], "M1a": S["M1a"]["pass"], "cmp": cmp, "pred": P}, indent=1))
