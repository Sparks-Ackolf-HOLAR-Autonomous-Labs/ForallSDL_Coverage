#!/usr/bin/env python3
"""Large-radius bracket for the SiC contrast (not pre-registered; reporting
check after ledger seq 105).

The MD checks measured the effective SiC contrast lam2 only for R = 1-3 nm,
where it falls about as R^-4. P1's numbers for R >= 10 nm still use one
constant contrast (the density edge adopted at seq 61). This script asks
whether any P1 claim about large SiC depends on that choice, by evaluating
each claim under the two limits: (a) the constant density-edge contrast,
and (b) the measured curve carried past 3 nm at its last local slope (the
rule of r30/p1_numbers.py). If the claim holds at both limits, the untested
range cannot overturn it. Output: large_radius_bracket.json.
"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
NS = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, NS)
os.environ.setdefault("SIGEC_M1_SUMMARY", os.path.join(NS, "m1c", "m1c_summary.json"))
import analyze_sigec as A  # noqa: E402

r = json.load(open(os.path.join(HERE, "r30_results.json")))
xs = np.log([c[0] for c in r["measured_curve"]]); ys = np.log([c[1] for c in r["measured_curve"]])


def lam_meas(Rn):
    if Rn <= 3.0:
        return float(np.exp(np.interp(np.log(Rn), xs, ys)))
    s = min(0.0, (ys[-1] - ys[-2]) / (xs[-1] - xs[-2]))
    return float(np.exp(ys[-1] + s * (np.log(Rn) - xs[-1])))


lam_const = lambda Rn: A.contrast(0.2, "density")  # noqa: E731
m = A.SModes(0.2)
out = {"lam2_density_edge_x0.2": lam_const(1), "rows": [], "zhao": []}
for Rn in (2.5, 5, 10, 25, 50, 100, 150):
    k0 = m.kappa2(300, 0.3e-6)
    row = {"R_nm": Rn}
    for lab, lf in (("constant", lam_const), ("measured_extrapolated", lam_meas)):
        k1 = m.kappa2(300, 0.3e-6, 0.0, np_R=Rn * 1e-9, np_phi=A.phi(0.2, 0.01), lam2=lf(Rn))
        row[lab + "_pct"] = 100 * (k1 / k0 - 1)
    out["rows"].append(row)
for Rn in (50, 100, 150):
    for d in A.DEFF:
        row = {"R_nm": Rn, "d_eff_m": d}
        k0 = A.lnT_interp(m.kappa2(900, d), m.kappa2(1100, d))
        for lab, lf in (("constant", lam_const), ("measured_extrapolated", lam_meas)):
            kk = [m.kappa2(T, d, 0.0, np_R=Rn * 1e-9, np_phi=A.phi(0.2, 0.03), lam2=lf(Rn)) for T in (900, 1100)]
            row[lab + "_pct"] = 100 * (A.lnT_interp(*kk) / k0 - 1)
        out["zhao"].append(row)
out["verdict"] = {
    "large_SiC_ge100nm_below_2pct_both_limits": all(abs(r[k]) < 2 for r in out["rows"] if r["R_nm"] >= 100 for k in r if k.endswith("pct")),
    "zhao_below_10pct_both_limits": all(abs(r[k]) < 10 for r in out["zhao"] for k in r if k.endswith("pct")),
    "zhao_largest_pct": {lab: min(r[lab + "_pct"] for r in out["zhao"]) for lab in ("constant", "measured_extrapolated")}}
json.dump(out, open(os.path.join(HERE, "large_radius_bracket.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
