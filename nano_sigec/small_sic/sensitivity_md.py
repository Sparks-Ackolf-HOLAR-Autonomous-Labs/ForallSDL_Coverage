#!/usr/bin/env python3
"""Post-hoc sensitivity (NOT pre-registered) for the small-SiC MD check.

The registered model side matches the boundary length d to the shell-only
control at each radius. PrS3 failed: the Si-rich pocket raises kappa by
13 to 64%, so d absorbs a composition effect the model host does not have.
Here d is instead matched once, to the bare host run (x = 0.2, no pocket),
and the particle reduction is also measured against the host. Output:
md_sensitivity.json. Reported beside, never in place of, md_results.json."""
import json, os, runpy, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
A = runpy.run_path(os.path.join(HERE, "analyze_md.py"), run_name="sens")
md, runs, SEEDS, r_model, bisect, base, T = A["md"], A["runs"], A["SEEDS"], A["r_model"], A["bisect"], A["base"], A["T"]
kh = float(np.mean([runs[("host", s)]["kappa"] for s in SEEDS]))
d_host = bisect(lambda dd: base.kappa2(T, dd) - kh, 1e-11, 1e-2)
out = {"d_host_m": d_host, "kappa_host_mean": kh, "rows": {}}
for R in (10, 15, 20):
    kp = float(np.mean([runs[(f"p{R}", s)]["kappa"] for s in SEEDS]))
    ks = md[R]["kappa_shell_mean"]
    ph = md[R]["phi_eff"]
    out["rows"][f"R{R}A"] = {
        "r_md_vs_shell": md[R]["r_mean"], "r_md_vs_host": 1 - kp / kh, "shell_gain_vs_host": ks / kh - 1,
        "r_model_density_dhost": r_model(R, d_host, ph, A["LAM_DEN"]),
        "r_model_stiffness_dhost": r_model(R, d_host, ph, A["LAM_STF"]),
        "r_model_cap_dhost": r_model(R, d_host, ph, A["LAM_CAP"])}
json.dump(out, open(os.path.join(HERE, "md_sensitivity.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
