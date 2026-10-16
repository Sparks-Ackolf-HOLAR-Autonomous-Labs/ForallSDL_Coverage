#!/usr/bin/env python3
"""Figure for the nano campaign from results/nano_results.json."""
import json, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
r = json.load(open(os.path.join(HERE, "results", "nano_results.json")))
DS = np.array(r["grain_sizes_m"]) * 1e9
fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
cols = plt.cm.viridis(np.linspace(0, 0.9, 7))
for c, x in zip(cols, [0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 0.85]):
    g = r["D1_undoped"][f"{x:g}_300"]
    ax[0].plot(DS, g["d"], color=c, label=f"x = {x:g}")
ax[0].set_xscale("log"); ax[0].set_yscale("log")
ax[0].set_xlabel("grain size d (nm)"); ax[0].set_ylabel("lattice kappa (W/mK)")
ax[0].set_title("(a) kappa_L vs grain size, 300 K (gray boundary model)")
ax[0].legend(fontsize=8)
for key, ls in (("T300_k2.5", "-"), ("T300_k4.0", "--"), ("T1100_k2.5", "-"), ("T1100_k4.0", "--")):
    row = [q for q in r["D2_frontier"][key] if q["d_m"] and q["x_min"] is not None]
    ax[1].plot([q["d_m"] * 1e9 for q in row], [q["x_min"] for q in row], ls,
               color="C0" if key.startswith("T300") else "C3", label=key.replace("_", ", ").replace("k", "kappa_L <= "))
ax[1].axhline(0.2, color="gray", lw=0.8); ax[1].text(6, 0.205, "RTG alloy x = 0.2", fontsize=8, color="gray")
ax[1].set_xscale("log"); ax[1].set_xlabel("grain size d (nm)"); ax[1].set_ylabel("minimum Ge fraction x")
ax[1].set_title("(b) Ge-savings frontier (model, uncalibrated)"); ax[1].legend(fontsize=8)
rows = r["pred"]["Pr5"]["rows"]
ax[2].scatter([q["kappa_total"] for q in rows], [q["d_eff_m"] * 1e9 for q in rows], color="C2")
ax[2].axhspan(10, 300, color="C1", alpha=0.15, label="pre-registered Pr5 band")
ax[2].set_yscale("log"); ax[2].set_xlabel("measured total kappa, 300 K (Starrydata, x = 0.2)")
ax[2].set_ylabel("effective grain size the model needs (nm)")
ax[2].set_title("(c) sintered Si0.8Ge0.2: model needs large d"); ax[2].legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "results", "fig_nano.png"), dpi=140)
print("ok")
