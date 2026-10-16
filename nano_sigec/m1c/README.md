# Stage M1c: converged MACE rerun (where carbon goes)

Run 2026-09-30 UTC. Pre-registered in `m1c_design.json`
(campaign record seq 53, pushed in commit 21426d4f before submission).

- Job: CHPC array **15898834** (51 tasks, notchpeak-guest / owner-guest, 16 cores, 8 h limit, requeue + 25-step checkpoints)
- Run directory on CHPC: `nano_sigec_m1c/` (out/, ckpt/, logs/)
- Script: `chpc/mace_relax_m1c.py` (FIRE to fmax 0.05 or 800 steps, then LBFGS to 0.005; cap 5000 total)
- Smoke test before submission: bulk Si converged in 8 steps, a = 5.4554 A (identical to M1)

Processing order: raw `out/` files hashed onto the campaign record by `commit_m1c.py` (batch_result, before analysis), then `summarize_m1c.py`, then the decision rule in the design.
