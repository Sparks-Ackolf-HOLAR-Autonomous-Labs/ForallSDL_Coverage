# Fixed-volume-fraction SiC radius check at R = 5 nm

Pre-registered (`r50_design.json`, campaign-record entry in `PREREGISTERED_SHA256`)
and pushed with every script before any run existed. Run 2026-10-07.

Why: after campaign-record entries 105 the measured SiC contrast curve stops at R = 3 nm and
P1 extends it to 5 nm with the 2-to-3 nm slope. The 5 nm Ge saving (9% at
300 K, 7% at 1100 K) is the only P1 number resting on that extrapolation.

Runs (4, x = 0.2, 300 K, protocol unchanged from campaign records 58 and 103):
R = 5.0 nm in a 30 x 30 x 100 cell box (720000 sites, phi_eff 10%, same as
R = 3.0 nm in campaign-record entries 104); pocket-only control and particle, two seeds.

Model predictions written before any run: reduction r = 0.18 if the
contrast keeps falling with the measured slope, 0.46 if it stops falling
(ledger-104 seed scatter at 3 nm: 0.247 / 0.254), so one run pair decides.

Decision: ratio_R5 = lam2(5 nm) / lam2(3 nm) at fixed phi against 0.369.
All below: strengthen (label becomes "measured to 5 nm"); all above:
correct the 5 nm saving and the pick; straddling: retain.

Workflow
1. `python3 preregister_r50.py`, commit + push
2. rsync this folder to CHPC; grn061 if no sparks-partition job runs there
   (`sbatch -M granite -p sparks-gpu-grn -A sparks -w grn061 --array=0-3 sbatch_r50.sh`),
   otherwise `sbatch --array=0-3 sbatch_r50.sh` (notchpeak-guest)
3. fetch `results/`, `python3 commit_r50.py <jobid>`, commit + push
4. `python3 analyze_r50.py` (MACE environment), then `python3 apply_r50_decision.py`

Submitted 2026-10-07:
CHPC array **16114853** (`--array=0-3`, notchpeak-guest / owner-guest),
design sha256 checked against campaign-record entries 106 on the CHPC copy before `sbatch`. grn061 not used (guest nodes idle, four
runs in parallel instead of serial on one node).

Collected 2026-10-08: all four runs exit 0 (one particle run
preempted once and requeued on notch171). Raw committed at campaign-record entries 108 before
analysis. Gate G1 missed in one run (p50_s20261007 ended at 5.62 A against
5.7 A, campaign-record entries 109); the branch without that run is the same, so the decision
was taken. ratio_R5 = 0.144-0.163 over the four seed pairings (threshold
0.369; R^-4 expectation 0.13); r(5 nm) = 0.190 / 0.198. Branch
**keeps_falling** -> **strengthen** (campaign-record entries 110). Measured-curve claims: Ge
saving 16.2% / 13.7% at 2.5 nm and 9.3% / 7.8% at 5 nm (300 / 1100 K); pick
unchanged (x = 0.3 at 300 K, 3 at% C, 2.5 nm). PrF1-PrF3 all pass. P1 label
"extrapolated beyond 3 nm" replaced by "measured to 5 nm".
