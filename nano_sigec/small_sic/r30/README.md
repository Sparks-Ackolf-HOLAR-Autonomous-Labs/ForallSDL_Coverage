# Fixed-volume-fraction SiC radius check (R = 3 nm)

Pre-registered (`r30_design.json`, campaign-record entry in `PREREGISTERED_SHA256`)
and pushed with every script before any run existed. Run 2026-10-06.

Why: the post-hoc fit in `../../radius_contrast/` found that the effective
SiC scattering contrast falls about 15x from R = 1 to 2 nm, and two forms fit
equally well but disagree beyond 2 nm. The P1 Ge saving at equal kappa_L is
27-33% if the fall stops and 7-16% if it continues. The earlier runs also
raised the volume fraction with the radius.

Runs (8, x = 0.2, 300 K, protocol unchanged from campaign-record entries 58):
- R = 3.0 nm in an 18 x 18 x 60 cell box (phi_eff 10%, same as R = 2.0 nm before);
- R = 2.0 nm in a 16 x 16 x 53 cell box (phi_eff 4.3%, same as R = 1.5 nm before);
- for each: pocket-only control and particle, two seeds.

Decision: ratio_R = lam2(3 nm) / lam2(2 nm) at fixed phi, against the
threshold 0.435 (log-midpoint of "constant" and the fitted power law).
All ratios above: strengthen the about-30% saving; all below: correct the
saving and preferred radius to the measured curve; straddling: retain the
qualified range.

Workflow
1. `python3 preregister_r30.py`, commit + push
2. rsync this folder to CHPC, `sbatch --array=0-7 sbatch_r30.sh`
3. fetch `results/`, `python3 commit_r30.py <jobid>`, commit + push
4. `python3 analyze_r30.py` (mace venv), then `python3 apply_r30_decision.py`

## Result (2026-10-06)

- Array 16081201 (notchpeak-guest, 64 ranks): all 8 runs finished in one
  attempt each, 49-86 min. Raw outputs at campaign record seq 104 before analysis.
- Gates G1-G4 pass (closest Ge-C gap 6.0 A, C coordination 3.80-3.87,
  pocket seeds within 3.4%). Predictions PrR1-PrR3 all pass.
- kappa (W/mK): R = 3 nm box: pocket 1.896 / 1.906, particle 1.428 / 1.421
  (r = 0.247, 0.254); R = 2 nm box at phi 4.3%: pocket 1.541 / 1.490,
  particle 1.272 / 1.273 (r = 0.175, 0.146).
- ratio_R = 0.18-0.23 (threshold 0.435; power law predicts 0.19);
  ratio_phi = 0.91-1.48; ratio_Rlo = 0.24-0.43. Local slope of the
  contrast between 2 and 3 nm: -3.9.
- Decision (seq 105): branch keeps_falling, `correct`. Measured-curve Ge
  saving at equal kappa_L (x = 0.2, 0.3 um, 1 at% C): 16% / 14% at 2.5 nm
  and 9% / 7% at 5 nm (300 / 1100 K). Pick: x = 0.3, 3 at% C, 2.5 nm at
  300 K (2.30 W/mK); at 1100 K x = 0.5 is 1.2% below x = 0.3 (near tie).
  `p1_numbers.py` derives the other numbers quoted in P1.
