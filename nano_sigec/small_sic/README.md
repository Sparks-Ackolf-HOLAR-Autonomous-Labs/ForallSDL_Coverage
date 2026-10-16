# Small-SiC scattering strength (R <= 5 nm)

The SiGeC campaign ([../RESULTS.md](../RESULTS.md)) models phonon scattering
by SiC particles with the Rayleigh-to-geometric form (Kim and Majumdar 2006;
Mingo et al. 2009). The single contrast parameter is unknown, so it was
bracketed by the density-only and stiffness-only contrast. At R <= 5 nm the
bracket spans a factor of about 2 in the kappa_L reduction.

## Step 1: does the uncertainty change a decision? (no compute)

`decision_robustness.py` sweeps the contrast over every value the model form
allows, from the density-only bracket to the limit where the cross-section
saturates at the geometric value 2 pi R^2. Output: `decision_robustness.json`.

| Output (x = 0.2, 0.3 um grains, 1 at% C) | density bracket | stiffness bracket | geometric cap |
|---|---|---|---|
| kappa_L reduction, R = 2.5 nm, 300 K | 12.5% | 27.2% | 27.7% |
| kappa_L reduction, R = 5 nm, 300 K | 14.7% | 18.2% | 18.2% |
| Ge saving at equal kappa_L, R = 2.5 nm, 300 K | 30% | 51% | 51% |
| Ge saving at equal kappa_L, R = 2.5 nm, 1100 K | 27% | 50% | 50% |
| Screening pick, 300 K and 1100 K | x = 0.3, 3 at% C, R = 5 nm | x = 0.3, 3 at% C, R = 2.5 nm | x = 0.3, 3 at% C, R = 2.5 nm |

Findings:

1. The stiffness bracket already sits at the model's ceiling. Within this
   model form the band is one-sided: the only open question is how far below
   the ceiling the true contrast lies.
2. Composition and carbon content of the best recipe do not depend on the
   contrast. The "carbon replaces germanium" claim holds across the band
   (at least 27% Ge saving). What depends on it: the size of the saving
   (27% to 51%) and whether 2.5 nm or 5 nm particles are preferred.
3. So the remaining experiment is worth running only to size the saving and
   to test the model form itself, not to choose the composition.

## Step 2: atomistic check (pre-registered, campaign record seq 58; runs need CHPC)

The model form ignores the relaxed, 20%-lattice-mismatched SiC/SiGe interface
and its strain field, which is the physics most likely to matter at a few nm.
[`md_design.json`](md_design.json) registers a Tersoff-MD check: Si0.8Ge0.2
with two embedded 3C-SiC spheres (R = 1.0, 1.5, 2.0 nm), against seed-paired
controls with the same Si-rich pocket but no SiC, and the bare host
(14 runs, notchpeak-guest, about 1 to 2 h each). The decision rule maps the
measured reduction onto the model band: inside it narrows the Ge-saving
interval; above the ceiling corrects the model form (effective radius); below
the density bracket downweights the few-nm claims.

Tersoff has no Ge-C parameters. `md/make_pot.py` adds the 12 missing Ge-C
triplets by the file's own mixing rules only so LAMMPS will start; gate G1
asserts that no Ge atom ever comes within 5.7 A of a C atom, which exceeds
the Si-C plus Si-Ge cutoffs, so no term involving both is ever evaluated.
A 1200-step smoke test on a 10k-atom box (DGX CPU, 2 threads, seconds)
checked the input; no production run has been made.

Submit (after the design commit is pushed):

```
cd nano_sigec/small_sic/md && sbatch --array=0-13 sbatch_sic.sh
```

## Step 3: results (CHPC array 15904258; raw campaign record seq 59, gate record seq 60, decision seq 61)

All 14 runs completed on notchpeak-guest (one preemption, requeued), about
40 to 70 minutes each on 32 ranks. Analysis: `analyze_md.py` (committed before
any result existed) -> `md_results.json`; decision: `apply_md_decision.py`.
A first analysis pass had a bug (the Tersoff rate-scaling bisection never
applied its trial factor); it is kept as `md_results_firstpass_scalebug.json`
and disclosed in the output's `analysis_time_choices`.

| R (A) | kappa host | kappa shell-only | kappa particle | r_MD (vs shell) | model band (density to cap) | registered class |
|---|---|---|---|---|---|---|
| 10 | 1.094 | 1.232 | 1.137 | 7.6 +/- 3.3% | 1.5% to 46% | inside |
| 15 | 1.094 | 1.429 | 1.164 | 18.6 +/- 1.4% | 12% to 63% | inside |
| 20 | 1.094 | 1.793 | 1.249 | 30.3 +/- 1.5% | 40% to 75% | below density edge |

Gates: G2, G3 (C coordination 3.73 to 3.80) and G4 (host seeds 0.8% apart)
pass. **G1 fails for one run**: p10 seed 20260930 ends with a 5.63 A
geometric gap against the 5.7 A limit (the other five particle runs end at
5.80 to 6.51 A). The geometric gap is a lower bound on the true Ge-C distance,
and the other R = 10 A seed gives the same kappa within 0.3%, but the gate is
reported as failed. Predictions: PrS1 pass (30% at R = 20 A), PrS2 fail,
PrS3 fail (the Si-rich pocket alone raises kappa by 13% at R = 10 A).

Registered decision (branch `below_density`, seq 61): **downweight** the
few-nm SiC claims and report the saving at the density edge of the band:

| At x = 0.2, 0.3 um grains, 1 at% C | before (density to stiffness) | after |
|---|---|---|
| kappa_L reduction, R = 2.5 nm, 300 K | 12.5% to 27.2% | 12.5% |
| Ge saving at equal kappa_L, R = 2.5 nm, 300 / 1100 K | 30 to 51% / 27 to 50% | 30% / 27% |
| Ge saving at equal kappa_L, R = 5 nm, 300 / 1100 K | 33 to 39% / 30 to 36% | 33% / 30% |
| Screening pick | R = 2.5 or 5 nm | x = 0.3, 3 at% C, R = 5 nm (density-edge pick) |

Post-hoc sensitivity (`sensitivity_md.py`, not registered): PrS3 shows the
registered model side is confounded, because matching d to the shell-only
control makes d absorb the pocket's composition effect (d grows from 0.44 to
7.3 um with R). Matching d once to the bare host instead puts all three radii
inside the band, at its lower end (R = 20 A: MD 30% against 28% to 60%). Both
readings agree: the small-particle scattering sits at or near the density
edge, not the stiffness edge.

Materials reading:

1. The explicit, relaxed SiC/SiGe interface does not scatter more than the
   continuum model allows; the stiffness-contrast upper estimate is not
   supported at R = 1 to 2 nm.
2. The Si-rich zone around each particle raises kappa by 13% to 64% relative
   to the plain alloy, more than the particle removes: every particle-plus-
   pocket run conducts 4% to 14% better than the plain Si0.8Ge0.2 host. MACE
   says carbon prefers Si neighbours, so a Ge-depleted zone around SiC
   precipitates is expected; if it forms, it works against the particle.
   Caveat: the 7 A pocket width was imposed so that no Ge-C term is ever
   evaluated, not relaxed from physics, and the volume fractions (1% to 10%)
   are higher than in the campaign.
3. The "carbon replaces germanium" claim holds at the lower edge, about 30%
   less Ge; the 51% upper figure is withdrawn.

## Step 4: composition check at density match, x = 0.3 (`x03/`; prereg seq 78, raw seq 79, decision seq 80)

Trigger: the density-only contrast vanishes near x = 0.3 (host density equals
SiC), so the recipe composition was scope-limited to x = 0.3 to 0.5 (seq 77).
CHPC array 15961370 (notchpeak-guest, 10 runs, 50 to 77 min each, all exit 0)
repeated the embedded-SiC MD in a Si0.7Ge0.3 host. All four gates and all three
predictions passed.

| R | reduction at x = 0.2 | reduction at x = 0.3 |
|---|---|---|
| 1.5 nm | 18.6% | 14.1% |
| 2.0 nm | 30.3% | 24.6% |

rho = 0.784 (registered threshold 0.75), branch `insensitive`: SiC scatters
mainly through stiffness and interface strain, so the x = 0.2 calibration
carries to other compositions. Correct operation (seq 80): pick restored to
x = 0.3, 3 at% C, R = 5 nm (kappa_L 1.91 / 1.81 W/mK at 300 / 1100 K), seq-77
scope limit lifted. Margin is narrow: recombining single seeds gives
rho 0.64 to 0.95, 9 of 16 combinations above 0.75. Details in
`x03/x03_results.json`.
