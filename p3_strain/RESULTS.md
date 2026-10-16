# P3 strain seed batch: results (2026-09-29)

CHPC array 15890855, 50 of 50 runs COMPLETED on notchpeak-guest, zero
preemptions. Raw hashes on the campaign record (seq 33) before this analysis;
analysis exactly as pre-registered (`analyze_strain.py`, output
`analysis_strain.json`); figure `fig_strain.png`.

## Headline: the strain response of SiGe films is anisotropic

Fitted model (ln kappa; rms residual 0.049, i.e. about 5%, the seed
scatter level): the strain sensitivity d ln kappa / d eps is

- **cross-plane (heat along the film normal): -(3.9 + 1.0 x + 3.2 T/1000)**,
  about -6 at x = 0.5, 700 K;
- **in-plane: that plus 7.0 [95% CI 5.3 to 8.8]**, so about +1, with a
  confidence interval that includes zero.

In words: compressing a SiGe film in-plane raises its cross-plane
conductivity substantially and leaves its in-plane conductivity
essentially unchanged; tension does the opposite for cross-plane.

## The materials call (pre-registered), relative to relaxed bulk

| film | T | cross-plane | in-plane |
|---|---|---|---|
| Si0.5Ge0.5 on Si(001) (eps = -2.1%) | 300 K | **+12%** [8, 16] | -3% [-7, +1] |
| Si0.5Ge0.5 on Si(001) | 1100 K | **+18%** [13, 23] | +2% [-2, +5] |
| Si0.9Ge0.1 on Si(001) (eps = -0.4%) | 300 K | +2% [1, 3] | -1% [-2, 0] |
| Si0.5Ge0.5 on relaxed Ge (eps = +2.0%) | 300 K | **-10%** [-13, -7] | +3% [-1, +7] |
| Si0.5Ge0.5 on relaxed Ge | 1100 K | **-15%** [-18, -11] | -2% [-5, +2] |

(Other compositions in `analysis_strain.json`; |eps| > 3% is flagged as
extrapolation, e.g. Si0.9Ge0.1 on Ge.)

**Why it matters for materials design.** A cross-plane thermoelectric or
a thermal barrier made from a pseudomorphic SiGe film on Si carries up to
about 18% more heat across the film than bulk data predicts, while
in-plane heat spreading is unaffected. Tabulated bulk kappa is therefore
the wrong input for cross-plane device models of strained films, and the
choice of substrate (Si vs relaxed Ge or a relaxed SiGe buffer) is a
roughly 30% lever on cross-plane kappa at mid-composition. We found no
anisotropic strain map for SiGe alloys in the literature searches (Asta
q33, q34); the nearest work is strain engineering of pure Si kappa
(Extreme Mech. Lett. 2026) and 2D materials.

## Predictions (stated before outcomes)

- Compressive in-plane strain raises kappa (cross-plane, 300 K, x = 0.5):
  **held**.
- Strain sensitivity smaller in the alloy than in dilute Si: **did not
  hold** (the x term is not significant; CI -4.2 to +2.9). The prediction
  assumed mass-disorder scattering would swamp strain effects; it does
  not in the cross-plane direction.

## Erasure analysis (P3's design claim)

The guest queue produced zero preemptions, so the physical erasure
channel was silent in this batch; both C1 (25 runs) and C2 (50 runs)
kept full pair coverage (81 of 81 pairs). The design comparison therefore
comes from exhaustive/random deletion on the committed design:

| runs lost | C1 (1-fold, 25 runs): pairs lost, P(any pair lost) | C2 (2-fold, 50 runs) |
|---|---|---|
| 1 | 2.0%, 100% | 0%, 0% |
| 2 | 4.2%, 100% | 0%, 2% |
| 3 | 6.4%, 100% | 0.1%, 7% |

C2 costs exactly twice C1 (the lower bound for this design) and removes
the single-loss failure mode entirely. The live erasure demonstration
still needs a batch that actually suffers preemptions (a larger or
longer batch, or a busier queue window).

## Caveats (method off-data, carried as scope labels)

1. **Tersoff over-relaxes the film normal**: effective biaxial Poisson
   ratio 0.88 to 1.05 vs about 0.77 for experimental Si (GATE.md). The
   cross-plane geometry is exactly where this bias acts, so the magnitude
   of the cross-plane effect is Tersoff-parented. A different method
   family (MACE on CHPC guest GPUs, or strained-lattice BTE) must check it
   before the paper states the number as physics.
2. **Box-length artifact, bounded**: at fixed cell count, the strained
   box length along the heat flow changes by the normal strain (about
   +1.8% cross-plane, -2% in-plane at eps = -2%), which through the
   finite-size slope shifts kappa by roughly +1% and -1%. That explains at
   most about 2 of the roughly 15 percentage-point anisotropy.
3. Single size (ncz = 48) and one seed per point; the absolute kappa
   carries the known Tersoff alloy suppression (P1), the relative strain
   response is the claim.
