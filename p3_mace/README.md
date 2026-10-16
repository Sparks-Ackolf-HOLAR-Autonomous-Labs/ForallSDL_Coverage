# P3: scoring the batch designs on MACE-MD evidence the planner never sees

Run 2026-09-30. Design: `p3_mace_design.json` (campaign record seq 64). Hidden
set: `hidden_set.json` (campaign record seq 65), drawn from the hash of the seq-64
entry, so it was fixed by the record, not by the planner or the author.

## Why

The P3 design study (`p3_sigec/`, campaign record seq 56-57) scored covering, random
and one-factor-at-a-time batches against the same physics model the planner
was built around. That shows the designs behave as intended on a known
surface, not that they help when the measurement is independent. Here the
measurement is a different method family: molecular dynamics with the
MACE-MP foundation potential.

## What runs

- Space the MD box can hold: x_Ge (5 levels) x C content (0, 1, 3 at%) x
  carbon form (dispersed on lattice sites / clustered as a sub-nanometre
  3C-SiC cluster) x dopant (none, B, P at 0.2 at%) x T (300, 700, 1100 K):
  225 distinct points.
- Planner batches from factor levels only (`design_p3_mace.py`): covering
  CA1 = 15 runs (equal to the lower bound 5 x 3), 1-loss-tolerant CA2 = 30,
  one-factor-at-a-time 12, random pool 30.
- 16 hidden points plus 2 seed repeats. 89 runs in total.
- Each run: 1152 atoms, relaxed, 8 ps Langevin, 50 ps Muller-Plathe NVE
  (`mace_rnemd.py`), CHPC notchpeak guest GPUs, array 15922934, about 80 min
  on an A6000, checkpointed every 5000 steps.

## Results (2026-09-30)

89 of 89 runs collected (19 rerun as array 15950354 after an uncorrectable GPU
ECC error on notch458 killed them at start-up). Raw outputs on the campaign record at
seq 72 before analysis; registered analysis at seq 73
(`results/p3_mace_results.json`). Post-hoc context, not used for any verdict:
`sensitivity_p3_mace.py` -> `results/p3_mace_posthoc.json`.

| Batch | Runs | Hidden log RMSE | Spearman | Regret of its pick | Pair coverage |
|---|---|---|---|---|---|
| Covering CA1 | 15 | 0.119 | 0.96 | 19% | 100% |
| Loss-tolerant covering CA2 | 30 | 0.090 | 0.96 | 0% | 100% |
| One factor at a time | 12 | 0.117 | 0.90 | 40% | 54% |
| Random pool | 30 | 0.076 | 0.98 | 0% | 93% |
| Random 15 (median of 200) | 15 | 0.087 | 0.95 | 0% | |
| Random 12 (median of 200) | 12 | 0.095 | 0.94 | 0% | |

Replicate seeds differ by 0.08 and 0.12 in log kappa; hidden log kappa spread
0.24 (0.32 to 0.77 W/mK).

Registered predictions: Q1 **refuted** (covering 15 sits at the 93rd percentile
of random 15-run batches); Q2 holds (one-factor-at-a-time worse than random,
80th percentile, and its pick is 40% off); Q3 **refuted** (ranking x_Ge >> T >
C content > dopant > carbon form, but the last three are all below the seed
noise, so their order carries no information); Q4 holds nominally
(dispersed/clustered 0.99 over 15 matched pairs, inside noise); Q5 holds
(carbon-free minimum at x = 0.5).

Reading:
- In a 1152-atom MACE box Ge content carries the response and temperature
  moves it more weakly (permutation importance 0.148 and 0.012, ln kappa
  shifts 0.1 to 0.32 from 300 to 1100 K, above the 0.08 to 0.12 seed noise;
  corrected 2026-10-11, earlier text said one factor); carbon
  at 1-3 at% (dispersed or as a sub-nanometre cluster) and 0.2 at% B or P change
  kappa by less than the seed-to-seed scatter. This agrees with the SiGeC
  response-model campaign: dopants are thermally inert and carbon only matters as few-nm
  SiC, which this box cannot hold.
- With three of five factors inert, covering buys no prediction accuracy over
  random: its balanced spread across inert factors gives the small-batch
  surrogate splits to overfit. Post hoc, dropping the inert factors brings
  covering 15 from 0.119 to 0.081 and random 30 from 0.076 to 0.064.
- What survives on independent evidence is the guarantee, not the accuracy:
  both covering designs sample every pair of levels, one-factor-at-a-time
  samples 54% and picks a recipe 40% off, and a random batch of the same size
  as the loss-tolerant design misses 7% of pairs without saying which.
