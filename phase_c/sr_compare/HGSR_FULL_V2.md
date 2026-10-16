# HG-SR full v2: the five improvement items, tested

Run 2026-10-02. The items are the five next steps at the end of
`HGSR_V2.md`. DGX CPU only; no GPU, no CHPC, no campaign-record entry
(software comparison on committed evidence).

Files: `hgsr2.py` (all items as options; with every option off it
reproduces the v2 results of `HGSR_V2.md` exactly, checked row by row), `propose.py` (I2),
`full_v2.py` (ladder, three tasks), `demo_full_props.py` (properties with
PySR in the loop), `adopted_v2.py` (adopted configuration, beam check,
scaling, properties). Results in `results/full_v2_*.json`,
`results/demo_full_props*.json`, `results/adopted_v2.json`; the first
I2 version is kept as `*_unbalanced.json`.

## What each item became

| Item | Implementation | Verdict |
|---|---|---|
| I1 physics constraints | (a) term-wise Matthiessen signs by NNLS with a required x(1-x) term; (b) a shape filter on the whole law over a grid of the declared domain: r > 0, a bigger box never conducts worse, Nordheim peak near x = 0.5; (c) the same filter also demanding a physical L -> infinity limit | (b) **adopted**; (a) and (c) rejected (hold-out 4.5% -> 20.6% and 18.6%) |
| I2 basis proposals | PySR (sqrt allowed) on the training fold, hall-of-fame terms expanded and stripped of coefficients, kept if finite on the domain and not collinear with the basis; PySR sees one cell-mean row per design cell | **optional, off by default**: improves in-campaign CV, worsens the external hold-out |
| I3 sparse search | deterministic forward beam search (width 200) | **adopted**: same law as exhaustive on all three tasks, 3 to 11x fewer fits |
| I4 noise weights | per-cell weight 1 / shrunk variance of ln kappa over distinct values (exact copies are not seeds) | **adopted**: small gain, invariance kept |
| I5 second system | SiGeC + dopant MACE-MD (P3), 71 planner points -> 16 hidden points fixed by campaign record seq 65 | **done**: v2 best on log error |

Adopted configuration ("v2.1"): F1-F4 + shape filter + noise weights +
beam search; PySR proposals on request.

## Ladder, SiGe Tersoff MD (MedAPE % / RMSE W/mK)

Trained on the alloy domain (18 runs); hold-out = 111 unseen production
runs at 40 new compositions and temperatures.

| Configuration | cv5 | LOCO | external hold-out |
|---|---|---|---|
| v2 (`HGSR_V2.md`) | 4.8 / 0.27 | 10.5 / 0.23 | 4.5 / 0.151 |
| + I1 shape filter | 5.0 / 0.12 | 8.2 / 0.16 | 4.5 / 0.151 |
| + I4 noise weights (= adopted, beam gives the same law) | 4.9 / 0.11 | 9.3 / 0.16 | **4.3 / 0.150** |
| + I2 PySR proposals | 2.2 / 0.09 | 4.0 / 0.13 | 11.6 / 0.436 |
| PySR (`HGSR_V2.md`) | 3.5 / 0.19 | 9.6 / 0.21 | 17.2 / 0.51 |
| SISSO (`HGSR_V2.md`) | 3.6 / 0.07 | 3.9 / 0.10 | 16.0 / 0.54 |
| rejected I1a, term-wise signs | | | 20.6 / 0.546 |
| rejected I1c, demanded bulk limit | | | 18.6 / 0.541 |

Trained with the ballistic endpoints (md task, 24 runs): v2 10.9%, +I1
10.9%, +I4 10.9%, +I2 10.8% (RMSE 0.31 -> 0.25), with LOCO worsening
from 3.7% to 10.7% under I2.

## Second system (I5): SiGeC + dopant, MACE-MP MD, 16 hidden points

| Engine | log RMSE | MedAPE | Spearman | pick regret |
|---|---|---|---|---|
| HG-SR v2 / adopted (same law) | **0.063** | 6.0% | **0.99** | 0 |
| HG-SR + PySR proposals (beam) | 0.065 | 4.2% | 0.99 | 0 |
| ExtraTrees (P3's registered surrogate, all 71 points) | 0.068 | 5.1% | 0.97 | 0 |
| PySR | 0.075 | 5.5% | 0.95 | 0 |
| SISSO | 0.111 | 10.2% | 0.97 | 19% |

Seed-to-seed scatter of the replicates is 0.08 to 0.12 in ln kappa, so
the HG-SR, ExtraTrees and PySR rows are within noise of each other; SISSO
is not, and its pick is 19% off. The selected law,
`r = 0.674 + 0.419 x + x(1-x) (8.56 - 2.58 t)`, contains Ge content and
temperature only: carbon (form and content) and B/P dopants were offered
and not selected, the same reading P3's registered permutation importance
gave from a black-box model. A symbolic engine states it as a law.

## Properties (md-alloy, adopted configuration)

Deterministic; replicating one design cell 5x leaves the law unchanged
(coefficient drift 2e-15); a +40% biased party in the high-T cells
creates one spurious term (1/t) and its blast radius is exactly that
term; unlearning the party restores the clean law bit for bit.

With PySR proposals in the loop (PySR rerun from scratch each time):
deterministic and exactly unlearnable, but the first version, which
showed PySR the raw rows, **broke replication invariance** (a copied
cell changed PySR's proposals and the selected law) and widened the
blast radius from one term to three. Showing PySR one cell-mean row per
design cell restores invariance (drift 1e-14); the blast radius is two
terms. Open search inherits the vote-counting problem HG-SR removes
unless it is fed the hyperedge-balanced view.

## Scaling (I3)

| Task | basis | exhaustive fits / s | beam fits / s | same law |
|---|---|---|---|---|
| md | 24 | 12 950 / 2.0 | 4 214 / 0.6 | yes |
| md-alloy | 24 | 12 950 / 0.6 | 4 733 / 0.2 | yes |
| MACE | 40 | 102 090 / 14.6 | 9 121 / 1.3 | yes |
| md-alloy, extended basis (squares, triples) | 43 | 136 697 / 6.3 | 12 190 / 0.6 | yes |
| same, 5 terms | 43 | (not run) | 18 479 / 0.9 | n/a |
| MACE + 11 PySR proposals | 51 | 272 051 / 72 | 13 178 / 2.0 | **no** |

The beam misses the exhaustive optimum once at 51 terms, so beam results
past about 40 terms are reported as beam results, not as optima. The
extended basis reads 3.5% on the hold-out, but it was not part of the
cross-validated ladder, so that number chose nothing and claims nothing.

## Reading

1. **The earlier v2 was already close to what this evidence supports.**
   The adopted items change the external hold-out from 4.5% to 4.3%; the
   gains are in reliability, not accuracy: every selected law is now
   certified physically admissible inside its declared domain (the shape
   filter rejects inadmissible fold laws, which is why LOCO improves from
   10.5% to 8.2%), search scales, and the same engine works unchanged on
   a second system and method family.
2. **Physics as a filter on the law, not a sign on each coefficient.**
   Term-wise Matthiessen signs assume the basis functions are separate
   mechanisms; they are not (x, 1-x and x(1-x) overlap), and the
   constraint removes the alloy term altogether. Shape conditions on the
   law's output are the physically stated requirement.
3. **The bulk limit is not identified by finite boxes.** Every law
   accurate inside 130 to 800 A has an unphysical or undetermined
   L -> infinity limit, and demanding a physical one costs a factor of
   four in hold-out error. The one size-free Tersoff value (GPUMD HNEMD,
   x = 0.5, 300 K, 2.23 +/- 0.22 W/mK) sits 41% above the P1 law's 1.58,
   so the 1/L two-stage law is a lower bound at mid-composition. Settling
   the form needs size-free values across x and T (an evidence purchase,
   not a constraint).
4. **Open search improves interpolation and hurts transfer.** PySR
   proposals cut in-campaign CV error by half and raised the unseen-run
   error from 4.3% to 11.6%, with terms such as
   `sqrt(5.49 sqrt(il) x + ... - 0.0235)` that carry no physical reading.
   The composition rule stands: open search proposes, HG-SR decides, and
   a proposal is adopted only when it passes the unseen-campaign test.
5. **On the second system the engine finds the factor ranking itself**:
   a two-variable law at the accuracy of the registered surrogate.

## Limits

- The shape filter's domain grid and the noise prior (3 pseudo-dof) are
  cold-start defaults set before scoring and not tuned.
- I2 was designed after seeing the v2 results; the external
  hold-out chose nothing, but one property, one protocol per system.
- The MACE hidden set is 16 points with 0.08 to 0.12 seed scatter; the
  top three engines are not separable on it.
