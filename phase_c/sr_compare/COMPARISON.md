# SR engine comparison + HG-SR (2026-09-27)

Aim: compare different symbolic learners (PySR /
SymbolicRegression.jl, SISSO, AI Feynman, QLattice/Feyn), then build a new one that is hypergraph-aware, since none of the
existing engines is.

## Protocol

Three benchmark tasks from the project's real evidence (`datasets.py`),
all fitting resistivity r = 1/kappa (Matthiessen additivity), scored on
kappa:

| Task | Rows | Features | What it probes |
|---|---|---|---|
| `literature` | 17 | x, T | the heterogeneous literature party |
| `md` | 24 | x, T, L | CHPC-MD production points incl. finite size |
| `combined` | 20 | x, T | two parties (literature + MD-extrapolated) |

Per engine and task (`run_compare.py`): full fit (expression,
complexity, runtime), fit-twice determinism check, deterministic 5-fold
CV. Errors on kappa: RMSE (endpoint-dominated) and median absolute
percentage error (scale-free). Raw outputs in `results/`.

## Availability on this stack

| Engine | Status |
|---|---|
| minisr (pilot) | runs; doubles as the HG-SR pooled ablation |
| HG-SR (new, `hgsr.py`) | runs |
| PySR 2.5.0 | runs on the DGX (aarch64), deterministic mode |
| SISSO (TorchSisso 0.1.8) | runs on the DGX after `pip install TorchSisso` |
| AI Feynman 2.0.7 | cannot build on the DGX (needs gfortran; none, no sudo); **built and ran on CHPC** (job 15855552, notchpeak-guest), which produced no solutions (see below) |
| QLattice / Feyn 3.x | closed-source wheels are x86_64/win/mac only, so no native DGX run; **but it runs on the DGX inside the amd64-emulated udocker container** (see `../../containers/README.md`): Feyn 3.5.0 installs, imports, and fits under qemu-x86_64 at roughly 1 min/epoch on these tables. Adapter in `engines.py`; results are written to `results/results_feyn_*.json` |

## Results (5-fold CV)

| Engine | Task | CV MedAPE | CV RMSE (W/mK) | Complexity | Fit time | Deterministic |
|---|---|---|---|---|---|---|
| minisr-pooled | literature | 100.0% | 122.3 | 3 | <0.1 s | yes |
| hgsr | literature | 100.0% | 123.5 | 3 | <0.1 s | yes |
| pysr | literature | 87.2% | 125.8 | 5 | 7.7 s | yes |
| sisso | literature | 96.0% | 131.7 | 4 | <0.1 s | yes |
| minisr-pooled | md | 8.3% | 23.2 | 6 | <0.1 s | yes |
| hgsr | md | 11.8% | 17.7 | 5 | <0.1 s | yes |
| pysr | md | 8.8% | 10.4 | 13 | 1.8 s | yes |
| sisso | md | 50.8% | 10.0 | 8 | 0.1 s | yes |
| feyn | md | 6.1% | 38.7 | 10 | 459 s | yes (seed-pinned) |
| feyn | literature | 80.3% | 117.2 | 2 | 405 s | yes (seed-pinned) |
| minisr-pooled | combined | 100.0% | 126.8 | 3 | <0.1 s | yes |
| hgsr | combined | 100.0% | 562.2 | 3 | <0.1 s | yes |
| pysr | combined | 100.0% | 119.9 | 7 | 1.8 s | yes |
| sisso | combined | 97.2% | 123.4 | 7 | <0.1 s | yes |
| feyn | combined | 81.1% | 109.5 | 4 | 415 s | yes (seed-pinned) |

Selected expressions (full fits, r = 1/kappa):

- minisr-pooled / md: `r = -0.148 + 3.19 x - 3.18 x^2 + 94.2/L`
- hgsr / md: `r = 3.11 x - 3.19 x^2 + 65.7/L`
- pysr / md: `r = 0.091 + (2.045 + 431.6/L) * x(1-x)`
- sisso / md: `r = -0.128 L/T + 225.3 x/L - 0.482 x + 0.757`

## Reading the results

1. **The md task is learnable; the literature task is not, and that is
   a shard finding, not an engine failure.** Every engine lands at
   87-100% CV MedAPE on the literature and combined tasks: the 17
   literature points mix bulk crystals, thin films, and nanostructured
   samples whose kappa at the same (x, T) differs by two orders of
   magnitude. No expression in (x, T) can fit that; the party needs
   scope-splitting (film vs bulk vs nano as separate scopes or an
   explicit microstructure feature) before a symbolic law is meaningful.
   The engines agree on the diagnosis by all failing the same way.

2. **PySR found real physics the fixed-basis engines missed.** Its md
   expression factors the size term as `(2.045 + 431.6/L) * x(1-x)`:
   the boundary correction rides on the alloy-scattering term
   (alloy-limited mean free path times a 1/L ballistic term), and it
   halves the CV RMSE vs the fixed-basis fits. We adopted `x(1-x)/L`
   into the HG-SR basis with a comment crediting the discovery: the
   open-search engine acts as a basis-discovery step for the
   deterministic engine. This engine-to-engine handoff is itself a
   meta-learning lesson for the campaign record.

3. **SISSO is fast and RMSE-competitive but its md expression
   (`L/T`, `x/L` terms) extrapolates worse mid-grid** (51% CV MedAPE):
   with n_expansion=2 its feature space contains many size-composition
   ratios that fit the endpoints but not the alloy plateau.

4. **Determinism: every engine that runs here passes the fit-twice
   check** (PySR needs `deterministic=True, parallelism='serial',
   random_state=0`). But only the exhaustive-basis engines (minisr,
   HG-SR) make the *clean retrain* exact by construction; PySR's
   determinism is seed-pinned, so a retrain after unlearning is
   reproducible but the epsilon=delta=0 argument holds only for the
   pinned seed, and SISSO's screening is deterministic given the
   feature space.

5. **hgsr vs minisr-pooled on md** shows the fairness trade: pooled
   optimizes the densely sampled regimes (better MedAPE, 8.3 vs 11.8),
   balanced optimizes cross-regime error (better RMSE, 17.7 vs 23.2)
   and refuses to let dense regimes outvote thin ones.

## What HG-SR adds that no compared engine has (`hgsr.py`, demos in `demo_bias.py`)

Flat engines consume (X, y); HG-SR consumes the evidence hypergraph
(rows tagged with regime/party/source hyperedges):

- **Hyperedge-weighted fitting** (weight 1/|regime edge|, BIC on the
  hyperedge count): provably invariant to replicating a regime's rows.
  Demo 1: replicating the x=0.5/300 K regime 5x drifts the pooled fit's
  coefficients by 26% and changes nothing in HG-SR (drift exactly 0).
- **Per-term hyperedge provenance and scope labels**: leave-one-
  hyperedge-out refits tag every term with the evidence groups it
  depends on (`multi-supported` / `regime-sensitive` / `scope:<party>`).
- **Blast radius**: inverting the provenance gives, per hyperedge,
  which terms fall if it is unlearned. Demo 2: a +40% biased party
  dominating the thin high-T regime manufactures a spurious `T` term;
  the blast radius of `party:BIASED-MD` is exactly `["T"]`, and
  unlearning the party restores the clean model **bit-for-bit**
  (epsilon = delta = 0, checked with atol=0).
- **Coverage-scoped validity**: the domain is the set of covered regime
  cells, not a bounding box.

None of PySR, SISSO, AI Feynman, or QLattice can express these: they
have no concept of evidence groups, so they can neither localize a
biased party nor certify its removal. HG-SR's search core is still
exhaustive-basis OLS, so it deliberately trades expression discovery
(where PySR wins) for exact revision semantics; the two compose (point
2 above).

## AI Feynman verdict (CHPC job 15855552, artifacts in `aif_chpc/`)

The build succeeded on CHPC (recipe in `aif_chpc/README.md`: python
3.11 for numpy.distutils, `--no-build-isolation`, the deprecated
`sklearn` override, then scipy/scikit-learn). The run is the verdict:

- `literature` and `combined`: the full pipeline (brute force over all
  its variable transforms, polyfit, compositionality and symmetry
  checks) completed in ~40 s each and wrote **zero-byte solution
  files** — nothing met its acceptance threshold on 17-20 noisy,
  heterogeneous real points.
- `md`: crashed inside the pipeline with
  `LinAlgError('Array must not contain infs or NaNs')` at 54 s.

This matches AI Feynman's design envelope: it targets large,
near-noiseless synthetic datasets (the Feynman equations corpus, 1e4+
rows). On this project's small real-world evidence tables it is not a
usable engine, while PySR/SISSO/HG-SR all return laws on the same
files. Comparison row: *ran, no output* — which is itself a result.

## Feyn results (emulated container on the DGX)

On the literature task, Feyn's
QLattice reaches **80.3% CV MedAPE / 117.2 RMSE** with expression
`r = 0.979 - 0.688(1.171 - 2.335 x)^2` (complexity 2, 405 s/fit under
qemu emulation). Best MedAPE of any engine on that task, and still
unusable: Feyn agrees with all four other engines that the raw
literature party admits no law in (x, T), reinforcing the
scope-splitting diagnosis (point 1 above) with a fifth independent
architecture (simulated-annealing QLattice vs GP evolution vs
screening vs brute force vs exhaustive OLS).

The md task: **6.1% CV MedAPE, the best
interpolation of any engine**, but the worst CV RMSE (38.7 W/mK) on
the endpoint-dominated folds. Its expression is a sum of Gaussian
bumps in (x, L) plus a linear T term (complexity 10): QLattice's
kernel-like basis interpolates the alloy plateau beautifully and
extrapolates poorly exactly where phonon transport leaves the sampled
regime, and none of its terms map onto scattering physics the way
PySR's x(1-x)/L does. Same division-of-labor conclusion from a new
direction: Feyn is a strong interpolating surrogate, not a
theory-term discoverer, and it cannot carry revision semantics. The
`combined` row (collected 2026-09-27) closes the run: 81.1% CV MedAPE,
best RMSE of the field on that task but the same order-of-magnitude
failure as every other engine, because the combined task inherits the
unpooled-literature pathology. All three Feyn rows are now committed;
the benchmark matrix is complete for every available engine.

The literature-party scope-splitting finding is taken up by the CHD
stage-0 pass (`../chd/STAGE0_RESULTS.md`).

## Update 2026-10-02: HG-SR v2 (see `HGSR_V2.md`)

The HG-SR row above is v1. Its trailing accuracy came from four search
and selection settings (BIC with an effective sample size of 5, regime
hyperedges that ignore box length, absolute loss on r, a hand-priced
fixed basis), not from the hyperedge weighting. v2 (`hgsr2.py`) fixes
all four and keeps determinism, replication invariance, blast radius and
bit-exact unlearning. On the md task its five-fold CV error drops from
11.8% / 17.7 W/mK to 4.2% / 2.2 W/mK. On 111 unseen production runs,
trained on the alloy domain, it reads 4.5% against 17.2% (PySR) and
16.0% (SISSO), while PySR and SISSO still win the row-CV interpolation
contest on that small set.

## Update 2026-10-02: full v2 (see `HGSR_FULL_V2.md`)

The five improvement items were implemented and tested. Adopted: a
physics shape filter on the whole law, per-cell noise weights, and beam
search (same law as exhaustive on all three tasks, 3 to 11x fewer fits).
Rejected: term-wise Matthiessen signs and a demanded bulk limit (hold-out
4.5% -> 20.6% / 18.6%). PySR basis proposals stay optional: they halve
in-campaign CV error but raise the unseen-run error from 4.3% to 11.6%.
On a second system (SiGeC + dopant MACE-MD, 16 hidden points fixed by the
campaign record) HG-SR reads log RMSE 0.063 against 0.068 (ExtraTrees), 0.075
(PySR), 0.111 (SISSO), selecting a law in Ge content and temperature only.
