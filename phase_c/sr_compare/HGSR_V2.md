# HG-SR v2: why v1 trailed PySR and SISSO, and what fixed it

Run 2026-10-02. Question: why does
HG-SR only reach results comparable to PySR or SISSO, and not better, and
how can it be improved?

Files: `hgsr2.py` (engine), `improve_hgsr.py` (ablation ladder and
external hold-out), `demo_v2_props.py` (property checks); raw numbers in
`results/improve_hgsr.json` and `results/demo_v2_props.json`. `hgsr.py`
(v1) is unchanged, so every earlier demo and campaign-record entry still
reproduces.

## Diagnosis: four causes, none of them the hypergraph idea itself

| | Cause in v1 | Effect on the md task |
|---|---|---|
| D1 | Model chosen by BIC with effective sample size = number of regime hyperedges (5 on md) | the complexity penalty dominates; selection stops at `x + x^2 + 1/L` although better 4-term laws sit on its own front |
| D2 | Regime hyperedge = (x band, T band), box length ignored | a 12-run finite-size series at x = 0.5, 300 K counts as one unit of evidence |
| D3 | Least squares on absolute r = 1/kappa | pure Si/Ge rows (r about 0.03) carry almost no weight, while the kappa score is dominated by them |
| D4 | Hand-written basis with hand-set term costs | the PySR-discovered `x(1-x)/L` was priced at cost 7 and never selected; no term carries a slow L^-1/2 size convergence |

So the gap came from the search and selection settings wrapped around
the hypergraph loss, not from the hyperedge weighting.

## Fixes (v2)

- F1: select along the complexity front by leave-one-hyperedge-out
  cross-validation on log kappa (deterministic, so unlearning stays exact).
- F2: hyperedges are design cells (x, T, size); seed replicates of a cell
  stay together.
- F3: relative (percentage) loss, weights 1/r^2.
- F4: generated basis: constant, seven physical primitives (x, 1-x,
  x(1-x), T, 1/T, 1/L, L^-1/2), and products of two primitives; every
  primitive factor costs 1. 24 candidate terms, exhaustive up to 4 terms,
  0.5 s per fit.

## Results

Three yardsticks: the five-fold row CV of the original comparison
(`COMPARISON.md`, 2026-09-27; cv5), leave-one-design-
cell-out CV (LOCO, stricter: an unseen composition-size cell), and an
external hold-out the engines never see: 111 standard-protocol runs from
the 2026-09-28 production campaign at 40 new compositions and
temperatures (x 0.05 to 0.90, T below the near-melting probes). Median
absolute percentage error on kappa; RMSE in W/mK.

Trained on the original md task (24 runs, includes the ballistic pure
Si and Ge endpoints):

| Engine | cv5 | LOCO | external hold-out |
|---|---|---|---|
| HG-SR v1 | 11.8% / 17.7 | 13.8% / 18.8 | 9.3% / 0.20 |
| + F1 selection | 6.4% / 3.8 | 8.9% / 3.3 | 11.4% / 0.26 |
| + F2 cells | 6.5% / 10.1 | 9.2% / 3.4 | 11.4% / 0.27 |
| + F3 relative loss | 7.4% / 1.1 | 6.8% / 1.1 | 12.5% / 0.28 |
| + F4 basis (= v2) | **4.2% / 2.2** | **3.9% / 2.0** | 10.9% / 0.31 |
| PySR | 8.8% / 10.4 | 12.4% / 14.1 | 13.1% / 0.24 |
| SISSO | 50.8% / 10.0 | 74.7% / 10.3 | 21.0% / 0.45 |

Trained on the alloy domain only (x in [0.05, 0.95], the scope-limit
ledgered on 2026-09-27; 18 runs):

| Engine | cv5 | LOCO | external hold-out |
|---|---|---|---|
| HG-SR v1 | 6.0% / 0.28 | 9.8% / 0.30 | 22.3% / 0.57 |
| + F1 selection | 5.2% / 0.12 | 7.3% / 0.14 | 14.5% / 0.44 |
| + F2 cells | 7.1% / 0.11 | 11.9% / 0.17 | 27.7% / 0.59 |
| + F3 relative loss | 9.3% / 0.11 | 12.1% / 0.16 | 11.5% / 0.32 |
| + F4 basis (= v2) | 4.8% / 0.27 | 10.5% / 0.23 | **4.5% / 0.15** |
| PySR | 3.5% / 0.19 | 9.6% / 0.21 | 17.2% / 0.51 |
| SISSO | 3.6% / 0.07 | 3.9% / 0.10 | 16.0% / 0.54 |

v2 law (alloy domain), l = L / 100 A:
`r = 0.319 (1-x) - 0.120 / l + 5.90 x(1-x) l^-1/2`.

## Reading

1. **The original comparison measured interpolation, not transfer.** On
   18 to 24 training points, five-fold row CV mostly asks an engine to
   fill gaps between seeds of a cell it has already seen. PySR and SISSO
   win that contest on the alloy set (3.5% and 3.6%), then fall to 16 to
   17% on the 111 unseen production runs. v2 reads 4.8% on row CV and
   4.5% on the unseen runs: three to four times lower error than PySR and
   SISSO where it matters for a campaign, i.e. predicting compositions and
   temperatures nobody has simulated yet.
2. **Scope before fitting.** Trained with the ballistic pure endpoints
   included, every engine does worse on the alloy hold-out (v2 10.9%,
   PySR 13.1%), because terms bend to fit rows the campaign already
   ruled out of scope. The pipeline order (screen and split, then fit) is
   what turns v2's basis into an accuracy gain.
3. **No single fix explains the gain.** F2 alone hurts (finer cells with
   the old basis overfit), F3 alone helps the hold-out, and F4 is only
   useful once selection (F1) can afford it. The four fixes act together.
4. **Materials reading of the selected law.** The size term that wins is
   `x(1-x) L^-1/2`, not `1/L`: the size effect scales with alloy disorder
   and converges slowly, as expected when long-wavelength phonons are
   only weakly (Rayleigh) scattered by mass disorder. The Phase B bulk
   extrapolation assumed 1/L; with an L^-1/2 tail the bulk alloy value
   would be higher, which bears on the Tersoff-vs-BTE gap. This is a
   hypothesis for a size series to test, not a result: the v2 law's
   infinite-size limit is not physical (it puts x = 0.15 below x = 0.5),
   so the law is valid only inside the sampled 130 to 800 A box range.
5. **The HG-SR properties are kept** (`demo_v2_props.py`): fit-twice
   determinism; replicating one design cell 5x leaves the law unchanged
   (coefficient drift 1e-15, machine precision); a +40% biased party
   added in the high-T cells manufactures a spurious `t * l^-1/2` term,
   the blast radius of that party is exactly that term, and unlearning
   the party restores the clean law bit for bit.

## Limits of this comparison

- v2's four fixes were designed after looking at the md benchmark; the
  external hold-out was added afterwards and chose nothing. It is still
  one protocol (Tersoff rNEMD, 4x4 cross-section) and one property.
- PySR and SISSO ran with the settings of the 2026-09-27 comparison (no unary operators
  for PySR, n_expansion = 2 for SISSO). Giving PySR `sqrt` would likely
  let it find the L^-1/2 term too; the point is that HG-SR can now hold
  such terms with exact revision semantics, not that PySR cannot.
- Exhaustive search over 24 terms is fast; past about 60 terms it is not
  (number of 4-term subsets grows as the fourth power).

## Next improvements, in order of expected value

1. Physics constraints on the fit (Matthiessen: every scattering term
   adds non-negative resistance; require a bulk x(1-x) term), via
   non-negative least squares. Fixes the unphysical infinite-size limit.
2. Basis proposals from PySR or SISSO, harvested automatically from their
   hall of fame into v2's library (the hand adoption of `x(1-x)/L` made
   systematic), so open search proposes and HG-SR decides.
3. Sparse search beyond exhaustive enumeration (L0 mixed-integer or beam
   search) to scale the basis past about 60 terms, keeping determinism.
4. Per-cell noise weights from seed scatter instead of uniform cell
   weights.
5. Repeat the hold-out test on a second property or system before any
   paper claims v2 generalises.
