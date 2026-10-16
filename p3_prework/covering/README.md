# Party-regime assignment as cost-aware covering: SiGe worked example

Pre-work for P3 claim 4 (`../../papers/paper3_crypto_design_outline.md`).
Script: `covering.py` (numpy, scipy 1.18 `milp`, matplotlib; runs in under
a second). Outputs: `results.json` and `covering_grid.png`.

## Problem statement

The campaign has several evidence parties, and each one can measure some
regimes and not others at a known cost. The design question is which
party should measure which regime. We want every regime cell measured
by **at least two parties from at least two method families**, so a
cross-party agreement check (the input the shard machinery acts on) is
possible in every cell. We want the cheapest assignment that achieves
this.

Formally this is weighted set multicover with side constraints. Let C be
the regime cells and P the askable parties. Binary y[p,c] = 1 buys party p
in cell c, allowed only where p is valid. Each cell must satisfy
(parties covering c, including fixed ones) >= 2 and (distinct families
among them) >= 2. The objective is to minimise the sum of cost[p] times
y[p,c]. BTE-DFT is bought per composition run: an integer count n[col]
with at most 3 temperature bands per run. Fixed parties (literature and
DB-CITRINE) count wherever they already have data, at zero cost, and
cannot be asked for new regimes. Where no combination of parties can
meet the requirement, the demand for that cell is capped at what is
achievable, and the cell is reported as infeasible.

**Grid (35 cells).** The x columns are the hypergraph's 0.2-wide bands
(`../../pilot/kg/extract_kg.py`, `regime()`), with pure Si (x < 0.05) and
pure Ge (x > 0.95) split out because that is where the rNEMD validity
boundary lies. The temperature bands are the hypergraph's own T0 to T4
(300-500, 500-700, 700-900, 900-1100, >= 1100 K). The data come from
`../../pilot/kg/kg_hypergraph.json` (136 hyperedges). Points with unknown
x or T are dropped. Literature points labelled model-fit or unclear do
not count as coverage.

**Parties, families, costs.** Each cost is in wall-clock minutes per
point on the party's own resource, taken as the midpoint of the measured
range.

| Party | Family | Cost | Validity (scenario A) | Source |
|---|---|---|---|---|
| CHPC-MD (Tersoff rNEMD, 16 ranks) | Tersoff-MD | 23 min (6-40) | not at x < 0.05 or x > 0.95 | `../../chpc/RESULTS.md` |
| DGX-MD (Tersoff rNEMD, 4 threads) | Tersoff-MD | 12.5 min (10-15) | as CHPC | `../../pilot/PILOT_REPORT.md` |
| GPUMD-HNEMD (Tersoff, GB10) | Tersoff-MD | 20 min | everywhere | `../../gpumd/production/RESULTS.md` |
| MACE-MD (rNEMD, GB10) | MLIP-MD | 246 min (4.1 h) | as CHPC (it is rNEMD too) | `../../party4_mace/RESULTS.md` |
| BTE-DFT (almaBTE, 16 ranks) | DFT-BTE | 75 min per composition, 3 T | not at >= 1100 K | `../../party4_bte/RESULTS.md` |
| literature (fixed) | per point: experiment / other-MD / DFT-BTE | 0 | existing cells only | hypergraph |
| DB-CITRINE (fixed) | experiment | 0 | existing cells only | hypergraph |

All three Tersoff parties are one family, so pairing any two of them
does not satisfy the family requirement. MACE is treated as its own
family. Its PBE training data give it partial DFT parentage (see
`../../party4_mace/README.md`); merging it with BTE-DFT would make the
alloy cells at >= 1100 K infeasible.

**Scenarios.**

- A (base): the capabilities and midpoint costs above.
- B (fidelity-weighted): the cost of one *usable* estimate. rNEMD needs 3
  box lengths for a 1/L extrapolation, and HNEMD needs 2 driving forces
  times 2 seeds.
- C (measured GPUMD validity): GPUMD endpoint points count only at
  >= 1100 K. This is the only temperature where the pre-registered
  production campaign passed its linearity gate.
- D: BTE accepted at >= 1100 K.
- E: as A, but restricted to the three bands the campaign actually
  targeted (T0, T2, T4; 21 cells).

## Results

Rerun 2026-10-11 on the current evidence graph. The
earlier table (optimum 41.4 h) predates the hypergraph rebuilds that
added the open-database experimental values; those now cover the
Si-rich columns (Si, x0-x2) at most temperatures, which removes them
from the purchase. The numbered readings below were written for the
earlier table; the directions hold, the figures are superseded by the
table: greedy gap 4-7 %, BTE acceptance 24.5 -> 16.0 h.

| Scenario | Infeasible cells | Optimal (ILP) | Greedy | Greedy gap | Campaign cells met | Optimal for the campaign's cells | Optimal top-up from what was run |
|---|---|---|---|---|---|---|---|
| A base | pure Ge >= 1100 K | **24.5 h** | 26.2 h | 6.8 % | 24 / 35 | 8.7 h | 14.4 h |
| B fidelity | pure Ge >= 1100 K | 56.4 h | 58.7 h | 4.1 % | 24 / 35 | 17.8 h | 34.2 h |
| C measured GPUMD | 6 (pure Ge at every T, pure Si 500-700 K) | 22.8 h | 24.2 h | 5.8 % | 22 / 35 | 6.8 h | 13.4 h |
| D BTE ok at >= 1100 K | none | 16.0 h | 16.6 h | 4.2 % | 27 / 35 | 9.5 h | 6.2 h |
| E campaign T bands | pure Ge >= 1100 K | 16.4 h | 17.1 h | 4.1 % | 18 / 21 | n/a | 8.2 h |

What the campaign actually ran is priced at the scenario-A per-point
costs. It comes to about 40.3 h: CHPC 29 points (11.1 h), DGX 4 (0.8 h),
MACE 2 (8.2 h), BTE 6 compositions (7.5 h), and GPUMD 38 runs, meaning
5 validation runs plus 33 production runs (12.7 h). The CHPC figure uses
the midpoint cost, not each job's measured elapsed time.

**Optimal assignment (scenario A, figure right panel).**

- DGX-MD covers 23 alloy cells, because it is the cheapest valid party
  per point.
- BTE-DFT covers 11 composition runs, taking the second family in every
  alloy and endpoint cell below 1100 K.
- GPUMD covers the pure-Ge column and pure Si at 500-700 K and
  >= 1100 K.
- MACE-MD covers only the five alloy cells at >= 1100 K. It is bought
  there because it is the only second family available once BTE is
  excluded at high temperature.
- CHPC-MD is never chosen, because DGX is cheaper per point within the
  same family (see caveats).

Of the 24.5 h, BTE accounts for 8.75 h, MACE 8.2 h, DGX 5.2 h and
GPUMD 2.3 h.

**Readings.**

1. **Greedy is close but not optimal: a 6 to 9 % excess in every
   scenario.** On this instance, which has 35 cells and 5 askable
   parties, the exact ILP solves in milliseconds, so the heuristic is
   not needed. The gap does show that plain ratio-greedy is not free
   once there are family constraints and a batched party (BTE's
   3-temperature runs).
2. **The cost is driven by validity judgments more than by the
   solver.** Accepting BTE at >= 1100 K (scenario D) roughly halves the
   optimum, from 41.4 h to 21.8 h, because it removes the need for
   MACE. That single validity call is worth about 20 h, compared with
   2 to 4 h for the choice between greedy and exact. For the paper,
   the practical point is that party validity domains are
   first-class planner inputs and should be stated and ledgered, not
   left implicit.
3. **One cell is structurally uncoverable: pure Ge at >= 1100 K.** Only
   GPUMD is valid there, and no fixed party has data. With GPUMD's
   *measured* validity (scenario C), the whole pure-Ge column plus
   pure Si at 500-700 K become uncoverable. The formulation reports
   these gaps explicitly instead of silently skipping them, and they
   point directly at the next party to recruit (an experimental high-T
   Ge source, or BTE with a high-T validity argument).
4. **What the campaign missed.** Within its own temperature bands
   (scenario E), the campaign met the two-party, two-family rule in 9
   of 21 cells. The gaps are:
   - the whole x = 0.6-0.8 column (never measured by anyone);
   - the Si-rich and x = 0.5 alloys at 700 K (BTE only, no MD
     cross-check);
   - Ge-rich x = 0.8-0.95 at 300 and 700 K (BTE only);
   - every alloy cell at >= 1100 K (two Tersoff parties, or BTE, which
     does not count there);
   - pure Ge at >= 1100 K.

   Six of its points fall outside the rule's validity domains (CHPC
   rNEMD at the pure endpoints), and so do its six BTE points at
   1100 K under scenario A. The optimal top-up to complete the
   campaign's bands costs 23.8 h. MACE accounts for 20.5 h of that, so
   the cheap part of the fix is only about 3.3 h: 10 DGX points and one
   BTE run at x = 0.6-0.8.
5. **What the campaign "would have saved"** needs careful reading. The
   same 9 cells were achievable for 4.3 h (A) or 9.0 h (B), against
   about 40 h spent. Most of that spend did not go to coverage. It
   bought seed statistics, 1/L extrapolation lengths, the cross-machine
   reproducibility check, and P2's pre-registered endpoint campaign,
   which needs replicate runs by design. The fair statement is that
   **coverage and convergence are different budgets**. A covering
   planner allocates the first, and the campaign spent mostly on the
   second while leaving 12 of its own 21 cells without a two-family
   check.

## Strength of the evidence

This is a worked example, not strong evidence. The instance is small
enough that exact optimisation is trivial and the greedy gap (6-9 %) is
modest. Several choices move the numbers substantially:

- the costs mix CPU wall-clock (CHPC, BTE, DGX) with GPU wall-clock
  (GPUMD, MACE) in one unit, with no exchange rate;
- CHPC guest compute is free and runs in parallel, so a money or
  throughput cost model would favour CHPC over DGX, and the "CHPC never
  chosen" result is an artefact of the wall-clock unit;
- the cell grid and the 2-party, 2-family rule are design choices;
- the validity domains are judgment calls, and scenarios C and D show
  how much they matter.

What the example does demonstrate is qualitative and robust across all
five scenarios:

- the formulation is well posed on real campaign data;
- it exposes structurally uncoverable cells;
- it separates coverage spend from convergence spend;
- it shows validity-domain declarations dominating cost.

A stronger P3 result would need a larger instance (the ternary or
strain campaign, several method families per cell), a principled cost
unit, and a prospective test: plan with the ILP, run the plan, and
compare against an unplanned baseline.

## Files

- `covering.py`: data loading, ILP (`scipy.optimize.milp`), greedy,
  campaign comparison, figure.
- `results.json`: fixed-party coverage per cell, and per scenario the
  optimal and greedy assignments and costs, infeasible cells, the
  campaign's satisfied and missed cells, invalid points, and the top-up
  plan.
- `covering_grid.png`: left panel, what the campaign ran (valid
  measurements only); right panel, the scenario-A optimum. Pink marks
  cells where the requirement is not met.

## Reviewer note

The "what the campaign ran" panel reads coverage from
`pilot/kg/kg_hypergraph.json`, which has not been rebuilt since the
Phase C production campaign (132 CHPC runs over the 40-point coverage
batch, campaign record seq 3, 9-13) or the Ge-corner shard action (19 runs,
seq 19, 26). It therefore undercounts real coverage, most visibly the
x 0.6 to 0.8 column and the mid-temperature alloy cells, which the
production batch did sample. Its GPUMD pure-Si cell at 300-500 K also
predates the endpoint campaign's failed 300 K linearity gate (seq 18).
Rebuild the hypergraph with those rows before quoting the
campaign-versus-optimum comparison; the optimal-assignment panel and the
scenario table do not depend on it.
