# C-corner validation gate for the P3 Si-Ge-C campaign

Pre-work for P3's materials experiment (the "C-corner validation gate"
listed as a dependency in `../../papers/paper3_crypto_design_outline.md`).
Run on the DGX Spark, 2026-09-29, serial LAMMPS (22 Jul 2025 update 4,
pip wheel, aarch64), one thread, 38 s wall-clock for everything.

## Verdict

| Part | Verdict |
|---|---|
| **Si-Ge-C ternary campaign as planned** | **FAIL: cannot run.** The stock `SiCGe.tersoff` has no Ge-C parameters; LAMMPS refuses the ternary before step 0 |
| Supported parts (diamond C, 3C-SiC, Si and Ge controls, substitutional C in Si, Si:C MD stability) | **PASS** against the pre-stated tolerances |

The file header states it directly: the Si(D), C and Ge entries "can be
used [for] pure silicon, pure carbon, pure germanium, binary SiC, and
binary SiGe, but not binary GeC or ternary SiGeC", because Tersoff (PRB
39, 5566, 1989) published no Ge-C mixing parameter. The planned
Si0.8Ge0.2 + 1 % C cell (`logs/in.ternary_attempt`, 512 atoms, 102 Ge,
5 C) fails with, verbatim:

```
ERROR: Potential file is missing an entry for: Si(D) Ge C (src/src/MANYBODY/pair_tersoff.cpp:554)
```

The project documents (`../../papers/00_scoping.md`,
`../../05_findings_and_lessons.md`) describe the file as "already validated"
and covering the ternary. It is validated only for SiGe (the CHPC gate in
`../../chpc/in.validate` used `Si(D) Ge`), and it does not cover Si-Ge-C.
In a dilute-C SiGe alloy every C atom has Ge neighbours with probability
x per bond, so the campaign cannot avoid Ge-C interactions. No Ge-C
parameters were constructed or used anywhere in this folder.

## Pre-stated tolerances

Fixed in `gate.py` (constants `TOL_A`, `TOL_E`, `MD_CRIT`) before any
number was computed:

- lattice constant within 1.5 % of the experimental reference;
- cohesive energy within 3 % of the experimental reference;
- MD stable: mean temperature over the last 10 ps of NVT within 5 % of
  target; total-energy drift over a following 5 ps NVE segment below
  1 meV/atom; no new C-C pairs (< 2.0 A); every C keeps 4 neighbours
  within 2.6 A; RMS displacement below 0.5 A (no melting or
  amorphisation).

References (gate criteria are the experimental column): Si 5.431 A,
-4.63 eV/atom; Ge 5.658 A, -3.85 eV/atom; diamond 3.567 A, -7.37 eV/atom;
3C-SiC 4.36 A, -6.34 eV/atom (values as given in the task brief and in
the Phase B gate). The substitutional-C slope and the MD pressure are
reported, not gated.

## Results (`results.json`)

Method: box/relax conjugate-gradient minimisation on 3x3x3 conventional
cells (216 atoms), the same style as `../../chpc/in.validate`, cross-checked
by a +-1.5 % energy-volume scan with a third-order Birch-Murnaghan fit
(the two agree to 1e-5 A in every case).

| Crystal | a (A) | dev. a | Ecoh (eV/atom) | dev. Ecoh | B0 (GPa) | Gate |
|---|---|---|---|---|---|---|
| Si (control) | 5.4320 | +0.02 % | -4.630 | 0.01 % | 97.7 | PASS |
| Ge (control) | 5.6567 | -0.02 % | -3.851 | 0.02 % | 75.7 | PASS |
| diamond C | 3.5656 | -0.04 % | -7.371 | 0.01 % | 425.7 | PASS |
| 3C-SiC | 4.3211 | -0.89 % | -6.165 | 2.77 % (underbound) | 224.1 | PASS (Ecoh margin thin) |

The Si and Ge controls reproduce the Phase B gate values (5.432 A,
-4.63 eV; 5.657 A, -3.85 eV), confirming this LAMMPS build matches the
CHPC one. Tersoff-C reproduces diamond essentially exactly (it was fitted
there). 3C-SiC reproduces Tersoff's own 4.32 A; its cohesive energy is
2.8 % too shallow, inside the 3 % tolerance with little room.

**SiC cutoff kink (own finding).** In zincblende SiC the Si-Si second
neighbours sit at a/sqrt(2); they enter the Tersoff Si-Si cutoff shell
(S = 3.0 A) below a = 4.243 A, only 1.8 % below equilibrium. The wide
E(a) scan in `results.json` (`SiC_wide_EV_scan`) shows the resulting
kink. A first E-V fit over +-4 % straddled it and gave a meaningless
bulk modulus, so the fit window was narrowed to +-1.5 %. Consequence:
anything that compresses Si-C-rich regions by about 2 % (high pressure,
C clusters, SiC precipitates) is in an artificial region of the
potential. Isolated substitutional C in Si does not reach it (Si-Si
second-neighbour distances stay near 3.8 A).

**Substitutional C in Si** (4x4x4 cell, 512 atoms, exact C counts, three
random seeds each, full box and ion relaxation):

| y (C fraction) | C atoms | a relaxed (A), mean of 3 seeds | seed std (A) | a Vegard (A) |
|---|---|---|---|---|
| 0 | 0 | 5.43201 | (1 cell) | 5.43201 |
| 0.98 % | 5 | 5.40877 | 0.00011 | 5.41378 |
| 1.95 % | 10 | 5.38579 | 0.00009 | 5.39555 |
| 3.91 % | 20 | 5.34327 | 0.00310 | 5.35910 |

Slope da/dy = -2.25 A against a Vegard slope (Tersoff's own a_C - a_Si)
of -1.87 A: **Tersoff contracts Si 1.21 times as strongly as Vegard's
rule.** The task brief quoted an experimental contraction of about 0.65
of the Vegard slope; I could not find a source for that number and do
not use it. The verified experimental reference is De Salvador et al.,
PRB 61, 13005 (2000), doi:10.1103/PhysRevB.61.13005 (co-authored by
Berti and Kelires), which reports "a strong deviation from Vegard's
rule" and a Ge:C compensation ratio of nu = 12 for SiGeC lattice-matched
to Si. Vegard's rule gives nu of about 8.3. A nu above the Vegard value
means each C contracts the lattice *more* than Vegard predicts, which is
the same direction Tersoff shows. Combining the two supported binaries
(dilute SiGe, da/dx = 0.217 A from x = 0.1 and 0.2, three seeds each)
gives an additive Tersoff estimate of nu = 10.4, between Vegard (8.3)
and experiment (12). This estimate is not a ternary result: it has no
Ge-C bonds and no Ge-C correlation, which De Salvador et al. discuss as
relevant. Qualitatively, then, the Si-C part of the potential behaves
correctly for dilute substitutional C.

**MD stability, Si + 1 % C** (5 C in 512 atoms, supported binary; the
planned Si0.8Ge0.2 + 1 % C cell cannot be built). 20 ps NVT at the 0 K
relaxed volume (Nose-Hoover, 0.1 ps damping, 0.5 fs step), then 5 ps
NVE:

| T (K) | <T> last 10 ps | T std | NVE drift (meV/atom, 5 ps) | C-C pairs start/end | C coordination | RMS disp. (A) | P (kbar) | Gate |
|---|---|---|---|---|---|---|---|---|
| 300 | 300.3 | 11.6 | +0.001 | 0 / 0 | all 4 | 0.11 | 5.5 | PASS |
| 1000 | 1000.1 | 38.7 | -0.023 | 0 / 0 | all 4 | 0.22 | 17.5 | PASS |

No clustering, no loss of fourfold C coordination, energy conserved to
well under the 1 meV/atom bound. The positive pressure is thermal
pressure from running NVT at the 0 K volume; production runs should use
NPT-equilibrated boxes.

## Known caveats of Tersoff carbon for thermal transport

The original Tersoff carbon parameters misrepresent phonon transport in
sp2 carbon: Lindsay and Broido, PRB 81, 205441 (2010),
doi:10.1103/PhysRevB.81.205441, found poor acoustic-phonon velocities
and re-optimised the parameters, obtaining considerably improved lattice
thermal conductivity for graphene and nanotubes. Whether the original
set is equally off for sp3 diamond or for C-related phonon modes in Si
was not tested here. Separately, Green-Kubo or HNEMD with any Tersoff
potential requires a heat-current expression derived for many-body
potentials (Fan et al., PRB 92, 094301 (2015),
doi:10.1103/PhysRevB.92.094301); the momentum-swap rNEMD protocol does
not use a heat-current formula and is unaffected. Finally, the dilute
C regime sits near the Si endpoint, where the project's rNEMD is already
measured to be ballistic (`../../chpc/RESULTS.md`), so dilute-C-in-Si
kappa from rNEMD would be subject to the same size limitation.

## Options for P3's third axis

**(a) A published Si-Ge-C Tersoff-type parameterisation.** OpenAlex
searches ("Tersoff potential SiGeC", "Tersoff potential germanium
carbon", "Ge-C Tersoff parameters", "germanium carbide empirical
potential") returned no paper that publishes Ge-C Tersoff parameters.
The one lead is Kelires, PRL 75, 1114 (1995),
doi:10.1103/PhysRevLett.75.1114, a Monte Carlo study of Si1-x-yGexCy
alloys with an empirical potential (the abstract reports repulsive C-C
interaction and no Ge-C bond formation). Whether that paper states
usable Ge-C parameters could not be verified from the abstract; this
needs a full-text check. Requirement if it pans out: transcribe the
parameters with provenance, then run this gate again plus a ternary
check against De Salvador 2000 (nu = 12).

**(b) Ge-C cross terms from Tersoff's mixing rules with an assumed
chi_GeC.** Mechanically easy (the Si-C and Si-Ge entries follow the rules
with chi = 0.9776 and 1.00061). It would be an unvalidated invention:
the result would depend on a free parameter chosen by us, and a "first
kappa map of the C-doped SiGe triangle" would really be a map of that
choice. Requirement: its own gate against Ge-C and SiGeC data (at
minimum the nu = 12 compensation ratio, plus a Ge-C bond length or
zincblende-GeC reference from a verified source), stated in the paper
as a new parameterisation. Not recommended for this paper.

**(c) Change the third axis to one the validated file supports.**
Candidates: dilute Si1-yCy (supported and gated here, but it is a
binary near the Si endpoint, where rNEMD is ballistic); biaxial strain
on Si1-xGex (fully supported, alloy-dominated so rNEMD stays valid over
x in 0.05 to 0.95, and strained SiGe layers are the epitaxy context that
motivated carbon in the first place); isotope mass disorder on SiGe
(zero potential changes, but its effect on kappa is largest near the
pure endpoints, again where rNEMD is invalid unless the HNEMD party is
used); vacancies (Tersoff defect energetics would need their own gate).

**Recommendation: (c), with biaxial strain on Si1-xGex as the third
axis** (design space x, strain, T). Reasons: no new parameters, it
inherits the Phase B SiGe gate unchanged, the existing rNEMD protocol
is valid across the alloy interior where the design lives, the space is
genuinely three-dimensional for the covering-array and transversal
claims, and it keeps an epitaxy-relevant materials story. Keep (a) as a
low-cost literature check (one full-text read of Kelires 1995); if it
yields published Ge-C parameters, the ternary can return as a second
campaign with its own gate. The paper text and `00_scoping.md` should
stop describing `SiCGe.tersoff` as covering the ternary.

## Files

- `gate.py`: all calculations, tolerances fixed at the top.
- `run_gate.sh`: sets the MPICH library path the DGX venv LAMMPS needs
  and runs `gate.py` (serial, `OMP_NUM_THREADS=1`).
- `results.json`: every number above, per-seed points, MD traces summary,
  the captured LAMMPS error, and the wide SiC E(a) scan.
- `logs/`: LAMMPS logs, MD temperature/energy traces (`md_*.dat`), and
  `in.ternary_attempt` (the failing ternary input).
