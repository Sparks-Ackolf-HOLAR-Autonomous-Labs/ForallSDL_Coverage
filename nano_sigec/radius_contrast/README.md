# Radius-dependent effective SiC contrast (post hoc, 2026-10-04)

**Status: post hoc.** Not pre-registered, no campaign-record entry, no new MD. The
analysis reuses committed results only and was designed after those results
were known. Three MD radii constrain at most two parameters; every
two-parameter form below is fitted with one degree of freedom (x = 0.2) or
none (x = 0.3, two radii), so the fits describe the data, they do not test a
mechanism.

## Question

P1 sets the small-SiC scattering strength (effective contrast lambda^2 in the
Rayleigh-to-geometric cross-section) from Tersoff MD at R = 1.0, 1.5 and
2.0 nm and applies one value, the density edge lambda^2 = 9.3e-3, at the
claimed radii 2.5 and 5 nm. The bond-angle analysis (`../../bond_angles`)
shows that a nominal R = 2 nm particle is crystalline only to 1.45 nm. Does
any plausible radius dependence of lambda^2 change the recipe
(x = 0.3, 3 at% C, R = 5 nm) or the 27 to 33% Ge saving?

## Method (`fit_radius_contrast.py`, about 2 min on 2 cores)

1. Each MD reduction r = 1 - kappa(particle)/kappa(pocket-only control),
   per seed, is inverted to an effective lambda^2 with the model side of
   `../small_sic/analyze_md.py` (Tersoff-scaled almaBTE modes, MD volume
   fraction, boundary length d matched to the MD box). Two frames:
   *registered* (d matched to the pocket-only control at each radius, as in
   campaign record seq 59 to 61) and *host-matched* (d matched once to the bare host,
   as in `../small_sic/sensitivity_md.py`). Inputs: x = 0.2 at R = 1.0, 1.5,
   2.0 nm; x = 0.3 at R = 1.5, 2.0 nm (`../small_sic/x03/results`).
2. Five forms for lambda^2(R) are fitted by least squares in r over all seeds,
   and refitted for every seed recombination (2^3 = 8 at x = 0.2):
   - constant (the current assumption, refitted to MD);
   - (a) lambda^2 = L f_core, f_core = ((R - s)/R)^3, and (a2) L f_core^2,
     with interface-shell thickness s = 0.55 nm fixed from the bond angles
     (angle s.d. below 0.1 deg to r = 1.45 nm at nominal R = 2 nm);
   - (b) core-shell: lambda^2 = (a_c f_core + a_s (1 - f_core))^2, the
     volume-averaged Rayleigh amplitude with a signed shell contrast;
   - (c) power law: lambda^2 = A (R / 2 nm)^p.
   s = 0.45 and 0.65 nm are run as a sensitivity check.
3. lambda^2(2.5 nm), lambda^2(5 nm) and lambda^2(R) on the full radius grid
   are passed through the decision analysis of
   `../small_sic/decision_robustness.py` (unscaled almaBTE modes): Ge saving at
   equal kappa_L (x = 0.2, 0.3 um grains, 1 at% C, 300 and 1100 K) and the
   screening pick (argmin kappa_L over x, C content, R).

Outputs: `radius_contrast_results.json` (all numbers), `fit_log.txt`,
`fig_radius_contrast.{pdf,png}` (3.3 in wide, 8 pt; `make_figure.py`).

## Findings

**1. The effective contrast falls with radius, the opposite of the
core-fraction picture.** Per-seed lambda^2:

| R (nm) | x = 0.2, registered | x = 0.2, host-matched | x = 0.3, registered |
|---|---|---|---|
| 1.0 | 0.035, 0.136 | 0.039, 0.165 | (not run) |
| 1.5 | 0.016, 0.021 | 0.023, 0.032 | 0.010, 0.019 |
| 2.0 | 0.0046, 0.0056 | 0.0099, 0.013 | 0.0037, 0.0049 |

A crystalline-core fraction that grows with R predicts the reverse trend.

**2. Only forms that fall with R fit at the seed-scatter level.** RMS
residual in r (x = 0.2, registered; seed floor 0.023): constant 0.073,
(a) 0.086, (a2) 0.095, (b) 0.025, (c) 0.024. Host-matched: 0.052, 0.067,
0.079, 0.025, 0.023. The power-law exponent is p = -4.1 (registered) and
-3.0 (host-matched); x = 0.3 gives -4.2 and -2.9. With p = -4 the Rayleigh
cross-section scales as R^2, i.e. with particle surface rather than volume,
which is what scattering by the disordered interface band would give. The
fitted core-shell form has core and shell amplitudes of opposite sign
(+0.27 and -0.28), so its lambda^2 passes through zero near 2.7 nm; that
cancellation is a property of the form, not an observation. Forms (b) and
(c) fit equally well and disagree beyond 2 nm: three radii cannot choose
between them.

**3. Ge saving at equal kappa_L (x = 0.2 calibration; 300 K / 1100 K;
brackets: range over seed recombinations at 300 K).**

| Form | frame | R = 2.5 nm | R = 5 nm |
|---|---|---|---|
| P1 (density edge) | | 30 / 27 | 33 / 30 |
| constant | registered | 28 / 25 [27, 29] | 32 / 29 [31, 33] |
| (a) core fraction | registered | 29 / 26 [28, 30] | 34 / 31 [34, 34] |
| (b) core-shell | registered | 5 / 4 [4, 44] | 34 / 31 [29, 39] |
| (c) power law | registered | 16 / 14 [13, 19] | 8 / 7 [4, 14] |
| constant | host-matched | 35 / 31 | 35 / 32 |
| (b) core-shell | host-matched | 20 / 17 | 31 / 28 |
| (c) power law | host-matched | 27 / 23 | 20 / 17 |

At x = 0.2 the shell thickness (0.45 to 0.65 nm) moves forms (a) and (b)
by at most 3 points at 5 nm. The x = 0.3 calibration gives the same pattern
(power law, registered: 14 / 12 at 2.5 nm, 7 / 6 at 5 nm).

**4. Recipe pick.** Composition and carbon content are unchanged: x = 0.3
with 3 at% C is the 300 K pick in all 20 fits. The radius is not: every
power-law fit prefers 2.5 nm, and at 1100 K the power-law fits pick x = 0.5
at 2.5 nm, ahead of x = 0.3 by 1.2% in kappa_L (2.11 against 2.14 W/mK).
x = 0.3 and x = 0.5 were already within 0.3% at 1100 K in the current P1
numbers, so this is a near tie rather than a reversal. Under the power law
the P1 recipe at 5 nm conducts 8% more than the same composition at 2.5 nm
(2.48 against 2.30 W/mK at 300 K).

## Reading

- The 27 to 33% Ge saving and the 5 nm radius hold under a constant or a
  core-fraction contrast, but those forms fit the MD poorly. The two forms
  that fit the MD (core-shell, power law) both make lambda^2 fall through
  2.5 nm; under the power law the saving at 5 nm drops to 7 to 8%, and the
  best of the two radii (2.5 nm) gives 14 to 16% (registered) or 23 to 27%
  (host-matched). A radius-robust statement is therefore "at least about
  14% at 2.5 nm, up to about 33%", with the radius choice open.
- Carbon replacing germanium remains a positive saving under every form;
  only its size and the preferred particle radius depend on the
  extrapolation.
- Confounders: in the MD design the volume fraction rises with R (1.3%,
  4.2%, 10%), so a fall of lambda^2 with R cannot be separated from a fall
  with volume fraction (particle interaction or overlap of the Si-rich
  pockets); the pocket width (0.7 nm) was imposed, not relaxed; one R = 1 nm
  run failed gate G1 and the two R = 1 nm seeds differ by a factor of 2.5;
  the shell thickness is from one 0 K minimisation that stopped at its
  iteration cap. A fourth radius (R = 3 nm) at fixed volume fraction would
  separate the forms; that is not run here.
