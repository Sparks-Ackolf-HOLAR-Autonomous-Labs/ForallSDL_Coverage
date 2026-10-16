# Light dopants H, D, He, B, P in nanostructured SiGeC (2026-10-01)

Scope (2026-10-01): the dopant set is H, D (deuterium), He, B and P, not only B and P.

## Order of steps
1. `design_light_dopants.py` wrote `light_dopants_design.json`; `preregister_light_dopants.py`
   put its hash on the campaign record (seq 74). Pushed (`f3ac6e3`) before any dopant number existed.
2. `analyze_light_dopants.py` ran once on committed evidence only (almaBTE host modes, campaign-record entries 45;
   converged MACE summary, campaign-record entries 54-55). Seconds on one CPU core, no new simulation.
3. `apply_light_dopants_decision.py` recorded the analysis (seq 75), the registered decision
   (seq 76, strengthen) and a scope-limit on the earlier density-edge recipe pick (seq 77).

## How each element enters
| Element | Physical form in sintered SiGeC | Model |
|---|---|---|
| B, P | substitutional, 1e19 to 3e20 cm^-3 | Tamura mass + size term (unchanged from the nano campaign) |
| H, D | isolated interstitial (bond-centre) | added mass at a site; upper bound adds a generous local-strain term |
| He | isolated interstitial at low dose | same as H/D |
| He (or H) after implantation and anneal | nanocavities, R 2.5 to 25 nm, 0.1 to 3 vol% | void inclusions (same cross-section as SiC, contrast 1 to 4) plus the Maxwell-Eucken porosity factor |

## Results (reference: Si0.8Ge0.2, 0.3 um effective grains)
| Prediction | Result | Verdict |
|---|---|---|
| PrL1 isolated H, D (mass only) < 1% | at most 0.002% | held |
| PrL2 isolated H, D, He < 5% at the strain bound | at most 2.0% | held |
| PrL3 H minus D < 0.5 points | at most 0.28 points | held |
| PrL4 cavities beat SiC at 1 vol%, 5 nm | 19.4% vs 14.7% (300 K); 17.2% vs 12.6% (1100 K) | held |
| PrL5 cavities save >= 25% Ge | 40% (300 K), 38% (1100 K); SiC reference 33% / 30% | held |
| PrL6 SiC + cavities sub-additive | 27.0% together vs 34.0% summed (300 K) | held |
| PrL7 B, P on top of SiC < 1% | at most 0.43% | held |

Classes under the registered rule (>= 5% lever, < 1% inert): B, P inert; isolated H, D, He minor
(only under the generous strain bound); cavities a lever. No dopant form moves the screening pick.

Cavity numbers barely depend on the contrast bracket (lam2 1 vs 4 moves them by < 1 point at
R >= 5 nm), because a void a few nanometres across already scatters near its geometric limit.
About 1.5 points of the 1 vol% reduction is the porosity factor; the rest is phonon scattering.

## The flaw the run exposed (campaign record seq 77)
`small_sic/decision_robustness.py` held the SiC contrast at its x = 0.2 value for all x. The
registered model recomputes it per composition. Near x = 0.3 the host density (3.23 g/cm3)
equals that of SiC (3.21), so a density-only contrast is about zero and SiC does nothing in the
model. Recomputed per x, the density-edge pick is x = 0.5, 3 at% C, R = 2.5 nm (1.84 / 1.73 W/mK
at 300 / 1100 K), not x = 0.3, R = 5 nm. The MD check (campaign-record entries 58-61) calibrated the contrast only
at x = 0.2, so neither pick is supported away from x = 0.2: the pick is scope-limited to
x = 0.3 to 0.5. The Ge-saving claim (27-33%) is unaffected, because its matches sit at x <= 0.15.
Next acquisition: one embedded-SiC MD pair at x = 0.3 (same protocol as campaign-record entries 58).

## Not modelled (stated in the design before the run)
Electrical effects (B-H passivation; carrier scattering by cavities), which decide zT;
whether H, D, He stay in the solid at service temperature; cavity coarsening;
Si-H / Si-D local vibrational modes (far above the Si phonon band, no resonance with
heat-carrying phonons).
