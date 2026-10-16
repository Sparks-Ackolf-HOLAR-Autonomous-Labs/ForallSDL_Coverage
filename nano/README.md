# Nanostructured SiGe + light doping (P1 system, carbon-free host stage)

> **Scope note, 2026-09-29:** the P1 system is nanostructured
> **SiGeC** with light doping, not SiGe only. This folder is kept as the
> carbon-free host stage (its almaBTE per-mode data are reused); the SiGeC
> campaign, where carbon is a composition axis with substitutional and SiC
> nanoinclusion forms, is `../nano_sigec/` (campaign record seq 46 onward).

Run 2026-09-29. System choice rationale: `../papers/notes/05_p1_system_choice_and_experimental_party_2026-09-29.md`.

## The materials question

Which grain size and composition window minimises the lattice thermal
conductivity of a SiGe thermoelectric leg (the radioisotope-generator
material), how much germanium can nanostructuring save, and is light
doping (B, P, C) a useful extra thermal lever on top of alloying and grain
boundaries?

## Method

almaBTE 1.3.2 on CHPC (first-principles Si and Ge force constants from the
almaBTE database, virtual-crystal alloy, 24^3 grid; the same inputs as
evidence party 4a). `phononinfo` writes every phonon mode's heat capacity,
lifetime and group velocity; `analyze_nano.py` adds two scattering channels
mode by mode (Matthiessen's rule) and re-sums the conductivity:

- grain boundaries: gray, fully diffuse Casimir rate |v|/d for grain size d;
- dopants: Tamura point-defect rate with the Abeles mass-plus-size
  scattering strength; the size (strain) term is bracketed between
  eps = 0 (mass only) and eps = 39 (Abeles 1963), not computed.

## Order of events (campaign record `../phase_c/ledger_campaign.jsonl`)

| seq | entry |
|---|---|
| 40 | pre-registration of `nano_design.json` (grid, constants, gates G1-G3, predictions Pr1-Pr5, deliverables D1-D4), pushed in commit `d8184f9` before any production run; one exploratory smoke run (x=0.2, 300 K, bulk only) was disclosed in the entry |
| 45 | batch_result: Merkle root over all 45 raw phononinfo + 7 Tsweep hashes and the 9 npz hashes, before `analyze_nano.py` was run |

## Files

- `design_nano.py` -> `nano_design.json`; `preregister_nano.py`; `PREREGISTERED_SHA256`
- `chpc/run_nano_x.sh`, `chpc/sbatch_nano.sh`: CHPC notchpeak-guest array (one task per composition)
- `reduce_modes.py`: raw phononinfo (7 MB each, kept on CHPC home) -> committed `modes/x*.npz` + `modes/raw_sha256.json`
- `analyze_nano.py`: gates, grid, frontier, dopant lever, screening, prediction checks -> `results/`

Results: `RESULTS.md` (gates pass; Pr1 and Pr5 fail and point at partly
transparent grain boundaries; light B/P doping is not a thermal lever;
dilute C depends on the uncomputed strain term).
