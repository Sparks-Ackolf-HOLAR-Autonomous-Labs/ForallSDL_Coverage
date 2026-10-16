# Covering and loss-tolerant batch designs guarantee what autonomous campaigns sample in SiGeC semiconductors

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23413465.svg)](https://doi.org/10.5281/zenodo.23413465)

Code, data and figure scripts for the paper "Covering and loss-tolerant batch designs guarantee what autonomous campaigns sample in SiGeC semiconductors" (Molecular Systems Design & Engineering (RSC)). The manuscript itself is not part of this folder.

Licences: code (Python, shell and SLURM scripts) under the MIT License (`LICENSE`); data, result files, figures and text under CC BY 4.0 (`LICENSE-CC-BY-4.0.txt`). Archived third-party datasets keep their own terms and are cited in the paper: the Starrydata slice is CC BY 4.0 (Katsura et al., Sci. Technol. Adv. Mater.: Methods 5, 2506976 (2025), doi:10.1080/27660400.2025.2506976; cite it and the source paper of each curve used), and the Citrine release keeps its published terms.

Every path below is relative to this folder, and every script runs from it. The folder depends on no other folder or repository.

Analysis, figure and checker scripts run on any machine with the Python packages in `requirements.txt`. Acquisition scripts (SLURM job files and run generators) name the cluster paths and simulation binaries (LAMMPS, GPUMD, almaBTE, MACE) used at the time, so re-running a simulation needs those paths edited to the local installation. The recorded outputs of those runs are included, so no simulation has to be repeated to check a result.

Files whose SHA-256 the campaign record stores (designs, pre-registration digests, the scripts and raw outputs they name, and some notes) are kept byte-for-byte as they ran, because any edit would break verification. They keep their original wording, which includes operator decisions with dates, internal step numbers and the run paths of the time. Raw simulation logs are also kept as written.

Terms used in folder names and notes: a *party* is one evidence source of the campaign (for example Tersoff MD, MACE-MP rNEMD, almaBTE first-principles transport or a curated database); a *gate* is the acceptance check a run must pass before its result enters the campaign record; `phase_c/` holds the main campaign and `pilot/` the campaign-record library and the first trial campaign.

## Campaign material

| path | what it holds |
|---|---|
| `p3_sigec` | SiGeC + dopant design study (pre-registered, campaign-record entries 56-57) |
| `p3_mace` | hidden MACE-MD subset the planner could not see |
| `p3_strain` | biaxial-strain covering campaign (Tersoff rNEMD on CHPC) |
| `p3_prework` | carbon-corner gate and party-assignment covering example |
| `p3_limits` | run-cost, strength-three and recorded-loss studies behind the limitations |
| `p3_reserve` | replication reserve study (pre-registered; ESI Section S7) |
| `scalability` | scaling check of covering-batch size with the number of factors (ESI Section S9) |
| `nano` | grain-boundary calibration of the response model (carbon-free nanostructured SiGe) |
| `nano_sigec/m1c` | converged MACE carbon-site energies reused by the response model |
| `nano_sigec/small_sic` | embedded-SiC molecular dynamics behind the SiC contrast |
| `nano_sigec/radius_contrast` | radius-dependent SiC contrast |
| `nano_sigec/light_dopants` | light dopants and cavities |
| `phase_c/sr_compare` | resistivity-law fit and symbolic-regression engines |
| `phase_c/ledger_campaign.jsonl` | the campaign record |
| `pilot/src` | campaign-record library |

Files outside these folders that a script reads (traced):

- `nano_sigec/PREREGISTERED_SHA256`
- `nano_sigec/analyze_sigec.py`
- `nano_sigec/sigec_design.json`

## Figures

Each figure is drawn by one script at its print size and written to `paper/figure_code/figures/` (plus a 600 dpi PNG). Run from this folder with a Python that has numpy, scipy, matplotlib, scikit-learn and ase:

| figure | command |
|---|---|
| `fig_assignment` | `python3 paper/figure_code/regen_figs/fig_assignment.py` |
| `fig_hidden_mace` | `python3 paper/figure_code/regen_figs/fig_hidden_mace.py` |
| `fig_ranking` | `python3 paper/figure_code/regen_figs/fig_ranking.py` |
| `fig_reserve` | `python3 paper/figure_code/regen_figs/fig_reserve.py` |
| `fig_scaling` | `python3 paper/figure_code/regen_figs/fig_scaling.py` |
| `fig_sigec_design` | `python3 paper/figure_code/regen_figs/fig_sigec_design.py` |
| `fig_slice` | `python3 paper/figure_code/regen_figs/fig_slice.py` |
| `fig_strain` | `python3 paper/figure_code/regen_figs/fig_strain.py` |
| `fig_transversal` | `python3 paper/figure_code/regen_figs/fig_transversal.py` |

## Checkers and analyses

| command |
|---|
| `python3 p3_sigec/run_p3_sigec.py ` |
| `python3 p3_sigec/main_effects_measured_contrast.py ` |
| `python3 p3_mace/analyze_p3_mace.py ` |
| `python3 p3_strain/analyze_strain.py ` |

`traced_inputs.json` lists every file each command opened when the folder was built; `FILES.sha256` lists the digest of every file in the folder.
