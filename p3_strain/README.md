# P3 materials experiment: biaxially strained SiGe (2026-09-29)

Scope (2026-09-29): the Si-Ge-C
campaign was impossible (`SiCGe.tersoff` has no Ge-C parameters), so P3's
own materials experiment is **biaxial strain on SiGe**, runnable on the
validated Tersoff-1989 inventory.

**Materials question.** How much does the substrate (Si, relaxed SiGe,
Ge) change the in-plane and cross-plane lattice thermal conductivity of a
pseudomorphic SiGe film, and is the effect large enough to matter for
thin-film thermoelectrics or heat spreading?

| file | what |
|---|---|
| `in.mp_strain` | MP rNEMD deck: strain relative to the unstrained box at T, normal axis relaxed at zero stress, heat along z; `dir` 0 cross-plane, 1 in-plane |
| `GATE.md`, `gate_results/` | validation gate, PASSED (6 runs) |
| `design_strain_batch.py` -> `strain_batch_design.json`, `points_strain.txt` | 2-fold strength-2 covering seed batch (50 runs = lower bound 2 x 25), C1 (25-run 1-fold covering) embedded |
| `preregister_strain.py` | campaign record seq 32, pushed in commit `44af6e3` before submission |
| `sbatch_strain.sh` | notchpeak-guest array, `--requeue` (preemptions = the erasure channel, logged per run in `results/*.attempts`) |
| `commit_strain.py` | hashes the raw outputs and the sacct record onto the campaign record before analysis |

Seed batch: CHPC array **15890855** (50 runs, notchpeak-guest, about 15 to
40 min each). Pre-registered analysis and predictions are inside the
design JSON.
