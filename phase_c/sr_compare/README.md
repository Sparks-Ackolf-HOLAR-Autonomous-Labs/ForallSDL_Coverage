# sr_compare — multi-engine symbolic learning + HG-SR

Aim (2026-09-27): the symbolic (meta) learning part tries different
engines and compares them (PySR /
SymbolicRegression.jl, SISSO, AI Feynman, QLattice/Feyn), then adds a
new one that is hypergraph-aware, since none of the existing ones is.

| File | Purpose |
|---|---|
| `datasets.py` | the three benchmark tasks built from the project's real evidence (literature party, CHPC-MD finite-size, combined two-party) |
| `engines.py` | uniform adapters (fit/predict/describe/complexity) with availability probing per engine |
| `hgsr.py` | **HG-SR, the new engine**: hyperedge-weighted deterministic SR with per-term provenance, scope labels, blast radius, exact unlearning |
| `run_compare.py` | one protocol for all engines: full fit, fit-twice determinism check, deterministic 5-fold CV; writes `results/results_<engine>_<task>.json` |
| `demo_bias.py` | the two demos of what hypergraph-awareness buys (oversampling invariance; biased-party localization + bitwise unlearning); writes `results/demo_bias.json` |
| `COMPARISON.md` | results, expressions, reading, and availability table |

Run everything (DGX pilot venv has pysr + TorchSisso):

```
~/mat_theory_pilot/venv/bin/python run_compare.py
~/mat_theory_pilot/venv/bin/python demo_bias.py
```

Engine roles that fell out of the comparison (see COMPARISON.md):
open-search engines (PySR) discover basis terms; deterministic
exhaustive engines (HG-SR) own the revision semantics (exact retrain,
provenance); the two compose as a pipeline stage, not competitors.
