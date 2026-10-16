# Extended basis through the full comparison (2026-10-02)

Aim: run the 43-term extended basis (squares and
three-factor products) through the same cross-validated comparison as the
adopted HG-SR v2. Script: `extended_compare.py` (adoption rule in its
docstring, committed in `26dcf40` before it ran); numbers:
`results/extended_compare.json`.

| Task | Basis | 5-fold CV | LOCO | Unseen runs (111) | kappa at L -> inf, 300 K, x = 0.15 / 0.5 / 0.85 |
|---|---|---|---|---|---|
| md (with endpoints) | standard | 4.2% | 3.7% | 10.9% | undefined (law diverges) |
| md (with endpoints) | extended | same law selected | | | |
| md-alloy | standard | 4.9% | 9.3% | 4.3% | 3.68 / 6.26 / 20.9 W/mK |
| md-alloy | extended | **3.1%** | **8.4%** | **3.5%** | 3.57 / 10.3 / 115 W/mK |

| SiGeC + dopant (MACE-MD, 16 hidden points) | log RMSE | rank corr. | lowest-kappa pick |
|---|---|---|---|
| standard | 0.063 | 0.99 | correct |
| extended | same law selected | | |

Properties of the extended engine on md-alloy: deterministic; copying one
cell five times moves no coefficient (3e-15); forgetting the biased source
restores the clean law bit for bit. **Blast radius is no longer local**: a
+40% biased high-temperature source changes all three terms of the
selected law (standard engine: exactly one term, `adopted_v2.json` C).

## Verdict: not adopted

The stated rule required all four properties. The extended basis is more
accurate on md-alloy but it lets one bad source rewrite the whole law,
which is the property P1 relies on to make corrections targeted. Its
infinite-box readout is also further from physical (115 W/mK at
x = 0.85, against an alloy value of a few W/mK).

Disclosure: the first run of the script omitted the blast-radius check
from its code and printed `adopt_extended: True`; the check was added
(as the docstring already required) and the script rerun before this
note was written. The md-alloy hold-out number (3.5%) had already been
seen in the full v2 study (`HGSR_FULL_V2.md`, 2026-10-02).

## Reading for autonomous discovery

More expressive bases buy accuracy inside the sampled design and spend
it on two things a campaign needs when it revises beliefs: correction
locality and a sane extrapolation. The extended basis stays available as
a proposal channel (like PySR), judged by the same rule.
