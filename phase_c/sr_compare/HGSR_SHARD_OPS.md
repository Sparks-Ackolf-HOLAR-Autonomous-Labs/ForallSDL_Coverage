# How HG-SR does the SISA-like process, and what it does besides forget

Question addressed (2026-09-26): how does the hypergraph unlearning
part carry out a SISA-like process, and which operations besides
unlearning can it perform?

Code: `shard_ops.py`. Executable evidence: `demo_shard_ops.py` →
`results/demo_shard_ops.json` (run on the real md task evidence plus the
demo_bias biased party).

## The SISA correspondence

SISA (Bourtoule et al., IEEE S&P 2021) makes unlearning cheap through
structure: the training set is cut into disjoint shards, one isolated
sub-model is trained per shard, predictions are aggregated, and deleting
a point retrains only the shard that contained it. Two properties do the
work there: **isolation** (each shard influences only its own sub-model)
and **bounded recompute** (one sub-model per deletion).

HG-SR keeps both properties but relocates them, because our cost
structure is inverted: SISA's model is expensive to retrain and cheap to
aggregate; our model (exhaustive deterministic OLS over a fixed basis)
retrains in seconds, so there is nothing to save on compute. What still
needs bounding is *belief change* — which parts of the fitted law an
evidence group is allowed to move, and what falls if that group goes.

| SISA | HG-SR |
|---|---|
| shards = disjoint random partition of rows | hyperedges = semantic, overlapping groups (every row is in a regime edge AND a party edge AND a source edge) |
| isolation by architecture: shard i only touches sub-model i | isolation by loss: hyperedge-balanced weights cap every edge's vote at 1/(edge size) per row, so a big or oversampled shard cannot outvote the rest (demo_bias demo 1: replication drift exactly 0.0) |
| influence map is trivial: point → its shard's sub-model | influence map is computed: leave-one-hyperedge-out provenance names which terms each edge supports; inverting it gives the blast radius B(e) |
| unlearning = retrain one sub-model, aggregate again | any operation = one clean refit of the single global law; determinism makes it epsilon = delta = 0 exact, verified bitwise |
| aggregation costs accuracy (ensemble of weak sub-models) | no aggregation: one law, full data, the balancing lives in the weights |

So the one-line answer: HG-SR is SISA with the shard map made
semantic and overlapping, the isolation moved from compute to influence,
and the retrain-one-shard trick made unnecessary because the clean
retrain — the gold standard SISA approximates — is our native operator.

## Beyond forget: the operator family as row/weight primitives

Every operator in U reduces to one of four primitives on the evidence
campaign record — remove rows, reweight rows, edit rows, re-tag rows — followed by
the same deterministic refit and clean-retrain verification. Measured
results from `demo_shard_ops.py` (clean law: r = 3.115x − 3.193x² +
65.68/L; tainted with a +40% biased party in the thin high-T regime, the
law grows a spurious −0.00031·T term):

| Operator | Primitive | Demo result |
|---|---|---|
| **forget** | remove the edge's rows | restores the clean model bitwise (atol=0) |
| **quarantine** | weight 0, rows kept | model identical to forget's, but **restore** brings back the tainted fit bitwise — nothing was destroyed, so a wrong quarantine costs nothing |
| **downweight** | weight gamma in (0,1) | a continuous trust dial: the spurious T term persists at gamma 1.0 and 0.5, disappears from gamma 0.2 down — retain (1) and quarantine (0) are its endpoints |
| **correct** | edit the rows (kappa/1.4 recalibration), keep them | terms return to the clean law and the corrected rows now *add* evidence (coefficients tighten slightly vs clean) |
| **split** | re-tag rows into finer edges (biased rows become their own sub-regime) | the balanced loss caps the biased edge's vote: the T coefficient shrinks (−0.00031 → −0.00026) but does **not** vanish — split contains bias, it does not remove it; the right split use is real heterogeneity (literature bulk/film/nano), not bias handling |
| **replace** | forget + insert replacement rows, one transaction | corrected re-measurement lands, terms match clean |
| **scope-limit** | no refit: restrict validity to covered regimes / attach party scope label | `hgsr.provenance` already emits the labels; the campaign record stores the restriction |
| **retain / strengthen** | no refit: campaign record status change | witnessed by the provenance labels — a term moving from single-edge provenance to multi-supported after new evidence lands IS the strengthen event |

Two design points worth keeping:

1. **Reversibility ordering.** quarantine/downweight are weight
   operations (reversible, evidence intact); correct/split/replace are
   campaign record edits with the old rows retained in history; forget is the
   only destructive operator. That matches the project's two-bar rule: a
   statistical verdict suffices to suspend influence (quarantine), but
   removal needs an identified cause or independent confirmation.
2. **Split is not a forgetting substitute.** The demo shows split
   containing a biased edge's influence without deleting any belief —
   the correct response when the anomaly is real heterogeneity
   (world-is-richer), and measurably insufficient when the rows are
   actually wrong. The triage between split and forget/correct is
   exactly the expand-vs-unlearn decision the design doc assigns to
   stage 4.
