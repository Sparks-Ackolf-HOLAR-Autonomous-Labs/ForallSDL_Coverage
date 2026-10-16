#!/usr/bin/env python3
"""Apply the registered decision rule of light_dopants_design.json (ledger 74) to
light_dopants_results.json, and record the screening-pick inconsistency the run
exposed. Run once, after analyze_light_dopants.py. Payload keys are text only
(lesson of ledger seq 71)."""
import hashlib, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
MT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(MT, "pilot", "src"))
from ledger import Ledger  # noqa: E402
D = json.load(open(os.path.join(HERE, "light_dopants_design.json")))
R = json.load(open(os.path.join(HERE, "light_dopants_results.json")))
sha = lambda f: hashlib.sha256(open(os.path.join(HERE, f), "rb").read()).hexdigest()
led = Ledger(os.path.join(MT, "phase_c", "ledger_campaign.jsonl"))
if any(e["payload"].get("op_id") == "analysis-" + D["campaign"] for e in led.entries):
    sys.exit("already applied")
P = R["pred"]
a = led.append("analysis", {
    "op_id": "analysis-" + D["campaign"], "campaign": D["campaign"], "prereg_seq": 74,
    "results_file": "nano_sigec/light_dopants/light_dopants_results.json",
    "results_sha256": sha("light_dopants_results.json"),
    "predictions": {k: bool(v["pass"]) for k, v in P.items()},
    "classes": {k: v["class"] for k, v in R["decision"].items()},
    "statement": "registered analysis run once on committed almaBTE host modes; no new simulation"})
s = led.append("shard_operation", {
    "op_id": "decision-" + D["campaign"], "campaign": D["campaign"], "op": "strengthen",
    "target": "P1 claim 'dopants are not a thermal control variable' and the SiGeC recipe",
    "cause": ("registered rule: isolated H, D, He change kappa_L by at most 2.0% (generous strain bound), "
              "B and P by at most 0.43%; H minus D at most 0.28 percentage points; no dopant form moves "
              "x, carbon content or SiC radius of the screening pick"),
    "effect": {"lever": [k for k, v in R["decision"].items() if v["class"] == "lever"],
               "inert": [k for k, v in R["decision"].items() if v["class"] == "inert"],
               "minor": [k for k, v in R["decision"].items() if v["class"] == "minor"],
               "new_factor": "He (or H) implantation cavities enter the factor list as a microstructure lever, not a dopant"},
    "evidence": {"analysis_seq": str(a["seq"]), "results_sha256": sha("light_dopants_results.json")},
    "decision": D["decision_rule"]["recipe"]})
sp = R["screening"]
c = led.append("shard_operation", {
    "op_id": "scope-density-edge-pick-2026-10-01", "campaign": "sigec-small-sic-md-2026-09-30", "op": "scope_limit",
    "target": "density-edge screening pick 'x = 0.3, 3 at% C, R = 5 nm' (small_sic/decision_robustness.json, cited in P1)",
    "cause": ("identified: decision_robustness.py held the SiC contrast at its x = 0.2 value for every x, while the "
              "registered SiGeC model (analyze_sigec.contrast) recomputes it per x. At x = 0.3 the host density "
              "(3.23 g/cm3) matches SiC (3.21), so the density-only contrast is about 0 and SiC is invisible to the "
              "model; the registered model then puts the pick at x = 0.5, R = 2.5 nm (1.84 / 1.73 W/mK at 300 / 1100 K). "
              "The density edge was checked by MD only at x = 0.2, so neither pick is supported away from x = 0.2."),
    "effect": {"withdrawn": "x = 0.3 as the density-edge pick",
               "claims_now": "recipe composition x = 0.3 to 0.5 with 3 at% C as few-nm SiC; the Ge saving claim at the "
                             "x = 0.2 reference (27-33%) is unchanged because its matches sit at x <= 0.15",
               "registered_model_pick": {"T300": sp["T300_no_cavities"], "T1100": sp["T1100_no_cavities"]},
               "next_acquisition": "one embedded-SiC MD pair at x = 0.3 (same protocol as ledger 58) decides whether a "
                                   "density-matched particle still scatters through its stiffness contrast"},
    "evidence": {"analysis_seq": str(a["seq"])},
    "decision": "scope-limit (cause identified, but the correct value needs new evidence, so no correct operation)"})
assert led.verify_chain()
print("seq", a["seq"], s["seq"], c["seq"])
