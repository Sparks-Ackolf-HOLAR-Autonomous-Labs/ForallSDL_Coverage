#!/usr/bin/env python3
"""Apply the registered decision rule of md_design.json to md_results.json:
one gate_failure entry for the failing gates and one shard_operation entry
for the branch the rule selects. Run once, after analyze_md.py."""
import hashlib, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
MT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(MT, "pilot", "src"))
from ledger import Ledger  # noqa: E402
D = json.load(open(os.path.join(HERE, "md_design.json")))
R = json.load(open(os.path.join(HERE, "md_results.json")))
sha = lambda f: hashlib.sha256(open(os.path.join(HERE, f), "rb").read()).hexdigest()
led = Ledger(os.path.join(MT, "phase_c", "ledger_campaign.jsonl"))
if any(e["payload"].get("op_id") == "decision-" + D["campaign"] for e in led.entries):
    sys.exit("already applied")
fail = [g for g, v in R["gates"].items() if not v["pass"]]
if fail:
    led.append("gate_failure", {
        "op_id": "gates-" + D["campaign"], "campaign": D["campaign"], "failing_gates": fail,
        "passing": [g for g in R["gates"] if g not in fail], "analysis_file": "nano_sigec/small_sic/md_results.json",
        "analysis_sha256": sha("md_results.json"),
        "detail": {"G1": R["gates"]["G1"]["rows"]} if "G1" in fail else {},
        "decision": "reported as registered; the design names no consequence for a gate failure, so the "
                    "decision rule is still applied and the failing run is named beside it"})
dec = R["decision"]
op = {"inside_band": "strengthen", "above_cap": "correct", "below_density": "downweight"}[dec["branch"]]
e = led.append("shard_operation", {
    "op_id": "decision-" + D["campaign"], "campaign": D["campaign"], "op": op,
    "target": "P1 few-nm SiC claims (reduction and Ge saving at R = 2.5 and 5 nm)",
    "cause": f"registered rule branch '{dec['branch']}': classifications {dec['classifications']}",
    "effect": {"claims_now": dec["claims"], "prior_band": dec["prior_band"]},
    "evidence": {"md_results.json": sha("md_results.json"), "prereg_seq": 58, "raw_seq": 59},
    "predictions": {k: v["pass"] for k, v in R["pred"].items()},
    "decision": D["decision_rule"][dec["branch"]]})
assert led.verify_chain()
print("seq", e["seq"], op)
