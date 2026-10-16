#!/usr/bin/env python3
"""Apply the registered decision rule of x03_design.json to x03_results.json. Run once."""
import hashlib, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
MT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, os.path.join(MT, "pilot", "src"))
from ledger import Ledger  # noqa: E402
D = json.load(open(os.path.join(HERE, "x03_design.json")))
R = json.load(open(os.path.join(HERE, "x03_results.json")))
led = Ledger(os.path.join(MT, "phase_c", "ledger_campaign.jsonl"))
if any(e["payload"].get("op_id") == "decision-" + D["campaign"] for e in led.entries):
    sys.exit("already applied")
raw = [e["seq"] for e in led.entries if e["payload"].get("op_id") == "batch-" + D["campaign"]]
fail = [g for g, v in R["gates"].items() if not v["pass"]]
if fail:
    led.append("gate_failure", {"op_id": "gates-" + D["campaign"], "campaign": D["campaign"],
                                "failing_gates": fail, "detail": R.get("gates_fail", {}),
                                "decision": D["decision_rule"]["gates_fail"]})
b = R["branch"] if R["decision_taken"] else "intermediate"
op, target = {"density_tracks": ("correct", "recipe pick -> x = 0.5, 3 at% C, R = 2.5 nm; seq-77 scope limit lifted"),
              "insensitive": ("correct", "recipe pick -> x = 0.3, 3 at% C, R = 5 nm; seq-77 scope limit lifted"),
              "intermediate": ("retain", "seq-77 scope limit x = 0.3 to 0.5 kept")}[b]
e = led.append("shard_operation", {
    "op_id": "decision-" + D["campaign"], "campaign": D["campaign"], "op": op, "target": target,
    "cause": f"registered rule branch '{b}', rho = {R['rho']:.3f} (decision_taken={R['decision_taken']})",
    "evidence": {"x03_results.json": hashlib.sha256(open(os.path.join(HERE, "x03_results.json"), "rb").read()).hexdigest(),
                 "prereg_seq": 78, "raw_seq": raw},
    "predictions": {k: v["pass"] for k, v in R["pred"].items()},
    "not_affected": D["decision_rule"]["not_affected"], "decision": D["decision_rule"][b]})
assert led.verify_chain()
print("seq", e["seq"], op, b)
