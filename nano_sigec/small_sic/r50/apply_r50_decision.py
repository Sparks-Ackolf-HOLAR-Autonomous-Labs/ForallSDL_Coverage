#!/usr/bin/env python3
"""Apply the registered decision rule of r50_design.json to r50_results.json. Run once."""
import hashlib, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
MT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, os.path.join(MT, "pilot", "src"))
from ledger import Ledger  # noqa: E402
D = json.load(open(os.path.join(HERE, "r50_design.json")))
R = json.load(open(os.path.join(HERE, "r50_results.json")))
led = Ledger(os.path.join(MT, "phase_c", "ledger_campaign.jsonl"))
if any(e["payload"].get("op_id") == "decision-" + D["campaign"] for e in led.entries):
    sys.exit("already applied")
pre = [e["seq"] for e in led.entries if e["payload"].get("op_id") == "prereg-" + D["campaign"]]
raw = [e["seq"] for e in led.entries if e["payload"].get("op_id") == "batch-" + D["campaign"]]
fail = [g for g, v in R["gates"].items() if not v["pass"]]
if fail:
    led.append("gate_failure", {"op_id": "gates-" + D["campaign"], "campaign": D["campaign"],
                                "failing_gates": fail, "detail": R.get("gates_fail", {}),
                                "decision": D["decision_rule"]["gates_fail"]})
b = R["branch"] if R["decision_taken"] else "undecided"
cl = R.get("claims_measured_curve", {})
summ = {k: cl.get(k) for k in ("ge_saving_pct_R2.5_T300", "ge_saving_pct_R5_T300",
                               "ge_saving_pct_R2.5_T1100", "ge_saving_pct_R5_T1100")}
op, target = {"keeps_falling": ("strengthen", "P1 measured-curve claims (ledger 105); 5 nm saving and pick recomputed with the measured 5 nm point, label 'measured to 5 nm'"),
              "stops_falling": ("correct", "P1 5 nm Ge saving and screening pick -> measured-curve values with the 5 nm point; extrapolated 7-9% at 5 nm withdrawn"),
              "undecided": ("retain", "P1 ledger-105 values kept with the 'extrapolated beyond 3 nm' label")}[b]
e = led.append("shard_operation", {
    "op_id": "decision-" + D["campaign"], "campaign": D["campaign"], "op": op, "target": target,
    "cause": f"registered rule branch '{b}' (decision_taken={R['decision_taken']}); ratio_R5 = {R['ratios']['ratio_R5']}",
    "evidence": {"r50_results.json": hashlib.sha256(open(os.path.join(HERE, "r50_results.json"), "rb").read()).hexdigest(),
                 "prereg_seq": pre, "raw_seq": raw},
    "measured_curve_claims": summ, "measured_curve_pick": R.get("pick_measured_curve"),
    "predictions": {k: v["pass"] for k, v in R["pred"].items()},
    "not_affected": D["decision_rule"]["not_affected"], "decision": D["decision_rule"][b]})
assert led.verify_chain()
print("seq", e["seq"], op, b)
