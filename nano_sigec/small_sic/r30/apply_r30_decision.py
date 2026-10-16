#!/usr/bin/env python3
"""Apply the registered decision rule of r30_design.json to r30_results.json. Run once."""
import hashlib, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
MT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, os.path.join(MT, "pilot", "src"))
from ledger import Ledger  # noqa: E402
D = json.load(open(os.path.join(HERE, "r30_design.json")))
R = json.load(open(os.path.join(HERE, "r30_results.json")))
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
op, target = {"stops_falling": ("strengthen", "P1 Ge-saving claim -> measured-curve values; radius qualifier lifted"),
              "keeps_falling": ("correct", "P1 Ge-saving claim and preferred radius -> measured-curve values; 27-33% withdrawn"),
              "undecided": ("retain", "P1 qualified Ge-saving range (14-33%) kept")}[b]
e = led.append("shard_operation", {
    "op_id": "decision-" + D["campaign"], "campaign": D["campaign"], "op": op, "target": target,
    "cause": f"registered rule branch '{b}' (decision_taken={R['decision_taken']}); ratio_R = {R['ratios']['ratio_R']}",
    "evidence": {"r30_results.json": hashlib.sha256(open(os.path.join(HERE, "r30_results.json"), "rb").read()).hexdigest(),
                 "prereg_seq": pre, "raw_seq": raw},
    "measured_curve_claims": summ, "measured_curve_pick": R.get("pick_measured_curve"),
    "ratio_phi": R["ratios"]["ratio_phi"], "ratio_Rlo": R["ratios"]["ratio_Rlo"],
    "predictions": {k: v["pass"] for k, v in R["pred"].items()},
    "not_affected": D["decision_rule"]["not_affected"], "decision": D["decision_rule"][b]})
assert led.verify_chain()
print("seq", e["seq"], op, b)
