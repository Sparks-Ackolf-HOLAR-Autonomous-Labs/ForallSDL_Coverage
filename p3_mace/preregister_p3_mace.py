#!/usr/bin/env python3
"""Register the P3 hidden-subset study before the hidden set is drawn and before any production run."""
import hashlib, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
MT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(MT, "pilot", "src"))
from ledger import Ledger  # noqa: E402
files = ["p3_mace_design.json", "design_p3_mace.py", "planner_designs.json", "chpc/points.json", "mace_rnemd.py", "reveal_hidden.py"]
sha = {f: hashlib.sha256(open(os.path.join(HERE, f), "rb").read()).hexdigest() for f in files}
d = json.load(open(os.path.join(HERE, "p3_mace_design.json")))
led = Ledger(os.path.join(MT, "phase_c", "ledger_campaign.jsonl"))
if any(e["payload"].get("design_sha256") == sha["p3_mace_design.json"] for e in led.entries):
    sys.exit("already registered")
e = led.append("preregistration", {"op_id": "prereg-" + d["campaign"], "campaign": d["campaign"], "paper": "P3",
    "design_file": "p3_mace/p3_mace_design.json", "design_sha256": sha["p3_mace_design.json"], "files_sha256": sha,
    "predictions": sorted(d["predictions"]),
    "statement": "space, planner designs, measurement protocol, surrogate, metrics and predictions fixed before the hidden set is drawn and before any production run; only two 300-step timing runs were made"})
assert led.verify_chain()
open(os.path.join(HERE, "PREREGISTERED_SHA256"), "w").write(f"{sha['p3_mace_design.json']}  p3_mace_design.json  ledger_seq={e['seq']}  entry_hash={e['entry_hash']}\n")
print("seq", e["seq"])
