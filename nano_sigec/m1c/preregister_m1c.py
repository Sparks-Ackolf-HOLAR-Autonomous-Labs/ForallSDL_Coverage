#!/usr/bin/env python3
"""Append the M1c converged-MACE pre-registration to the campaign ledger.
Idempotent on the design hash. Commit and push BEFORE submitting the job."""
import hashlib, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
MT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(MT, "pilot", "src"))
from ledger import Ledger  # noqa: E402

blob = open(os.path.join(HERE, "m1c_design.json"), "rb").read()
sha = hashlib.sha256(blob).hexdigest()
d = json.loads(blob)
led = Ledger(os.path.join(MT, "phase_c", "ledger_campaign.jsonl"))
if any(e["payload"].get("design_sha256") == sha for e in led.entries):
    sys.exit("already registered")
e = led.append("preregistration", {
    "op_id": "prereg-" + d["campaign"], "campaign": d["campaign"], "paper": d["paper"],
    "parent_campaign": d["parent_campaign"],
    "design_file": "nano_sigec/m1c/m1c_design.json", "design_sha256": sha,
    "gates": sorted(d["gates"]), "predictions": sorted(d["predictions"]),
    "statement": "structures, thresholds, gates, predictions and decision rule fixed before the M1c job was submitted"})
assert led.verify_chain()
open(os.path.join(HERE, "PREREGISTERED_SHA256"), "w").write(
    f"{sha}  m1c_design.json  ledger_seq={e['seq']}  entry_hash={e['entry_hash']}\n")
print("seq", e["seq"], sha)
