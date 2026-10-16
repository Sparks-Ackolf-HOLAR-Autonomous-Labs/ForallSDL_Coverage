#!/usr/bin/env python3
"""Append the nano campaign pre-registration to the campaign ledger.
Idempotent on the design hash. Commit and push before any production run."""
import hashlib, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
MT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(MT, "pilot", "src"))
from ledger import Ledger  # noqa: E402

LEDGER = os.path.join(MT, "phase_c", "ledger_campaign.jsonl")
blob = open(os.path.join(HERE, "nano_design.json"), "rb").read()
sha = hashlib.sha256(blob).hexdigest()
d = json.loads(blob)
led = Ledger(LEDGER)
if any(e["payload"].get("design_sha256") == sha for e in led.entries):
    sys.exit("already registered")
e = led.append("preregistration", {
    "op_id": "prereg-" + d["campaign"], "campaign": d["campaign"], "paper": d["paper"],
    "design_file": "nano/nano_design.json", "design_sha256": sha,
    "gates": sorted(d["gates_N0"]), "predictions": sorted(d["predictions"]),
    "statement": ("grid, constants, gates, predictions and deliverables fixed before any "
                  "production phononinfo run; the only prior output is one exploratory "
                  "smoke run (x=0.2, 300 K, bulk) used to confirm the phononinfo column "
                  "format and that the per-mode sum reproduces kappa_Tsweep RTA (7.3780 "
                  "vs 7.37796 W/mK); no boundary or dopant term was evaluated"),
})
assert led.verify_chain()
open(os.path.join(HERE, "PREREGISTERED_SHA256"), "w").write(
    f"{sha}  nano_design.json  ledger_seq={e['seq']}  entry_hash={e['entry_hash']}\n")
print("seq", e["seq"], sha)
