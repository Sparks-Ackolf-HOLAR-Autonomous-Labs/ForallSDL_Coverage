#!/usr/bin/env python3
"""Append the P3 SiGeC design-study pre-registration to the campaign ledger.
Idempotent on the design hash. Commit and push before run_p3_sigec.py runs."""
import hashlib, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
MT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(MT, "pilot", "src"))
from ledger import Ledger  # noqa: E402

LEDGER = os.path.join(MT, "phase_c", "ledger_campaign.jsonl")
blob = open(os.path.join(HERE, "p3_sigec_design.json"), "rb").read()
sha = hashlib.sha256(blob).hexdigest()
d = json.loads(blob)
led = Ledger(LEDGER)
if any(e["payload"].get("design_sha256") == sha for e in led.entries):
    sys.exit("already registered")
e = led.append("preregistration", {
    "op_id": "prereg-" + d["campaign"], "campaign": d["campaign"], "paper": d["paper"],
    "design_file": "p3_sigec/p3_sigec_design.json", "design_sha256": sha,
    "predictions": sorted(d["predictions"]),
    "statement": ("factors, designs, surrogate, metrics and predictions fixed before any design was "
                  "built or evaluated on the SiGeC response model; no exploratory run exists"),
})
assert led.verify_chain()
open(os.path.join(HERE, "PREREGISTERED_SHA256"), "w").write(
    f"{sha}  p3_sigec_design.json  ledger_seq={e['seq']}  entry_hash={e['entry_hash']}\n")
print("seq", e["seq"], sha)
