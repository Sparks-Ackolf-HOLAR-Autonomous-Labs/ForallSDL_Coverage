#!/usr/bin/env python3
"""Append the fixed-volume-fraction SiC radius check (R = 3 nm) pre-registration to the campaign ledger.
Idempotent on the design hash. Commit and push BEFORE submitting the job."""
import hashlib, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
MT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, os.path.join(MT, "pilot", "src"))
from ledger import Ledger  # noqa: E402

files = ["r30_design.json", "in.sic_particle", "SiCGe_GeC.tersoff", "analyze_mp.py", "points.txt",
         "sbatch_r30.sh", "analyze_r30.py", "commit_r30.py", "apply_r30_decision.py"]
sha = {f: hashlib.sha256(open(os.path.join(HERE, f), "rb").read()).hexdigest() for f in files}
d = json.load(open(os.path.join(HERE, "r30_design.json")))
led = Ledger(os.path.join(MT, "phase_c", "ledger_campaign.jsonl"))
if any(e["payload"].get("design_sha256") == sha["r30_design.json"] for e in led.entries):
    sys.exit("already registered")
e = led.append("preregistration", {
    "op_id": "prereg-" + d["campaign"], "campaign": d["campaign"], "paper": d["paper"],
    "parent_campaign": d["parent_campaign"], "design_file": "nano_sigec/small_sic/r30/r30_design.json",
    "design_sha256": sha["r30_design.json"], "files_sha256": sha,
    "gates": sorted(d["gates"]), "predictions": sorted(d["predictions"]),
    "statement": "design, inputs, predictions, decision rule and analysis script fixed before any R = 3 nm or fixed-phi run"})
assert led.verify_chain()
open(os.path.join(HERE, "PREREGISTERED_SHA256"), "w").write(
    f"{sha['r30_design.json']}  r30_design.json  ledger_seq={e['seq']}  entry_hash={e['entry_hash']}\n")
print("seq", e["seq"], sha["r30_design.json"])
