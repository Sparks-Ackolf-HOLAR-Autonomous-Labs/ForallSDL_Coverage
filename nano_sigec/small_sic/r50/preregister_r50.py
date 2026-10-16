#!/usr/bin/env python3
"""Append the fixed-volume-fraction SiC radius check (R = 5 nm) pre-registration to the campaign ledger.
Idempotent on the design hash. Commit and push BEFORE submitting the job."""
import hashlib, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
MT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, os.path.join(MT, "pilot", "src"))
from ledger import Ledger  # noqa: E402

files = ["r50_design.json", "in.sic_particle", "SiCGe_GeC.tersoff", "analyze_mp.py", "points.txt",
         "sbatch_r50.sh", "analyze_r50.py", "commit_r50.py", "apply_r50_decision.py"]
sha = {f: hashlib.sha256(open(os.path.join(HERE, f), "rb").read()).hexdigest() for f in files}
d = json.load(open(os.path.join(HERE, "r50_design.json")))
led = Ledger(os.path.join(MT, "phase_c", "ledger_campaign.jsonl"))
if any(e["payload"].get("design_sha256") == sha["r50_design.json"] for e in led.entries):
    sys.exit("already registered")
e = led.append("preregistration", {
    "op_id": "prereg-" + d["campaign"], "campaign": d["campaign"], "paper": d["paper"],
    "parent_campaign": d["parent_campaign"], "design_file": "nano_sigec/small_sic/r50/r50_design.json",
    "design_sha256": sha["r50_design.json"], "files_sha256": sha,
    "gates": sorted(d["gates"]), "predictions": sorted(d["predictions"]),
    "statement": "design, inputs, predictions, decision rule and analysis script fixed before any R = 5 nm run"})
assert led.verify_chain()
open(os.path.join(HERE, "PREREGISTERED_SHA256"), "w").write(
    f"{sha['r50_design.json']}  r50_design.json  ledger_seq={e['seq']}  entry_hash={e['entry_hash']}\n")
print("seq", e["seq"], sha["r50_design.json"])
