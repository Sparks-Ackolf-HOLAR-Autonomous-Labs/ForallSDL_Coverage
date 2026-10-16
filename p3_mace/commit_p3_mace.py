#!/usr/bin/env python3
"""Ledger the P3 hidden-subset MACE-MD outputs BEFORE analyze_p3_mace.py runs: one
batch_result entry with the Merkle root over every raw/out/<id>.json. Refuses if any of the
89 registered points is missing or a file would be dropped by .gitignore."""
import hashlib, json, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
MT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(MT, "pilot", "src"))
from ledger import Ledger, leaf_hash, merkle_root  # noqa: E402

D = json.load(open(os.path.join(HERE, "p3_mace_design.json")))
OUT = os.path.join(HERE, "raw", "out")
P = json.load(open(os.path.join(HERE, "chpc", "points.json")))["points"]
miss = [p for p in P if not os.path.exists(os.path.join(OUT, p + ".json"))]
if miss:
    sys.exit(f"{len(miss)} registered points have no output: {miss[:5]}")
files = {}
for f in sorted(os.listdir(OUT)):
    p = os.path.join(OUT, f)
    if subprocess.run(["git", "check-ignore", "-q", p]).returncode == 0:
        sys.exit(f"raw/out/{f} would be dropped by .gitignore")
    files[f"raw/out/{f}"] = hashlib.sha256(open(p, "rb").read()).hexdigest()
leaves = [leaf_hash({"file": k, "sha256": v}) for k, v in sorted(files.items())]
led = Ledger(os.path.join(MT, "phase_c", "ledger_campaign.jsonl"))
op = "batch-" + D["campaign"]
if any(e["payload"].get("op_id") == op for e in led.entries):
    sys.exit("already committed")
e = led.append("batch_result", {
    "op_id": op, "campaign": D["campaign"], "stage": "all 89 runs", "raw_merkle_root": merkle_root(leaves),
    "n_files": len(files), "files_sha256": files,
    "chpc_job": "15922934 (70 runs); 15950354 (19 runs resubmitted after an uncorrectable ECC GPU error "
                "on notch458 killed them at start-up; same script, points and seeds)",
    "statement": "committed before analyze_p3_mace.py was run on these files"})
assert led.verify_chain()
print("seq", e["seq"])
