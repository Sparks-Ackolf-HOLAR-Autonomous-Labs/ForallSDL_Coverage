#!/usr/bin/env python3
"""Ledger the small-SiC MD outputs BEFORE analyze_md.py runs: one
batch_result entry with the Merkle root over every md/results/ file.
Refuses if any of the 14 runs is missing, exited non-zero, or would be
dropped by .gitignore."""
import hashlib, json, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
MT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(MT, "pilot", "src"))
from ledger import Ledger, leaf_hash, merkle_root  # noqa: E402

D = json.load(open(os.path.join(HERE, "md_design.json")))
RES = os.path.join(HERE, "md", "results")
pts = [l.split() for l in open(os.path.join(HERE, "md", "points.txt")) if l.strip()]
if len(pts) != 14:
    sys.exit("points.txt must list 14 runs")
for tag, _, _, _, seed in pts:
    st = os.path.join(RES, f"{tag}_s{seed}.exit_status")
    if not os.path.exists(st) or open(st).read().strip() != "0":
        sys.exit(f"{tag}_s{seed} missing or did not exit 0")
files = {}
for f in sorted(os.listdir(RES)):
    p = os.path.join(RES, f)
    if subprocess.run(["git", "check-ignore", "-q", p]).returncode == 0:
        sys.exit(f"{f} would be dropped by .gitignore")
    files[f] = hashlib.sha256(open(p, "rb").read()).hexdigest()
leaves = [leaf_hash({"file": k, "sha256": v}) for k, v in sorted(files.items())]
led = Ledger(os.path.join(MT, "phase_c", "ledger_campaign.jsonl"))
op = "batch-" + D["campaign"]
if any(e["payload"].get("op_id") == op for e in led.entries):
    sys.exit("already committed")
e = led.append("batch_result", {
    "op_id": op, "campaign": D["campaign"], "stage": "embedded-SiC-MD", "raw_merkle_root": merkle_root(leaves),
    "n_files": len(files), "files_sha256": files, "chpc_job": "15904258 (notchpeak-guest)",
    "statement": "committed before analyze_md.py was run on these files"})
assert led.verify_chain()
print("seq", e["seq"])
