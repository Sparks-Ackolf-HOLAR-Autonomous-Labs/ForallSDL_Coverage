#!/usr/bin/env python3
"""Ledger the r30 MD outputs BEFORE analyze_r30.py runs (Merkle root over results/).
Refuses if a run is missing, exited non-zero, or would be dropped by .gitignore."""
import hashlib, json, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
MT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, os.path.join(MT, "pilot", "src"))
from ledger import Ledger, leaf_hash, merkle_root  # noqa: E402

D = json.load(open(os.path.join(HERE, "r30_design.json")))
RES = os.path.join(HERE, "results")
for tag, _, _, _, _, _, seed in (l.split() for l in open(os.path.join(HERE, "points.txt")) if l.strip()):
    st = os.path.join(RES, f"{tag}_s{seed}.exit_status")
    if not os.path.exists(st) or open(st).read().strip() != "0":
        sys.exit(f"{tag}_s{seed} missing or did not exit 0")
files = {}
for f in sorted(os.listdir(RES)):
    p = os.path.join(RES, f)
    if subprocess.run(["git", "check-ignore", "-q", p]).returncode == 0:
        sys.exit(f"{f} would be dropped by .gitignore")
    files[f] = hashlib.sha256(open(p, "rb").read()).hexdigest()
led = Ledger(os.path.join(MT, "phase_c", "ledger_campaign.jsonl"))
op = "batch-" + D["campaign"]
if any(e["payload"].get("op_id") == op for e in led.entries):
    sys.exit("already committed")
e = led.append("batch_result", {
    "op_id": op, "campaign": D["campaign"], "stage": "embedded-SiC-MD-r30",
    "raw_merkle_root": merkle_root([leaf_hash({"file": k, "sha256": v}) for k, v in sorted(files.items())]),
    "n_files": len(files), "files_sha256": files, "chpc_job": sys.argv[1] if len(sys.argv) > 1 else "",
    "statement": "committed before analyze_r30.py was run on these files"})
assert led.verify_chain()
print("seq", e["seq"])
