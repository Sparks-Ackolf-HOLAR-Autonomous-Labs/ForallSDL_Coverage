#!/usr/bin/env python3
"""Ledger the M1c outputs BEFORE summarize_m1c.py runs: one batch_result
entry with the Merkle root over every out/ file. Refuses if any of the 51
tasks is missing, exited non-zero, or would be dropped by .gitignore."""
import hashlib, json, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
MT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(MT, "pilot", "src"))
sys.path.insert(0, os.path.join(HERE, "chpc"))
from ledger import Ledger, leaf_hash, merkle_root  # noqa: E402
D = json.load(open(os.path.join(HERE, "m1c_design.json")))
OUT = os.path.join(HERE, "out")
tasks = [l for l in (json.load(open(os.path.join(OUT, f)))["task"] for f in sorted(os.listdir(OUT)) if f.endswith(".json"))]
if len(tasks) != 51:
    sys.exit(f"expected 51 task outputs, found {len(tasks)}")
files = {}
for f in sorted(os.listdir(OUT)):
    p = os.path.join(OUT, f)
    if subprocess.run(["git", "check-ignore", "-q", p]).returncode == 0:
        sys.exit(f"{f} would be dropped by .gitignore")
    files[f] = hashlib.sha256(open(p, "rb").read()).hexdigest()
for i in range(51):
    if open(os.path.join(OUT, f"task{i}.exit_status")).read().strip() != "0":
        sys.exit(f"task {i} did not exit 0")
leaves = [leaf_hash({"file": k, "sha256": v}) for k, v in sorted(files.items())]
led = Ledger(os.path.join(MT, "phase_c", "ledger_campaign.jsonl"))
op = "batch-" + D["campaign"]
if any(e["payload"].get("op_id") == op for e in led.entries):
    sys.exit("already committed")
e = led.append("batch_result", {
    "op_id": op, "campaign": D["campaign"], "stage": "M1c-converged", "raw_merkle_root": merkle_root(leaves),
    "n_files": len(files), "files_sha256": files, "chpc_job": "15898834 (notchpeak-guest)",
    "statement": "committed before summarize_m1c.py was run on these files"})
assert led.verify_chain()
print("seq", e["seq"])
