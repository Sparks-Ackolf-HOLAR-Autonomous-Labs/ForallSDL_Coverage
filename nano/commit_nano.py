#!/usr/bin/env python3
"""Ledger the nano per-mode evidence BEFORE analysis: one batch_result entry
with the Merkle root over (raw file name, sha256) leaves and the sha256 of
every committed npz. Refuses if any pre-registered composition is missing."""
import hashlib, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
MT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(MT, "pilot", "src"))
from ledger import Ledger, leaf_hash, merkle_root  # noqa: E402
D = json.load(open(os.path.join(HERE, "nano_design.json")))
raw = json.load(open(os.path.join(HERE, "modes", "raw_sha256.json")))
npz = {}
for x in D["grid"]["x_Ge"]:
    p = os.path.join(HERE, "modes", f"x{x:g}.npz")
    if not os.path.exists(p):
        sys.exit(f"missing composition x={x:g}")
    npz[f"x{x:g}.npz"] = hashlib.sha256(open(p, "rb").read()).hexdigest()
leaves = [leaf_hash({"file": k, "sha256": v}) for k, v in sorted(raw.items())]
led = Ledger(os.path.join(MT, "phase_c", "ledger_campaign.jsonl"))
if any(e["payload"].get("op_id") == "batch-" + D["campaign"] for e in led.entries):
    sys.exit("already committed")
e = led.append("batch_result", {
    "op_id": "batch-" + D["campaign"], "campaign": D["campaign"], "stage": "N0-N1-modes",
    "raw_merkle_root": merkle_root(leaves), "n_raw_files": len(raw), "npz_sha256": npz,
    "raw_location": "CHPC home ~/mat_theory/nano/modes/ (phononinfo, Tsweep); npz reductions in git",
    "statement": "committed before analyze_nano.py was run on any of these files"})
assert led.verify_chain()
print("seq", e["seq"])
