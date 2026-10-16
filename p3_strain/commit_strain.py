#!/usr/bin/env python3
"""Commit the P3 strain seed batch's raw outputs to the ledger BEFORE
analysis: per-run leaves (run_id + sha256 of json/log/profile/attempts)
-> Merkle root, appended as a batch_result. Completed runs only; the
list of runs never completed is recorded (erasures). Idempotent.
Usage: python3 commit_strain.py"""
import hashlib
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
MT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(MT, "pilot", "src"))
from ledger import Ledger, leaf_hash, merkle_root  # noqa: E402

LEDGER = os.path.join(MT, "phase_c", "ledger_campaign.jsonl")


def main():
    d = json.load(open(os.path.join(HERE, "strain_batch_design.json")))
    led = Ledger(LEDGER)
    op = "commit-" + d["campaign"]
    if any(e["payload"].get("op_id") == op for e in led.entries):
        sys.exit("already committed")
    leaves, done, missing = [], [], []
    for p in d["points"]:
        files = {}
        for ext in ("json", "log", "profile", "attempts"):
            f = os.path.join(HERE, "results", f"{p['run_id']}.{ext}")
            if os.path.exists(f):
                if subprocess.run(["git", "check-ignore", "-q", f], cwd=HERE).returncode == 0:
                    sys.exit(f"REFUSING: {f} is git-ignored")
                files[ext] = hashlib.sha256(open(f, "rb").read()).hexdigest()
        if "json" in files:
            done.append(p["run_id"])
            leaves.append(leaf_hash({"run_id": p["run_id"], "files": files}))
        else:
            missing.append(p["run_id"])
    e = led.append("batch_result", {
        "op_id": op, "campaign": d["campaign"], "n_completed": len(done),
        "never_completed": missing, "merkle_root": merkle_root(leaves),
        "statement": "raw outputs committed before analysis; requeue record "
                     "in results/sacct_erasures.txt"})
    assert led.verify_chain()
    print("committed seq", e["seq"], len(done), "done,", len(missing), "missing")


if __name__ == "__main__":
    main()
