#!/usr/bin/env python3
"""Pre-register the P3 strain seed batch on the campaign ledger: SHA-256
of strain_batch_design.json + Merkle root over its points. Pushed to git
before submission. Idempotent. Usage: python3 preregister_strain.py"""
import hashlib
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
MT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(MT, "pilot", "src"))
from ledger import Ledger, leaf_hash, merkle_root  # noqa: E402

LEDGER = os.path.join(MT, "phase_c", "ledger_campaign.jsonl")


def main():
    blob = open(os.path.join(HERE, "strain_batch_design.json"), "rb").read()
    sha = hashlib.sha256(blob).hexdigest()
    d = json.loads(blob)
    led = Ledger(LEDGER)
    for e in led.entries:
        if e["payload"].get("design_sha256") == sha:
            print("already pre-registered at seq", e["seq"])
            return
    e = led.append("preregistration", {
        "op_id": "prereg-" + d["campaign"], "campaign": d["campaign"],
        "paper": d["paper"], "design_file": "p3_strain/strain_batch_design.json",
        "design_sha256": sha,
        "runs_merkle_root": merkle_root([leaf_hash(p) for p in d["points"]]),
        "n_runs": len(d["points"]), "C1_size": d["C1_size"],
        "statement": "coverage-designed seed batch (2-fold strength-2 covering, "
                     "C1 embedded), analysis and predictions fixed before "
                     "submission; gate in p3_strain/GATE.md passed first"})
    assert led.verify_chain()
    open(os.path.join(HERE, "PREREGISTERED_SHA256"), "w").write(
        f"{sha}  strain_batch_design.json  ledger_seq={e['seq']}\n")
    print("pre-registered seq", e["seq"], sha)


if __name__ == "__main__":
    main()
