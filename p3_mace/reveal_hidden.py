#!/usr/bin/env python3
"""Draw the hidden test set from the hash of the registration ledger entry
(commit-reveal: the planner's designs were fixed before the seed existed).
Appends hidden points to chpc/points.json and a ledger entry. Idempotent."""
import hashlib, json, os, random, sys
HERE = os.path.dirname(os.path.abspath(__file__))
MT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(MT, "pilot", "src")); sys.path.insert(0, HERE)
from ledger import Ledger  # noqa: E402
import design_p3_mace as D  # noqa: E402
reg = open(os.path.join(HERE, "PREREGISTERED_SHA256")).read().split()
entry_hash = [t for t in reg if t.startswith("entry_hash=")][0].split("=")[1]
led = Ledger(os.path.join(MT, "phase_c", "ledger_campaign.jsonl"))
if any(e["payload"].get("op_id") == "reveal-hidden-p3-mace-2026-09-30" for e in led.entries):
    sys.exit("already revealed")
des = json.load(open(os.path.join(HERE, "planner_designs.json")))
used = {tuple(r) for v in des["designs"].values() for r in v}
rng = random.Random(int(hashlib.sha256(entry_hash.encode()).hexdigest(), 16))
hidden = rng.sample([r for r in D.FULL if r not in used], 16)
reps = rng.sample(hidden, 2)
P = json.load(open(os.path.join(HERE, "chpc", "points.json")))
hid = {}
for i, r in enumerate(hidden):
    hid["H_" + D.pid(r)] = D.point(r, 500 + i)
for i, r in enumerate(reps):
    hid["R_" + D.pid(r)] = D.point(r, 900 + i)
P["points"].update(hid)
json.dump(P, open(os.path.join(HERE, "chpc", "points.json"), "w"), indent=1)
json.dump({"seed_from_entry_hash": entry_hash, "hidden": [list(r) for r in hidden], "replicates": [list(r) for r in reps],
           "ids": list(hid)}, open(os.path.join(HERE, "hidden_set.json"), "w"), indent=1)
sha = hashlib.sha256(open(os.path.join(HERE, "chpc", "points.json"), "rb").read()).hexdigest()
e = led.append("design", {"op_id": "reveal-hidden-p3-mace-2026-09-30", "campaign": "p3-mace-hidden-subset-2026-09-30",
                          "paper": "P3", "seed_from_entry_hash": entry_hash, "n_hidden": 16, "n_replicates": 2,
                          "points_json_sha256": sha,
                          "statement": "hidden set drawn from the registration entry hash; all run points fixed before any production run"})
assert led.verify_chain()
print("seq", e["seq"], len(P["points"]), "points")
