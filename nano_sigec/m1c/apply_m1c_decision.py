#!/usr/bin/env python3
"""Apply the M1c pre-registered decision rule (m1c_design.json): gates pass
and PM2-PM4 hold -> ledger 'strengthen' on the where-carbon-goes claims
(the M1 'unvalidated' flag is lifted) and record the analyze_sigec re-run on
the converged summary."""
import hashlib, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
MT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(MT, "pilot", "src"))
from ledger import Ledger  # noqa: E402
S = json.load(open(os.path.join(HERE, "m1c_summary.json")))
P = S["predictions"]
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
res = os.path.join(HERE, "..", "results", "sigec_results_m1c.json")
led = Ledger(os.path.join(MT, "phase_c", "ledger_campaign.jsonl"))
op_id = "m1c-strengthen-where-carbon-goes-2026-09-30"
if any(e["payload"].get("op_id") == op_id for e in led.entries):
    sys.exit("already applied")
if not (S["gates_pass"] and all(P[k]["pass"] for k in ("PM2", "PM3", "PM4"))):
    sys.exit("decision rule does not prescribe strengthen; apply the gate_failure / correct branch by hand")
e = led.append("shard_operation", {
    "op_id": op_id, "op": "strengthen", "campaign": "nano-sigec-m1c-converged-2026-09-30",
    "target": {"claims": ["dE_prec (carbon leaves the lattice to form SiC)", "dE_site (carbon prefers Si-rich shells)",
                          "delta_C (substitutional size misfit)"],
               "supersedes_flag": "ledger seq 48 (M1b gate failure: energies unvalidated)"},
    "cause": "pre-registered rule (seq 53): gates M1c-a/b/c pass (51/51 relaxations, max force "
             f"{S['M1b']['max_fmax_eV_A']:.4f} eV/A, max |p| {S['M1b']['max_abs_pressure_GPa']:.3f} GPa) "
             "and PM2-PM4 hold (shifts <= 0.02 eV, delta_C <= 1.9%)",
    "decision": {"by": "rule (pre-registered), applied by agent", "second_bar": "independent converged rerun"},
    "evidence": {"summary_file": "nano_sigec/m1c/m1c_summary.json", "summary_sha256": sha(os.path.join(HERE, "m1c_summary.json")),
                 "rerun_results_file": "nano_sigec/results/sigec_results_m1c.json", "rerun_results_sha256": sha(res),
                 "raw_batch": "ledger seq 54"},
    "effect": "MACE delta_C now used per x (analyze_sigec pre-registered branch); only the substitutional strain "
              "upper bound moves (x=0.2, 300 K, 0.3 um, 1 at% C: 2.533 -> 2.497 W/mK); predictions PrC1-PrC6 verdicts unchanged"})
assert led.verify_chain()
print("seq", e["seq"])
