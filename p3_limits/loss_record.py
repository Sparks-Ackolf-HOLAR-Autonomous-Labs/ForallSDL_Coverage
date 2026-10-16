#!/usr/bin/env python3
"""Limitation 3: parse the per-run attempt records the project kept on CHPC guest partitions.

Sources with one record per run (an `attempts` file: one line per start with
timestamp, SLURM job id, requeue counter and host, plus annotations such as
`cuda_init_fail`):
  gpumd/production_v2/raw/*/attempts   P2 v2 endpoint campaign, 63 runs
  gpumd/sizefree/raw/*/attempts        size-free HNEMD campaign, 66 runs
  p3_strain/results/*.attempts         P3 strained-film campaign, 50 runs
Aggregate statements only (no per-run attempt record committed):
  p3_mace/README.md                    19 of 89 runs killed at start-up on notch458, rerun
  chpc/RESULTS.md                      Phase B, 19 + 10 tasks, zero preemptions
  party4_bte/RESULTS.md                6 tasks, zero preemptions
Writes results/loss_record.json.
"""
import glob, os, re
import common as C

SRC = {"gpumd_endpoint_v2": "gpumd/production_v2/raw/*/attempts",
       "gpumd_sizefree": "gpumd/sizefree/raw/*/attempts",
       "p3_strain": "p3_strain/results/*.attempts"}
LINE = re.compile(r"^(\S+Z) job=(\S+) restart=(\d+) host=(\S+)")


def parse(path):
    att = []
    for ln in open(path):
        m = LINE.match(ln)
        if m:
            att.append({"t": m.group(1), "job": m.group(2), "restart": int(m.group(3)), "host": m.group(4), "notes": []})
        elif ln.strip() and att:
            att[-1]["notes"].append(ln.strip())
    first = att[0]
    if len(att) == 1:
        cause = "completed"
    elif any("cuda_init_fail" in n for n in first["notes"]):
        cause = "node_fault_at_startup"
    elif att[1]["restart"] > first["restart"]:
        cause = "preempted_requeued"
    else:
        cause = "cancelled_by_operator_and_resubmitted"
    return {"run": os.path.basename(os.path.dirname(path)) if path.endswith("/attempts") else os.path.basename(path)[:-9],
            "first_start": first["t"], "first_job": first["job"], "first_host": first["host"],
            "n_starts": len(att), "first_attempt_outcome": cause,
            "midrun_preemptions": sum(1 for a in att if a["restart"] > 0),
            "fault_hosts": sorted({re.search(r"host=(\S+)", n).group(1) for a in att for n in a["notes"]
                                   if "cuda_init_fail" in n})}


def main():
    out = {"campaigns": {}}
    for name, pat in SRC.items():
        runs = sorted((parse(p) for p in glob.glob(os.path.join(C.MT, pat))), key=lambda r: (r["first_start"], r["first_job"]))
        by = {}
        for r in runs:
            by[r["first_attempt_outcome"]] = by.get(r["first_attempt_outcome"], 0) + 1
        hosts = {}
        for r in runs:
            if r["first_attempt_outcome"] != "completed":
                hosts[r["first_host"]] = hosts.get(r["first_host"], 0) + 1
        out["campaigns"][name] = {"source": pat, "runs": len(runs), "starts": sum(r["n_starts"] for r in runs),
                                  "first_attempt_outcomes": by,
                                  "first_attempt_loss_rate": 1 - by.get("completed", 0) / len(runs),
                                  "node_fault_rate": by.get("node_fault_at_startup", 0) / len(runs),
                                  "midrun_preemptions": sum(r["midrun_preemptions"] for r in runs),
                                  "first_attempt_losses_by_host": hosts,
                                  "sequence_first_attempt_lost": [int(r["first_attempt_outcome"] != "completed") for r in runs],
                                  "sequence_node_fault": [int(r["first_attempt_outcome"] == "node_fault_at_startup") for r in runs],
                                  "per_run": runs}
    out["aggregate_only"] = {
        "p3_mace": {"runs": 89, "node_fault_at_startup": 19, "host": "notch458",
                    "source": "p3_mace/README.md (uncorrectable GPU ECC error; rerun as array 15950354); which 19 runs is not recorded"},
        "chpc_phase_b": {"runs": 29, "losses": 0, "source": "chpc/RESULTS.md, chpc/JOBS.md"},
        "party4_bte": {"runs": 6, "losses": 0, "source": "party4_bte/RESULTS.md"}}
    out["readme_statement_not_reproduced"] = (
        "gpumd/production_v2/README.md states 205 attempts for the 63 endpoint runs; the committed per-run files "
        "record %d starts. The difference is not in any committed per-run record." % out["campaigns"]["gpumd_endpoint_v2"]["starts"])
    C.jdump(out, "loss_record.json")
    for k, v in out["campaigns"].items():
        print(k, v["runs"], v["starts"], v["first_attempt_outcomes"], v["first_attempt_losses_by_host"], v["midrun_preemptions"])


if __name__ == "__main__":
    main()
