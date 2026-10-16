#!/usr/bin/env python3
"""Limitation 3: replay the losses recorded on CHPC guest partitions onto the P3 covering batches,
with no requeue.

Designs: the registered uniform-cost single (CA1) and doubled (CA2) covering batches
(20 construction seeds each; the reference batches are CA1_ref = 33 runs, CA2_ref = 62 runs).
Loss patterns (from results/loss_record.json, written by loss_record.py):
  replay     the observed first-attempt outcome sequence of a campaign, in submission order,
             laid over the batch in its submission order, at every cyclic offset (keeps the
             clustering of losses on one faulty node within consecutive array tasks)
  bernoulli  independent losses at a recorded rate (5000 trials), including the aggregate-only
             MACE campaign and the one mid-run preemption among the per-run records
Batches are submitted either in construction order or in a random order (50 orders per batch).
Reports the fraction of trials that still sample every two-factor combination and the median
fraction kept; for the reference batches also the screening regret of the remaining runs.
Writes results/loss_replay.json.
"""
import json, os, time
import numpy as np
import common as C

SEEDS = range(20)


def kept(rows, lost):
    return [r for r, l in zip(rows, lost) if not l]


def stats(rows_list, masks):
    cov = [C.coverage_t(kept(r, m), 2) for r, m in zip(rows_list, masks)]
    return {"complete_frac": float(np.mean(np.array(cov) == 1)), "coverage_median": float(np.median(cov)),
            "coverage_q10": float(np.quantile(cov, 0.1)), "lost_mean": float(np.mean([sum(m) for m in masks]))}


def replay_masks(seq, n, offsets):
    seq = np.asarray(seq)
    return [seq[(o + np.arange(n)) % len(seq)] for o in offsets]


def main():
    t0 = time.time()
    L = json.load(open(os.path.join(C.OUT, "loss_record.json")))
    y = C.load_response()
    U1 = [C.R.covering(1, s) for s in SEEDS]
    U2 = [C.R.covering(2, 1000 + s, U1[s]) for s in SEEDS]
    ref1, ref2 = U1[16], U2[16]  # registered reference batches (33 and 62 runs)
    D = json.load(open(os.path.join(C.MT, "p3_sigec", "results", "designs.json")))
    assert [list(r) for r in ref1] == D["CA1_ref"] and [list(r) for r in ref2] == D["CA2_ref"]
    rng = np.random.default_rng(29)
    camp = L["campaigns"]
    seqs = {}
    for k in ("gpumd_endpoint_v2", "gpumd_sizefree", "p3_strain"):
        seqs[f"{k}:first_attempt"] = camp[k]["sequence_first_attempt_lost"]
        seqs[f"{k}:node_fault_only"] = camp[k]["sequence_node_fault"]
    del seqs["p3_strain:node_fault_only"]
    per_run_runs = sum(camp[k]["runs"] for k in camp)
    pre = sum(camp[k]["midrun_preemptions"] for k in camp)
    agg = L["aggregate_only"]
    rates = {f"{k}": float(np.mean(v)) for k, v in seqs.items()}
    rates["p3_mace:node_fault (aggregate only)"] = agg["p3_mace"]["node_fault_at_startup"] / agg["p3_mace"]["runs"]
    rates["midrun_preemption (per-run records + zero-loss aggregates)"] = pre / (
        per_run_runs + agg["chpc_phase_b"]["runs"] + agg["party4_bte"]["runs"])
    lost_all = sum(sum(camp[k]["sequence_first_attempt_lost"]) for k in camp) + agg["p3_mace"]["node_fault_at_startup"]
    rates["pooled first-attempt (all campaigns)"] = lost_all / (per_run_runs + agg["p3_mace"]["runs"]
                                                               + agg["chpc_phase_b"]["runs"] + agg["party4_bte"]["runs"])
    res = {"rates": rates, "replay": {}, "bernoulli": {}, "reference_regret": {}}
    for sname, seq in seqs.items():
        res["replay"][sname] = {}
        for dname, B in (("single", U1), ("doubled", U2)):
            n = len(seq)
            rl, ml = [], []
            for b in B:  # construction order
                for m in replay_masks(seq, len(b), range(n)):
                    rl.append(b); ml.append(m)
            rr, mr = [], []
            for b in B:  # random submission order
                for _ in range(50):
                    p = rng.permutation(len(b))
                    bb = [b[i] for i in p]
                    m = replay_masks(seq, len(b), [int(rng.integers(n))])[0]
                    rr.append(bb); mr.append(m)
            res["replay"][sname][dname] = {"construction_order": stats(rl, ml), "random_order": stats(rr, mr)}
    for rname, p in rates.items():
        res["bernoulli"][rname] = {}
        for dname, B in (("single", U1), ("doubled", U2)):
            rl, ml = [], []
            for b in B:
                for _ in range(250):
                    rl.append(b); ml.append(rng.random(len(b)) < p)
            res["bernoulli"][rname][dname] = stats(rl, ml)
    # screening score of the reference batches after replayed losses (construction order, every offset)
    for sname in ("gpumd_endpoint_v2:first_attempt", "gpumd_sizefree:first_attempt", "gpumd_endpoint_v2:node_fault_only"):
        seq = seqs[sname]
        res["reference_regret"][sname] = {}
        for dname, b in (("CA1_ref", ref1), ("CA2_ref", ref2)):
            reg = []
            for m in replay_masks(seq, len(b), range(len(seq))):
                kk = kept(b, m)
                reg.append(C.R.evaluate(kk, y)[0]["regret"] if len(kk) >= 2 else float("nan"))
            reg = np.array(reg)
            res["reference_regret"][sname][dname] = {"regret_median": float(np.nanmedian(reg)),
                                                    "found_best_frac": float(np.nanmean(reg == 0))}
    res["no_loss_reference_regret"] = {"CA1_ref": C.R.evaluate(ref1, y)[0]["regret"], "CA2_ref": C.R.evaluate(ref2, y)[0]["regret"]}
    # how many runs a requeue-free batch would need to be submitted at the recorded rates:
    # smallest multiple of covering (greedy, multiplicity m) whose median coverage stays complete in >= 95% of trials
    need = {}
    for rname in ("midrun_preemption (per-run records + zero-loss aggregates)", "p3_mace:node_fault (aggregate only)",
                  "pooled first-attempt (all campaigns)"):
        p = rates[rname]
        need[rname] = {}
        prev = None
        for mult in (1, 2, 3, 4, 6, 8):
            b = C.R.covering(mult, 2000 + mult, prev) if prev is not None else C.R.covering(1, 16)
            prev = b
            ok = np.mean([C.coverage_t(kept(b, rng.random(len(b)) < p), 2) == 1 for _ in range(300)])
            need[rname][mult] = {"runs": len(b), "complete_frac": float(ok)}
            if ok >= 0.95:
                break
    res["multiplicity_needed"] = need
    res["wall_seconds"] = time.time() - t0
    C.jdump(res, "loss_replay.json")
    print(json.dumps(rates, indent=0))
    for s, v in res["replay"].items():
        print(s, {d: (round(x["construction_order"]["complete_frac"], 3), round(x["construction_order"]["coverage_median"], 3),
                      round(x["random_order"]["complete_frac"], 3), round(x["random_order"]["coverage_median"], 3)) for d, x in v.items()})
    for s, v in res["bernoulli"].items():
        print(s, {d: (round(x["complete_frac"], 3), round(x["coverage_median"], 3)) for d, x in v.items()})
    print(res["reference_regret"], res["no_loss_reference_regret"])
    print(json.dumps(need))


if __name__ == "__main__":
    main()
