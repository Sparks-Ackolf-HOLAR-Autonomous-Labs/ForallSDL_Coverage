#!/usr/bin/env python3
"""Limitation 1: run cost is not uniform.

On the registered P3 SiGeC space and scoring (../p3_sigec), compares for three
cost models (cost_models.py):
  uniform        the registered uniform-cost greedy covering (seeds 0-19; doubled
                 from 1000-1019 on top of its own single batch)
  cost_weighted  greedy on (deficient pair cells touched) / (marginal cost of the run)
  pellet_greedy  lab_shared_pellet only: greedy whose unit is a pellet measured at a
                 set of temperatures
  ilp            exact minimum-cost covering (HiGHS via scipy.milp, 120 s limit),
                 reported with its proven lower bound
  random_at_*    random batches filled up to the median cost of a design (200 each)
Every batch is scored with the registered surrogate and screening score.
Also: pellet-level losses for the shared-pellet model (a cracked pellet removes
every measurement made on it). Writes results/cost_covering.json.
"""
import time
import numpy as np
import common as C
import cost_models as M

SEEDS = range(20)
NRND = 200


def summ(v):
    v = np.asarray(v, float)
    return {"median": float(np.median(v)), "q25": float(np.quantile(v, .25)), "q75": float(np.quantile(v, .75)),
            "min": float(v.min()), "max": float(v.max())}


def score(sets, y, yd=None):
    ev = [C.R.evaluate(s, y)[0] for s in sets]
    out = {"runs": summ([len(s) for s in sets]), "coverage2": summ([C.coverage_t(s, 2) for s in sets]),
           "complete_frac": float(np.mean([C.coverage_t(s, 2) == 1 for s in sets])),
           "regret": summ([e["regret"] for e in ev]), "log_rmse": summ([e["log_rmse"] for e in ev]),
           "found_best_frac": float(np.mean([e["regret"] == 0 for e in ev]))}
    if yd is not None:
        out["regret_density_contrast"] = summ([C.R.evaluate(s, yd)[0]["regret"] for s in sets])
    return out


def random_equal_cost(model, budget, rng):
    """Random rows in random order; a row is added if its marginal cost still fits the budget."""
    rows, st, spent = [], model.init([]), 0.0
    for i in rng.permutation(C.N):
        m = model.marg1(st, i)
        if spent + m > budget:
            continue
        rows.append(C.FULL[i]); model.add(st, i); spent += m
        if budget - spent < 1.0:
            break
    return rows


def levels(B):
    return {C.NAMES[j]: [float(np.median([sum(1 for r in b if r[j] == l) for b in B])) for l in range(C.NL[j])]
            for j in (2, 3, 5)}


def pellet_losses(rows, k, trials, rng):
    rec = np.array([M.RECIPE[C.R.IDX[tuple(r)]] for r in rows])
    u = np.unique(rec)
    ok = 0
    for _ in range(trials):
        lost = set(rng.choice(u, k, replace=False))
        ok += C.coverage_t([r for r, q in zip(rows, rec) if q not in lost], 2) == 1
    return ok / trials


def main():
    t0 = time.time()
    y = C.load_response()
    yd = C.R.response(C.R.D["model_settings"]["SiC_contrast_check"])
    res = {"cost_model_doc": M.__doc__, "models": {}}
    U = {1: [C.R.covering(1, s) for s in SEEDS]}
    U[2] = [C.R.covering(2, 1000 + s, U[1][s]) for s in SEEDS]
    SU = {m: score(U[m], y, yd) for m in (1, 2)}
    for model in M.MODELS:
        out = {"unit": model.unit}
        W = {1: [C.greedy_cover(2, 1, s, cost=model) for s in SEEDS]}
        W[2] = [C.greedy_cover(2, 2, 1000 + s, cost=model, base=W[1][s]) for s in SEEDS]
        designs = {"uniform": U, "cost_weighted": W}
        if model.name == "lab_shared_pellet":
            PG = {1: [M.greedy_pellet(1, s) for s in SEEDS]}
            PG[2] = [M.greedy_pellet(2, 1000 + s, PG[1][s]) for s in SEEDS]
            designs["pellet_greedy"] = PG
        for mult in (1, 2):
            tag = {1: "single", 2: "doubled"}[mult]
            d = {}
            rng = np.random.default_rng(17 + mult)
            for name, B in designs.items():
                if mult not in B:
                    continue
                cost = [model.total(r) for r in B[mult]]
                d[name] = {"cost": summ(cost), **(SU[mult] if name == "uniform" else score(B[mult], y, yd)),
                           "level_counts_median": levels(B[mult])}
                if model.name == "lab_shared_pellet":
                    d[name]["pellets"] = summ([len({M.RECIPE[C.R.IDX[tuple(r)]] for r in b}) for b in B[mult]])
                    if mult == 2:
                        d[name]["pellet_loss_complete_frac"] = {
                            f"k{k}": float(np.median([pellet_losses(b, k, 100, rng) for b in B[2]])) for k in (1, 2)}
            if model.name != "lab_shared_pellet":
                o = M.lp_bound(model, mult, integer=True, time_limit=120)
                d["ilp"] = {"cost": o.get("objective"), "lower_bound": o.get("dual_bound"), "status": o["message"],
                            "lp_relaxation": M.lp_bound(model, mult)["objective"]}
                if "rows" in o:
                    d["ilp"].update(score([o["rows"]], y, yd))
                    d["ilp"]["rows"] = [list(map(int, r)) for r in o["rows"]]
                    designs.setdefault("ilp", {})[mult] = [o["rows"]]
            for name in list(designs):
                budget = float(np.median([model.total(r) for r in designs[name][mult]]))
                RB = [random_equal_cost(model, budget, rng) for _ in range(NRND)]
                d[f"random_at_{name}_cost"] = {"budget": budget, **score(RB, y)}
            base = float(np.median([model.total(r) for r in U[mult]]))
            d["saving_vs_uniform_median"] = {n: 1 - float(np.median([model.total(r) for r in designs[n][mult]])) / base
                                             for n in designs if n != "uniform"}
            out[tag] = d
            print(model.name, tag, {k: round(v, 3) for k, v in d["saving_vs_uniform_median"].items()},
                  round(time.time() - t0), flush=True)
        res["models"][model.name] = out
    res["wall_seconds"] = time.time() - t0
    C.jdump(res, "cost_covering.json")


if __name__ == "__main__":
    main()
