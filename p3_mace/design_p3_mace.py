#!/usr/bin/env python3
"""P3 hidden-subset study, step 1: build the design space that MACE-MD can
represent, let the planner build its batches from factor levels alone (the same
greedy covering, random and one-factor-at-a-time rules as p3_sigec), and write
points.json with every design point. The hidden test set is NOT drawn here: it
is drawn by reveal_hidden.py from the hash of the ledger entry that registers
this file, so neither the planner nor the author can choose it.
Deterministic; writes p3_mace_design.json and chpc/points.json."""
import itertools, json, os, random

HERE = os.path.dirname(os.path.abspath(__file__))
F = {"x_Ge": [0.1, 0.2, 0.3, 0.5, 0.7],
     "y_C": [0.0, 0.01, 0.03],
     "carbon_form": ["dispersed", "clustered"],
     "dopant": ["none", "B", "P"],
     "T_K": [300, 700, 1100]}
NAMES = list(F)
NL = [len(F[k]) for k in NAMES]
# carbon form is meaningless without carbon: y_C = 0 rows carry form index 0
FULL = sorted({tuple(r[:2]) + ((0,) if r[1] == 0 else (r[2],)) + tuple(r[3:])
               for r in itertools.product(*[range(n) for n in NL])})
PAIRS = [(a, b) for a, b in itertools.combinations(range(5), 2)]


def keys():
    k = set()
    for r in FULL:
        for a, b in PAIRS:
            if 2 in (a, b) and r[1] == 0:
                continue
            k.add((a, b, r[a], r[b]))
    return k


KEYS = keys()


def counts(rows):
    c = {}
    for r in rows:
        for a, b in PAIRS:
            if 2 in (a, b) and r[1] == 0:
                continue
            c[(a, b, r[a], r[b])] = c.get((a, b, r[a], r[b]), 0) + 1
    return c


def covering(mult, seed, base=()):
    rng = random.Random(seed)
    rows, full = list(base), list(FULL)
    while True:
        c = counts(rows)
        deficit = {k for k in KEYS if c.get(k, 0) < mult}
        if not deficit:
            return rows
        rng.shuffle(full)
        best, bs = None, -1
        for cand in full:
            if cand in rows:
                continue
            s = sum(1 for a, b in PAIRS if (a, b, cand[a], cand[b]) in deficit and not (2 in (a, b) and cand[1] == 0))
            if s > bs:
                best, bs = cand, s
        rows.append(best)


def ofat():
    base = (F["x_Ge"].index(0.2), F["y_C"].index(0.01), F["carbon_form"].index("clustered"), 0, F["T_K"].index(700))
    rows = [base]
    for j in range(5):
        for l in range(NL[j]):
            if l != base[j]:
                r = list(base); r[j] = l
                if r[1] == 0:
                    r[2] = 0
                if tuple(r) not in rows:
                    rows.append(tuple(r))
    return rows


def pid(r):
    x, y, f, d, T = (F[NAMES[j]][r[j]] for j in range(5))
    return f"x{x}_C{y}_{f[:4] if y else 'none'}_{d}_T{T}"


def point(r, seed):
    x, y, f, d, T = (F[NAMES[j]][r[j]] for j in range(5))
    return {"x_Ge": x, "y_C": y, "carbon_form": f if y else "none", "dopant": d, "dopant_at_frac": 0.002,
            "T_K": T, "seed": seed}


def main():
    ca1s = [covering(1, s) for s in range(20)]
    ref = min(range(20), key=lambda s: len(ca1s[s]))
    ca1 = ca1s[ref]
    ca2 = covering(2, 1000 + ref, ca1)
    of = ofat()
    rng = random.Random(20260930)
    pool = [r for r in FULL if r not in set(ca2) | set(of)]
    rnd = rng.sample(pool, len(ca2))
    design = {
        "campaign": "p3-mace-hidden-subset-2026-09-30", "paper": "P3",
        "operator_go": "issue 55, 2026-09-30: 'score the designs against a MACE-MD subset that the planner cannot see'",
        "factors": F, "n_full": len(FULL), "n_pair_keys": len(KEYS),
        "designs": {"CA1": [list(r) for r in ca1], "CA2": [list(r) for r in ca2], "OFAT": [list(r) for r in of],
                    "RND_pool": [list(r) for r in rnd]},
        "sizes": {"CA1": len(ca1), "CA2": len(ca2), "OFAT": len(of), "RND_pool": len(rnd),
                  "lower_bound_CA1": sorted(NL)[-1] * sorted(NL)[-2]},
    }
    rows = list(dict.fromkeys(ca2 + of + rnd))
    pts = {pid(r): point(r, 100 + i) for i, r in enumerate(rows)}
    json.dump(design, open(os.path.join(HERE, "planner_designs.json"), "w"), indent=1)
    json.dump({"protocol": {"model": "small", "ncx": 3, "ncy": 3, "ncz": 16, "relax_steps": 150, "dt_fs": 1.0,
                            "equil": 8000, "prod": 50000, "swap_every": 50}, "points": pts},
              open(os.path.join(HERE, "chpc", "points.json"), "w"), indent=1)
    print(design["sizes"], "design points to run:", len(pts))


if __name__ == "__main__":
    main()
