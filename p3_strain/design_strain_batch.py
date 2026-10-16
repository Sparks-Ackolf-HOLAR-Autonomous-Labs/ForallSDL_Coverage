#!/usr/bin/env python3
"""Coverage-designed seed batch for the P3 biaxial-strain SiGe campaign
(scope set 2026-09-29).

Factors (levels):  x_Ge {0.1,0.3,0.5,0.7,0.9}; in-plane strain eps
{-0.03,-0.02,-0.01,0,+0.01} (relative to the alloy's own lattice at T;
pseudomorphic SiGe on Si(001) has eps ~ -0.042 x, on relaxed Ge ~ +0.04(1-x));
T {300,700,1100} K; transport direction {cross-plane, in-plane}.
Full factorial = 150 points.

Two designs are built deterministically and both are recorded:
  * C1: a strength-2 covering array (every pair of levels of every two
    factors appears at least once) of minimum size 25 = 5 x 5;
  * C2: a 2-fold strength-2 covering (every pair at least twice), the
    1-erasure-tolerant variant: any single lost run leaves every pair
    covered. Its extra size over C1 is the measured price of tolerance.
C2 is the batch that runs; C1 is its embedded subset (the first 25 rows),
so the paper can compare "what C1 alone would have seen" after real
guest-partition preemptions (the erasure channel) are counted.
Every point runs at ncz = 48 (13 nm long box, 6144 atoms), one seed each.
Writes strain_batch_design.json + points_strain.txt; the SHA-256 of the
json is appended to the ledger (preregister_strain.py) and pushed before
submission.
"""
import hashlib
import itertools
import json
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
X = [0.1, 0.3, 0.5, 0.7, 0.9]
EPS = [-0.03, -0.02, -0.01, 0.0, 0.01]
T = [300, 700, 1100]
DIR = [0, 1]
LEVELS = [X, EPS, T, DIR]
PAIRS = list(itertools.combinations(range(4), 2))


def counts(rows):
    c = {}
    for r in rows:
        for a, b in PAIRS:
            c[(a, b, r[a], r[b])] = c.get((a, b, r[a], r[b]), 0) + 1
    return c


def all_pair_keys():
    return [(a, b, u, v) for a, b in PAIRS for u in LEVELS[a] for v in LEVELS[b]]


def covering(mult, seed, base=None):
    rng = random.Random(seed)
    rows = list(base or [])
    full = list(itertools.product(*LEVELS))
    keys = all_pair_keys()
    while True:
        c = counts(rows)
        deficit = {k: mult - c.get(k, 0) for k in keys if c.get(k, 0) < mult}
        if not deficit:
            return rows
        best, bs = None, -1
        rng.shuffle(full)
        for cand in full:
            s = sum(1 for a, b in PAIRS if (a, b, cand[a], cand[b]) in deficit)
            if s > bs:
                best, bs = cand, s
        rows.append(best)


def c1_orthogonal():
    """25-row strength-2 covering: x and eps form a 5x5 grid; T and dir
    are assigned by a small deterministic search that covers the rest."""
    grid = list(itertools.product(range(5), range(5)))
    rng = random.Random(20260929)
    for _ in range(20000):
        rows = [(X[i], EPS[j], rng.choice(T), rng.choice(DIR)) for i, j in grid]
        c = counts(rows)
        if all(c.get(k, 0) >= 1 for k in all_pair_keys()):
            return rows
    raise RuntimeError("no 25-row covering found")


def main():
    c1 = c1_orthogonal()
    c2 = covering(2, 7, base=c1)
    pts, lines = [], []
    for i, (x, e, t, d) in enumerate(c2):
        rid = f"S{i:02d}_x{x}_e{e:+.2f}_T{t}_{'cp' if d == 0 else 'ip'}"
        row = {"run_id": rid, "x": x, "eps": e, "T": t, "dir": d, "ncz": 48,
               "seed": 29100 + i, "in_C1": i < len(c1)}
        pts.append(row)
        lines.append(f"{rid} {x} {t} {e} {d} 48 {row['seed']}")
    design = {
        "campaign": "p3-strain-seed-2026-09-29",
        "paper": "P3 (coverage-designed campaigns), own materials experiment",
        "factors": {"x": X, "eps": EPS, "T": T, "dir": {"0": "cross-plane", "1": "in-plane"}},
        "full_factorial": 150, "C1_size": len(c1), "C2_size": len(c2),
        "price_of_1_erasure_tolerance_runs": len(c2) - len(c1),
        "protocol": "p3_strain/in.mp_strain (Tersoff-1989 SiCGe.tersoff, MP rNEMD "
                    "along z, 500 ps production, 4x4 cross-section, swap 50), "
                    "notchpeak-guest 16 ranks, requeue on preemption (erasure "
                    "events counted from sacct + per-run attempts files)",
        "gate": "p3_strain/GATE.md must pass before submission",
        "analysis_preregistered": [
            "per point: kappa from analyze_mp.py; ratio to the unstrained "
            "reference kappa(x,T,eps=0,dir) from a fitted smooth model",
            "strain law: fit ln kappa = f(x,T) + (a + b x + c T/1000) eps + "
            "d eps * dir by weighted least squares over C2; report a, b, c, d "
            "with bootstrap CIs",
            "erasure analysis: for every run lost to preemption (never "
            "completed), check pair coverage of C1-only vs C2; report the "
            "fraction of pairs C1 would have lost",
            "materials call: kappa change of pseudomorphic SiGe on Si(001) "
            "(eps = -0.042 x) and on relaxed Ge (eps = +0.04 (1-x)) vs relaxed "
            "bulk, cross-plane and in-plane, at 300 and 1100 K",
        ],
        "predictions": [
            "compressive in-plane strain raises kappa and tensile lowers it "
            "(stiffer bonds, higher group velocities), in line with strained-Si "
            "literature (Parrish et al. PRB 2014)",
            "the strain sensitivity |d ln kappa / d eps| is smaller in the "
            "alloy than in pure Si because mass-disorder scattering dominates",
        ],
        "points": pts,
    }
    blob = json.dumps(design, indent=1, sort_keys=True).encode() + b"\n"
    open(os.path.join(HERE, "strain_batch_design.json"), "wb").write(blob)
    open(os.path.join(HERE, "points_strain.txt"), "w").write("\n".join(lines) + "\n")
    print("C1", len(c1), "C2", len(c2), "sha256", hashlib.sha256(blob).hexdigest())


if __name__ == "__main__":
    main()
