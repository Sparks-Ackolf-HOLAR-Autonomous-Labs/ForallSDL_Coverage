#!/usr/bin/env python3
"""Write SiCGe_GeC.tersoff: the LAMMPS SiCGe.tersoff file (Tersoff 1989 Si/C
and Si/Ge entries, unchanged) plus the 12 triplets that involve a Ge-C pair,
built with the same mixing rules the upstream file uses (R and S geometric
means, lambda1 and lambda2 arithmetic means, A and B geometric means, with
chi_GeC = 1 because Tersoff 1989 gives no Ge-C value).

These Ge-C entries exist only so LAMMPS accepts a three-element run. The
structure builder turns every Ge within R + 7 A of a particle centre into Si,
and the run asserts (gate G1 in ../md_design.json) that no Ge atom has a C
atom within 5.7 A. That distance exceeds the Si-C cutoff (2.51 A) plus the
Si-Ge cutoff (3.05 A), so no pair or triplet term involving both Ge and C is
ever evaluated, and the added parameters cannot affect the result."""
import itertools, math, os
H = os.path.dirname(os.path.abspath(__file__))
# element: c d h n beta lambda2 B R D lambda1 A   (m=3, gamma=1, lambda3=0)
P = {"Si(D)": (1.0039e5, 16.217, -0.59825, 0.78734, 1.1000e-6, 1.7322, 471.18, 2.85, 0.15, 2.4799, 1830.8),
     "Ge": (1.0643e5, 15.652, -0.43884, 0.75627, 9.0166e-7, 1.7047, 419.23, 2.95, 0.15, 2.4451, 1769.0),
     "C": (3.8049e4, 4.3484, -0.57058, 0.72751, 1.5724e-7, 2.2119, 346.74, 1.95, 0.15, 3.4879, 1393.6)}


def pair(i, k):
    Ri, Di = P[i][7], P[i][8]
    Rk, Dk = P[k][7], P[k][8]
    R = math.sqrt(Ri * Rk)
    return R, math.sqrt((Ri + Di) * (Rk + Dk)) - R


up = open(os.path.join(H, "SiCGe_lammps_upstream.tersoff")).read()
have = set()
for line in up.splitlines():
    t = line.split()
    if len(t) >= 3 and not line.startswith("#") and t[0] in P and t[1] in P and t[2] in P:
        have.add(tuple(t[:3]))
new = []
for i, j, k in itertools.product(P, repeat=3):
    if (i, j, k) in have or not ({"Ge", "C"} <= {i, j, k}):
        continue
    c, d, h, n, beta = P[i][:5]
    R, D = pair(i, k)
    if j == k:
        lam2 = (P[i][5] + P[j][5]) / 2
        B = math.sqrt(P[i][6] * P[j][6])  # chi = 1
        lam1 = (P[i][9] + P[j][9]) / 2
        A = math.sqrt(P[i][10] * P[j][10])
    else:
        n = beta = lam2 = B = lam1 = A = 0.0
    new.append(f"{i:6s} {j:6s} {k:6s} 3.0 1.0 0.0 {c:.4e} {d:.4f} {h:.5f} {n:.5f} {beta:.4e} "
               f"{lam2:.5f} {B:.4f} {R:.4f} {D:.4f} {lam1:.4f} {A:.4f}")
assert len(new) == 12, len(new)
open(os.path.join(H, "SiCGe_GeC.tersoff"), "w").write(
    up + "\n# ---- added by make_pot.py: Ge-C triplets by the upstream mixing rules, chi_GeC = 1.\n"
         "# Never evaluated in the embedded-SiC runs (Si shell + gate G1 asserts zero Ge-C pairs within 2.551 A).\n"
    + "\n".join(new) + "\n")
print("\n".join(new))
