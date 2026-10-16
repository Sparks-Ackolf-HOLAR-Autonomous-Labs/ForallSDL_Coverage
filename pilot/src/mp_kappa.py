"""Muller-Plathe reverse-NEMD kappa point for Si1-xGex with SiCGe.tersoff.

Usage: python mp_kappa.py <x_ge> <T> <outjson> [threads] [ncell_z]

Why not Green-Kubo: this LAMMPS build (pip wheel 22 Jul 2025) reports
centroid/stress/atom unsupported for pair tersoff (many-body pair styles
are CENTROID_NOTAVAIL in mainline LAMMPS), and plain stress/atom heat
flux is quantitatively unreliable for many-body potentials. Muller-Plathe
(fix thermal/conductivity) imposes the flux by kinetic-energy swaps and
measures the gradient, so no per-atom heat flux is needed at all, and
OpenMP threading is fine.

Geometry: 4 x 4 x ncell_z diamond cells (default 24 -> 3072 atoms).
Protocol: NPT equil 50 ps, MP warmup 100 ps (discarded), production
500 ps with running-averaged temperature profile. kappa from
E_swapped / (2 A t) / <|dT/dz|>, slopes fit on both halves of the
profile excluding the swap slabs. Finite-size caveat: at these box
lengths alloy kappa is underestimated (long-wavelength phonons cut off);
run two ncell_z values to see the trend. Production-grade sizes are
CHPC work.
"""

from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
from lammps import lammps

POT = os.path.join(os.path.dirname(__import__("lammps").__file__),
                   "share", "lammps", "potentials", "SiCGe.tersoff")

NBIN = 20
SWAP_EVERY = 50
EQUIL, WARM, PROD = 50_000, 100_000, 500_000
EV_PS_A_TO_W_MK = 1.602176634e-19 / (1e-12 * 1e-10)   # eV/(ps A K) -> W/(m K)


def main():
    x, T, out = float(sys.argv[1]), float(sys.argv[2]), sys.argv[3]
    threads = int(sys.argv[4]) if len(sys.argv) > 4 else 5
    ncz = int(sys.argv[5]) if len(sys.argv) > 5 else 24
    t0 = time.time()
    a0 = 5.431 + (5.658 - 5.431) * x
    seed = 20260925
    lmp = lammps(cmdargs=["-log", out + ".log", "-screen", "none",
                          "-sf", "omp", "-pk", "omp", str(threads)])
    lmp.commands_string(f"""
units metal
boundary p p p
lattice diamond {a0}
region box block 0 4 0 4 0 {ncz}
create_box 2 box
create_atoms 1 box
{"set type 1 type/fraction 2 " + str(x) + " " + str(seed) if x > 0 else ""}
mass 1 28.0855
mass 2 72.63
pair_style tersoff
pair_coeff * * {POT} Si(D) Ge
timestep 0.001
velocity all create {2*T} {seed} dist gaussian
fix npt all npt temp {T} {T} 0.1 aniso 0.0 0.0 1.0
run {EQUIL}
unfix npt
fix nve all nve
fix flux all thermal/conductivity {SWAP_EVERY} z {NBIN}
run {WARM}
""")
    e_warm = lmp.extract_fix("flux", 0, 0)  # cumulative swapped KE, eV
    lmp.commands_string(f"""
compute keatom all ke/atom
compute bins all chunk/atom bin/1d z lower {1.0/NBIN} units reduced
fix prof all ave/chunk 10 100 1000 bins temp ave running overwrite file {out}.profile
run {PROD}
""")
    e_prod = lmp.extract_fix("flux", 0, 0)
    lx = lmp.extract_box()[1][0] - lmp.extract_box()[0][0]
    ly = lmp.extract_box()[1][1] - lmp.extract_box()[0][1]
    lz = lmp.extract_box()[1][2] - lmp.extract_box()[0][2]
    natoms = lmp.get_natoms()
    lmp.close()

    # last running-averaged profile block: rows "bin coord ncount temp"
    rows = []
    for line in open(out + ".profile"):
        p = line.split()
        if len(p) == 4 and not line.startswith("#"):
            rows.append((float(p[1]), float(p[3])))
    prof = np.array(rows[-NBIN:])          # (coord in reduced z, T)
    z = prof[:, 0] * lz
    tt = prof[:, 1]
    h = NBIN // 2
    # exclude the two swap slabs (bin 0 and bin NBIN/2) and their neighbors
    s1 = np.polyfit(z[2:h-1], tt[2:h-1], 1)[0]
    s2 = np.polyfit(z[h+2:NBIN-1], tt[h+2:NBIN-1], 1)[0]
    slope = 0.5 * (abs(s1) + abs(s2))      # K per Angstrom
    de = e_prod - e_warm                   # eV over production
    j = de / (2 * lx * ly * PROD * 0.001)  # eV / (ps A^2), 2 flux directions
    kappa = j / slope * EV_PS_A_TO_W_MK
    grad_K = slope * lz / 2

    result = {"x": x, "T": T, "method": "Muller-Plathe rNEMD",
              "natoms": natoms, "box_A": [round(lx, 2), round(ly, 2), round(lz, 2)],
              "kappa_W_mK": round(float(kappa), 3),
              "slopes_K_per_A": [round(float(s1), 5), round(float(s2), 5)],
              "deltaT_across_half_box_K": round(float(grad_K), 1),
              "swapped_eV": round(float(de), 1),
              "equil_warm_prod_steps": [EQUIL, WARM, PROD],
              "threads": threads, "wall_min": round((time.time() - t0) / 60, 1)}
    with open(out, "w") as f:
        json.dump(result, f, indent=2)
    print(json.dumps(result))


if __name__ == "__main__":
    main()
