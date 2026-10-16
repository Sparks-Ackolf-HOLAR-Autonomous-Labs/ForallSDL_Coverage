"""Green-Kubo kappa point for Si1-xGex with SiCGe.tersoff.

Usage: python gk_kappa.py <x_ge> <T> <outjson> [threads]

Follows the standard LAMMPS GK recipe (fix ave/correlate + trap), with
centroid/stress/atom for the many-body virial (required for Tersoff;
plain stress/atom heat flux is wrong for many-body potentials).

Pilot settings: 1728 atoms, 1 fs step, 100 ps NPT equil, 1 ns NVE
sampling, 20 ps correlation window. Adequate for alloy compositions
(alloy scattering shortens phonon lifetimes); NOT converged for pure
Si/Ge, which need bigger cells and longer correlations (CHPC scale).

NOTE: runs SERIAL on purpose. The OPENMP variant tersoff/omp does not
implement centroid/stress/atom, which the many-body heat flux needs
(plain stress/atom heat flux is quantitatively wrong for Tersoff).
On CHPC, use MPI ranks instead of OpenMP threads for GK points.
"""

from __future__ import annotations

import json
import os
import sys
import time

from lammps import lammps

POT = os.path.join(os.path.dirname(__import__("lammps").__file__),
                   "share", "lammps", "potentials", "SiCGe.tersoff")

S, P = 10, 2000          # sample every 10 steps; 2000 windows -> 20 ps
D = S * P                # dump correlation every 20 ps
EQUIL, SAMPLE = 100_000, 1_000_000
NCELL = 6                # 6x6x6 diamond cells = 1728 atoms


def main():
    x, T, out = float(sys.argv[1]), float(sys.argv[2]), sys.argv[3]
    threads = 1   # serial: see NOTE in module docstring
    t0 = time.time()
    lmp = lammps(cmdargs=["-log", out + ".log", "-screen", "none"])
    a0 = 5.431 + (5.658 - 5.431) * x   # Vegard start; NPT settles the rest
    seed = 20260925
    lmp.commands_string(f"""
units metal
boundary p p p
lattice diamond {a0}
region box block 0 {NCELL} 0 {NCELL} 0 {NCELL}
create_box 2 box
create_atoms 1 box
{"set type 1 type/fraction 2 " + str(x) + " " + str(seed) if x > 0 else ""}
mass 1 28.0855
mass 2 72.63
pair_style tersoff
pair_coeff * * {POT} Si(D) Ge
timestep 0.001
velocity all create {2*T} {seed} dist gaussian
fix npt all npt temp {T} {T} 0.1 iso 0.0 0.0 1.0
run {EQUIL}
unfix npt
fix nve all nve
compute myKE all ke/atom
compute myPE all pe/atom
compute myStress all centroid/stress/atom NULL virial
compute flux all heat/flux myKE myPE myStress
fix JJ all ave/correlate {S} {P} {D} c_flux[1] c_flux[2] c_flux[3] type auto file {out}.J0Jt.dat ave running
variable scale equal 1.60218e-19*1.60218e-19/1.3806504e-23/{T}/{T}/vol/1e-12/1e-10*{S}*0.001
variable k11 equal trap(f_JJ[3])*v_scale
variable k22 equal trap(f_JJ[4])*v_scale
variable k33 equal trap(f_JJ[5])*v_scale
thermo {D}
thermo_style custom step temp v_k11 v_k22 v_k33
run {SAMPLE}
""")
    k = [lmp.extract_variable(f"k{i}{i}", None, 0) for i in (1, 2, 3)]
    vol = lmp.get_thermo("vol")
    lmp.close()
    result = {"x": x, "T": T, "natoms": 8 * NCELL**3, "vol_A3": vol,
              "kappa_W_mK": round(sum(k) / 3, 3),
              "kappa_xyz": [round(v, 3) for v in k],
              "equil_steps": EQUIL, "sample_steps": SAMPLE,
              "corr_window_ps": S * P * 0.001, "threads": threads,
              "wall_min": round((time.time() - t0) / 60, 1)}
    with open(out, "w") as f:
        json.dump(result, f, indent=2)
    print(json.dumps(result))


if __name__ == "__main__":
    main()
