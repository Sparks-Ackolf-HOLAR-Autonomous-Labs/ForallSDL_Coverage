"""Pilot B: real LAMMPS/Tersoff bring-up on the DGX (pip wheel, aarch64).

1. Lattice-constant + cohesive-energy validation for Si, Ge, Si0.5Ge0.5
   with the bundled SiCGe.tersoff (Erhart-Albe 2005; elements Si(D), Ge).
   Phase B gate rehearsal: values must land near literature before any
   CHPC campaign.
2. MD throughput benchmark (NVT, 4096 atoms) to price a Green-Kubo
   kappa point on this box and on CHPC.

Run:  python lammps_pilot.py <outdir> [threads]
"""

from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

from lammps import lammps

POT = os.path.join(os.path.dirname(__import__("lammps").__file__),
                   "share", "lammps", "potentials", "SiCGe.tersoff")

REF = {  # experimental lattice constants, Angstrom
    "Si": 5.431, "Ge": 5.658, "Si0.5Ge0.5_vegard": 5.545,
}


def new_lmp(threads: int):
    args = ["-log", "none", "-screen", "none"]
    if threads > 1:
        args += ["-sf", "omp", "-pk", "omp", str(threads)]
    return lammps(cmdargs=args)


def relax(x_ge: float, threads: int, a0: float = 5.5, n: int = 3) -> dict:
    lmp = new_lmp(threads)
    cmds = f"""
units metal
boundary p p p
lattice diamond {a0}
region box block 0 {n} 0 {n} 0 {n}
create_box 2 box
create_atoms 1 box
"""
    if x_ge > 0:
        cmds += f"set type 1 type/fraction 2 {x_ge} 20260925\n"
    cmds += f"""
mass 1 28.0855
mass 2 72.63
pair_style tersoff
pair_coeff * * {POT} Si(D) Ge
fix rel all box/relax iso 0.0
minimize 1e-12 1e-12 5000 50000
"""
    lmp.commands_string(cmds)
    natoms = lmp.get_natoms()
    lx = lmp.extract_box()[1][0] - lmp.extract_box()[0][0]
    pe = lmp.get_thermo("pe")
    lmp.close()
    return {"a": round(lx / n, 4), "E_coh_eV_per_atom": round(pe / natoms, 4),
            "natoms": natoms}


def md_benchmark(threads: int, n: int = 8, steps: int = 5000) -> dict:
    lmp = new_lmp(threads)
    lmp.commands_string(f"""
units metal
boundary p p p
lattice diamond 5.545
region box block 0 {n} 0 {n} 0 {n}
create_box 2 box
create_atoms 1 box
set type 1 type/fraction 2 0.5 20260925
mass 1 28.0855
mass 2 72.63
pair_style tersoff
pair_coeff * * {POT} Si(D) Ge
velocity all create 600 4928459 dist gaussian
timestep 0.001
fix nvt all nvt temp 300 300 0.1
run 200
""")
    natoms = lmp.get_natoms()
    t0 = time.time()
    lmp.command(f"run {steps}")
    dt = time.time() - t0
    lmp.close()
    atom_steps_per_s = natoms * steps / dt
    # A conservative Green-Kubo kappa point: 4096 atoms, 2 M steps
    # (0.5 ns equil + 1.5 ns sampling at 1 fs)
    gk_cost_h = 4096 * 2_000_000 / atom_steps_per_s / 3600
    return {"natoms": natoms, "steps": steps, "threads": threads,
            "wall_s": round(dt, 2),
            "atom_steps_per_s": int(atom_steps_per_s),
            "est_greenkubo_point_hours": round(gk_cost_h, 2)}


def main():
    outdir = Path(sys.argv[1] if len(sys.argv) > 1 else "../results")
    threads = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    outdir.mkdir(parents=True, exist_ok=True)
    res = {"potential": "SiCGe.tersoff (Erhart-Albe 2005)", "reference_a": REF}
    res["Si"] = relax(0.0, threads)
    res["Ge"] = relax(1.0, threads, a0=5.65)
    res["Si0.5Ge0.5"] = relax(0.5, threads)
    res["benchmark"] = md_benchmark(threads)
    (outdir / "lammps_pilot.json").write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
