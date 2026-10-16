#!/usr/bin/env python3
"""Stage M1c (m1c_design.json): one converged MACE-MP-0 relaxation per task,
checkpointed so a preempted guest job resumes where it stopped.
Usage: mace_relax_m1c.py <task_index>; writes out/<task>.json."""
import json, os, sys, time
import numpy as np
from ase.build import bulk
from ase.filters import FrechetCellFilter
from ase.io import read, write
from ase.neighborlist import neighbor_list
from ase.optimize import FIRE, LBFGS
from ase.units import GPa
from mace.calculators import mace_mp

NREP, CAP, S1CAP, F1, FGATE, PGATE = 3, 5000, 800, 0.05, 0.01, 0.1


def task_list():
    t = ["bulk_Si", "bulk_Ge", "bulk_SiC"]
    for x in ("0", "0.2", "0.5"):
        for s in (1, 2, 3):
            labs = ["nC0", "nC1", "nC2", "nC4"] + (["maxGe", "minGe"] if x != "0" else [])
            t += [f"{x}_{s}_{l}" for l in labs]
    return t


def build(task):
    if task.startswith("bulk_"):
        name = task[5:]
        a0, kind = {"Si": (5.431, "diamond"), "Ge": (5.658, "diamond"), "SiC": (4.360, "zincblende")}[name]
        return bulk(name, kind, a=a0, cubic=True), {}
    xs, ss, lab = task.split("_")
    x, seed = float(xs), int(ss)
    a_host = 5.431 * (1 - x) + 5.658 * x
    host = bulk("Si", "diamond", a=a_host, cubic=True).repeat((NREP,) * 3)
    n = len(host)
    rng = np.random.default_rng(seed)
    ge = rng.choice(n, int(round(x * n)), replace=False)
    sym = np.array(["Si"] * n, dtype=object)
    sym[ge] = "Ge"
    host.set_chemical_symbols(list(sym))
    si_sites = np.where(sym == "Si")[0]
    meta = {"x": x, "seed": seed, "n_sites": n}
    s = list(sym)
    if lab.startswith("nC"):
        nc = int(lab[2:])
        for i in np.random.default_rng(1000 + seed).permutation(si_sites)[:nc]:
            s[i] = "C"
        meta["n_C"] = nc
    else:  # site preference, identical selection to M1
        i_, j_ = neighbor_list("ij", host, 0.5 * np.sqrt(3) / 2 * a_host + 0.3)
        nge_nn = np.zeros(n, int)
        for a, b in zip(i_, j_):
            if sym[b] == "Ge":
                nge_nn[a] += 1
        cand = si_sites
        if lab == "maxGe":
            idx = int(cand[np.argmax(nge_nn[cand])])
        else:
            lo_c = cand[nge_nn[cand] == 0]
            idx = int(lo_c[0]) if len(lo_c) else int(cand[np.argmin(nge_nn[cand])])
        s[idx] = "C"
        meta["n_Ge_nn"] = int(nge_nn[idx])
    host.set_chemical_symbols(s)
    return host, meta


def main(i):
    task = task_list()[i]
    os.makedirs("out", exist_ok=True); os.makedirs("ckpt", exist_ok=True)
    ck, st = f"ckpt/{task}.extxyz", f"ckpt/{task}.state.json"
    atoms, meta = build(task)
    state = {"steps": 0, "stage": 1, "wall_s": 0.0, "restarts": 0}
    if os.path.exists(ck) and os.path.exists(st):
        prev = read(ck)
        assert list(prev.get_chemical_symbols()) == list(atoms.get_chemical_symbols())
        atoms = prev
        state = json.load(open(st)); state["restarts"] += 1
    atoms.calc = mace_mp(model="medium", default_dtype="float64", device="cpu")
    filt = FrechetCellFilter(atoms, hydrostatic_strain=True)
    t0 = time.time()

    def save():
        write(ck, atoms)
        json.dump({**state, "wall_s": state["wall_s"] + time.time() - t0}, open(st, "w"))

    def run(opt, fmax, cap):
        while state["steps"] < cap:
            n0 = opt.nsteps
            done = opt.run(fmax=fmax, steps=min(25, cap - state["steps"]))
            state["steps"] += opt.nsteps - n0
            save()
            if done or opt.nsteps == n0:
                return bool(done)
        return False

    if state["stage"] == 1:
        run(FIRE(filt, logfile=None), F1, S1CAP)
        state["stage"] = 2; save()
    conv = run(LBFGS(filt, maxstep=0.1, logfile=None), 0.005, CAP)
    f = float(np.linalg.norm(atoms.get_forces(), axis=1).max())
    fcomp = float(np.abs(atoms.get_forces()).max())
    p = float(-np.trace(atoms.get_stress(voigt=False)) / 3 / GPa)
    n = len(atoms)
    res = {"task": task, "model": "MACE-MP-0 medium float64", **meta,
           "E_eV": float(atoms.get_potential_energy()), "V_A3": float(atoms.get_volume()),
           "a_A": float((atoms.get_volume() / (NREP ** 3 if n > 8 else 1)) ** (1 / 3)),
           "n_atoms": n, "fmax_norm": f, "fmax_component": fcomp, "pressure_GPa": p,
           "optimizer_converged": conv, "steps": state["steps"], "restarts": state["restarts"],
           "gate_pass": bool(f <= FGATE and abs(p) <= PGATE),
           "wall_s": state["wall_s"] + time.time() - t0}
    json.dump(res, open(f"out/{task}.json", "w"), indent=1)
    print("done", task, res["gate_pass"], f, p, res["steps"])


if __name__ == "__main__":
    main(int(sys.argv[1]))
