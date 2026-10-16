#!/usr/bin/env python3
"""P3 hidden-subset study: one MACE-MP Muller-Plathe rNEMD point on a SiGeC(+dopant)
random alloy (p3_mace_design.json). Same protocol and unit handling as
party4_mace/mp_rnemd_ase.py (the corrected version, ledger seq 1), extended to:
- carbon placed either dispersed (isolated substitutional C, no C-C neighbours)
  or clustered (C on one fcc sublattice around a centre, its neighbours kept Si,
  i.e. a sub-nanometre 3C-SiC cluster), then Ge and dopants on remaining sites;
- a short MACE hydrostatic cell + position relaxation before MD;
- species-matched velocity swaps for every species present (Si, Ge, C, B, P);
- checkpoints every CKPT steps, so a preempted guest job resumes.

Usage: mace_rnemd.py <point_json_string_or_id> --points points.json --out out/ID.json
"""
import argparse, json, os, random, time
from pathlib import Path
import numpy as np
from ase.build import bulk
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase import units

NSLAB, CKPT = 20, 5000


def build(p, ncx, ncy, ncz, seed):
    x, y = p["x_Ge"], p["y_C"]
    a = 5.431 * (1 - x) + 5.658 * x
    at = bulk("Si", "diamond", a=a, cubic=True) * (ncx, ncy, ncz)
    n = len(at)
    rng = random.Random(seed)
    pos = at.get_positions()
    syms = np.array(["Si"] * n, dtype=object)
    nC = round(y * n)
    fixed_si = set()
    if nC:
        d = at.get_all_distances(mic=True)
        nn = [set(np.where((row > 0.1) & (row < a * 0.45))[0]) for row in d]
        if p["carbon_form"] == "dispersed":
            idx = list(range(n)); rng.shuffle(idx)
            chosen = []
            for i in idx:
                if len(chosen) == nC:
                    break
                if not (nn[i] & set(chosen)):
                    chosen.append(i)
        else:  # clustered: sublattice-A sites nearest a random centre
            c = pos[rng.randrange(n)]
            sub = np.array([round((q / (a / 4)).sum()) % 2 == 0 for q in (pos - pos[0])])
            dc = np.linalg.norm(((pos - c + at.cell.diagonal() / 2) % at.cell.diagonal()) - at.cell.diagonal() / 2, axis=1)
            cand = [i for i in np.argsort(dc) if sub[i] == sub[np.argmin(dc)]]
            chosen = cand[:nC]
            for i in chosen:
                fixed_si |= nn[i]
        syms[chosen] = "C"
    free = [i for i in range(n) if syms[i] == "Si" and i not in fixed_si]
    rng.shuffle(free)
    nge = round(x * n)
    syms[free[:nge]] = "Ge"
    rest = free[nge:]
    if p["dopant"] != "none":
        nd = max(1, round(p["dopant_at_frac"] * n))
        syms[rest[:nd]] = p["dopant"]
    at.set_chemical_symbols(syms.tolist())
    return at


def slab_of(at):
    z = at.get_positions()[:, 2] % at.cell[2, 2]
    return (z / at.cell[2, 2] * NSLAB).astype(int) % NSLAB


def slab_temps(at, slabs):
    v, m = at.get_velocities(), at.get_masses()
    ke = 0.5 * m * (v ** 2).sum(axis=1)
    return np.array([ke[slabs == s].mean() * 2 / (3 * units.kB) for s in range(NSLAB)])


def mp_swap(at, slabs, species):
    v, m = at.get_velocities(), at.get_masses()
    ke = 0.5 * m * (v ** 2).sum(axis=1)
    syms = np.array(at.get_chemical_symbols())
    cold, hot = slabs == 0, slabs == NSLAB // 2
    best = None
    for sp in species:
        ci = np.where(cold & (syms == sp))[0]
        hi = np.where(hot & (syms == sp))[0]
        if not len(ci) or not len(hi):
            continue
        i, j = ci[np.argmax(ke[ci])], hi[np.argmin(ke[hi])]
        g = ke[i] - ke[j]
        if g > 0 and (best is None or g > best[2]):
            best = (i, j, g)
    if best is None:
        return 0.0
    i, j, g = best
    v[[i, j]] = v[[j, i]]
    at.set_velocities(v)
    return float(g)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pid")
    ap.add_argument("--points", default="points.json")
    ap.add_argument("--out", required=True)
    ap.add_argument("--ckpt-dir", default="ckpt")
    ap.add_argument("--bench", type=int, default=0, help="only time this many MD steps")
    a = ap.parse_args()
    P = json.load(open(a.points))
    S = P["protocol"]
    p = P["points"][a.pid]
    import torch
    from mace.calculators import mace_mp
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    if dev == "cuda":
        torch.zeros(1).cuda()  # fail fast on unusable guest GPUs
    calc = mace_mp(model=S["model"], device=dev, default_dtype="float32")
    ck = Path(a.ckpt_dir) / f"{a.pid}.npz"
    t0 = time.time()
    at = build(p, S["ncx"], S["ncy"], S["ncz"], p["seed"])
    at.calc = calc
    species = sorted(set(at.get_chemical_symbols()))
    state = {"phase": "relax", "done": 0, "transferred": 0.0, "profiles": np.zeros(NSLAB), "nprof": 0}
    if ck.exists():
        z = np.load(ck, allow_pickle=True)
        at.set_cell(z["cell"]); at.set_positions(z["pos"]); at.set_velocities(z["vel"])
        state = z["state"].item()
        print("resumed", state["phase"], state["done"], flush=True)

    def save():
        ck.parent.mkdir(parents=True, exist_ok=True)
        np.savez(str(ck) + ".tmp.npz", cell=at.cell[:], pos=at.get_positions(), vel=at.get_velocities(),
                 state=np.array(state, dtype=object))
        os.replace(str(ck) + ".tmp.npz", ck)

    if state["phase"] == "relax":
        from ase.filters import FrechetCellFilter
        from ase.optimize import FIRE
        FIRE(FrechetCellFilter(at, hydrostatic_strain=True), logfile=None).run(fmax=0.05, steps=S["relax_steps"])
        MaxwellBoltzmannDistribution(at, temperature_K=p["T_K"], rng=np.random.default_rng(p["seed"]))
        state["phase"] = "equil"; save()
    if a.bench:
        dyn = VelocityVerlet(at, S["dt_fs"] * units.fs)
        dyn.run(20); tb = time.time(); dyn.run(a.bench)
        print(json.dumps({"ms_per_step": (time.time() - tb) / a.bench * 1e3, "natoms": len(at),
                          "gpu": torch.cuda.get_device_name(0) if dev == "cuda" else "cpu"}))
        return
    slabs = slab_of(at)
    area, lz = at.cell[0, 0] * at.cell[1, 1], at.cell[2, 2]
    if state["phase"] == "equil":
        dyn = Langevin(at, S["dt_fs"] * units.fs, temperature_K=p["T_K"], friction=0.02)
        while state["done"] < S["equil"]:
            dyn.run(CKPT); state["done"] += CKPT; save()
        state.update(phase="prod", done=0); save()
    nve = VelocityVerlet(at, S["dt_fs"] * units.fs)
    discard = int(S["prod"] * 0.2)
    while state["done"] < S["prod"]:
        nve.run(S["swap_every"]); state["done"] += S["swap_every"]
        g = mp_swap(at, slabs, species)
        if state["done"] > discard:
            state["transferred"] += g; state["profiles"] += slab_temps(at, slabs); state["nprof"] += 1
        if state["done"] % CKPT == 0:
            save(); print(f"prod {state['done']}/{S['prod']} {time.time()-t0:.0f}s", flush=True)
    prof = state["profiles"] / max(state["nprof"], 1)
    zp = (np.arange(NSLAB) + 0.5) * lz / NSLAB
    h1, h2 = slice(1, NSLAB // 2), slice(NSLAB // 2 + 1, NSLAB)
    f1, f2 = np.polyfit(zp[h1], prof[h1], 1), np.polyfit(zp[h2], prof[h2], 1)
    r2 = []
    for h, f in ((h1, f1), (h2, f2)):
        res = prof[h] - np.polyval(f, zp[h]); r2.append(1 - res.var() / prof[h].var())
    grad = (abs(f1[0]) + abs(f2[0])) / 2
    t_prod = (S["prod"] - discard) * S["dt_fs"] * units.fs
    kappa = state["transferred"] / (2 * area * t_prod) / grad * units.fs * 1.602176634e-19 / 1e-10 / 1e-15
    out = {"pid": a.pid, **p, "kappa_W_mK": float(kappa), "grad_K_per_A": float(grad),
           "profile_R2": [float(v) for v in r2], "dT_half_K": float(grad * lz / 2),
           "natoms": len(at), "composition": {s: int((np.array(at.get_chemical_symbols()) == s).sum()) for s in species},
           "cell_A": [float(v) for v in at.cell.diagonal()], "profile_K": [float(v) for v in prof],
           "engine": f"ASE MP-rNEMD + mace-mp {S['model']} float32 ({torch.cuda.get_device_name(0) if dev == 'cuda' else 'cpu'})",
           "wall_s_this_attempt": time.time() - t0}
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(out, indent=1))
    print(json.dumps({k: out[k] for k in ("pid", "kappa_W_mK", "profile_R2")}))


if __name__ == "__main__":
    main()
