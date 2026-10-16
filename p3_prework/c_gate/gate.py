"""C-corner validation gate for the P3 Si-Ge-C ternary campaign.

Run via ./run_gate.sh (sets the MPICH library path the DGX venv LAMMPS needs).
Serial LAMMPS, one process, one OpenMP thread (well inside the 4-thread cap).
Writes results.json and logs/ in this directory. Uses only the stock
SiCGe.tersoff; no Ge-C parameters are constructed anywhere.

Tolerances are fixed HERE, before any number is computed:
  lattice constant within 1.5 % of the experimental reference,
  cohesive energy within 3 % of the experimental reference,
  MD stable (criteria in MD_CRIT below).
"""
import json, os, time, math, sys
import numpy as np
from lammps import lammps

HERE = os.path.dirname(os.path.abspath(__file__))
STOCK = os.path.expanduser(
    "~/mat_theory_pilot/venv/lib/python3.12/site-packages/lammps/share/lammps/potentials/SiCGe.tersoff")
LOG = os.path.join(HERE, "logs"); os.makedirs(LOG, exist_ok=True)

TOL_A, TOL_E = 0.015, 0.03
MD_CRIT = dict(T_mean_rel=0.05,       # |<T>-T0|/T0 over last 10 ps of NVT
               nve_drift_meV_per_atom=1.0,  # |E_end-E_start| over 5 ps NVE
               new_CC_pairs=0,        # no C-C bonds (< 2.0 A) formed during MD
               C_coord_change=0,      # every C keeps 4 neighbours within 2.6 A
               rms_disp_max_A=0.5)    # no melting / amorphisation
# references: (a in A, Ec in eV/atom). exp = gate reference; lit = Tersoff-model
# values quoted in the task brief (for context only, not gate criteria).
REF = {
 "Si_diamond": dict(exp=(5.431, -4.63), lit=(5.432, -4.63), elems="Si(D)", masses=[28.0855]),
 "Ge_diamond": dict(exp=(5.658, -3.85), lit=(5.657, -3.85), elems="Ge", masses=[72.63]),
 "C_diamond":  dict(exp=(3.567, -7.37), lit=(3.566, -7.37), elems="C", masses=[12.011]),
 "SiC_3C":     dict(exp=(4.36, -6.34), lit=(4.32, None), elems="Si(D) C", masses=[28.0855, 12.011]),
}
MASS = {"Si(D)": 28.0855, "Ge": 72.63, "C": 12.011}


def L(tag):
    return lammps(cmdargs=["-log", os.path.join(LOG, tag + ".log"), "-screen", "none"])


def setup_crystal(lmp, a, n, elems, zincblende=False):
    ntype = len(elems.split())
    lmp.commands_string(f"""
units metal
boundary p p p
atom_modify map array
lattice diamond {a}
region box block 0 {n} 0 {n} 0 {n}
create_box {ntype} box
""")
    if zincblende:
        lmp.command("create_atoms 1 box basis 5 2 basis 6 2 basis 7 2 basis 8 2")
    else:
        lmp.command("create_atoms 1 box")
    for t, e in enumerate(elems.split(), 1):
        lmp.command(f"mass {t} {MASS[e]}")


def crystal(name):
    r = REF[name]; elems = r["elems"]; zb = name == "SiC_3C"
    a0 = r["exp"][0]
    # (i) box/relax minimisation, 3x3x3 conventional cells (216 atoms), like chpc/in.validate
    lmp = L(name + "_relax")
    setup_crystal(lmp, a0 * 1.01, 3, elems, zb)
    lmp.commands_string(f"""
pair_style tersoff
pair_coeff * * {STOCK} {elems}
fix relax all box/relax iso 0.0
min_style cg
minimize 1e-14 1e-14 10000 100000
""")
    n = lmp.get_natoms()
    a_relax = lmp.extract_box()[1][0] / 3 - lmp.extract_box()[0][0] / 3
    e_relax = lmp.get_thermo("pe") / n
    lmp.close()
    # (ii) energy-volume scan, +-1.5 %, Birch-Murnaghan fit
    aa = np.linspace(0.985, 1.015, 13) * a_relax   # +-1.5 %: stays above the SiC cutoff kink
    ee = []
    for a in aa:
        lmp = L(name + "_ev")
        setup_crystal(lmp, a, 1, elems, zb)
        lmp.commands_string(f"pair_style tersoff\npair_coeff * * {STOCK} {elems}\nrun 0")
        ee.append(lmp.get_thermo("pe") / lmp.get_natoms()); lmp.close()
    ee = np.array(ee); V = aa ** 3 / 8.0  # volume per atom
    # 3rd-order BM: E = E0 + 9V0B0/16 [ (eta-1)^3 B' + (eta-1)^2 (6-4 eta) ], eta=(V0/V)^(2/3)
    from scipy.optimize import curve_fit
    def bm(V, E0, V0, B0, Bp):
        eta = (V0 / V) ** (2.0 / 3.0)
        return E0 + 9 * V0 * B0 / 16 * ((eta - 1) ** 3 * Bp + (eta - 1) ** 2 * (6 - 4 * eta))
    p, _ = curve_fit(bm, V, ee, p0=[ee.min(), V[np.argmin(ee)], 0.5, 4.0])
    a_ev = (8 * p[1]) ** (1 / 3)
    ea, ee_ = r["exp"]
    out = dict(a_relax=a_relax, Ecoh_relax=e_relax, a_EV=a_ev, Ecoh_EV=p[0],
               B0_GPa=p[2] * 160.21766, Bprime=p[3],
               ref_exp=dict(a=ea, Ecoh=ee_), ref_tersoff_lit=dict(a=r["lit"][0], Ecoh=r["lit"][1]),
               dev_a_rel=(a_relax - ea) / ea, dev_E_rel=(e_relax - ee_) / abs(ee_))
    out["pass_a"] = abs(out["dev_a_rel"]) <= TOL_A
    out["pass_E"] = abs(out["dev_E_rel"]) <= TOL_E
    return out


def supercell(xge, nC, seed, pot, relax=True, n=4, tag=None):
    """4x4x4 supercell (512 atoms) of Si1-yCy or Si1-xGex (stock file supports
    only these binaries), exact counts via type/subset."""
    assert not (xge and nC), "Si-Ge-C is not supported by SiCGe.tersoff"
    elems = "Si(D) Ge" if xge else "Si(D) C"
    tag = tag or f"sc_x{xge}_nC{nC}_s{seed}"
    lmp = L(tag)
    a_guess = 5.432 + (5.657 - 5.432) * xge
    setup_crystal(lmp, a_guess, n, elems)
    N = lmp.get_natoms()
    nGe = int(round(xge * N))
    ctype = 2
    if nGe:
        lmp.command(f"set type 1 type/subset 2 {nGe} {seed}")
    if nC:
        lmp.command(f"set type 1 type/subset {ctype} {nC} {seed + 7}")
    lmp.commands_string(f"pair_style tersoff\npair_coeff * * {pot} {elems}")
    if relax:
        lmp.commands_string("""
fix relax all box/relax iso 0.0 vmax 0.001
min_style cg
minimize 1e-14 1e-12 20000 200000
unfix relax
minimize 1e-14 1e-12 20000 200000
fix relax all box/relax iso 0.0 vmax 0.001
minimize 1e-14 1e-12 20000 200000
""")
    vol = lmp.get_thermo("vol")
    a = (vol ** (1 / 3)) / n
    res = dict(N=N, nGe=nGe, nC=nC, x=nGe / N, y=nC / N, a=a,
               pe_per_atom=lmp.get_thermo("pe") / N, press_bar=lmp.get_thermo("press"))
    return lmp, res, ctype


def c_env(lmp, ctype):
    """C-C pairs < 2.0 A and C coordination (neighbours within 2.6 A)."""
    x = np.array(lmp.numpy.extract_atom("x"))[:lmp.get_natoms()].copy()
    t = np.array(lmp.numpy.extract_atom("type"))[:lmp.get_natoms()].copy()
    lo, hi = lmp.extract_box()[0], lmp.extract_box()[1]
    Lb = np.array(hi) - np.array(lo)
    ci = np.where(t == ctype)[0]
    cc, coord, mind = 0, [], []
    for k, i in enumerate(ci):
        d = x - x[i]; d -= Lb * np.round(d / Lb)
        r = np.sqrt((d ** 2).sum(1)); r[i] = 99
        coord.append(int((r < 2.6).sum()))
        rc = r[ci]; mind.append(float(rc.min()) if len(ci) > 1 else None)
        cc += int((rc < 2.0).sum())
    return dict(CC_pairs=cc // 2, C_coord=coord,
                min_CC_dist=min([m for m in mind if m is not None], default=None)), x


def md(xge, nC, seed, pot, T, label):
    lmp, res, ctype = supercell(xge, nC, seed, pot, tag=f"md_{label}_T{T}")
    env0, x0 = c_env(lmp, ctype)
    lmp.commands_string(f"""
reset_timestep 0
timestep 0.0005
velocity all create {2*T} {seed} mom yes rot yes dist gaussian
fix nvt all nvt temp {T} {T} 0.1
thermo_style custom step temp pe etotal press
thermo 200
compute msd all msd
variable T equal temp
variable E equal etotal
variable P equal press
variable M equal c_msd[4]
fix av all print 200 "$(step) ${{T}} ${{E}} ${{P}} ${{M}}" file {LOG}/md_{label}_T{T}.dat screen no
run 40000
""")
    # 40000 x 0.5 fs = 20 ps NVT
    dat = np.loadtxt(f"{LOG}/md_{label}_T{T}.dat", comments="#")
    last = dat[dat[:, 0] >= 20000]
    Tm, Ts = last[:, 1].mean(), last[:, 1].std()
    Pm = last[:, 3].mean()
    rms_nvt = math.sqrt(dat[-1, 4])
    lmp.commands_string(f"""
unfix nvt
unfix av
fix nve all nve
fix av2 all print 100 "$(step) ${{T}} ${{E}}" file {LOG}/md_{label}_T{T}_nve.dat screen no
run 10000
""")
    nve = np.loadtxt(f"{LOG}/md_{label}_T{T}_nve.dat", comments="#")
    N = lmp.get_natoms()
    drift = (nve[-1, 2] - nve[0, 2]) / N * 1000.0
    fluct = nve[:, 2].std() / N * 1000.0
    env1, x1 = c_env(lmp, ctype)
    lmp.close()
    out = dict(T_target=T, NVT_ps=20.0, dt_fs=0.5, T_mean_last10ps=Tm, T_std_last10ps=Ts,
               P_mean_bar_last10ps=Pm, rms_disp_A=rms_nvt,
               NVE_ps=5.0, NVE_drift_meV_per_atom=drift, NVE_Etot_std_meV_per_atom=fluct,
               CC_pairs_start=env0["CC_pairs"], CC_pairs_end=env1["CC_pairs"],
               min_CC_dist_start=env0["min_CC_dist"], min_CC_dist_end=env1["min_CC_dist"],
               C_coord_start=env0["C_coord"], C_coord_end=env1["C_coord"], cell=res)
    out["stable"] = bool(abs(Tm - T) / T <= MD_CRIT["T_mean_rel"]
                         and abs(drift) <= MD_CRIT["nve_drift_meV_per_atom"]
                         and env1["CC_pairs"] - env0["CC_pairs"] <= MD_CRIT["new_CC_pairs"]
                         and env1["C_coord"] == env0["C_coord"]
                         and rms_nvt <= MD_CRIT["rms_disp_max_A"])
    return out


TERNARY_IN = """units metal
boundary p p p
lattice diamond 5.475
region box block 0 4 0 4 0 4
create_box 3 box
create_atoms 1 box
set type 1 type/subset 2 102 101
set type 1 type/subset 3 5 108
mass 1 28.0855
mass 2 72.63
mass 3 12.011
pair_style tersoff
pair_coeff * * {pot} Si(D) Ge C
run 0
"""


def ternary_attempt():
    """Try the planned Si0.8Ge0.2 + 1 % C cell with the stock file; capture LAMMPS's own error."""
    import subprocess
    inp = os.path.join(LOG, "in.ternary_attempt")
    open(inp, "w").write(TERNARY_IN.format(pot=STOCK))
    lmp_bin = os.path.expanduser("~/mat_theory_pilot/venv/bin/lmp")
    p = subprocess.run([lmp_bin, "-in", inp, "-log", os.path.join(LOG, "ternary_attempt.log")],
                       capture_output=True, text=True)
    err = [l for l in (p.stdout + p.stderr).splitlines() if "ERROR" in l]
    return dict(input=os.path.relpath(inp, HERE), exit_code=p.returncode,
                lammps_error=err[0] if err else None, ran=p.returncode == 0)


def main():
    t0 = time.time()
    R = dict(meta=dict(date=time.strftime("%Y-%m-%d %H:%M"), potential=STOCK,
                       lammps="LAMMPS 22 Jul 2025 update 4 (pip wheel, aarch64), serial, 1 thread",
                       tolerances=dict(a_rel=TOL_A, Ecoh_rel=TOL_E, md=MD_CRIT)))
    R["ternary_attempt"] = ternary_attempt()
    print("ternary attempt:", R["ternary_attempt"], flush=True)

    R["crystals"] = {k: crystal(k) for k in REF}
    for k, v in R["crystals"].items():
        print(k, {q: v[q] for q in ("a_relax", "Ecoh_relax", "a_EV", "B0_GPa")}, flush=True)
    # SiC E(a) on a wide grid: documents the Si-Si second-neighbour cutoff kink
    scan = []
    for a in np.round(np.arange(4.15, 4.451, 0.01), 3):
        lmp = L("SiC_wide_scan"); setup_crystal(lmp, a, 1, "Si(D) C", True)
        lmp.commands_string(f"pair_style tersoff\npair_coeff * * {STOCK} Si(D) C\nrun 0")
        scan.append([float(a), lmp.get_thermo("pe") / 8]); lmp.close()
    R["SiC_wide_EV_scan"] = dict(a_E=scan, kink_lattice_constant_A=3.0 * math.sqrt(2),
        note="Si-Si second neighbours sit at a/sqrt(2); the Tersoff Si-Si cutoff S=3.0 A is "
             "entered below a=4.243 A, producing a kink 1.8 % below equilibrium")

    aSi, aC = R["crystals"]["Si_diamond"]["a_relax"], R["crystals"]["C_diamond"]["a_relax"]
    aGe = R["crystals"]["Ge_diamond"]["a_relax"]
    # substitutional C in Si, 3 seeds per concentration
    sub = []
    for nC in (0, 5, 10, 20):
        for seed in ((101,) if nC == 0 else (101, 202, 303)):
            lmp, res, _ = supercell(0.0, nC, seed, STOCK); lmp.close()
            res["seed"] = seed; res["a_vegard"] = (1 - res["y"]) * aSi + res["y"] * aC
            sub.append(res); print("Si:C", nC, seed, res["a"], flush=True)
    ys = np.array([s["y"] for s in sub]); As = np.array([s["a"] for s in sub])
    slope = float(np.polyfit(ys, As, 1)[0])
    # dilute SiGe for da/dx (needed for the Ge:C compensation ratio)
    sige = []
    for x in (0.1, 0.2):
        for seed in (101, 202, 303):
            lmp, res, _ = supercell(x, 0, seed, STOCK); lmp.close()
            res["seed"] = seed; sige.append(res); print("SiGe", x, seed, res["a"], flush=True)
    xs = np.array([0.0] + [s["x"] for s in sige]); Ax = np.array([aSi] + [s["a"] for s in sige])
    dadx = float(np.polyfit(xs, Ax, 1)[0])
    R["SiC_substitutional"] = dict(points=sub, slope_dady_A=slope, vegard_slope_A=aC - aSi,
                                   ratio_to_vegard=slope / (aC - aSi))
    R["SiGe_dilute"] = dict(points=sige, slope_dadx_A=dadx, vegard_slope_A=aGe - aSi)
    R["compensation_ratio"] = dict(
        nu_tersoff_additive=-slope / dadx, nu_vegard_tersoff=-(aC - aSi) / (aGe - aSi),
        nu_experiment_DeSalvador2000=12.0,
        note="additive estimate from the two supported binaries (no Ge-C bonds in the model); "
             "experiment: De Salvador et al., PRB 61, 13005 (2000)")
    # MD stability, Si:C 1 % (stock file, supported binary)
    R["MD"] = {}
    for T in (300, 1000):
        R["MD"][f"Si_C1pct_T{T}"] = md(0.0, 5, 101, STOCK, T, "SiC1")
        v = R["MD"][f"Si_C1pct_T{T}"]
        print("MD", T, v["T_mean_last10ps"], v["NVE_drift_meV_per_atom"], v["stable"], flush=True)

    c = R["crystals"]
    supported_pass = all(c[k]["pass_a"] and c[k]["pass_E"] for k in c) and all(
        v["stable"] for v in R["MD"].values())
    R["verdict"] = dict(supported_parts="PASS" if supported_pass else "FAIL",
                        ternary_as_planned="FAIL" if not R["ternary_attempt"]["ran"] else "RUNS")
    R["meta"]["wall_s"] = time.time() - t0
    json.dump(R, open(os.path.join(HERE, "results.json"), "w"), indent=1, default=float)
    print("verdict", R["verdict"], "wall", R["meta"]["wall_s"])


if __name__ == "__main__":
    main()
