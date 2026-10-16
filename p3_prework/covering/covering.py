"""Party-regime assignment as cost-aware covering: worked example on the SiGe
evidence hypergraph (P3 claim 4).

Problem (weighted set multicover with capability and diversity side
constraints). Regime cells c = (x column, T band). Every cell must end up
covered by >= 2 distinct evidence parties drawn from >= 2 distinct method
families, so a cross-party agreement check is possible in every cell.
Fixed parties (literature, DB-CITRINE, DB-STARRY) contribute whatever they already
cover at zero cost and cannot be asked for more. Askable parties (CHPC-MD,
DGX-MD, MACE-MD, BTE-DFT, GPUMD-HNEMD) can be assigned to any cell inside
their validity domain at a measured per-point cost. BTE-DFT is bought per
composition run (one run = one x column, up to 3 temperature bands).
Minimise total cost. Solved exactly with scipy.optimize.milp and with a
greedy ratio heuristic; compared with what the campaign actually ran.

Run: python covering.py   (writes results.json and covering_grid.png here)
"""
from __future__ import annotations

import json
import os
from collections import defaultdict

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp

HERE = os.path.dirname(os.path.abspath(__file__))
MT = os.path.abspath(os.path.join(HERE, "..", ".."))
KG = os.path.join(MT, "pilot", "kg", "kg_hypergraph.json")

# ---------------------------------------------------------------- the grid
# x columns: the hypergraph's 0.2-wide bands, with the pure endpoints split
# out because that is where the validity boundaries of the parties sit.
XCOLS = ["Si", "x0", "x1", "x2", "x3", "x4", "Ge"]
XLABEL = {"Si": "pure Si\nx<0.05", "x0": "0.05-0.2", "x1": "0.2-0.4", "x2": "0.4-0.6",
          "x3": "0.6-0.8", "x4": "0.8-0.95", "Ge": "pure Ge\nx>0.95"}
# T bands: the hypergraph's own banding, int((T-300)/200) clipped to 0..4
TBANDS = ["T0", "T1", "T2", "T3", "T4"]
TLABEL = {"T0": "300-500 K", "T1": "500-700 K", "T2": "700-900 K", "T3": "900-1100 K",
          "T4": ">=1100 K"}
CELLS = [(xc, tb) for tb in TBANDS for xc in XCOLS]


def xcol(x):
    if x < 0.05:
        return "Si"
    if x > 0.95:
        return "Ge"
    return f"x{min(int(x * 5), 4)}"


def tband(T):
    return f"T{min(max(int((T - 300) / 200), 0), 4)}"


# ---------------------------------------------------------------- parties
FAMILY = {  # method family of each askable party
    "CHPC-MD": "Tersoff-MD", "DGX-MD": "Tersoff-MD", "GPUMD-HNEMD": "Tersoff-MD",
    "MACE-MD": "MLIP-MD", "BTE-DFT": "DFT-BTE",
}
LIT_FAMILY = {"experiment": "experiment", "MD": "other-MD", "BTE/DFT": "DFT-BTE"}
# model-fit and unclear literature points are not counted as coverage (a fit
# or an unattributed number cannot anchor an agreement check).

# Cost per point in wall-clock minutes of the party's own resource (midpoint
# of the measured range). Sources are listed in README.md.
COST_BASE = {"CHPC-MD": 23.0, "DGX-MD": 12.5, "MACE-MD": 246.0,
             "BTE-DFT": 75.0, "GPUMD-HNEMD": 20.0}
# Fidelity-weighted: what one *usable* cell estimate costs in practice.
# rNEMD needs >= 3 box lengths for a 1/L extrapolation; HNEMD needs >= 2
# driving forces x 2 seeds for the linearity gate (production protocol);
# BTE's 75 min is already a converged bulk value.
COST_FIDELITY = {"CHPC-MD": 3 * 23.0, "DGX-MD": 3 * 12.5, "MACE-MD": 3 * 246.0,
                 "BTE-DFT": 75.0, "GPUMD-HNEMD": 4 * 20.0}
ASKABLE = list(COST_BASE)


def valid(party, cell, scen):
    xc, tb = cell
    endpoint = xc in ("Si", "Ge")
    if party in ("CHPC-MD", "DGX-MD", "MACE-MD"):   # all three are rNEMD
        return not endpoint
    if party == "BTE-DFT":
        return scen.get("bte_T4", False) or tb != "T4"
    if party == "GPUMD-HNEMD":
        if scen.get("gpumd_measured", False) and endpoint:
            return tb == "T4"      # production campaign: linearity gate passed only at 1100 K
        return True
    raise KeyError(party)


SCENARIOS = {
    "A_base": dict(desc="task-brief capabilities, per-point midpoint costs", cost=COST_BASE),
    "B_fidelity": dict(desc="as A, costs per usable estimate (rNEMD x3 lengths, HNEMD x4 runs)",
                       cost=COST_FIDELITY),
    "C_measured_gpumd": dict(desc="as A, but GPUMD endpoint validity as measured by the "
                                  "pre-registered campaign (linearity gate passed only at 1100 K)",
                             cost=COST_BASE, gpumd_measured=True),
    "D_bte_T4_ok": dict(desc="as A, but BTE-DFT accepted at >=1100 K", cost=COST_BASE, bte_T4=True),
}


# ---------------------------------------------------------------- data
def load_fixed_and_actual():
    kg = json.load(open(KG))
    fixed = defaultdict(lambda: defaultdict(set))   # cell -> party -> families
    actual_pts = defaultdict(lambda: defaultdict(int))  # party -> cell -> n points
    for h in kg["hyperedges"]:
        m = dict(s.split(":", 1) for s in h["members"])
        x, T = h.get("x_ge"), h.get("T_K")
        if x is None or T is None:
            continue
        cell = (xcol(float(x)), tband(float(T)))
        paper = m["paper"]
        corpus = str(h.get("provenance", {}).get("corpusId", ""))
        if paper in ASKABLE:
            actual_pts[paper][cell] += 1
        elif paper == "DB-CITRINE":
            fixed[cell]["DB-CITRINE"].add("experiment")
        elif corpus == "DB-STARRY":
            # Starrydata rows name their upstream paper as the paper member;
            # the source is the database slice (one method family plus one
            # provenance chain, the counting rule of all three papers).
            fixed[cell]["DB-STARRY"].add("experiment")
        else:
            fam = LIT_FAMILY.get(m["method"])
            if fam:
                fixed[cell]["LIT"].add(fam)
    # GPUMD HNEMD (not in the hypergraph): validation gate (5 runs, 300 K:
    # Si x3, Ge, SiGe x=0.5) + pre-registered production campaign (33 runs)
    gp = actual_pts["GPUMD-HNEMD"]
    gp[("Si", "T0")] += 3; gp[("Ge", "T0")] += 1; gp[("x2", "T0")] += 1
    design = json.load(open(os.path.join(MT, "gpumd", "production", "endpoint_campaign_design.json")))
    for r in design["runs"]:
        xg = float(r["x_ge"]) if "x_ge" in r else {"Si": 0.0, "Ge": 1.0}[r.get("element") or r["system"]]
        gp[(xcol(xg), tband(r["T_K"]))] += 1
    return fixed, actual_pts


# ---------------------------------------------------------------- model
def build_options(scen):
    """Options = purchasable units. Each option: (party, cells covered, cost).
    Point parties: one option per valid cell. BTE: handled by milp with
    per-column run counts; for greedy, one option per (column, set of up to 3
    valid T bands)."""
    cost = scen["cost"]
    opts = []
    for p in ASKABLE:
        if p == "BTE-DFT":
            continue
        for c in CELLS:
            if valid(p, c, scen):
                opts.append((p, (c,), cost[p]))
    return opts


def coverage_state(fixed, chosen):
    """chosen: iterable of (party, cell). Returns per-cell parties and families."""
    parties = defaultdict(set); fams = defaultdict(set)
    for c in CELLS:
        for p, f in fixed.get(c, {}).items():
            parties[c].add(p); fams[c] |= f
    for p, c in chosen:
        parties[c].add(p); fams[c].add(FAMILY[p])
    return parties, fams


def requirements(fixed, scen):
    """Per-cell demand (parties, families), capped at what is achievable."""
    everything = [(p, c) for p in ASKABLE for c in CELLS if valid(p, c, scen)]
    parties, fams = coverage_state(fixed, everything)
    req = {c: (min(2, len(parties[c])), min(2, len(fams[c]))) for c in CELLS}
    infeasible = [c for c in CELLS if req[c] != (2, 2)]
    return req, infeasible


def solve_milp(fixed, scen, restrict_cells=None):
    cost = scen["cost"]
    req, infeasible = requirements(fixed, scen)
    target = CELLS if restrict_cells is None else restrict_cells
    # variables: y[p,c] for point parties and BTE per-cell; n[col] BTE runs
    var = []
    for p in ASKABLE:
        for c in CELLS:
            if valid(p, c, scen):
                var.append(("y", p, c))
    for xc in XCOLS:
        var.append(("n", "BTE-DFT", xc))
    fam_names = sorted(set(FAMILY.values()) | set(LIT_FAMILY.values()) | {"experiment"})
    for c in target:
        for F in fam_names:
            var.append(("f", F, c))
    idx = {v: i for i, v in enumerate(var)}
    nv = len(var)
    cvec = np.zeros(nv)
    for v, i in idx.items():
        if v[0] == "y" and v[1] != "BTE-DFT":
            cvec[i] = cost[v[1]]
        elif v[0] == "n":
            cvec[i] = cost["BTE-DFT"]
    rows, lo, hi = [], [], []

    def add(coefs, l, h):
        r = np.zeros(nv)
        for k, a in coefs.items():
            r[idx[k]] += a
        rows.append(r); lo.append(l); hi.append(h)

    # BTE: covered T bands in a column <= 3 * runs
    for xc in XCOLS:
        co = {("y", "BTE-DFT", (xc, tb)): 1.0 for tb in TBANDS if ("y", "BTE-DFT", (xc, tb)) in idx}
        co[("n", "BTE-DFT", xc)] = -3.0
        add(co, -np.inf, 0.0)
    for c in target:
        fx_parties = set(fixed.get(c, {}))
        fx_fams = set().union(*fixed[c].values()) if c in fixed and fixed[c] else set()
        rp, rf = req[c]
        # parties: sum y + |fixed| >= rp
        add({("y", p, c): 1.0 for p in ASKABLE if ("y", p, c) in idx}, rp - len(fx_parties), np.inf)
        # families: f[F,c] <= fixed_F + sum_{p in F} y[p,c]; sum_F f >= rf
        for F in fam_names:
            co = {("f", F, c): 1.0}
            for p in ASKABLE:
                if FAMILY[p] == F and ("y", p, c) in idx:
                    co[("y", p, c)] = -1.0
            add(co, -np.inf, 1.0 if F in fx_fams else 0.0)
        add({("f", F, c): 1.0 for F in fam_names}, rf, np.inf)
    ub = np.ones(nv); ub[[idx[("n", "BTE-DFT", xc)] for xc in XCOLS]] = 2
    res = milp(cvec, constraints=LinearConstraint(np.array(rows), lo, hi),
               integrality=np.ones(nv), bounds=Bounds(np.zeros(nv), ub))
    assert res.success, res.message
    x = np.round(res.x).astype(int)
    chosen = [(v[1], v[2]) for v, i in idx.items() if v[0] == "y" and x[i]]
    runs = {v[2]: int(x[i]) for v, i in idx.items() if v[0] == "n" and x[i]}
    return dict(cost_min=float(res.fun), chosen=chosen, bte_runs=runs, infeasible=infeasible)


def deficit(fixed, chosen, req, cells):
    parties, fams = coverage_state(fixed, chosen)
    return sum(max(0, req[c][0] - len(parties[c])) + max(0, req[c][1] - len(fams[c]))
               for c in cells)


def solve_greedy(fixed, scen, restrict_cells=None):
    """Classic ratio greedy: repeatedly buy the option with the largest
    (deficit reduction) / cost until every target cell meets its demand.
    BTE options: one composition run covering the (up to) 3 valid T bands of
    a column that reduce deficit most."""
    import itertools
    cost = scen["cost"]
    req, infeasible = requirements(fixed, scen)
    target = CELLS if restrict_cells is None else restrict_cells
    opts = build_options(scen)
    for xc in XCOLS:
        vb = [tb for tb in TBANDS if valid("BTE-DFT", (xc, tb), scen)]
        for k in range(1, 4):
            for combo in itertools.combinations(vb, k):
                opts.append(("BTE-DFT", tuple((xc, tb) for tb in combo), cost["BTE-DFT"]))
    chosen, spent, bte_runs = [], 0.0, defaultdict(int)
    while True:
        d0 = deficit(fixed, chosen, req, target)
        if d0 == 0:
            break
        best = None
        for p, cells, cst in opts:
            new = [(p, c) for c in cells if (p, c) not in chosen]
            if not new:
                continue
            gain = d0 - deficit(fixed, chosen + new, req, target)
            if gain > 0 and (best is None or gain / cst > best[0] + 1e-12):
                best = (gain / cst, p, new, cst)
        assert best is not None, "greedy stuck"
        _, p, new, cst = best
        chosen += new; spent += cst
        if p == "BTE-DFT":
            bte_runs[new[0][1][0]] += 1
    return dict(cost=spent, chosen=chosen, bte_runs=dict(bte_runs))


def summarize_assignment(chosen):
    by_party = defaultdict(list)
    for p, c in chosen:
        by_party[p].append(f"{c[0]}/{c[1]}")
    return {p: sorted(v) for p, v in by_party.items()}


def actual_campaign(fixed, actual_pts, scen):
    cost = scen["cost"]
    # cost of what was actually run, at base per-point costs (BTE per composition)
    n_bte_comps = 6  # six compositions x 3 temperatures (party4_bte/RESULTS.md)
    spent = sum(COST_BASE[p] * sum(v.values()) for p, v in actual_pts.items() if p != "BTE-DFT")
    spent += COST_BASE["BTE-DFT"] * n_bte_comps
    npts = {p: int(sum(v.values())) for p, v in actual_pts.items()}
    # coverage the campaign earned under the rule (only valid measurements count)
    chosen = [(p, c) for p, v in actual_pts.items() for c in v if valid(p, c, scen)]
    invalid = sorted({(p, f"{c[0]}/{c[1]}") for p, v in actual_pts.items() for c in v
                      if not valid(p, c, scen)})
    parties, fams = coverage_state(fixed, chosen)
    satisfied = [c for c in CELLS if len(parties[c]) >= 2 and len(fams[c]) >= 2]
    return dict(spent_min=spent, points=npts, chosen=chosen, invalid_points=invalid,
                satisfied=satisfied)


# ---------------------------------------------------------------- figure
def figure(fixed, sol, act, scen, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    colors = {"CHPC-MD": "#1f77b4", "DGX-MD": "#9ecae1", "GPUMD-HNEMD": "#08306b",
              "MACE-MD": "#d95f02", "BTE-DFT": "#1b9e77", "LIT": "#7f7f7f", "DB-CITRINE": "#bdbdbd",
              "DB-STARRY": "#d9d9d9"}
    short = {"CHPC-MD": "CHPC", "DGX-MD": "DGX", "GPUMD-HNEMD": "GPUMD", "MACE-MD": "MACE",
             "BTE-DFT": "BTE", "LIT": "LIT", "DB-CITRINE": "CIT",
             "DB-STARRY": "STAR"}
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.6), constrained_layout=True)
    panels = [("What the campaign ran (valid measurements only)", act["chosen"], act["satisfied"]),
              (f"Optimal assignment, scenario A ({sol['cost_min']/60:.1f} h)", sol["chosen"], None)]
    for ax, (title, chosen, satisfied) in zip(axes, panels):
        parties, fams = coverage_state(fixed, chosen)
        for j, tb in enumerate(TBANDS):
            for i, xc in enumerate(XCOLS):
                c = (xc, tb)
                ok = len(parties[c]) >= 2 and len(fams[c]) >= 2
                ax.add_patch(Rectangle((i, j), 1, 1, facecolor="#f7f7f7" if ok else "#fde0dd",
                                       edgecolor="white", lw=2))
                plist = sorted(parties[c], key=lambda p: list(colors).index(p))
                w = 0.84 / max(len(plist), 1)
                for k, p in enumerate(plist):
                    ax.add_patch(Rectangle((i + 0.08 + k * w, j + 0.42), w * 0.92, 0.3,
                                           facecolor=colors[p], edgecolor="none"))
                    ax.text(i + 0.08 + (k + 0.5) * w, j + 0.25, short[p], ha="center",
                            va="center", fontsize=6.2, rotation=0 if len(plist) < 4 else 90,
                            color="#333333")
                if not ok:
                    lab = "infeasible" if c in sol["infeasible"] else "not met"
                    ax.text(i + 0.5, j + 0.86, lab, ha="center", va="center",
                            fontsize=7, color="#b2182b")
        ax.set_xlim(0, len(XCOLS)); ax.set_ylim(0, len(TBANDS))
        ax.set_xticks(np.arange(len(XCOLS)) + 0.5, [XLABEL[x] for x in XCOLS], fontsize=8)
        ax.set_yticks(np.arange(len(TBANDS)) + 0.5, [TLABEL[t] for t in TBANDS], fontsize=8)
        ax.set_xlabel("Ge fraction x"); ax.set_title(title, fontsize=11)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.tick_params(length=0)
    handles = [Rectangle((0, 0), 1, 1, color=colors[p]) for p in colors]
    fig.legend(handles, [f"{short[p]} ({p})" if short[p] != p else p for p in colors],
               loc="outside lower center", ncol=7, fontsize=8, frameon=False)
    fig.suptitle("Party-regime covering on the SiGe evidence grid: every cell needs >= 2 parties "
                 "from >= 2 method families (pink = requirement not met)", fontsize=11)
    fig.savefig(path, dpi=160)


# ---------------------------------------------------------------- main
def main():
    fixed, actual_pts = load_fixed_and_actual()
    out = dict(grid=dict(x_columns=XCOLS, T_bands=TBANDS, n_cells=len(CELLS)),
               families=FAMILY, lit_families=LIT_FAMILY, cost_base_min=COST_BASE,
               cost_fidelity_min=COST_FIDELITY,
               fixed_coverage={f"{c[0]}/{c[1]}": {p: sorted(f) for p, f in v.items()}
                               for c, v in fixed.items() if v},
               scenarios={})
    for name, scen in SCENARIOS.items():
        opt = solve_milp(fixed, scen)
        gr = solve_greedy(fixed, scen)
        act = actual_campaign(fixed, actual_pts, scen)
        # optimal / greedy cost to earn exactly the cells the campaign earned
        opt_same = solve_milp(fixed, scen, restrict_cells=act["satisfied"])
        gr_same = solve_greedy(fixed, scen, restrict_cells=act["satisfied"])
        # cost to complete the campaign: optimal top-up given what it ran
        req, _ = requirements(fixed, scen)
        fixed_plus = defaultdict(lambda: defaultdict(set))
        for c, v in fixed.items():
            for p, f in v.items():
                fixed_plus[c][p] |= f
        for p, c in act["chosen"]:
            fixed_plus[c][p].add(FAMILY[p])
        topup = solve_milp(fixed_plus, scen)
        missed = [c for c in CELLS if c not in act["satisfied"] and c not in opt["infeasible"]]
        out["scenarios"][name] = dict(
            description=scen["desc"],
            infeasible_cells=[f"{c[0]}/{c[1]}" for c in opt["infeasible"]],
            optimal=dict(cost_min=opt["cost_min"], cost_h=opt["cost_min"] / 60,
                         n_assignments=len(opt["chosen"]), bte_runs=opt["bte_runs"],
                         assignment=summarize_assignment(opt["chosen"])),
            greedy=dict(cost_min=gr["cost"], cost_h=gr["cost"] / 60,
                        gap_rel=gr["cost"] / opt["cost_min"] - 1 if opt["cost_min"] else 0.0,
                        assignment=summarize_assignment(gr["chosen"])),
            actual=dict(spent_min_at_base_costs=act["spent_min"], spent_h=act["spent_min"] / 60,
                        points=act["points"], invalid_points=act["invalid_points"],
                        cells_satisfied=[f"{c[0]}/{c[1]}" for c in act["satisfied"]],
                        n_cells_satisfied=len(act["satisfied"]),
                        cells_missed_but_feasible=[f"{c[0]}/{c[1]}" for c in missed]),
            optimal_for_same_cells_as_campaign=dict(cost_min=opt_same["cost_min"],
                                                    greedy_cost_min=gr_same["cost"],
                                                    assignment=summarize_assignment(opt_same["chosen"])),
            topup_to_full_grid=dict(cost_min=topup["cost_min"],
                                    assignment=summarize_assignment(topup["chosen"]),
                                    bte_runs=topup["bte_runs"]),
        )
        s = out["scenarios"][name]
        print(f"{name}: opt {s['optimal']['cost_h']:.2f} h, greedy {s['greedy']['cost_h']:.2f} h "
              f"(gap {100*s['greedy']['gap_rel']:.1f} %), actual {s['actual']['spent_h']:.1f} h "
              f"satisfying {s['actual']['n_cells_satisfied']} cells; opt for same cells "
              f"{opt_same['cost_min']/60:.2f} h; top-up {topup['cost_min']/60:.2f} h; "
              f"infeasible {s['infeasible_cells']}")
        if name == "A_base":
            figure(fixed, opt, act, scen, os.path.join(HERE, "covering_grid.png"))
    # Scenario E: only the campaign's own temperature bands (300, 700, 1100 K
    # land in T0, T2, T4), base costs. Fairer comparison: the campaign never
    # aimed at the 500-700 K and 900-1100 K bands.
    scen = SCENARIOS["A_base"]
    tgt = [c for c in CELLS if c[1] in ("T0", "T2", "T4")]
    opt = solve_milp(fixed, scen, restrict_cells=tgt)
    gr = solve_greedy(fixed, scen, restrict_cells=tgt)
    act = actual_campaign(fixed, actual_pts, scen)
    fixed_plus = defaultdict(lambda: defaultdict(set))
    for c, v in fixed.items():
        for p, f in v.items():
            fixed_plus[c][p] |= f
    for p, c in act["chosen"]:
        fixed_plus[c][p].add(FAMILY[p])
    topup = solve_milp(fixed_plus, scen, restrict_cells=tgt)
    sat = [c for c in act["satisfied"] if c in tgt]
    out["scenarios"]["E_campaign_T_bands"] = dict(
        description="as A, but only the T0/T2/T4 bands the campaign targeted (21 cells)",
        infeasible_cells=[f"{c[0]}/{c[1]}" for c in opt["infeasible"] if c in tgt],
        optimal=dict(cost_min=opt["cost_min"], cost_h=opt["cost_min"] / 60,
                     bte_runs=opt["bte_runs"], assignment=summarize_assignment(opt["chosen"])),
        greedy=dict(cost_min=gr["cost"], cost_h=gr["cost"] / 60,
                    gap_rel=gr["cost"] / opt["cost_min"] - 1),
        actual=dict(spent_h=act["spent_min"] / 60, n_cells_satisfied=len(sat),
                    cells_satisfied=[f"{c[0]}/{c[1]}" for c in sat]),
        topup_to_full_grid=dict(cost_min=topup["cost_min"],
                                assignment=summarize_assignment(topup["chosen"]),
                                bte_runs=topup["bte_runs"]))
    e = out["scenarios"]["E_campaign_T_bands"]
    print(f"E: opt {e['optimal']['cost_h']:.2f} h, greedy {e['greedy']['cost_h']:.2f} h "
          f"(gap {100*e['greedy']['gap_rel']:.1f} %), actual satisfies {len(sat)}/21, "
          f"top-up {topup['cost_min']/60:.2f} h, infeasible {e['infeasible_cells']}")
    json.dump(out, open(os.path.join(HERE, "results.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
