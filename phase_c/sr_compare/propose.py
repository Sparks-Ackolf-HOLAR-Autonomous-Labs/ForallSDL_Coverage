#!/usr/bin/env python3
"""I2: basis proposals harvested from an open-search engine (PySR) into
HG-SR v2's term library, so open search proposes and HG-SR decides
(2026-10-02).

This makes systematic an earlier hand step (COMPARISON.md: adopting the
PySR-found x(1-x)/L into HG-SR's basis). Procedure, deterministic given
the training rows:

1. Run PySR (fixed seed, deterministic mode, `sqrt` allowed so slow
   size convergence is expressible) on the TRAINING rows only, in the
   same scaled variables HG-SR uses (x, t = T/1000, il = 100 A / L and
   any task flags), target r = 1/kappa.
2. Expand every hall-of-fame equation into additive terms and strip
   each term's leading numeric coefficient; the remaining expression
   is a candidate basis function (constants inside a term, as in
   sqrt(x*(3.3*il - 0.16)), stay fixed).
3. Keep a candidate only if it is finite on the whole declared domain
   (L -> infinity included) and not collinear (|corr| > 0.9999 on the
   training rows) with an existing basis term or an earlier candidate.
4. Price it like a product of primitives (cost 2, or 3 when it carries
   a fitted constant inside), so it competes on equal terms.

Leakage rule: in cross-validation the proposals are harvested inside
each training fold; a hold-out is never shown to PySR.
"""
from __future__ import annotations

import numpy as np
import sympy

import hgsr2

MAX_PROPOSALS = 12


def scaled(rows, var_map):
    """var_map: name -> fn(row) giving the scaled variable."""
    return np.column_stack([[float(f(r)) for r in rows] for f in var_map.values()])


def run_pysr(rows, var_map, seed=0, niterations=60):
    from pysr import PySRRegressor
    X = scaled(rows, var_map)
    y = 1.0 / np.array([r["kappa"] for r in rows])
    m = PySRRegressor(
        niterations=niterations, populations=15, maxsize=22,
        binary_operators=["+", "-", "*", "/"], unary_operators=["sqrt"],
        model_selection="best", random_state=seed, deterministic=True,
        procs=0, multithreading=False, progress=False, verbosity=0,
        temp_equation_file=True, delete_tempfiles=True)
    m.fit(X, y, variable_names=list(var_map))
    return [str(e) for e in m.equations_["sympy_format"]]


def harvest(equations, var_map, basis, rows, domain_rows):
    syms = {v: sympy.Symbol(v) for v in var_map}
    cand, seen = [], set()
    for eq in equations:
        expr = sympy.expand(sympy.sympify(eq, locals=syms))
        for term in sympy.Add.make_args(expr):
            coef, rest = term.as_coeff_Mul()
            if rest.is_number or rest == 1:
                continue
            key = sympy.srepr(rest)
            if key in seen:
                continue
            seen.add(key)
            cand.append(rest)

    def colfn(rest):
        f = sympy.lambdify([syms[v] for v in var_map], rest, "numpy")
        fns = list(var_map.values())

        def fn(d):
            n = len(d["x"])
            rs = [{k: d[k][i] for k in d} for i in range(n)]
            args = [np.array([float(g(r)) for r in rs]) for g in fns]
            with np.errstate(all="ignore"):
                return np.broadcast_to(np.asarray(f(*args), float), (n,)).astype(float)
        return fn

    d_tr = hgsr2._cols(rows)
    d_dom = hgsr2._cols([dict(r, kappa=1.0) for r in domain_rows])
    cols = [b[2](d_tr) for b in basis]
    out = []
    for rest in cand:
        fn = colfn(rest)
        dom = fn(d_dom)
        tr = fn(d_tr)
        if not (np.all(np.isfinite(dom)) and np.all(np.isfinite(tr))):
            continue
        if np.std(tr) < 1e-12:
            continue
        dup = False
        for c in cols:
            if np.std(c) > 1e-12 and abs(np.corrcoef(c, tr)[0, 1]) > 0.9999:
                dup = True
                break
        if dup:
            continue
        has_float = any(a.is_Float for a in sympy.preorder_traversal(rest))
        name = "pysr:" + str(rest)
        out.append((name, 3 if has_float else 2, fn))
        cols.append(tr)
        if len(out) >= MAX_PROPOSALS:
            break
    return out


def cell_means(rows):
    """One row per design-cell hyperedge (geometric mean of kappa over the
    cell's DISTINCT values, in a fixed cell order). PySR then sees each
    cell once, so replicating a batch cannot change what it proposes:
    the replication invariance of HG-SR extends to the open search."""
    g = {}
    for r in rows:
        g.setdefault(hgsr2.cell(r), []).append(r)
    out = []
    for k in sorted(g):
        rs = g[k]
        lk = np.unique(np.round(np.log([r["kappa"] for r in rs]), 12))
        rep = dict(rs[0])
        for f in ("x", "T", "L"):
            if f in rep:
                rep[f] = float(np.mean([r[f] for r in rs]))
        rep["kappa"] = float(np.exp(lk.mean()))
        out.append(rep)
    return out


def proposals_for(rows, features, var_map, domain_rows, basis_kw=None,
                  cache=None, balanced=True):
    """Harvested proposals for one training set (memoised on the set of
    training sources so repeated fits of the same fold reuse them).
    balanced=True (default since 2026-10-02) shows PySR cell means only;
    balanced=False is the first version, which broke replication
    invariance (results/demo_full_props_unbalanced.json)."""
    key = (balanced,) + tuple(sorted(r["source"] for r in rows))
    if cache is not None and key in cache:
        return cache[key]
    eqs = run_pysr(cell_means(rows) if balanced else rows, var_map)
    basis = hgsr2.make_basis(features, **(basis_kw or {}))
    out = harvest(eqs, var_map, basis, rows, domain_rows)
    if cache is not None:
        cache[key] = out
    return out
