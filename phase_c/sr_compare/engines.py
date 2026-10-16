#!/usr/bin/env python3
"""Uniform adapters for the SR engine comparison.

Adapter contract (all engines fit resistivity r = 1/kappa and predict
kappa back):
    fit(rows, features) -> model
    predict(model, rows) -> np.array of kappa predictions
    describe(model) -> str        (the selected expression)
    complexity(model) -> int      (engine's own complexity measure)
Availability is probed at import; unavailable engines carry the reason
so the comparison table reports why an engine has no numbers
on this machine (aarch64 DGX).
"""
from __future__ import annotations

import numpy as np

import hgsr

REGISTRY = {}


def register(name, available=True, reason=""):
    def deco(cls):
        REGISTRY[name] = {"cls": cls, "available": available, "reason": reason}
        return cls
    return deco


def _xmat(rows, features):
    return np.column_stack([[float(r[f]) for r in rows] for f in features])


# --- our engines -----------------------------------------------------------

@register("minisr-pooled")
class MiniSR:
    """The pilot's exhaustive-basis SR (pooled row MSE). This is also the
    ablation of HG-SR's hyperedge-balanced loss."""
    pooled = True

    def fit(self, rows, features):
        return hgsr.fit(rows, features, pooled=self.pooled)

    def predict(self, model, rows):
        return hgsr.predict(model, rows)

    def describe(self, model):
        s = model["selected"]
        return "1/r, r = " + " + ".join(
            f"{c:.3g}*{n}" for n, c in zip(s["names"], s["coefs"]))

    def complexity(self, model):
        return model["selected"]["complexity"]


@register("hgsr")
class HGSR(MiniSR):
    """Hypergraph-provenance SR: hyperedge-balanced loss + per-term
    provenance/blast radius (reported separately)."""
    pooled = False


# --- PySR -------------------------------------------------------------------

try:
    from pysr import PySRRegressor
    _PYSR_OK, _PYSR_WHY = True, ""
except Exception as e:                                     # pragma: no cover
    _PYSR_OK, _PYSR_WHY = False, f"import failed: {e}"


@register("pysr", available=_PYSR_OK, reason=_PYSR_WHY)
class PySR:
    def fit(self, rows, features):
        X = _xmat(rows, features)
        y = 1.0 / np.array([r["kappa"] for r in rows])
        m = PySRRegressor(
            niterations=60, populations=15, maxsize=22,
            binary_operators=["+", "-", "*", "/"], unary_operators=[],
            model_selection="best", random_state=0, deterministic=True,
            procs=0, multithreading=False, progress=False, verbosity=0,
            temp_equation_file=True, delete_tempfiles=True)
        m.fit(X, y, variable_names=list(features))
        return {"m": m, "features": features}

    def predict(self, model, rows):
        r = model["m"].predict(_xmat(rows, model["features"]))
        return 1.0 / r

    def describe(self, model):
        return "1/r, r = " + str(model["m"].sympy())

    def complexity(self, model):
        return int(model["m"].get_best()["complexity"])


# --- SISSO (TorchSisso implementation) ---------------------------------------

try:
    import pandas as pd
    import sympy
    from TorchSisso import SissoModel  # noqa: F401
    _SISSO_OK, _SISSO_WHY = True, ""
except Exception as e:                                     # pragma: no cover
    _SISSO_OK, _SISSO_WHY = False, f"import failed: {e}"


@register("sisso", available=_SISSO_OK, reason=_SISSO_WHY)
class Sisso:
    """SISSO (sure-independence screening + sparsifying operator), via
    the TorchSisso package (Muthyala et al. 2024, Digital Discovery)."""

    def fit(self, rows, features):
        df = pd.DataFrame({"target": [1.0 / r["kappa"] for r in rows]})
        for f in features:
            df[f] = [float(r[f]) for r in rows]
        sm = SissoModel(data=df, operators=["+", "-", "*", "/"],
                        n_expansion=2, n_term=3, k=10)
        out = sm.fit()          # (rmse, equation, r2, ...)
        eq = str(out[1])
        expr = sympy.sympify(eq, locals={f: sympy.Symbol(f) for f in features})
        fn = sympy.lambdify([sympy.Symbol(f) for f in features], expr, "numpy")
        return {"expr": expr, "fn": fn, "features": features}

    def predict(self, model, rows):
        cols = [np.array([float(r[f]) for r in rows]) for f in model["features"]]
        r = np.broadcast_to(np.asarray(model["fn"](*cols), float), (len(rows),))
        return 1.0 / r

    def describe(self, model):
        return "1/r, r = " + str(model["expr"])

    def complexity(self, model):
        return int(sympy.count_ops(model["expr"]))


# --- AI Feynman ---------------------------------------------------------------

@register("aifeynman", available=False,
          reason="needs gfortran, so it cannot build on the DGX (aarch64); "
                 "built and ran on CHPC instead (job 15855552, see "
                 "aif_chpc/): zero-byte solution files on literature/"
                 "combined, LinAlgError on md. Designed for large "
                 "near-noiseless synthetic data; unusable on these tables.")
class AIFeynman:
    pass


# --- QLattice / Feyn ----------------------------------------------------------

try:
    import pandas as _pd
    import sympy as _sympy
    import feyn as _feyn
    _FEYN_OK, _FEYN_WHY = True, ""
except Exception as e:                                     # pragma: no cover
    _FEYN_OK, _FEYN_WHY = (False,
        "closed-source wheels are x86_64/win/mac only (no linux-aarch64), "
        "so it cannot import natively on the DGX; run it inside the "
        "amd64-emulated udocker container instead, see "
        "../../containers/README.md (import failed: %s)" % e)


@register("feyn", available=_FEYN_OK, reason=_FEYN_WHY)
class Feyn:
    """QLattice (Abzu's Feyn 3.x, closed source). Runs on the DGX only
    under qemu-x86_64 emulation inside the udocker container; wall times
    are emulated-CPU times, not native."""

    def _df(self, rows, features):
        df = _pd.DataFrame({f: [float(r[f]) for r in rows] for f in features})
        df["r"] = [1.0 / r["kappa"] for r in rows]
        return df

    def fit(self, rows, features):
        ql = _feyn.QLattice(random_seed=1)
        # n_epochs=5 (not the default 10): under qemu emulation each epoch
        # costs ~1 min and 5 already recovered the planted smoke-test law;
        # threads=1 for the determinism check.
        models = ql.auto_run(self._df(rows, features), output_name="r",
                             n_epochs=5, max_complexity=10, threads=1)
        return {"m": models[0], "features": features}

    def predict(self, model, rows):
        r = model["m"].predict(self._df(rows, model["features"]))
        return 1.0 / np.asarray(r, dtype=float)

    def describe(self, model):
        return "1/r, r = " + str(model["m"].sympify(signif=6))

    def complexity(self, model):
        return int(model["m"].edge_count)
