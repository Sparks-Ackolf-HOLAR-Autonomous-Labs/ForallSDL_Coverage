#!/usr/bin/env python3
"""Benchmark datasets for the SR engine comparison.

Two tasks, both regressing thermal resistivity r = 1/kappa (Matthiessen:
scattering mechanisms add in r, so the additive model class every engine
searches is physically right for r, not for kappa):

  Task A  "literature"  kappa(x, T) on the literature party's fit-usable
          rows from ../evidence_table.json (same rows fit_party.py uses).
  Task B  "md"          kappa(x, T, L) on the CHPC-MD production points
          (Muller-Plathe, 4x4 cross-section, swap every 50). L is the
          box length Lz in Angstrom; the known ballistic suppression
          means an engine SHOULD discover a 1/L term. Protocol-variant
          probes (6x6/8x8, swap 100/200) are excluded: they vary
          knobs no engine is given as features.

Every row carries its hyperedge tags (party, regime, source) so the
hypergraph-aware engine can use them; plain engines see only x, T, L, r.
"""
import glob
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
PHASE_C = os.path.dirname(HERE)
CHPC_RESULTS = os.path.join(PHASE_C, "..", "chpc", "results")

X_BANDS = [0.0, 0.05, 0.25, 0.45, 0.65, 0.85, 1.01]   # pure / dilute / mid bands
T_BANDS = [0.0, 450.0, 750.0, 1050.0, 1e9]


def regime(x, T):
    xi = next(i for i in range(len(X_BANDS) - 1) if X_BANDS[i] <= x < X_BANDS[i + 1])
    ti = next(i for i in range(len(T_BANDS) - 1) if T_BANDS[i] <= T < T_BANDS[i + 1])
    return f"x{xi}T{ti}"


def literature_task():
    rows = json.load(open(os.path.join(PHASE_C, "evidence_table.json")))
    out = []
    for r in rows:
        if r["party"] == "literature" and r["fit_usable"] and r["T"] is not None:
            out.append({"x": float(r["x"]), "T": float(r["T"]), "L": 0.0,
                        "kappa": float(r["kappa"]),
                        "party": "literature",
                        "regime": regime(float(r["x"]), float(r["T"])),
                        "source": str(r.get("source", "?"))})
    return {"name": "literature", "features": ["x", "T"], "rows": out}


def md_task():
    out = []
    for f in sorted(glob.glob(os.path.join(CHPC_RESULTS, "x*.json"))):
        base = os.path.basename(f)
        # production protocol only: no variant suffix, or explicit a4_w50
        if ("_a" in base) and ("_a4_w50" not in base):
            continue
        d = json.load(open(f))
        out.append({"x": float(d["x"]), "T": float(d["T"]),
                    "L": float(d["box_A"][2]), "kappa": float(d["kappa_W_mK"]),
                    "party": d["party"],
                    "regime": regime(float(d["x"]), float(d["T"])),
                    "source": f"{base[:-5]}"})
    return {"name": "md", "features": ["x", "T", "L"], "rows": out}


def combined_task():
    """Literature + the 1/L-extrapolated MD points, two parties over
    (x, T). Small, but the only task where party-edge provenance can
    fire (single-party terms get scope labels, not confirmed status)."""
    lit = literature_task()["rows"]
    rows = json.load(open(os.path.join(PHASE_C, "evidence_table.json")))
    for r in rows:
        if r["party"] == "MD-extrapolated" and r["fit_usable"]:
            lit = lit + [{"x": float(r["x"]), "T": float(r["T"]), "L": 0.0,
                          "kappa": float(r["kappa"]),
                          "party": "MD-extrapolated",
                          "regime": regime(float(r["x"]), float(r["T"])),
                          "source": str(r.get("source", "?"))}]
    return {"name": "combined", "features": ["x", "T"], "rows": lit}


TASKS = {"literature": literature_task, "md": md_task,
         "combined": combined_task}

if __name__ == "__main__":
    for name, fn in TASKS.items():
        t = fn()
        print(f"{name}: {len(t['rows'])} rows, features {t['features']}")
