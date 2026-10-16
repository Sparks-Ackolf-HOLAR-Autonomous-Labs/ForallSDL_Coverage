"""Tally of the evidence sources in the method-to-regime assignment, from committed files only.

Counting rule shared by the project's papers: one evidence source is one
method family plus one provenance chain. Tersoff MD on the DGX Spark, on
CHPC and the GPUMD HNEMD runs are three sources of one method family; the
literature extraction and the two curated compilations (Citrine,
Starrydata2) are three sources of the measurement family.

The sources are read from the assignment solver's own output
(p3_prework/covering/results.json: askable sources and every fixed source
that covers a cell), so the counts quoted in the draft cannot drift from
the problem that was solved.

Output: audit/source_tally.json. Run from any directory:
    python3 <path to this folder>/source_tally.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROJ = HERE.parents[2]  # project root
RES = PROJ / "p3_prework" / "covering" / "results.json"

# solver party id -> (reader-facing name, method family under the shared rule)
NAMES = {
    "DGX-MD": ("Tersoff MD (DGX Spark)", "Tersoff MD"),
    "CHPC-MD": ("Tersoff MD (CHPC)", "Tersoff MD"),
    "GPUMD-HNEMD": ("Tersoff HNEMD (GPUMD)", "Tersoff MD"),
    "MACE-MD": ("MACE-MP-0 MD", "machine-learned MD"),
    "BTE-DFT": ("first-principles phonon BTE", "first-principles BTE"),
    "LIT": ("literature extraction", "measurement"),
    "DB-CITRINE": ("Citrine compilation", "measurement"),
    "DB-STARRY": ("Starrydata2 compilation", "measurement"),
}


def main() -> None:
    r = json.loads(RES.read_text())
    askable = sorted(r["families"])
    fixed = sorted({p for cell in r["fixed_coverage"].values() for p in cell})
    srcs = askable + [p for p in fixed if p not in askable]
    unknown = [p for p in srcs if p not in NAMES]
    assert not unknown, unknown
    out = {
        "counting_rule": "one method family plus one provenance chain",
        "n_sources": len(srcs),
        "n_askable": len(askable),
        "n_fixed": len(srcs) - len(askable),
        "n_method_families": len({NAMES[p][1] for p in srcs}),
        "sources": {p: {"name": NAMES[p][0], "method_family": NAMES[p][1],
                        "role": "askable" if p in askable else "fixed"} for p in srcs},
    }
    (HERE / "source_tally.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
