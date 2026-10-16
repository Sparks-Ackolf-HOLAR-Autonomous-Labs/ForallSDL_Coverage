"""Figure style for this paper (Molecular Systems Design & Engineering, RSC two-column layout: column 8.3 cm, full width 17.1 cm).

Every data figure is drawn by one script in this folder from committed result
files. Each figure is created at its true print size (column or full journal
width), so a 7 pt tick label is 7 pt on the page; the draft includes the file
at natural size (no width= scaling).

Fonts: ticks 7 pt, axis labels 8 pt, panel titles 9 pt, legends 7 pt.
Colours: Okabe-Ito (colour-vision-deficiency safe); one colour per evidence
source, the same in every figure.
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

COL = 3.27      # one-column width, inches
FULL = 6.73    # two-column (full page) width, inches
TOC_W, TOC_H = 8 / 2.54, 4 / 2.54  # table-of-contents graphic, 8 x 4 cm

HERE = Path(__file__).resolve().parent
PAPER = HERE.parent                 # the paper folder (draft.tex lives here)
FIGDIR = PAPER / "figures"
PROJECT = PAPER.parents[1]          # project root

INK = "#333333"
GREY = "#8C8C8C"
LIGHT = "#EBEBEB"

# One colour per evidence source (Okabe-Ito), the same in every figure.
SOURCE_COLOURS = {
    "literature": "#0072B2",        # machine-read published values
    "tersoff_md": "#E69F00",        # Tersoff molecular dynamics
    "mace_md": "#CC79A7",           # MACE-MP molecular dynamics
    "bte_dft": "#009E73",           # first-principles phonon BTE
    "database": "#56B4E9",          # curated experimental databases
    "starrydata": "#D55E00",        # digitized curves (Starrydata2)
    "law": "#000000",               # fitted law
}
SOURCE_NAMES = {
    "literature": "published values",
    "tersoff_md": "Tersoff MD",
    "mace_md": "MACE-MP MD",
    "bte_dft": "DFT phonon BTE",
    "database": "experimental database",
    "starrydata": "Starrydata2 curves",
    "law": "fitted law",
}


def apply_style():
    plt.rcParams.update({
        "figure.dpi": 100,
        "savefig.dpi": 600,
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans"],
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": INK,
        "axes.labelcolor": INK,
        "text.color": INK,
        "xtick.color": INK,
        "ytick.color": INK,
        "axes.grid": True,
        "grid.color": "#DDDDDD",
        "grid.linewidth": 0.5,
        "axes.axisbelow": True,
        "lines.linewidth": 1.4,
        "legend.frameon": False,
        "figure.constrained_layout.use": True,
        "pdf.fonttype": 42,
    })


def save(fig, name):
    """Write figures/<name>.pdf (vector, used by the draft) and a 600 dpi PNG
    (for review). Never use bbox_inches='tight': it changes the print size."""
    FIGDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGDIR / f"{name}.pdf", facecolor="white", metadata={"CreationDate": None})
    fig.savefig(FIGDIR / f"{name}.png", dpi=600, facecolor="white")
    w, h = fig.get_size_inches()
    plt.close(fig)
    print(f"wrote figures/{name}.pdf ({w:.2f} x {h:.2f} in)")
