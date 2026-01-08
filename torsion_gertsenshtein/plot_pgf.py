from __future__ import annotations

from typing import Literal

import matplotlib as mpl

mpl.use("pgf")
from matplotlib import rcParams


def enable_pgf(
    texsystem: Literal["pdflatex", "xelatex", "lualatex"] = "pdflatex",
    *,
    serif: bool = True,
    base_size: int = 10,
) -> None:
    """Enable PGF backend for matplotlib with LaTeX rendering."""
    rcParams.update(
        {
            "text.usetex": True,
            "pgf.texsystem": texsystem,
            "pgf.rcfonts": False,
            "font.family": "serif" if serif else "sans-serif",
            "font.size": base_size,
            "axes.labelsize": base_size,
            "axes.titlesize": base_size,
            "legend.fontsize": base_size - 2,
            "xtick.labelsize": base_size - 2,
            "ytick.labelsize": base_size - 2,
            "axes.unicode_minus": False,
        }
    )
