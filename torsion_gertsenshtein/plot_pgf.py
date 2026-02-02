from __future__ import annotations

import shutil
import warnings
from typing import Literal

import matplotlib as mpl

# Try to use PGF backend, but fall back to default if unavailable
try:
    mpl.use("pgf")
    _pgf_available = True
except ImportError:
    _pgf_available = False
    warnings.warn(
        "PGF backend not available. Using default backend instead.",
        UserWarning,
        stacklevel=2,
    )

from matplotlib import rcParams


def _check_tex_available() -> bool:
    """Check if LaTeX executables are available on the system."""
    return any(
        shutil.which(cmd) is not None
        for cmd in ["pdflatex", "xelatex", "lualatex"]
    )


def enable_pgf(
    texsystem: Literal["pdflatex", "xelatex", "lualatex"] = "pdflatex",
    *,
    serif: bool = True,
    base_size: int = 10,
    fallback_on_error: bool = True,
) -> None:
    """Enable PGF backend for matplotlib with LaTeX rendering.
    
    If LaTeX or PGF backend is not available, falls back to standard
    matplotlib settings with similar font configuration.
    
    Args:
        texsystem: LaTeX system to use (pdflatex, xelatex, or lualatex)
        serif: Whether to use serif fonts
        base_size: Base font size
        fallback_on_error: Whether to fall back to standard matplotlib on error
    """
    if not _pgf_available:
        if fallback_on_error:
            _setup_fallback_fonts(serif=serif, base_size=base_size)
            return
        else:
            raise RuntimeError("PGF backend not available and fallback disabled")
    
    if not _check_tex_available():
        if fallback_on_error:
            warnings.warn(
                f"LaTeX system '{texsystem}' not found. Falling back to standard matplotlib.",
                UserWarning,
                stacklevel=2,
            )
            _setup_fallback_fonts(serif=serif, base_size=base_size)
            return
        else:
            raise RuntimeError(f"LaTeX system '{texsystem}' not available and fallback disabled")
    
    try:
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
    except Exception as e:
        if fallback_on_error:
            warnings.warn(
                f"Failed to configure PGF backend: {e}. Falling back to standard matplotlib.",
                UserWarning,
                stacklevel=2,
            )
            _setup_fallback_fonts(serif=serif, base_size=base_size)
        else:
            raise


def _setup_fallback_fonts(*, serif: bool = True, base_size: int = 10) -> None:
    """Set up fallback font configuration when PGF/LaTeX is not available."""
    rcParams.update(
        {
            "text.usetex": False,
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
