from __future__ import annotations

import shutil
import warnings
from typing import Literal

import matplotlib as mpl
from matplotlib import rcParams

# Try to use PGF backend, but fall back to default if unavailable.
# Check for TeX availability FIRST — mpl.use("pgf") can succeed even when
# no LaTeX executable is installed, causing deferred crashes on render.
pgf_available = False
_pgf_backend_error = None

if any(shutil.which(cmd) is not None for cmd in ["pdflatex", "xelatex", "lualatex"]):
    try:
        mpl.use("pgf")
        pgf_available = True
    except RuntimeError as e:
        _pgf_backend_error = str(e)
        warnings.warn(
            f"PGF backend not available: {e}. Using default backend instead.",
            UserWarning,
            stacklevel=2,
        )
else:
    _pgf_backend_error = "No LaTeX executable found (pdflatex, xelatex, lualatex)"


def check_tex_available() -> bool:
    """Check if LaTeX executables are available on the system."""
    return any(
        shutil.which(cmd) is not None for cmd in ["pdflatex", "xelatex", "lualatex"]
    )


def _check_pdf_converter_available() -> bool:
    """Check if PDF to PNG converter is available (ghostscript or similar)."""
    return shutil.which("gs") is not None or shutil.which("convert") is not None


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

    Raises
    ------
    RuntimeError
        If PGF backend, LaTeX system, or PDF converter is not available and
        fallback_on_error is False.
    """
    if not pgf_available:
        if fallback_on_error:
            _setup_fallback_fonts(serif=serif, base_size=base_size)
            return
        msg = f"PGF backend not available ({_pgf_backend_error}) and fallback disabled"
        raise RuntimeError(msg)

    if not check_tex_available():
        if fallback_on_error:
            warnings.warn(
                f"LaTeX system '{texsystem}' not found. Falling back to standard matplotlib.",
                UserWarning,
                stacklevel=2,
            )
            _setup_fallback_fonts(serif=serif, base_size=base_size)
            return
        msg = f"LaTeX system '{texsystem}' not available and fallback disabled"
        raise RuntimeError(msg)

    if not _check_pdf_converter_available():
        if fallback_on_error:
            warnings.warn(
                "PDF to PNG converter not found (ghostscript). "
                "Falling back to standard matplotlib.",
                UserWarning,
                stacklevel=2,
            )
            _setup_fallback_fonts(serif=serif, base_size=base_size)
            return
        msg = "PDF to PNG converter (ghostscript) not available and fallback disabled"
        raise RuntimeError(msg)

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
    # Ensure we're not using the pgf backend
    if mpl.get_backend() == "pgf":
        mpl.use("Agg")  # Use Agg backend which is reliable and supports all formats

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
