"""TIDAL: Tensor Integration and Derivation for Any Lagrangian.

A symbolic physics pipeline for deriving PDEs from Lagrangians using
xAct/Mathematica and simulating them with py-pde.

This package provides:
- `kgsim`: Klein-Gordon equation simulations
- `symbolic`: Lagrangian-to-PDE pipeline (symbolic computation layer)
- `vectorfield`: Multi-component field simulations
"""

from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as _pkg_version

try:
    __version__ = _pkg_version("tidal")
except PackageNotFoundError:
    __version__ = "0.0.0+unknown"

try:
    from .plot_pgf import check_tex_available, pgf_available

    _tex_support_available = pgf_available and check_tex_available()
except ImportError:
    _tex_support_available = False


def has_tex_support() -> bool:
    """Check if TeX support is available for high-quality plotting."""
    return _tex_support_available


__all__ = ["__version__", "has_tex_support"]
