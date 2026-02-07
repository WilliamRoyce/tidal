"""Torsion-Gertsenshtein: Klein-Gordon PDE simulations for gravitational wave research.

This package provides:
- `kgsim`: Klein-Gordon equation simulations
- `symbolic`: Lagrangian-to-PDE pipeline (symbolic computation layer)
- `vectorfield`: Multi-component field simulations
"""

__version__ = "0.2.7"

try:
    from .plot_pgf import check_tex_available, pgf_available

    _tex_support_available = pgf_available and check_tex_available()
except ImportError:
    _tex_support_available = False


def has_tex_support() -> bool:
    """Check if TeX support is available for high-quality plotting."""
    return _tex_support_available


__all__ = ["__version__", "has_tex_support"]
