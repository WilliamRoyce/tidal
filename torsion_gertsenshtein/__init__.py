"""Torsion-Gertsenshtein: Klein-Gordon PDE simulations for gravitational wave research.

This package provides:
- `kgsim`: Klein-Gordon equation simulations
- `symbolic`: Lagrangian-to-PDE pipeline (symbolic computation layer)
- `vectorfield`: Multi-component field simulations
"""

__version__ = "0.1.0"


def has_tex_support() -> bool:
    """Check if TeX support is available for high-quality plotting."""
    try:
        from .plot_pgf import _check_tex_available, _pgf_available

        return _pgf_available and _check_tex_available()
    except ImportError:
        return False


__all__ = ["__version__", "has_tex_support"]
