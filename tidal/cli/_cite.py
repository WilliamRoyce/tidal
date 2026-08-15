"""``tidal --cite`` — Print citation information for TIDAL and dependencies."""

from __future__ import annotations

_CITATION = """\
If you use TIDAL in your research, please cite:

  Royce, W. "TIDAL: Tensor Integration and Derivation for Any Lagrangian"
  https://github.com/WilliamRoyce/torsion-gertsenshtein

Dependencies you may also wish to cite:

  SUNDIALS (IDA/CVODE solvers):
    Hindmarsh et al., "SUNDIALS: Suite of Nonlinear and Differential/
    Algebraic Equation Solvers", ACM Trans. Math. Softw. 31(3), 2005.
    DOI: 10.1145/1089014.1089020

  SciPy:
    Virtanen et al., "SciPy 1.0: Fundamental Algorithms for Scientific
    Computing in Python", Nature Methods 17, 2020.
    DOI: 10.1038/s41592-019-0686-2

  NumPy:
    Harris et al., "Array programming with NumPy", Nature 585, 2020.
    DOI: 10.1038/s41586-020-2649-2

  xAct (Wolfram symbolic tensor algebra):
    Martín-García, J.M., "xAct: Efficient tensor computer algebra for
    the Wolfram Language", http://www.xact.es/
"""


def print_citation() -> None:
    """Print citation information to stdout."""
    print(_CITATION)
