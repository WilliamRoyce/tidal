"""Bayesian inference for TIDAL parameter estimation.

Wraps the existing simulation + measurement pipeline as a likelihood
function for Monte Carlo and nested sampling.  Nested sampling uses
**PolyChord** (Handley et al. 2015) with **anesthetic** (Handley 2019)
for analysis and visualization.

Setup::

    pip install tidal[inference]            # anesthetic only
    bash scripts/install_polychord.sh       # PolyChord (requires gfortran)

References
----------
Skilling, J. (2004) "Nested Sampling", AIP Conference Proceedings 735.
Handley, W. et al. (2015) "PolyChord: next-generation nested sampling", MNRAS 453.
Handley, W. (2019) "anesthetic: nested sampling visualization", JOSS 4.
"""

from __future__ import annotations

from tidal.inference._constraints import ConstraintSet, parse_constraint
from tidal.inference._importance import ParameterImportanceResult
from tidal.inference._prior import Prior
from tidal.inference._results import InferenceResult

__all__ = [
    "ConstraintSet",
    "InferenceResult",
    "ParameterImportanceResult",
    "Prior",
    "parse_constraint",
]
