"""Bayesian inference for TIDAL parameter estimation.

Wraps the existing simulation + measurement pipeline as a likelihood
function for Monte Carlo and nested sampling.

Optional dependencies::

    pip install tidal[inference]   # dynesty + anesthetic

References
----------
Skilling, J. (2004) "Nested Sampling", AIP Conference Proceedings 735.
Speagle, J.S. (2020) "dynesty: a dynamic nested sampling package", MNRAS 493.
Handley, W. et al. (2015) "PolyChord: next-generation nested sampling", MNRAS 453.
Handley, W. (2019) "anesthetic: nested sampling visualization", JOSS 4.
"""

from __future__ import annotations

from tidal.inference._constraints import ConstraintSet, parse_constraint
from tidal.inference._prior import Prior
from tidal.inference._results import InferenceResult

__all__ = [
    "ConstraintSet",
    "InferenceResult",
    "Prior",
    "parse_constraint",
]
