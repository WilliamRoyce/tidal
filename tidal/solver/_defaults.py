"""Global solver defaults — single source of truth for tolerance and other parameters.

Leapfrog is fixed-step symplectic and does not use tolerances.
"""

from __future__ import annotations

# Relative/absolute tolerances for adaptive time integrators (CVODE, IDA, scipy).
DEFAULT_RTOL: float = 1e-8
DEFAULT_ATOL: float = 1e-10

# Time-derivative order threshold: equations with time_order >= SECOND_ORDER
# are wave-type (need velocity slots); lower orders are diffusion or constraint.
SECOND_ORDER: int = 2
