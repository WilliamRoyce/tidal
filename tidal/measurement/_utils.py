"""Shared private utilities for the measurement package.

Functions here are used by sibling modules (_conversion, _spectral_conversion,
_dispersion) within the measurement package.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Sequence

    from tidal.measurement._io import SimulationData


def normalize_group(
    fields: str | Sequence[str],
) -> tuple[str, ...]:
    """Normalize a string or sequence to a tuple of field names."""
    if isinstance(fields, str):
        return (fields,)
    return tuple(fields)


def check_no_position_dependent_terms(
    data: SimulationData,
    context: str,
) -> None:
    """Raise ValueError if any RHS term is position-dependent.

    FFT-based spectral methods (P(k,t), omega(k)) assume spatial translation
    invariance.  Position-dependent coefficients break this assumption and
    produce physically meaningless results.

    Raises
    ------
    ValueError
        If any equation has a position-dependent RHS term.
    """
    for eq in data.spec.equations:
        for term in eq.rhs_terms:
            if term.position_dependent:
                msg = (
                    f"{context} requires a spatially uniform system, but "
                    f"equation '{eq.field_name}' has position-dependent term: "
                    f"{term.coefficient_symbolic or term.coefficient} * "
                    f"{term.operator}({term.field}). "
                    f"Position-dependent coefficients break translation "
                    f"invariance — Fourier decomposition is not physical."
                )
                raise ValueError(msg)
