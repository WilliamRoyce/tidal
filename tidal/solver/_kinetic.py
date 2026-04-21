"""Diagonal inverse-kinetic helper for time-domain solvers.

Canonical form for equations with non-trivial mass matrix
``M · d²ₜ q = K(q)``: the JSON spec carries ``kinetic_coefficient_symbolic``
on the LHS and the RHS is left un-normalized. Modal solves the generalized
eigenvalue problem ``(K - λM) v = 0`` directly (`modal.py:641-669`, `1797-1820`).

Time-domain solvers (CVODE, IDA, leapfrog, scipy) need to apply ``M⁻¹`` once
at setup so the RHS evaluator can produce ``d²ₜ q = M⁻¹ · K(q)`` cleanly.
For the diagonal case (current scope of #301) this is a scalar per field.

See GitHub #302 (Bug B of #301). Extends to off-diagonal kinetic mixing
(kinetic_matrix_symbolic) under #305.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from tidal.symbolic._eval_utils import evaluate_coefficient

if TYPE_CHECKING:
    from tidal.symbolic.json_loader import EquationSystem


_UNIT_TOLERANCE = 1e-12


def build_inverse_kinetic_diag(
    spec: EquationSystem,
    params: dict[str, float],
) -> dict[str, float] | None:
    """Return per-field inverse kinetic coefficients, or ``None`` for the fast path.

    Iterates over every dynamical equation (``time_derivative_order > 0``) and
    evaluates its ``kinetic_coefficient_symbolic`` at ``params``. If every
    evaluated coefficient is within ``1e-12`` of ``1`` (including equations
    where the symbolic is ``None``, implicitly ``1``), returns ``None`` so the
    caller can skip the per-step multiply — this preserves the existing
    M = I fast path for theories unaffected by #301.

    Otherwise returns ``{field_name: 1 / M_ii}`` for every dynamical field.
    Fields with trivial M = 1 are omitted (the solver treats missing entries
    as unit). Zero kinetic coefficient raises — those fields should already
    have been demoted to constraints before reaching a time-domain solver,
    which the modal path (`modal.py:804`+) handles via Schur elimination.

    Parameters
    ----------
    spec
        Parsed equation system (post ``base_spec`` in perturbative flows).
    params
        Runtime parameter values, e.g. ``{"B0": 1.0, "rho": 0.01}``.

    Returns
    -------
    ``None`` if every dynamical field has ``M_ii ≈ 1``; else a dict
    ``{field_name: 1 / M_ii}`` with only the non-trivial entries.

    Raises
    ------
    ValueError
        If any kinetic coefficient evaluates to ``0`` (singular M).
        The caller should use modal, which demotes the field via Schur
        elimination rather than producing division-by-zero.
    """
    result: dict[str, float] = {}
    nontrivial = False

    for eq in spec.equations:
        if eq.time_derivative_order <= 0:
            continue

        kin_sym = eq.kinetic_coefficient_symbolic
        if kin_sym is None:
            continue

        value = evaluate_coefficient(kin_sym, params, spec.effective_coordinates)
        if not isinstance(value, float):
            value = float(value)  # type: ignore[arg-type]

        if value == 0.0:
            msg = (
                f"Kinetic coefficient for field '{eq.field_name}' evaluates to "
                f"zero at the given parameters ({kin_sym!r}). A time-domain "
                "solver (cvode/ida/leapfrog/scipy) cannot handle a singular "
                "mass matrix. Use modal, which demotes zero-kinetic fields to "
                "constraints via Schur elimination."
            )
            raise ValueError(msg)

        if abs(value - 1.0) > _UNIT_TOLERANCE:
            result[eq.field_name] = 1.0 / value
            nontrivial = True

    return result if nontrivial else None
