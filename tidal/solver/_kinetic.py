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

import numpy as np

from tidal.symbolic._eval_utils import evaluate_coefficient

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from tidal.solver.grid import GridInfo
    from tidal.symbolic.json_loader import EquationSystem


_UNIT_TOLERANCE = 1e-12


def build_inverse_kinetic_diag(
    spec: EquationSystem,
    params: dict[str, float],
    grid: GridInfo | None = None,
) -> dict[str, float | NDArray[np.float64]] | None:
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
    grid
        Optional grid. When provided, per-coordinate arrays are built from
        ``grid.coord_arrays()`` and passed to ``evaluate_coefficient`` so
        position-dependent kinetic coefficients (post-GH #380:
        kinetic_coefficient_symbolic containing ``x[]``, ``y[]``, …) evaluate
        to per-grid-point arrays. ``None`` (the default) preserves the
        constant-coefficient behavior. See GH #382.

    Returns
    -------
    ``None`` if every dynamical field has ``M_ii ≈ 1``; else a dict
    ``{field_name: 1 / M_ii}`` with only the non-trivial entries. Per-field
    values are scalars for constant kinetics and ``NDArray[float64]`` for
    position-dependent kinetics (when ``coord_arrays`` is supplied).

    Raises
    ------
    ValueError
        If any kinetic coefficient evaluates to ``0`` (singular M).
        The caller should use modal, which demotes the field via Schur
        elimination rather than producing division-by-zero.
    """
    result: dict[str, float | NDArray[np.float64]] = {}
    nontrivial = False

    # Build per-coordinate grid arrays once (only if a grid was provided).
    # The same loop is used by CoefficientEvaluator (coefficients.py:88-93).
    coord_arrays: dict[str, NDArray[np.float64]] | None
    if grid is not None:
        coord_arrays = {}
        grid_coords = grid.coord_arrays()
        spatial = [c for c in spec.effective_coordinates if c != "t"]
        for i, name in enumerate(spatial):
            if i < len(grid_coords):
                coord_arrays[name] = grid_coords[i]
    else:
        coord_arrays = None

    for eq in spec.equations:
        if eq.time_derivative_order <= 0:
            continue

        kin_sym = eq.kinetic_coefficient_symbolic
        if kin_sym is None:
            continue

        value = evaluate_coefficient(
            kin_sym,
            params,
            spec.effective_coordinates,
            coord_arrays,
        )

        if isinstance(value, np.ndarray):
            # Position-dependent kinetic. Zero-anywhere is a singular-M failure
            # at that grid point — same diagnostic as the scalar case.
            if (value == 0.0).any():
                msg = (
                    f"Kinetic coefficient for field '{eq.field_name}' evaluates "
                    f"to zero at at least one grid point ({kin_sym!r}). A "
                    "time-domain solver (cvode/ida/leapfrog/scipy) cannot "
                    "handle a singular mass matrix. Use modal, which demotes "
                    "zero-kinetic fields to constraints via Schur elimination."
                )
                raise ValueError(msg)
            inv = 1.0 / value
            # Only count as non-trivial if it actually departs from M=1
            # anywhere on the grid — otherwise treat as the fast-path identity.
            if not np.allclose(inv, 1.0, atol=_UNIT_TOLERANCE):
                result[eq.field_name] = inv.astype(np.float64)
                nontrivial = True
            continue

        scalar = float(value)
        if scalar == 0.0:
            msg = (
                f"Kinetic coefficient for field '{eq.field_name}' evaluates to "
                f"zero at the given parameters ({kin_sym!r}). A time-domain "
                "solver (cvode/ida/leapfrog/scipy) cannot handle a singular "
                "mass matrix. Use modal, which demotes zero-kinetic fields to "
                "constraints via Schur elimination."
            )
            raise ValueError(msg)

        if abs(scalar - 1.0) > _UNIT_TOLERANCE:
            result[eq.field_name] = 1.0 / scalar
            nontrivial = True

    return result if nontrivial else None


def velocity_row_scale(
    field_name: str,
    m_inv: dict[str, float | NDArray[np.float64]] | None,
) -> float | NDArray[np.float64]:
    """Return the ``M⁻¹`` scale for a 2nd-order field's velocity-row contributions.

    Single source of truth for the ``kinetic_coefficient_symbolic`` consumption
    pattern used by every matrix builder in :mod:`tidal.solver.modal` that
    emits entries into ``A[velocity_slot, target_slot]``. Builders MUST call
    this rather than inline the lookup, so the contract is reviewable from
    the function definition and so the regression guard in
    :class:`tests.test_solver_kinetic_consistency.TestAllModalPathsRespectKinetic`
    can pin all callers to the same behavior.

    ``m_inv`` is the dict returned by :func:`build_inverse_kinetic_diag`.
    ``None`` (every field has ``M ≈ 1``) and missing-from-dict (this field
    has ``M ≈ 1``) both yield ``1.0`` — zero-overhead for theories where
    the kinetic-coefficient class doesn't apply.

    Note: :func:`tidal.solver.modal._build_evolution_matrices` (the Schur /
    generalized-eig path) does NOT call this helper. That builder writes
    ``M`` directly onto the diagonal of ``M_mat`` for the generalized
    eigenvalue problem ``(K − λM) v = 0`` rather than post-multiplying a
    velocity-row contribution; it has a structurally different consumption
    pattern that's documented separately in :mod:`tidal.solver.modal`'s
    module docstring.

    Parameters
    ----------
    field_name : str
        Name of the dynamical (2nd-order) field whose velocity-row
        contribution is being assembled.
    m_inv : dict[str, float | NDArray[np.float64]] | None
        Output of :func:`build_inverse_kinetic_diag`. ``None`` denotes the
        fast path (no theory in the spec has non-trivial ``M``); a missing
        key denotes that this specific field has ``M ≈ 1``. Values can be
        scalars (constant kinetic) or arrays (position-dependent kinetic
        when GH #382's ``coord_arrays`` path is exercised).

    Returns
    -------
    float | NDArray[np.float64]
        Multiplicative factor to apply to every velocity-row contribution
        for ``field_name``. ``1.0`` for the no-op path.

        An array is only ever returned when ``grid`` is supplied, which
        builds the ``coord_arrays`` the expression needs (GH #382).  The
        four time-domain backends do this; the modal builders currently do
        not, so a position-dependent kinetic reaching them raises rather
        than being silently reduced — see GH #421.

    See Also
    --------
    build_inverse_kinetic_diag : Source of the ``m_inv`` dict.
    """
    return 1.0 if m_inv is None else m_inv.get(field_name, 1.0)
