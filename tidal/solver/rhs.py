"""Unified RHS evaluation for TIDAL solvers.

``RHSEvaluator`` applies spatial operators with resolved coefficients,
replacing the duplicated inner loops in ``ida.py`` and ``leapfrog.py``.

Depends on ``FieldSet`` (Phase 1) and ``CoefficientEvaluator`` (Phase 2).
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np

from tidal.solver.operators import (
    OPERATOR_REGISTRY,
    AxisBCSpec,
    BCSpec,
    _bc_from_grid,
    _normalize_bc,
    _resolve_axis_bc,
    apply_operator,
)

if TYPE_CHECKING:
    from tidal.solver.coefficients import CoefficientEvaluator
    from tidal.solver.fields import FieldSet
    from tidal.solver.grid import GridInfo
    from tidal.symbolic.json_loader import EquationSystem, OperatorTerm


class RHSEvaluator:
    """Evaluate RHS of field equations with resolved coefficients.

    Applies spatial operators to field data and multiplies by resolved
    coefficients (constant, parameter-overridden, position-dependent,
    or time-dependent).

    Parameters
    ----------
    spec : EquationSystem
        Parsed JSON equation specification.
    grid : GridInfo
        Spatial grid.
    coeff_eval : CoefficientEvaluator
        Coefficient resolver with caching.
    bc : str or tuple of str, optional
        Boundary conditions for spatial operators.
    """

    def __init__(
        self,
        spec: EquationSystem,
        grid: GridInfo,
        coeff_eval: CoefficientEvaluator,
        bc: BCSpec | None = None,
    ) -> None:
        self._spec = spec
        self._grid = grid
        self._coeff_eval = coeff_eval
        self._bc = bc
        self._eq_map: dict[str, int] = spec.equation_map
        self._result_buffer = np.zeros(grid.shape)
        self._term_buffer = np.zeros(grid.shape)

        # Pre-normalize BCs once (avoids per-call validation in operators)
        if bc is not None:
            self._normalized_bc: BCSpec = _normalize_bc(bc, grid)
        else:
            self._normalized_bc = _bc_from_grid(grid)

        # Pre-resolve per-axis BCs to AxisBCSpec (avoids isinstance +
        # _str_to_axis_bc object creation per operator call in hot path)
        self._resolved_bcs: tuple[AxisBCSpec, ...] = tuple(
            _resolve_axis_bc(bc_entry) for bc_entry in self._normalized_bc
        )

        # Pre-resolve operator functions (avoids dict lookup per call)
        self._resolved_ops: list[list[tuple[Any, str]]] = []
        for eq in spec.equations:
            ops: list[tuple[Any, str]] = []
            for term in eq.rhs_terms:
                if term.operator == "first_derivative_t":
                    ops.append((None, term.field))
                else:
                    fn = OPERATOR_REGISTRY.get(term.operator)
                    if fn is None:
                        msg = (
                            f"Unknown operator {term.operator!r}; "
                            f"known: {sorted(OPERATOR_REGISTRY)}"
                        )
                        raise ValueError(msg)
                    ops.append((fn, term.field))
            self._resolved_ops.append(ops)

    def begin_timestep(self, t: float) -> None:
        """Notify the coefficient evaluator of a new timestep."""
        self._coeff_eval.begin_timestep(t)

    def evaluate(
        self,
        eq_idx: int,
        fields: FieldSet,
        t: float = 0.0,
    ) -> np.ndarray:
        """Compute RHS for a single equation.

        Parameters
        ----------
        eq_idx : int
            Index of the equation in ``spec.equations``.
        fields : FieldSet
            Current field state.
        t : float
            Current simulation time.

        Returns
        -------
        np.ndarray
            Grid-shaped result array.  **Warning:** the returned array is
            an internal buffer and may be overwritten by the next call to
            ``evaluate()``.  Callers must copy if they need to persist it.
        """
        eq = self._spec.equations[eq_idx]
        result = self._result_buffer
        temp = self._term_buffer
        terms = eq.rhs_terms
        if not terms:
            result.fill(0.0)
            return result

        resolved = self._resolved_ops[eq_idx]

        # First term: write directly to result (eliminates fill(0))
        operated = self._apply_resolved(resolved[0], fields)
        coeff = self._coeff_eval.resolve(terms[0], t, eq_idx=eq_idx, term_idx=0)
        np.multiply(coeff, operated, out=result)

        # Remaining terms: accumulate
        for term_idx in range(1, len(terms)):
            operated = self._apply_resolved(resolved[term_idx], fields)
            coeff = self._coeff_eval.resolve(
                terms[term_idx], t, eq_idx=eq_idx, term_idx=term_idx
            )
            np.multiply(coeff, operated, out=temp)
            result += temp

        return result

    def evaluate_by_field(
        self,
        field_name: str,
        fields: FieldSet,
        t: float = 0.0,
    ) -> np.ndarray:
        """Compute RHS for the equation governing *field_name*.

        Raises
        ------
        KeyError
            If *field_name* has no associated equation.
        """
        eq_idx = self._eq_map.get(field_name)
        if eq_idx is None:
            msg = f"No equation for field '{field_name}'"
            raise KeyError(msg)
        return self.evaluate(eq_idx, fields, t)

    # ---- Internal ----

    def _apply_resolved(
        self,
        resolved: tuple[Any, str],
        fields: FieldSet,
    ) -> np.ndarray:
        """Apply a pre-resolved operator, returning the operated data.

        Uses pre-resolved function pointer and pre-normalized BCs to
        avoid per-call dict lookup and BC validation overhead.

        Raises
        ------
        ValueError
            If a ``first_derivative_t`` term references a field whose
            velocity slot is not present in the state.
        """
        op_fn, field_name = resolved
        if op_fn is None:  # first_derivative_t
            vel_name = f"v_{field_name}"
            if vel_name not in fields:
                msg = (
                    f"Cannot resolve first_derivative_t({field_name}): "
                    f"velocity slot '{vel_name}' not found. "
                    f"Available: {sorted(fields.slot_names)}"
                )
                raise ValueError(msg)
            return fields[vel_name]
        target = self._get_field_data(field_name, fields)
        return op_fn(target, self._grid, self._resolved_bcs)

    def _apply_operator(
        self,
        term: OperatorTerm,
        fields: FieldSet,
    ) -> np.ndarray:
        """Apply the spatial operator for a term, returning the operated data.

        For ``first_derivative_t``, resolves to the velocity slot.

        Raises
        ------
        ValueError
            If a ``first_derivative_t`` term references a field whose
            velocity slot is not present in the state.
        """
        if term.operator == "first_derivative_t":
            vel_name = f"v_{term.field}"
            if vel_name not in fields:
                msg = (
                    f"Cannot resolve first_derivative_t({term.field}): "
                    f"velocity slot '{vel_name}' not found. "
                    f"Available: {sorted(fields.slot_names)}"
                )
                raise ValueError(msg)
            return fields[vel_name]
        target = self._get_field_data(term.field, fields)
        return apply_operator(term.operator, target, self._grid, self._bc)

    def _evaluate_term(
        self,
        term: OperatorTerm,
        fields: FieldSet,
        t: float,
        *,
        eq_idx: int,
        term_idx: int,
    ) -> np.ndarray:
        """Evaluate a single operator term (allocating variant).

        Delegates to ``_apply_operator`` for the spatial operator, then
        multiplies by the resolved coefficient.  Returns a new array.
        """
        operated = self._apply_operator(term, fields)
        coeff = self._coeff_eval.resolve(term, t, eq_idx=eq_idx, term_idx=term_idx)
        return coeff * operated

    @staticmethod
    def _get_field_data(field_name: str, fields: FieldSet) -> np.ndarray:
        """Get field data.  Raises on unknown field references.

        Raises
        ------
        ValueError
            If *field_name* cannot be resolved to any known field.
        """
        if field_name in fields:
            return fields[field_name]
        msg = (
            f"Unknown field reference '{field_name}'. "
            f"Available: {sorted(fields.slot_names)}"
        )
        raise ValueError(msg)
