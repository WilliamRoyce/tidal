"""Unified RHS evaluation for TIDAL solvers.

``RHSEvaluator`` applies spatial operators with resolved coefficients,
replacing the duplicated inner loops in ``ida.py`` and ``leapfrog.py``.

Depends on ``FieldSet`` (Phase 1) and ``CoefficientEvaluator`` (Phase 2).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from tidal.solver.operators import BCSpec, apply_operator

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
        self._eq_map: dict[str, int] = {
            eq.field_name: i for i, eq in enumerate(spec.equations)
        }

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
            Grid-shaped result array.
        """
        eq = self._spec.equations[eq_idx]
        result = np.zeros(self._grid.shape)
        for term_idx, term in enumerate(eq.rhs_terms):
            result += self._evaluate_term(
                term, fields, t, eq_idx=eq_idx, term_idx=term_idx
            )
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

    def _evaluate_term(
        self,
        term: OperatorTerm,
        fields: FieldSet,
        t: float,
        *,
        eq_idx: int,
        term_idx: int,
    ) -> np.ndarray:
        """Evaluate a single operator term.

        For ``first_derivative_t`` operators (E-L velocity form), resolves
        the time derivative by reading the velocity slot ``v_{field}``
        directly from the state vector.

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
            operated = fields[vel_name]
            coeff = self._coeff_eval.resolve(term, t, eq_idx=eq_idx, term_idx=term_idx)
            return coeff * operated
        target = self._get_field_data(term.field, fields)
        operated = apply_operator(term.operator, target, self._grid, self._bc)
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
