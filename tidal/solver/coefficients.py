"""Coefficient evaluation with multi-level caching for TIDAL solvers.

``CoefficientEvaluator`` resolves operator term coefficients, handling:
- Constant numeric coefficients (no symbolic expression)
- Parameter-overridable symbolic coefficients (e.g. ``"m2"``, ``"-kappa"``)
- Position-dependent coefficients (evaluated on the spatial grid)
- Time-dependent coefficients (re-evaluated each timestep)

Caching levels:
    L0: Constants pre-resolved at init (no coordinate/time dependence)
    L1: Mathematica → Python string cache (computed once)
    L2: Spatial-only arrays (position-dependent, NOT time-dependent)
    L3: Per-timestep deduplication (time+position dependent)

Reuses the evaluation pipeline from ``tidal.symbolic._eval_utils`` — no
reimplementation of Mathematica→Python conversion or eval() logic.
"""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, cast

import numpy as np

from tidal.symbolic._eval_utils import (
    build_eval_namespace,
    evaluate_coefficient,
    mathematica_to_python,
)

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from tidal.solver.grid import GridInfo
    from tidal.symbolic.json_loader import EquationSystem, OperatorTerm


_MISSING = object()  # sentinel for dict.get() fast path


class CoefficientEvaluator:
    """Resolve operator term coefficients with multi-level caching.

    Parameters
    ----------
    spec : EquationSystem
        Parsed JSON equation specification.
    grid : GridInfo
        Spatial grid (for position-dependent coefficient evaluation).
    parameters : dict[str, float]
        Runtime parameter overrides (e.g. ``{"m2": 1.0}``).

    Raises
    ------
    ValueError
        If a symbolic coefficient references an unknown parameter and
        no default numeric value is available.
    """

    def __init__(
        self,
        spec: EquationSystem,
        grid: GridInfo,
        parameters: dict[str, float] | None = None,
    ) -> None:
        self._spec = spec
        self._grid = grid
        self._parameters = parameters or {}
        self._coordinates = spec.effective_coordinates
        self._spatial_coords = spec.spatial_coordinates

        # L1: Mathematica → Python expression string cache
        self._expr_cache: dict[str, str] = {}

        # Static eval namespace (math funcs + parameters)
        self._namespace = build_eval_namespace(self._parameters)

        # Build coordinate arrays for position-dependent evaluation
        self._coord_arrays: dict[str, NDArray[np.float64]] = {}
        grid_coords = grid.coord_arrays()
        for i, name in enumerate(self._spatial_coords):
            if i < len(grid_coords):
                self._coord_arrays[name] = grid_coords[i]

        # L0: Pre-resolve constant coefficients
        self._constants: dict[tuple[int, int], float] = {}
        for eq_idx, eq in enumerate(spec.equations):
            for term_idx, term in enumerate(eq.rhs_terms):
                if self._is_constant(term):
                    self._constants[eq_idx, term_idx] = self._resolve_constant(term)

        # L2: Spatial-only cache (position-dependent, NOT time-dependent)
        self._spatial_cache: dict[tuple[int, int], NDArray[np.float64]] = {}
        self._precompute_spatial()

        # L3: Per-timestep cache
        self._timestep_cache: dict[
            tuple[int, int, float], float | NDArray[np.float64]
        ] = {}

        # Fail-fast: validate all unresolved symbolic terms at init
        self._validate_unresolved()

        # Diagnostic: warn about mass sign changes
        self._check_mass_sign()

    def resolve(
        self,
        term: OperatorTerm,
        t: float = 0.0,
        *,
        eq_idx: int = -1,
        term_idx: int = -1,
    ) -> float | NDArray[np.float64]:
        """Resolve effective coefficient for an operator term.

        Parameters
        ----------
        term : OperatorTerm
            The operator term whose coefficient to resolve.
        t : float
            Current simulation time.
        eq_idx : int
            Equation index (for cache lookup). -1 disables caching.
        term_idx : int
            Term index within equation (for cache lookup). -1 disables caching.

        Returns
        -------
        float | NDArray[np.float64]
            Scalar for constant coefficients, grid-shaped array for
            position-dependent ones.
        """
        key = (eq_idx, term_idx)

        # L0: Pre-resolved constant (single dict.get avoids two-op in+[])
        c = self._constants.get(key, _MISSING)
        if c is not _MISSING:
            return cast("float | NDArray[np.float64]", c)

        # No symbolic → return numeric coefficient directly
        if term.coefficient_symbolic is None:
            return term.coefficient

        # Spatial-only cache (position-dependent, not time-dependent)
        c = self._spatial_cache.get(key, _MISSING)
        if c is not _MISSING:
            return cast("NDArray[np.float64]", c)

        # L3: Per-timestep cache
        ts_key = (eq_idx, term_idx, t)
        c = self._timestep_cache.get(ts_key, _MISSING)
        if c is not _MISSING:
            return cast("float | NDArray[np.float64]", c)

        # Full evaluation
        result = self._evaluate_symbolic(term, t)

        # Cache if time-dependent (L3)
        if term.time_dependent and eq_idx >= 0:
            self._timestep_cache[ts_key] = result

        return result

    def begin_timestep(self, t: float) -> None:  # noqa: ARG002
        """Clear per-timestep cache (L3).

        Call at the start of each timestep to ensure time-dependent
        coefficients are re-evaluated.
        """
        self._timestep_cache.clear()

    def all_constant(self) -> bool:
        """Check if every RHS term has a constant (scalar) coefficient.

        Returns True when all coefficients are pre-resolved at L0 (no
        position or time dependence).  This enables the analytical
        Jacobian optimization for constant-coefficient linear systems.
        """
        for eq_idx, eq in enumerate(self._spec.equations):
            for term_idx, _term in enumerate(eq.rhs_terms):
                if (eq_idx, term_idx) not in self._constants:
                    return False
        return True

    def all_time_independent(self) -> bool:
        """Check if every RHS term has a time-independent coefficient.

        Returns True when no term has ``time_dependent=True``.  Position-
        dependent (but time-independent) coefficients still produce a
        constant Jacobian because the spatial grid is fixed, so the
        analytical Jacobian optimization applies.
        """
        for eq in self._spec.equations:
            for term in eq.rhs_terms:
                if term.time_dependent:
                    return False
        return True

    # ---- Internal helpers ----

    def _is_constant(self, term: OperatorTerm) -> bool:
        """Check if a term has a constant (non-varying) coefficient."""
        sym = term.coefficient_symbolic
        if sym is None:
            return True
        if term.time_dependent or term.position_dependent:
            return False
        # Simple parameter name or negated parameter
        if sym in self._parameters:
            return True
        if sym.startswith("-") and sym[1:] in self._parameters:
            return True
        # Try to evaluate as compound expression (no coords needed)
        try:
            self._evaluate_symbolic(term, 0.0)
        except (ValueError, TypeError, NameError):
            return False
        return True

    def _resolve_constant(self, term: OperatorTerm) -> float:
        """Resolve a constant coefficient to a float.

        Raises
        ------
        TypeError
            If a compound expression produces an array instead of scalar.
        """
        sym = term.coefficient_symbolic
        if sym is None:
            return term.coefficient

        # Fast path: direct parameter lookup
        if sym in self._parameters:
            return self._parameters[sym]
        if sym.startswith("-") and sym[1:] in self._parameters:
            return -self._parameters[sym[1:]]

        # Compound expression (e.g. "-2*m2")
        result = self._evaluate_symbolic(term, 0.0)
        if isinstance(result, np.ndarray):
            msg = (
                f"Expected scalar for constant coefficient '{sym}', "
                f"got array. Check coordinate_dependent flags."
            )
            raise TypeError(msg)
        return float(result)

    def _precompute_spatial(self) -> None:
        """Pre-compute L2 spatial-only coefficients."""
        for eq_idx, eq in enumerate(self._spec.equations):
            for term_idx, term in enumerate(eq.rhs_terms):
                key = (eq_idx, term_idx)
                if key in self._constants:
                    continue
                if term.coefficient_symbolic is None:
                    continue
                if term.position_dependent and not term.time_dependent:
                    result = self._evaluate_symbolic(term, 0.0)
                    if isinstance(result, np.ndarray):
                        self._spatial_cache[key] = result
                    else:
                        # Scalar result despite position_dependent flag —
                        # treat as constant
                        self._constants[key] = float(result)

    def _validate_unresolved(self) -> None:
        """Fail fast on symbolic terms that can't be evaluated at init.

        Any term with coefficient_symbolic that is NOT in L0 (constants),
        L2 (spatial), or expected to vary with time must be evaluable now.
        """
        for eq_idx, eq in enumerate(self._spec.equations):
            for term_idx, term in enumerate(eq.rhs_terms):
                key = (eq_idx, term_idx)
                if term.coefficient_symbolic is None:
                    continue
                if key in self._constants or key in self._spatial_cache:
                    continue
                if term.time_dependent:
                    # Time-dependent: try evaluating at t=0 to validate
                    self._evaluate_symbolic(term, 0.0)
                    continue
                # Not in any cache and not time-dependent → should have
                # been resolved. Try again — will raise on failure.
                self._evaluate_symbolic(term, 0.0)

    def _evaluate_symbolic(
        self, term: OperatorTerm, t: float
    ) -> float | NDArray[np.float64]:
        """Evaluate a symbolic coefficient expression."""
        sym = term.coefficient_symbolic
        if sym is None:
            return term.coefficient

        # Use coord_arrays only if position-dependent
        coord_arrays = self._coord_arrays if term.position_dependent else None

        return evaluate_coefficient(
            sym, self._parameters, self._coordinates, coord_arrays, t
        )

    def _get_python_expr(self, symbolic_expr: str) -> str:
        """Convert Mathematica → Python with L1 caching."""
        if symbolic_expr not in self._expr_cache:
            self._expr_cache[symbolic_expr] = mathematica_to_python(
                symbolic_expr, self._coordinates
            )
        return self._expr_cache[symbolic_expr]

    def _check_mass_sign(self) -> None:
        """Warn if position-dependent mass-like terms change sign on grid."""
        for eq in self._spec.equations:
            for term in eq.rhs_terms:
                if (
                    term.operator != "identity"
                    or term.field != eq.field_name
                    or term.coefficient_symbolic is None
                    or not term.position_dependent
                    or term.time_dependent
                ):
                    continue

                # Find the spatial cache entry
                for (ei, ti), arr in self._spatial_cache.items():
                    spec_term = self._spec.equations[ei].rhs_terms[ti]
                    if spec_term is term:
                        if float(arr.min()) * float(arr.max()) < 0:
                            warnings.warn(
                                f"Position-dependent mass term "
                                f"'{term.coefficient_symbolic}' for field "
                                f"'{eq.field_name}' changes sign across "
                                f"the grid (min={float(arr.min()):.4g}, "
                                f"max={float(arr.max()):.4g}). This may "
                                f"cause tachyonic instability at locations "
                                f"where the effective mass² is negative.",
                                UserWarning,
                                stacklevel=2,
                            )
                        break
