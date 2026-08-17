"""Tests for tidal.solver.validation — extracted spec validation."""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import pytest

from tidal.solver.coefficients import CoefficientEvaluator
from tidal.solver.grid import GridInfo
from tidal.solver.operators import AxisBCSpec, SideBCSpec
from tidal.solver.validation import (
    check_cfl_stability,
    check_mass_sign,
    check_robin_stability,
    validate_field_references,
    validate_operator_dimensions,
)
from tidal.symbolic.json_loader import EquationSystem

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_spec(
    terms: list[dict[str, Any]],
    *,
    dim: int = 2,
    field_name: str = "phi_0",
    order: int = 2,
) -> EquationSystem:
    """Build a minimal spec."""
    signature = [-1] + [1] * (dim - 1)
    lhs_expr = "d2_t(phi_0)" if order == 2 else "0"
    data: dict[str, Any] = {
        "spacetime": {"dimension": dim, "signature": signature},
        "fields": [{"name": field_name, "index": 0}],
        "equations": [
            {
                "field": field_name,
                "lhs": {"expression": lhs_expr, "order": {"time": order}},
                "rhs": {"type": "linear_combination", "terms": terms},
            },
        ],
    }
    return EquationSystem.from_dict(data)


# ---------------------------------------------------------------------------
# Tests — operator dimensions
# ---------------------------------------------------------------------------


class TestOperatorDimensions:
    def test_laplacian_in_1d_passes(self) -> None:
        """Laplacian works in any dimension."""
        spec = _make_spec(
            [{"coefficient": 1.0, "operator": "laplacian", "field": "phi_0"}],
            dim=2,
        )
        validate_operator_dimensions(spec)  # Should not raise

    def test_gradient_z_in_2d_raises(self) -> None:
        """gradient_z requires 3D but spec is 2D."""
        spec = _make_spec(
            [{"coefficient": 1.0, "operator": "gradient_z", "field": "phi_0"}],
            dim=3,  # 2+1D → spatial_dim = 2
        )
        with pytest.raises(ValueError, match="gradient_z"):
            validate_operator_dimensions(spec)

    def test_gradient_y_in_2d_passes(self) -> None:
        """gradient_y requires 2D, spec has spatial_dim=2."""
        spec = _make_spec(
            [{"coefficient": 1.0, "operator": "gradient_y", "field": "phi_0"}],
            dim=3,
        )
        validate_operator_dimensions(spec)  # Should not raise

    def test_gradient_x_in_1d_passes(self) -> None:
        """gradient_x requires 1D, spec has spatial_dim=1."""
        spec = _make_spec(
            [{"coefficient": 1.0, "operator": "gradient_x", "field": "phi_0"}],
            dim=2,
        )
        validate_operator_dimensions(spec)  # Should not raise


# ---------------------------------------------------------------------------
# Tests — field references
# ---------------------------------------------------------------------------


class TestFieldReferences:
    def test_valid_self_reference(self) -> None:
        """Term referencing its own field passes."""
        spec = _make_spec(
            [{"coefficient": 1.0, "operator": "identity", "field": "phi_0"}],
        )
        validate_field_references(spec)  # Should not raise

    def test_valid_cross_reference(self) -> None:
        """Cross-field reference between known fields passes."""
        data: dict[str, Any] = {
            "spacetime": {"dimension": 2, "signature": [-1, 1]},
            "fields": [
                {"name": "phi_0", "index": 0},
                {"name": "chi_0", "index": 1},
            ],
            "equations": [
                {
                    "field": "phi_0",
                    "lhs": {
                        "expression": "d2_t(phi_0)",
                        "order": {"time": 2},
                    },
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "identity",
                                "field": "chi_0",
                            },
                        ],
                    },
                },
                {
                    "field": "chi_0",
                    "lhs": {
                        "expression": "d2_t(chi_0)",
                        "order": {"time": 2},
                    },
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "identity",
                                "field": "phi_0",
                            },
                        ],
                    },
                },
            ],
        }
        spec = EquationSystem.from_dict(data)
        validate_field_references(spec)  # Should not raise


# ---------------------------------------------------------------------------
# Tests — CFL stability
# ---------------------------------------------------------------------------


class TestCFLStability:
    def test_cfl_satisfied(self) -> None:
        """Small enough dt → no warnings."""
        spec = _make_spec(
            [{"coefficient": 1.0, "operator": "laplacian", "field": "phi_0"}],
        )
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(64,), periodic=(True,))
        dx = grid.dx[0]
        dt = dx * 0.5  # Well within CFL
        warnings = check_cfl_stability(spec, grid, dt)
        assert len(warnings) == 0

    def test_cfl_violated(self) -> None:
        """Dt too large → warning returned."""
        spec = _make_spec(
            [{"coefficient": 1.0, "operator": "laplacian", "field": "phi_0"}],
        )
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(64,), periodic=(True,))
        dt = 10.0  # Way too large
        warnings = check_cfl_stability(spec, grid, dt)
        assert len(warnings) == 1
        assert "CFL violation" in warnings[0]

    def test_constraint_ignored(self) -> None:
        """Constraint equations (order=0) are not CFL-checked."""
        spec = _make_spec(
            [{"coefficient": 1.0, "operator": "laplacian", "field": "phi_0"}],
            order=0,
        )
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(64,), periodic=(True,))
        warnings = check_cfl_stability(spec, grid, dt=100.0)
        assert len(warnings) == 0


# ---------------------------------------------------------------------------
# Tests — mass sign diagnostic
# ---------------------------------------------------------------------------


class TestMassSign:
    def test_sign_change_detected(self) -> None:
        """Sin[x] mass term changes sign → warning returned."""
        spec = _make_spec(
            [
                {
                    "coefficient": 1.0,
                    "operator": "identity",
                    "field": "phi_0",
                    "coefficient_symbolic": "Sin[x[]]",
                    "coordinate_dependent": ["x"],
                },
            ],
        )
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(32,), periodic=(True,))
        import warnings as w

        with w.catch_warnings():
            w.simplefilter("ignore")
            coeff_eval = CoefficientEvaluator(spec, grid)
        result = check_mass_sign(coeff_eval, spec)
        assert len(result) == 1
        assert "changes sign" in result[0]

    def test_summed_mass_not_individual_terms(self) -> None:
        """A term crossing zero does not warn when the total mass does not.

        The physical mass of a field is the sum of its ``identity`` self-terms,
        and a component can carry several. ``Sin[x]`` spans [-1, +1] and so
        crosses zero on its own, but added to ``Cos[x]^2 + 2`` (which spans
        [+2, +3]) the total spans [+1, +3.25] and never does. Testing terms
        one at a time reported a tachyonic mass that is not there.
        """
        spec = _make_spec(
            [
                {
                    "coefficient": 1.0,
                    "operator": "identity",
                    "field": "phi_0",
                    "coefficient_symbolic": "Sin[x[]]",
                    "coordinate_dependent": ["x"],
                },
                {
                    "coefficient": 1.0,
                    "operator": "identity",
                    "field": "phi_0",
                    "coefficient_symbolic": "Cos[x[]]^2 + 2",
                    "coordinate_dependent": ["x"],
                },
            ],
        )
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(64,), periodic=(True,))
        import warnings as w

        with w.catch_warnings():
            w.simplefilter("ignore")
            coeff_eval = CoefficientEvaluator(spec, grid)
        assert check_mass_sign(coeff_eval, spec) == []

    def test_summed_mass_detects_what_single_terms_miss(self) -> None:
        """Two single-signed terms warn when their sum crosses zero.

        The converse trap: ``Sin[x] + 1.5`` stays in [+0.5, +2.5] and
        ``Cos[x] - 2`` in [-3, -1], so neither crosses zero alone, but the
        total spans [-1.91, +0.91] and does. Per-term testing stayed silent on
        a genuinely sign-changing mass.
        """
        spec = _make_spec(
            [
                {
                    "coefficient": 1.0,
                    "operator": "identity",
                    "field": "phi_0",
                    "coefficient_symbolic": "Sin[x[]] + 1.5",
                    "coordinate_dependent": ["x"],
                },
                {
                    "coefficient": 1.0,
                    "operator": "identity",
                    "field": "phi_0",
                    "coefficient_symbolic": "Cos[x[]] - 2",
                    "coordinate_dependent": ["x"],
                },
            ],
        )
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(64,), periodic=(True,))
        import warnings as w

        with w.catch_warnings():
            w.simplefilter("ignore")
            coeff_eval = CoefficientEvaluator(spec, grid)
        result = check_mass_sign(coeff_eval, spec)
        assert len(result) == 1
        assert "changes sign" in result[0]

    def test_positive_definite_no_warning(self) -> None:
        """Always-positive mass term → no warning."""
        spec = _make_spec(
            [
                {
                    "coefficient": 1.0,
                    "operator": "identity",
                    "field": "phi_0",
                    "coefficient_symbolic": "Cos[x[]]^2 + 1",
                    "coordinate_dependent": ["x"],
                },
            ],
        )
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(32,), periodic=(True,))
        coeff_eval = CoefficientEvaluator(spec, grid)
        result = check_mass_sign(coeff_eval, spec)
        assert len(result) == 0


# ---------------------------------------------------------------------------
# check_robin_stability
# ---------------------------------------------------------------------------


class TestCheckRobinStability:
    """Tests for Robin BC ghost-cell stability check."""

    def test_no_axis_bcs_returns_empty(self) -> None:
        grid = GridInfo(bounds=((0, 10),), shape=(64,), periodic=(False,))
        assert check_robin_stability(grid) == []

    def test_stable_robin(self) -> None:
        """gamma*dx < 2 should produce no warnings."""
        side = SideBCSpec(kind="robin", gamma=1.0, value=0.0)
        abc = AxisBCSpec(periodic=False, low=side, high=side)
        grid = GridInfo(
            bounds=((0, 10),),
            shape=(64,),
            periodic=(False,),
            axis_bcs=(abc,),
        )
        result = check_robin_stability(grid)
        assert len(result) == 0

    def test_unstable_robin_warns(self) -> None:
        """gamma*dx >= 2 should produce a warning."""
        # dx = 10/4 = 2.5, gamma = 1.0 -> gamma*dx = 2.5 >= 2
        side = SideBCSpec(kind="robin", gamma=1.0, value=0.0)
        abc = AxisBCSpec(periodic=False, low=side, high=side)
        grid = GridInfo(
            bounds=((0, 10),),
            shape=(4,),
            periodic=(False,),
            axis_bcs=(abc,),
        )
        result = check_robin_stability(grid)
        assert len(result) > 0
        assert "unstable" in result[0].lower()

    def test_periodic_axis_skipped(self) -> None:
        """Periodic axes are not checked for Robin stability."""
        abc = AxisBCSpec(periodic=True)
        grid = GridInfo(
            bounds=((0, 10),),
            shape=(4,),
            periodic=(True,),
            axis_bcs=(abc,),
        )
        assert check_robin_stability(grid) == []

    def test_neumann_not_warned(self) -> None:
        """Neumann BC should not trigger Robin stability warning."""
        side = SideBCSpec(kind="neumann")
        abc = AxisBCSpec(periodic=False, low=side, high=side)
        grid = GridInfo(
            bounds=((0, 10),),
            shape=(4,),
            periodic=(False,),
            axis_bcs=(abc,),
        )
        assert check_robin_stability(grid) == []


# ---------------------------------------------------------------------------
# Tests — pointwise mass stability
# ---------------------------------------------------------------------------


# TestPointwiseMassStability removed: check_pointwise_mass_stability
# was deleted. The modal solver eigenvalue pre-check is authoritative.


# ---------------------------------------------------------------------------
# Tests — expression sandboxing (GH #406)
# ---------------------------------------------------------------------------


class TestCoefficientExpressionSandbox:
    """Coefficient strings come from JSON spec files, so they are file input.

    Spec JSONs are committed, shared, and pulled from HPC. Before GH #406 the
    coefficient path evaluated them against a namespace binding the live numpy
    module with no validation at all, so any ``np.<func>(...)`` was reachable.
    """

    @pytest.mark.parametrize(
        "expr",
        [
            "np.zeros(2)",
            "np.sum(np.zeros(3))",
            "np.float64(2.0)",
            "np.pi",
            "np.load('payload.npy', None, True)",
            "np.save('/tmp/out.npy', 1)",
            "np.fromfile('/etc/passwd')",
        ],
    )
    def test_module_access_is_rejected(self, expr: str) -> None:
        """No attribute access survives in a coefficient expression."""
        from tidal.symbolic._eval_utils import evaluate_coefficient

        with pytest.raises(ValueError, match=r"Cannot evaluate|attribute access"):
            evaluate_coefficient(expr, {}, ())

    @pytest.mark.parametrize(
        ("expr", "params", "expected"),
        [
            ("2*m2", {"m2": 3.0}, 6.0),
            ("Pi", {}, math.pi),
            ("Exp[-1]", {}, math.exp(-1)),
            ("Sqrt[4]", {}, 2.0),
            ("Abs[-3]", {}, 3.0),
        ],
    )
    def test_legitimate_expressions_still_evaluate(
        self,
        expr: str,
        params: dict[str, float],
        expected: float,
    ) -> None:
        """The allowlisted maths surface is untouched.

        ``Pi`` in particular: the converter now emits a bare ``pi`` rather than
        ``np.pi``, which is what allows attribute access to be banned outright.
        """
        from tidal.symbolic._eval_utils import evaluate_coefficient

        assert evaluate_coefficient(expr, params, ()) == pytest.approx(expected)


class TestFormulaSandbox:
    """User formulas arrive from TOML sweep configs as well as the CLI."""

    @staticmethod
    def _namespace() -> dict[str, object]:
        from tidal.cli._simulate import FORMULA_NAMESPACE

        return {**FORMULA_NAMESPACE, "x": np.linspace(0.0, 1.0, 5)}

    @pytest.mark.parametrize(
        "expr",
        [
            "np.load('payload.npy', None, True)",
            "np.save('/tmp/out.npy', x)",
            "np.fromfile('/etc/passwd')",
            "().__class__",
            "[c for c in ()]",
            "lambda: 1",
        ],
    )
    def test_escapes_are_rejected(self, expr: str) -> None:
        """Non-allowlisted numpy attributes and escape constructs are refused."""
        from tidal.cli._simulate import safe_formula_eval

        with pytest.raises((ValueError, TypeError)):
            safe_formula_eval(expr, self._namespace())

    @pytest.mark.parametrize(
        "expr",
        [
            "np.exp(-(x - 0.5)**2 / 0.1)",
            "np.sin(2 * np.pi * x)",
            "np.cos(x)",
            "exp(-x**2)",
            "sin(x) * pi",
        ],
    )
    def test_committed_formula_spellings_still_work(self, expr: str) -> None:
        """Every formula spelling used in examples/ keeps working.

        ``np.exp`` / ``np.sin`` / ``np.cos`` / ``np.pi`` are the only numpy
        attributes any committed formula uses, so the allowlist costs nothing.
        """
        from tidal.cli._simulate import safe_formula_eval

        result = np.asarray(safe_formula_eval(expr, self._namespace()), dtype=float)
        assert np.all(np.isfinite(result))

    def test_numpy_module_is_not_reachable(self) -> None:
        """The namespace holds a shim, not numpy — defence beyond validation."""
        import types

        from tidal.cli._simulate import FORMULA_NAMESPACE

        assert not [
            k for k, v in FORMULA_NAMESPACE.items() if isinstance(v, types.ModuleType)
        ]
        with pytest.raises(AttributeError, match="not available in formulas"):
            _ = FORMULA_NAMESPACE["np"].load  # type: ignore[union-attr]
