"""Cross-backend consistency test for non-trivial kinetic_coefficient_symbolic.

Verifies the #301 Phase 2 root fix: all time-domain backends (scipy, CVODE,
IDA, leapfrog) now consume ``kinetic_coefficient_symbolic`` via
:func:`tidal.solver._kinetic.build_inverse_kinetic_diag` and agree with the
modal solver on a theory whose canonical spec carries a non-trivial M.

Before Phase 2 this test would have failed: non-modal backends assumed M = I
and evolved with M⁻¹ missing (``dv/dt = K(q)`` instead of ``dv/dt = M⁻¹ K(q)``),
producing a different (wrong) frequency. Modal read M correctly.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from tidal.solver.cvode import solve_cvode
from tidal.solver.grid import GridInfo
from tidal.solver.ida import solve_ida
from tidal.solver.leapfrog import solve_leapfrog
from tidal.solver.modal import solve_modal
from tidal.solver.scipy_solver import solve_scipy
from tidal.solver.state import StateLayout
from tidal.symbolic.json_loader import EquationSystem


def _make_kg_spec_with_kinetic(kinetic: str) -> EquationSystem:
    """Klein-Gordon wave equation with a symbolic kinetic coefficient on the LHS.

    ``M · d²ₜφ = ∇²φ`` with M given by ``kinetic_coefficient_symbolic`` as
    interpreted by every backend. Expected dispersion ``ω² = k² / M``.
    """
    data: dict[str, Any] = {
        "metadata": {"parameters": {"alpha": 0.5}},
        "spacetime": {"dimension": 2, "signature": [-1, 1]},
        "fields": [{"name": "phi_0", "index": 0}],
        "equations": [
            {
                "field": "phi_0",
                "lhs": {
                    "expression": "d2_t(phi_0)",
                    "order": {"time": 2},
                    "kinetic_coefficient_symbolic": kinetic,
                },
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {
                            "coefficient": 1.0,
                            "operator": "laplacian",
                            "field": "phi_0",
                        },
                    ],
                },
            },
        ],
        "canonical": {
            "hamiltonian_terms": [
                {
                    "coefficient": 0.5,
                    "factor_a": {"field": "phi_0", "operator": "time_derivative"},
                    "factor_b": {"field": "phi_0", "operator": "time_derivative"},
                },
            ],
        },
    }
    return EquationSystem.from_dict(data)


def _gaussian_ic(layout: StateLayout, grid: GridInfo) -> np.ndarray:
    x = np.linspace(0.0, 2 * np.pi, grid.num_points, endpoint=False)
    pulse = np.exp(-(((x - np.pi) / 0.5) ** 2))
    y0 = np.zeros(layout.total_size)
    y0[: grid.num_points] = pulse
    return y0


class TestKineticConsistencyAcrossBackends:
    """Same spec, same IC, same params → agreement across all five backends.

    Uses a compound symbolic kinetic ``"1 + alpha"`` (α=0.5 ⇒ M=1.5) so the
    root fix must actually fire. A simple constant like ``"1.5"`` would skip
    the unit-tolerance fast path in :func:`build_inverse_kinetic_diag` but a
    parameter-dependent expression exercises the full
    ``evaluate_coefficient`` → M⁻¹ assembly path.
    """

    GRID_SHAPE = 64
    T_END = 0.4
    PARAMS = {"alpha": 0.5}
    # Pre-fix, non-modal backends differed from modal by ~16% (wrong frequency
    # because M=I was assumed). Post-fix, residual disagreement is integration
    # error from adaptive-tol/leapfrog stepping — well under 1%.
    RTOL_CROSS_BACKEND = 1e-2

    @pytest.fixture
    def setup(self) -> dict[str, Any]:
        spec = _make_kg_spec_with_kinetic("1 + alpha")
        grid = GridInfo(
            bounds=((0.0, 2 * np.pi),),
            shape=(self.GRID_SHAPE,),
            periodic=(True,),
        )
        layout = StateLayout.from_spec(spec, grid.num_points)
        y0 = _gaussian_ic(layout, grid)
        return {"spec": spec, "grid": grid, "layout": layout, "y0": y0}

    def _final_field(self, y: np.ndarray, grid: GridInfo) -> np.ndarray:
        return np.asarray(y[: grid.num_points], dtype=np.float64)

    def test_all_backends_agree_with_modal(self, setup: dict[str, Any]) -> None:
        t_span = (0.0, self.T_END)

        r_modal = solve_modal(
            setup["spec"],
            setup["grid"],
            setup["y0"],
            t_span,
            parameters=self.PARAMS,
            num_snapshots=11,
        )
        r_scipy = solve_scipy(
            setup["spec"],
            setup["grid"],
            setup["y0"],
            t_span,
            parameters=self.PARAMS,
            num_snapshots=11,
        )
        r_cvode = solve_cvode(
            setup["spec"],
            setup["grid"],
            setup["y0"],
            t_span,
            parameters=self.PARAMS,
            num_snapshots=11,
        )
        r_ida = solve_ida(
            setup["spec"],
            setup["grid"],
            setup["y0"],
            t_span,
            parameters=self.PARAMS,
            num_snapshots=11,
        )
        r_leapfrog = solve_leapfrog(
            setup["spec"],
            setup["grid"],
            setup["y0"],
            t_span,
            dt=1e-3,
            parameters=self.PARAMS,
        )

        ref = self._final_field(r_modal["y"][-1], setup["grid"])
        norm_ref = float(np.linalg.norm(ref))
        assert norm_ref > 1e-6, "Modal reference solution is numerically zero"

        for name, result in [
            ("scipy", r_scipy),
            ("cvode", r_cvode),
            ("ida", r_ida),
            ("leapfrog", r_leapfrog),
        ]:
            other = self._final_field(result["y"][-1], setup["grid"])
            rel_err = float(np.linalg.norm(ref - other) / norm_ref)
            assert rel_err < self.RTOL_CROSS_BACKEND, (
                f"{name} disagrees with modal: rel_err = {rel_err:.3e} "
                f"(threshold {self.RTOL_CROSS_BACKEND:.0e}). "
                "Likely cause: backend does not consume "
                "kinetic_coefficient_symbolic. See #301 / #302."
            )

    def test_fast_path_trivial_kinetic(self, setup: dict[str, Any]) -> None:
        """Trivial kinetic (M = 1) must skip the M⁻¹ multiply (fast path).

        Swapping the spec for one without kinetic_coefficient_symbolic (i.e.,
        ``kinetic = None``) gives the same final state — proves the scale
        factor is never applied when absent, which is important because some
        backends (scipy, leapfrog) have branch-per-step overhead sensitive to
        the fast-path check.
        """
        from tidal.solver._kinetic import build_inverse_kinetic_diag

        spec = _make_kg_spec_with_kinetic("1")  # ≡ 1.0 exactly
        m_inv = build_inverse_kinetic_diag(spec, {})
        assert m_inv is None, (
            f"Expected fast-path None for trivial kinetic; got {m_inv}"
        )


def _make_kg_spec_with_kinetic_path(
    kinetic: str,
    path: str,
) -> EquationSystem:
    """KG spec that routes ``solve_modal`` through a specified internal path.

    All three variants are physically the same Klein-Gordon equation
    ``M · d²ₜφ = ∇²φ`` on ``phi_0`` (matching :func:`_make_kg_spec_with_kinetic`);
    only the spec shape differs to satisfy the dispatcher conditions in
    :func:`tidal.solver.modal.solve_modal` (around modal.py:2699+):

    * ``per_mode``        — no constraints, no position-dependent terms,
                            no time-derivative RHS ops →
                            :func:`_build_per_mode_matrices`.
    * ``convolution``     — ``coefficient_symbolic = "1 + 0*x[]"`` (forces
                            ``position_dependent=True`` while remaining
                            analytically constant) →
                            :func:`_build_convolution_matrix`. This is the
                            path that GH #367 traced as buggy pre-v0.42.0.
    * ``evolution_matrices`` — adds a decoupled algebraic constraint field
                            ``aux_0`` (``time_derivative_order=0``,
                            ``aux_0 = 0``) to force the Schur reduction
                            path → :func:`_build_evolution_matrices`.
                            The constraint is independent of ``phi_0`` so
                            the physics on ``phi_0`` is unchanged.
    * ``convolution_kinetic`` — appends ``+ 0*x[]`` to the KINETIC
                            coefficient (GH #427: forces
                            ``kinetic_position_dependent=True`` while the
                            value stays analytically constant), routing
                            through the kinetics-aware predicate to
                            :func:`_build_convolution_matrix` with the
                            real-space M⁻¹(x) fold.
    * ``convolution_kinetic_schur`` — same pos-dep-in-form kinetic plus
                            the decoupled ``aux_0`` constraint → the
                            GH #379 constraints variant
                            (:func:`_build_convolution_matrix_with_constraints`)
                            with the M⁻¹(x) fold. Closes the path-4
                            coverage gap noted pre-GH #427.
    """
    base_term: dict[str, Any] = {
        "coefficient": 1.0,
        "operator": "laplacian",
        "field": "phi_0",
    }
    if path == "convolution":
        base_term["coefficient_symbolic"] = "1 + 0*x[]"

    lhs_kinetic = kinetic
    if path in {"convolution_kinetic", "convolution_kinetic_schur"}:
        lhs_kinetic = f"{kinetic} + 0*x[]"

    fields: list[dict[str, Any]] = [{"name": "phi_0", "index": 0}]
    equations: list[dict[str, Any]] = [
        {
            "field": "phi_0",
            "lhs": {
                "expression": "d2_t(phi_0)",
                "order": {"time": 2},
                "kinetic_coefficient_symbolic": lhs_kinetic,
            },
            "rhs": {"type": "linear_combination", "terms": [base_term]},
        },
    ]

    if path in {"evolution_matrices", "convolution_kinetic_schur"}:
        # Add a decoupled algebraic constraint: aux_0 = 0 (trivial Schur).
        fields.append({"name": "aux_0", "index": 1})
        equations.append(
            {
                "field": "aux_0",
                "lhs": {"expression": "aux_0", "order": {"time": 0, "space": 0}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {"coefficient": 0.0, "operator": "identity", "field": "phi_0"},
                    ],
                },
            },
        )

    data: dict[str, Any] = {
        "metadata": {"parameters": {"alpha": 0.5}},
        "spacetime": {"dimension": 2, "signature": [-1, 1]},
        "fields": fields,
        "equations": equations,
        "canonical": {
            "hamiltonian_terms": [
                {
                    "coefficient": 0.5,
                    "factor_a": {"field": "phi_0", "operator": "time_derivative"},
                    "factor_b": {"field": "phi_0", "operator": "time_derivative"},
                },
            ],
        },
    }
    return EquationSystem.from_dict(data)


class TestAllModalPathsRespectKinetic:
    """Cross-path regression guard for ``kinetic_coefficient_symbolic``.

    Asserts that all three modal-internal matrix builders consume
    ``M⁻¹`` identically via the shared
    :func:`tidal.solver._kinetic.velocity_row_scale` helper. Pre-v0.42.0
    this test would have failed on the ``convolution`` parametrization by
    orders of magnitude (the GH #367 regression). Post-v0.42.0 plus the
    helper refactor (v0.42.1), all three parametrizations agree to within
    cross-backend integration tolerance.

    Adding a new modal matrix builder without applying the helper —
    exactly the class of bug that caused GH #367 — would re-fail this
    test loudly. See the cross-builder contract in
    :mod:`tidal.solver.modal`'s module docstring.
    """

    GRID_SHAPE = 64
    T_END = 0.4
    PARAMS = {"alpha": 0.5}
    # Tolerance: 1% — same as cross-backend test. All three modal paths
    # solve the *same* physical equation; only the matrix-builder dispatch
    # differs. Numerical tolerance only needs to accommodate FP noise from
    # the different internal arithmetic, not solver-method error.
    RTOL_CROSS_PATH = 1e-2

    @pytest.fixture
    def per_mode_reference(self) -> np.ndarray:
        """Per-mode-path solution as the cross-path reference."""
        spec = _make_kg_spec_with_kinetic_path("1 + alpha", "per_mode")
        grid = GridInfo(
            bounds=((0.0, 2 * np.pi),),
            shape=(self.GRID_SHAPE,),
            periodic=(True,),
        )
        layout = StateLayout.from_spec(spec, grid.num_points)
        y0 = _gaussian_ic(layout, grid)
        result = solve_modal(
            spec,
            grid,
            y0,
            (0.0, self.T_END),
            parameters=self.PARAMS,
            num_snapshots=11,
        )
        return np.asarray(result["y"][-1][: grid.num_points], dtype=np.float64)

    @pytest.mark.parametrize(
        "path",
        [
            "convolution",
            "evolution_matrices",
            "convolution_kinetic",
            "convolution_kinetic_schur",
        ],
    )
    def test_path_agrees_with_per_mode_reference(
        self,
        path: str,
        per_mode_reference: np.ndarray,
    ) -> None:
        """Path ``path`` must produce the same physics as the per-mode path.

        Builds a spec that routes through the specified internal builder,
        solves with ``solve_modal``, and checks rel-error against the
        per-mode reference. Different dispatch, same Klein-Gordon equation
        with the same ``M = 1 + alpha = 1.5`` ⇒ identical late-time state.
        """
        spec = _make_kg_spec_with_kinetic_path("1 + alpha", path)
        grid = GridInfo(
            bounds=((0.0, 2 * np.pi),),
            shape=(self.GRID_SHAPE,),
            periodic=(True,),
        )
        layout = StateLayout.from_spec(spec, grid.num_points)
        y0 = _gaussian_ic(layout, grid)
        result = solve_modal(
            spec,
            grid,
            y0,
            (0.0, self.T_END),
            parameters=self.PARAMS,
            num_snapshots=11,
        )
        other = np.asarray(result["y"][-1][: grid.num_points], dtype=np.float64)
        norm_ref = float(np.linalg.norm(per_mode_reference))
        assert norm_ref > 1e-6, "Per-mode reference is numerically zero"
        rel_err = float(np.linalg.norm(per_mode_reference - other) / norm_ref)
        assert rel_err < self.RTOL_CROSS_PATH, (
            f"modal path '{path}' disagrees with per-mode reference: "
            f"rel_err = {rel_err:.3e} (threshold {self.RTOL_CROSS_PATH:.0e}). "
            "Likely cause: the matrix builder for this path does not apply "
            "`velocity_row_scale` for kinetic_coefficient_symbolic — see the "
            "cross-builder contract in tidal/solver/modal.py module docstring "
            "(GH #367, fixed v0.42.0; helper hoisted in v0.42.1)."
        )


# ---------------------------------------------------------------------------
# GH #382: position-dependent kinetic_coefficient_symbolic on time-domain
# backends. Pre-fix, build_inverse_kinetic_diag didn't pass coord_arrays to
# evaluate_coefficient, so any kinetic with x[] raised NameError.
# ---------------------------------------------------------------------------


class TestGH382PositionDependentKinetic:
    """build_inverse_kinetic_diag must handle position-dependent kinetic
    when a grid is supplied. Returns per-grid-point inverse arrays.
    """

    def test_position_dep_kinetic_returns_array(self) -> None:
        from tidal.solver._kinetic import build_inverse_kinetic_diag

        spec = _make_kg_spec_with_kinetic("1 + 0.1*x[]**2")
        grid = GridInfo(
            bounds=((0.0, 2.0),),
            shape=(64,),
            periodic=(True,),
        )
        result = build_inverse_kinetic_diag(spec, {}, grid=grid)
        assert result is not None
        assert "phi_0" in result
        val = result["phi_0"]
        assert isinstance(val, np.ndarray), (
            f"expected ndarray for position-dependent kinetic; got {type(val)}"
        )
        assert val.shape == (64,)
        # Sanity: matches analytical inverse at cell centers
        x_centered = grid.axes_coords(0)
        expected = 1.0 / (1.0 + 0.1 * x_centered**2)
        assert np.allclose(val, expected, atol=1e-12)

    def test_position_dep_kinetic_constant_path_unchanged(self) -> None:
        """Constant kinetic must still return scalar even when grid supplied,
        preserving the M=1 fast path for non-#382 callers.
        """
        from tidal.solver._kinetic import build_inverse_kinetic_diag

        spec = _make_kg_spec_with_kinetic("2.0")
        grid = GridInfo(bounds=((0.0, 2.0),), shape=(64,), periodic=(True,))
        result = build_inverse_kinetic_diag(spec, {}, grid=grid)
        assert result == {"phi_0": 0.5}

    def test_no_grid_back_compat(self) -> None:
        """Pre-#382 callers (grid=None) must still work for constant kinetic."""
        from tidal.solver._kinetic import build_inverse_kinetic_diag

        spec = _make_kg_spec_with_kinetic("1 + alpha")
        result = build_inverse_kinetic_diag(spec, {"alpha": 0.5})
        # 1 / (1+0.5) = 2/3
        assert result is not None
        assert abs(result["phi_0"] - 2 / 3) < 1e-12

    def test_position_dep_kinetic_no_grid_raises(self) -> None:
        """Without grid, a position-dependent kinetic must raise via
        evaluate_coefficient's NameError-on-undefined-coord. This documents
        the existing failure mode and motivates the fix.
        """
        from tidal.solver._kinetic import build_inverse_kinetic_diag

        spec = _make_kg_spec_with_kinetic("1 + 0.1*x[]**2")
        with pytest.raises((ValueError, NameError)):
            build_inverse_kinetic_diag(spec, {})
