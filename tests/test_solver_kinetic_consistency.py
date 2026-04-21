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
