"""Stage 4 end-to-end tests for the Pass 1 Duhamel modal solver.

Validates :func:`tidal.solver.modal.solve_modal_pass1` against closed-form
analytical solutions of driven linear oscillators: the base Pass 0
system (unperturbed Klein-Gordon) is solved by the existing modal
eigendecomposition, and a small correction term — here an additional
mass shift ``ε·φ`` — is solved by the Duhamel kernel. The sum of the
two passes is compared to the exact full solution at the same ε.

Reference: the iterative expansion
    q(t) = q⁽⁰⁾(t) + ε·q⁽¹⁾(t) + O(ε²)
with q⁽⁰⁾ driving the Pass 1 source. See Parker & Simon 1993
(literature/gr-qc_9211002/) and the v6 implementation plan.
"""

from __future__ import annotations

import copy
from typing import Any, cast

import numpy as np
import pytest

from tidal.solver.grid import GridInfo
from tidal.solver.modal import solve_modal, solve_modal_pass1
from tidal.solver.state import StateLayout
from tidal.symbolic.json_loader import EquationSystem

# ---------------------------------------------------------------------------
# Shared KG base spec (same shape as test_solver_modal._KG_1D_SPEC)
# ---------------------------------------------------------------------------

_KG_BASE: dict[str, object] = {
    "metadata": {"source": "inline-test", "parameters": {"m2": 1.0, "eps": 0.1}},
    "spacetime": {"dimension": 2, "signature": [-1, 1], "coordinates": ["t", "x"]},
    "fields": [{"name": "phi_0", "index": 0, "is_dynamical": True}],
    "equations": [
        {
            "field": "phi_0",
            "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2, "space": 0}},
            "rhs": {
                "type": "linear_combination",
                "terms": [
                    # Base mass: -m2 · φ  (order 0)
                    {
                        "coefficient": -1.0,
                        "operator": "identity",
                        "field": "phi_0",
                        "coefficient_symbolic": "-m2",
                    },
                    # Base laplacian: +∇²φ  (order 0)
                    {"coefficient": 1.0, "operator": "laplacian_x", "field": "phi_0"},
                    # Correction: -eps · φ  (order 1 — mass shift)
                    {
                        "coefficient": -1.0,
                        "operator": "identity",
                        "field": "phi_0",
                        "coefficient_symbolic": "-eps",
                        "order_in_eps": 1,
                    },
                ],
            },
        }
    ],
    "coupling": {"mass_matrix_symbolic": [["-m2 - eps"]]},
}


def _make_spec(data: dict[str, Any]) -> EquationSystem:
    return EquationSystem.from_dict(copy.deepcopy(data))


class TestSolveModalPass1AnalyticalMatch:
    """Driven Klein-Gordon: analytical answer via iterative expansion."""

    @pytest.fixture
    def setup(self) -> dict[str, Any]:
        spec = _make_spec(_KG_BASE)
        n = 64
        length = 2 * np.pi
        grid = GridInfo(shape=(n,), bounds=((0.0, length),), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        m2 = 1.0
        eps = 0.05
        k_mode = 1.0
        x = np.linspace(0.0, length, n, endpoint=False)

        # IC: phi = sin(x), v = 0 — single Fourier mode
        y0 = np.zeros(layout.num_slots * grid.num_points)
        y0[:n] = np.sin(k_mode * x)
        return {
            "spec": spec,
            "grid": grid,
            "layout": layout,
            "y0": y0,
            "m2": m2,
            "eps": eps,
            "k_mode": k_mode,
            "x": x,
            "n": n,
        }

    def test_pass1_combined_matches_exact_to_eps2(self, setup: dict[str, Any]) -> None:
        """q⁽⁰⁾ + ε·q⁽¹⁾ matches exact full-ε solution to O(ε²)."""
        spec = setup["spec"]
        grid = setup["grid"]
        y0 = setup["y0"]
        m2 = setup["m2"]
        eps = setup["eps"]
        setup["k_mode"]
        setup["x"]
        n = setup["n"]
        t_end = 2.0
        num_snap = 21

        base_spec = spec.filter_by_order(0)
        correction_spec = spec.filter_by_order(1)

        pass0 = cast(
            "dict[str, Any]",
            solve_modal(
                base_spec,
                grid,
                y0,
                t_span=(0.0, t_end),
                parameters={"m2": m2, "eps": eps},
                num_snapshots=num_snap,
                return_eigendata=True,
            ),
        )
        t_eval = pass0["t"]
        eigendata = pass0["eigendata"]

        pass1 = solve_modal_pass1(
            eigendata,
            correction_spec,
            grid,
            t_eval,
            parameters={"m2": m2, "eps": eps},
        )

        # Combined approximation: q_total = q⁽⁰⁾ + q⁽¹⁾
        # (ε is already baked into the correction coefficients by the
        # coefficient evaluator, so no extra ε multiplier is needed.)
        q_total = pass0["y"] + pass1["y"]

        # Exact: evolve the FULL spec (both order=0 and order=1 terms) to
        # ω_full = sqrt(k² + m² + ε). Use modal solver as ground truth.
        full = cast(
            "dict[str, Any]",
            solve_modal(
                spec,
                grid,
                y0,
                t_span=(0.0, t_end),
                parameters={"m2": m2, "eps": eps},
                num_snapshots=num_snap,
            ),
        )

        # The iterative approximation error is O(ε²). With ε=0.05 and
        # ω·t ~ (sqrt(1+1+0.05))·2 ≈ 2.9, the expected truncation error
        # is roughly ε² · (ω·t) · max|q| ≈ 2.5e-3 · 3 ≈ 7e-3. Use 1e-2.
        final_err = float(np.max(np.abs(q_total[-1, :n] - full["y"][-1, :n])))
        assert final_err < 1e-2, (
            f"Pass 0 + Pass 1 truncation error {final_err:.2e} exceeds O(ε²) bound"
        )

        # Sanity: the pure Pass 0 error is the omitted O(ε) mass-shift term,
        # much larger than the Pass 1-corrected error. Verify that the
        # correction substantially improves accuracy.
        pass0_only_err = float(np.max(np.abs(pass0["y"][-1, :n] - full["y"][-1, :n])))
        assert final_err < 0.25 * pass0_only_err, (
            f"Pass 1 failed to improve accuracy: "
            f"pass0_err={pass0_only_err:.3e}, combined_err={final_err:.3e}"
        )

    def test_pass1_zero_source_returns_zero(self, setup: dict[str, Any]) -> None:
        """With no order-1 terms, Pass 1 must return the zero solution."""
        spec = setup["spec"]
        grid = setup["grid"]
        y0 = setup["y0"]

        base_spec = spec.filter_by_order(0)
        empty_spec = spec.filter_by_order(99)  # no such order → all-empty rhs
        pass0 = cast(
            "dict[str, Any]",
            solve_modal(
                base_spec,
                grid,
                y0,
                t_span=(0.0, 1.0),
                parameters={"m2": 1.0, "eps": 0.05},
                num_snapshots=5,
                return_eigendata=True,
            ),
        )
        pass1 = solve_modal_pass1(
            pass0["eigendata"],
            empty_spec,
            grid,
            pass0["t"],
            parameters={"m2": 1.0, "eps": 0.05},
        )
        assert np.allclose(pass1["y"], 0.0, atol=1e-14)

    def test_pass1_initial_condition_is_zero(self, setup: dict[str, Any]) -> None:
        """Pass 1 by construction satisfies q⁽¹⁾(0) = 0."""
        spec = setup["spec"]
        grid = setup["grid"]
        y0 = setup["y0"]

        base_spec = spec.filter_by_order(0)
        correction_spec = spec.filter_by_order(1)
        pass0 = cast(
            "dict[str, Any]",
            solve_modal(
                base_spec,
                grid,
                y0,
                t_span=(0.0, 1.0),
                parameters={"m2": 1.0, "eps": 0.05},
                num_snapshots=3,
                return_eigendata=True,
            ),
        )
        pass1 = solve_modal_pass1(
            pass0["eigendata"],
            correction_spec,
            grid,
            pass0["t"],
            parameters={"m2": 1.0, "eps": 0.05},
        )
        # At t = t_eval[0], the Pass 1 correction is zero (IC).
        assert np.max(np.abs(pass1["y"][0])) < 1e-14
