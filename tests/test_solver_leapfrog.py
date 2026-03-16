"""Tests for tidal.solver.leapfrog — Störmer-Verlet and Yoshida symplectic integrators."""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from tidal.solver.grid import GridInfo
from tidal.solver.leapfrog import YOSHIDA_WEIGHTS, compute_force, solve_leapfrog
from tidal.solver.state import StateLayout
from tidal.symbolic.json_loader import EquationSystem


def _make_kg_spec() -> EquationSystem:
    """Klein-Gordon wave equation: d²φ/dt² = ∇²φ - m²φ."""
    data: dict[str, Any] = {
        "spacetime": {"dimension": 2, "signature": [-1, 1]},
        "fields": [{"name": "phi_0", "index": 0}],
        "equations": [
            {
                "field": "phi_0",
                "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2}},
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
            }
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


class TestLeapfrogBasic:
    def test_zero_state_stays_zero(self) -> None:
        """Zero state should remain zero (trivial equilibrium)."""
        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(32,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        y0 = np.zeros(layout.total_size)
        result = solve_leapfrog(
            spec,
            grid,
            y0,
            t_span=(0.0, 1.0),
            dt=0.01,
        )
        assert result["success"]
        np.testing.assert_allclose(result["y"][-1], 0, atol=1e-14)

    def test_standing_wave(self) -> None:
        """Standing wave sin(x)cos(t) should be preserved symplectically."""
        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(64,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        n = grid.num_points
        x = grid.axes_coords(0)

        # IC: phi = sin(x), pi = 0 (standing wave)
        y0 = np.zeros(layout.total_size)
        y0[0:n] = np.sin(x)

        result = solve_leapfrog(
            spec,
            grid,
            y0,
            t_span=(0.0, 1.0),
            dt=0.01,
            snapshot_interval=1.0,
        )
        assert result["success"]

        # At t=1, analytic: phi = sin(x)*cos(1) ≈ 0.5403*sin(x)
        phi_final = result["y"][-1][0:n]
        expected = np.sin(x) * np.cos(1.0)
        np.testing.assert_allclose(phi_final, expected, atol=0.01)


class TestLeapfrogEnergy:
    def test_energy_conservation(self) -> None:
        """Leapfrog should preserve the discrete Hamiltonian."""
        from tidal.solver.operators import laplacian

        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(64,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        n = grid.num_points
        x = grid.axes_coords(0)
        dx = grid.dx[0]

        # IC: phi = sin(x), pi = 0
        y0 = np.zeros(layout.total_size)
        y0[0:n] = np.sin(x)

        # Collect the discrete Hamiltonian using the SAME laplacian as the
        # solver. H = ½ ∫ [π² - φ·∇²φ] dx (integration by parts on the
        # discrete 3-point stencil, NOT the central-difference gradient).
        energies: list[float] = []

        def _energy_callback(_t: float, y: np.ndarray) -> None:
            phi = y[0:n]
            pi = y[n : 2 * n]
            lap_phi = laplacian(phi, grid)
            energy = 0.5 * dx * float(np.sum(pi**2 - phi * lap_phi))
            energies.append(energy)

        result = solve_leapfrog(
            spec,
            grid,
            y0,
            t_span=(0.0, 10.0),
            dt=0.01,
            snapshot_interval=0.1,
            snapshot_callback=_energy_callback,
        )
        assert result["success"]
        assert len(energies) > 10

        # Symplectic: shadow Hamiltonian preserved to O(dt²) ≈ 1e-4.
        e0 = energies[0]
        max_drift = max(abs(e - e0) / e0 for e in energies)
        assert max_drift < 1e-4, f"Energy drift {max_drift:.2e} exceeds tolerance"


class TestLeapfrogValidation:
    def test_warns_on_constraint_system(self) -> None:
        """Leapfrog should warn (not reject) systems with constraints."""
        data: dict[str, Any] = {
            "spacetime": {"dimension": 3, "signature": [-1, 1, 1]},
            "fields": [
                {"name": "A_0", "index": 0},
                {"name": "A_1", "index": 1},
            ],
            "equations": [
                {
                    "field": "A_0",
                    "lhs": {"expression": "0", "order": {"time": 0}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian",
                                "field": "A_0",
                            },
                        ],
                    },
                },
                {
                    "field": "A_1",
                    "lhs": {"expression": "d2_t(A_1)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian",
                                "field": "A_1",
                            },
                        ],
                    },
                },
            ],
        }
        spec = EquationSystem.from_dict(data)
        grid = GridInfo(bounds=((0, 10), (0, 10)), shape=(8, 8), periodic=(True, True))
        layout = StateLayout.from_spec(spec, grid.num_points)

        y0 = np.zeros(layout.total_size)
        with pytest.warns(UserWarning, match="frozen"):
            result = solve_leapfrog(spec, grid, y0, t_span=(0.0, 0.1), dt=0.01)
        assert result["success"]

    def test_rejects_first_order_system(self) -> None:
        """Leapfrog should reject systems with first-order equations."""
        data: dict[str, Any] = {
            "spacetime": {"dimension": 3, "signature": [-1, 1, 1]},
            "fields": [{"name": "T_0", "index": 0}],
            "equations": [
                {
                    "field": "T_0",
                    "lhs": {"expression": "d_t(T_0)", "order": {"time": 1}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian",
                                "field": "T_0",
                            },
                        ],
                    },
                },
            ],
        }
        spec = EquationSystem.from_dict(data)
        grid = GridInfo(bounds=((0, 10), (0, 10)), shape=(8, 8), periodic=(True, True))
        layout = StateLayout.from_spec(spec, grid.num_points)

        y0 = np.zeros(layout.total_size)
        with pytest.raises(ValueError, match="first-order"):
            solve_leapfrog(spec, grid, y0, t_span=(0.0, 1.0), dt=0.01)

    def test_snapshot_callback(self) -> None:
        """Snapshot callback should be called at correct intervals."""
        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(16,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        y0 = np.zeros(layout.total_size)
        callback_times: list[float] = []

        def _cb(t: float, _y: np.ndarray) -> None:
            callback_times.append(t)

        solve_leapfrog(
            spec,
            grid,
            y0,
            t_span=(0.0, 1.0),
            dt=0.1,
            snapshot_interval=0.5,
            snapshot_callback=_cb,
        )
        # Should have t=0.0, t≈0.5, t≈1.0 (3 snapshots)
        assert len(callback_times) >= 3


class TestLeapfrogSnapshotCount:
    """Verify integer-based snapshot indexing produces exact counts."""

    def test_snapshot_count_exact(self) -> None:
        """Leapfrog produces exactly the expected number of snapshots."""
        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(32,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        y0 = np.zeros(layout.total_size)
        dt = 0.01
        t_end = 1.0
        snapshot_interval = 0.1

        result = solve_leapfrog(
            spec,
            grid,
            y0,
            t_span=(0.0, t_end),
            dt=dt,
            snapshot_interval=snapshot_interval,
        )
        assert result["success"]

        # Expected: int(1.0 / 0.1) + 1 = 11 snapshots (t=0, 0.1, ..., 1.0)
        expected = int(t_end / snapshot_interval) + 1
        assert len(result["t"]) == expected, (
            f"Expected {expected} snapshots, got {len(result['t'])}"
        )

    def test_snapshot_count_when_dt_exceeds_interval(self) -> None:
        """When dt > snapshot_interval, one snapshot per step."""
        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(32,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        y0 = np.zeros(layout.total_size)
        dt = 0.15  # larger than snapshot_interval
        t_end = 1.0
        snapshot_interval = 0.1

        result = solve_leapfrog(
            spec,
            grid,
            y0,
            t_span=(0.0, t_end),
            dt=dt,
            snapshot_interval=snapshot_interval,
        )
        assert result["success"]

        # With dt=0.15 > interval=0.1, each step triggers a snapshot.
        # n_steps = int(1.0/0.15) = 6, plus initial = 7, plus final ≈ 8
        n_steps = int(t_end / dt)
        n_snapshots = len(result["t"])
        # Should be close to n_steps + 1 (+ possible final)
        assert n_snapshots >= n_steps + 1
        assert n_snapshots <= n_steps + 2  # at most +1 final save


# ===================================================================
# Yoshida 4th-order integrator tests
# ===================================================================
# Ref: Yoshida (1990), Physics Letters A 150(5-7), pp. 262-268.


class TestYoshidaConstants:
    """Verify Yoshida weight properties."""

    def test_weights_sum_to_one(self) -> None:
        """Yoshida weights must sum to 1 (consistency condition)."""
        assert abs(sum(YOSHIDA_WEIGHTS) - 1.0) < 1e-15

    def test_weight_symmetry(self) -> None:
        """w₁ = w₃ (time-reversibility)."""
        w1, _w2, w3 = YOSHIDA_WEIGHTS
        assert w1 == w3

    def test_middle_weight_negative(self) -> None:
        """w₂ < 0 (backward sub-step that cancels O(dt²) error)."""
        assert YOSHIDA_WEIGHTS[1] < 0


class TestYoshidaBasic:
    def test_zero_state_stays_zero(self) -> None:
        """Zero state should remain zero with Yoshida integrator."""
        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(32,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        y0 = np.zeros(layout.total_size)
        result = solve_leapfrog(
            spec,
            grid,
            y0,
            t_span=(0.0, 1.0),
            dt=0.01,
            order=4,
        )
        assert result["success"]
        assert "Yoshida" in result["message"]
        np.testing.assert_allclose(result["y"][-1], 0, atol=1e-14)

    def test_standing_wave(self) -> None:
        """Yoshida should preserve standing wave sin(x)cos(t) accurately."""
        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(64,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        n = grid.num_points
        x = grid.axes_coords(0)

        y0 = np.zeros(layout.total_size)
        y0[0:n] = np.sin(x)

        result = solve_leapfrog(
            spec,
            grid,
            y0,
            t_span=(0.0, 1.0),
            dt=0.01,
            snapshot_interval=1.0,
            order=4,
        )
        assert result["success"]

        # At t=1, analytic: phi = sin(x)*cos(1) ≈ 0.5403*sin(x)
        phi_final = result["y"][-1][0:n]
        expected = np.sin(x) * np.cos(1.0)
        np.testing.assert_allclose(phi_final, expected, atol=0.01)

    def test_rejects_invalid_order(self) -> None:
        """Order must be 2 or 4."""
        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(16,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        y0 = np.zeros(layout.total_size)

        with pytest.raises(ValueError, match="must be 2 or 4"):
            solve_leapfrog(spec, grid, y0, t_span=(0.0, 1.0), dt=0.01, order=3)


class TestYoshidaEnergy:
    def test_energy_conservation(self) -> None:
        """Yoshida 4th-order should conserve the discrete Hamiltonian."""
        from tidal.solver.operators import laplacian

        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(64,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        n = grid.num_points
        x = grid.axes_coords(0)
        dx = grid.dx[0]

        y0 = np.zeros(layout.total_size)
        y0[0:n] = np.sin(x)

        energies: list[float] = []

        def _energy_callback(_t: float, y: np.ndarray) -> None:
            phi = y[0:n]
            pi = y[n : 2 * n]
            lap_phi = laplacian(phi, grid)
            energy = 0.5 * dx * float(np.sum(pi**2 - phi * lap_phi))
            energies.append(energy)

        result = solve_leapfrog(
            spec,
            grid,
            y0,
            t_span=(0.0, 10.0),
            dt=0.01,
            snapshot_interval=0.1,
            snapshot_callback=_energy_callback,
            order=4,
        )
        assert result["success"]
        assert len(energies) > 10

        # Yoshida shadow Hamiltonian error is O(dt⁴) ≈ 1e-8 — much tighter
        # than Verlet's O(dt²) ≈ 1e-4.
        e0 = energies[0]
        max_drift = max(abs(e - e0) / e0 for e in energies)
        assert max_drift < 1e-6, f"Energy drift {max_drift:.2e} exceeds tolerance"


class TestYoshidaConvergenceOrder:
    """Verify that Yoshida achieves 4th-order temporal convergence.

    Measures temporal error by comparing against a fine-dt reference solution
    on the SAME spatial grid, eliminating spatial error entirely.

    Ref: Yoshida (1990), Physics Letters A 150(5-7), pp. 262-268.
    """

    @staticmethod
    def _make_convergence_setup(
        order: int,
    ) -> tuple[EquationSystem, GridInfo, np.ndarray, np.ndarray]:
        """Compute a reference solution at very fine dt."""
        spec = _make_kg_spec()
        # N=64 is enough — spatial error cancels via reference comparison
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(64,), periodic=(True,))
        n = grid.num_points
        x = grid.axes_coords(0)

        y0 = np.zeros(2 * n)
        y0[0:n] = np.sin(x)

        t_end = 2.0
        # Reference: same order integrator at very fine dt (temporal error negligible)
        dt_ref = 0.001
        if order == 4:
            dt_ref = 0.002  # Yoshida already very accurate at this dt
        ref = solve_leapfrog(
            spec,
            grid,
            y0,
            t_span=(0.0, t_end),
            dt=dt_ref,
            snapshot_interval=t_end,
            order=order,
        )
        ref_phi = ref["y"][-1][0:n]
        return spec, grid, y0, ref_phi

    def test_convergence_order_4(self) -> None:
        """Error should scale as dt⁴ for Yoshida integrator."""
        spec, grid, y0, ref_phi = self._make_convergence_setup(order=4)
        n = grid.num_points

        # CFL at N=64: dx ≈ 0.098.  Yoshida max |sub_dt| = dt * 1.70,
        # so outer dt must be < 0.098 / 1.70 ≈ 0.058.
        dts = [0.05, 0.025, 0.0125, 0.00625]
        errors: list[float] = []

        for dt in dts:
            result = solve_leapfrog(
                spec,
                grid,
                y0,
                t_span=(0.0, 2.0),
                dt=dt,
                snapshot_interval=2.0,
                order=4,
            )
            phi_final = result["y"][-1][0:n]
            errors.append(float(np.max(np.abs(phi_final - ref_phi))))

        # Compute log-log slopes between consecutive points
        slopes = []
        for i in range(len(errors) - 1):
            slope = np.log(errors[i] / errors[i + 1]) / np.log(dts[i] / dts[i + 1])
            slopes.append(slope)

        # All slopes should be close to 4 (Yoshida 4th-order)
        avg_slope = np.mean(slopes)
        assert avg_slope > 3.5, (
            f"Yoshida convergence order {avg_slope:.2f} < 3.5 "
            f"(slopes: {[f'{s:.2f}' for s in slopes]})"
        )

    def test_convergence_order_2_for_verlet(self) -> None:
        """Verify order 2 baseline: error should scale as dt²."""
        spec, grid, y0, ref_phi = self._make_convergence_setup(order=2)
        n = grid.num_points

        # CFL at N=64: dx ≈ 0.098.
        dts = [0.08, 0.04, 0.02, 0.01]
        errors: list[float] = []

        for dt in dts:
            result = solve_leapfrog(
                spec,
                grid,
                y0,
                t_span=(0.0, 2.0),
                dt=dt,
                snapshot_interval=2.0,
                order=2,
            )
            phi_final = result["y"][-1][0:n]
            errors.append(float(np.max(np.abs(phi_final - ref_phi))))

        slopes = []
        for i in range(len(errors) - 1):
            slope = np.log(errors[i] / errors[i + 1]) / np.log(dts[i] / dts[i + 1])
            slopes.append(slope)

        avg_slope = np.mean(slopes)
        assert 1.5 < avg_slope < 2.5, (
            f"Verlet convergence order {avg_slope:.2f} not near 2 "
            f"(slopes: {[f'{s:.2f}' for s in slopes]})"
        )


class TestYoshidaFusedKicks:
    """Verify the fused-kick optimization for Yoshida order 4.

    The Forest & Ruth (1990) fused-kick scheme merges adjacent half-kicks
    between sub-steps, reducing 6N → 3N+1 force evaluations. Tests verify
    both the force count and numerical equivalence with the unfused form.
    """

    def test_force_eval_count(self) -> None:
        """Yoshida should use exactly 3N+1 force evaluations for N steps."""
        from unittest.mock import patch

        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(16,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        y0 = np.zeros(layout.total_size)
        x = grid.axes_coords(0)
        y0[0:16] = np.sin(x)

        call_count = 0
        original_compute_force = compute_force

        def counting_force(*args, **kwargs) -> np.ndarray:  # noqa: ANN002, ANN003
            nonlocal call_count
            call_count += 1
            return original_compute_force(*args, **kwargs)

        with patch("tidal.solver.leapfrog.compute_force", side_effect=counting_force):
            n_steps = 10
            t_end = 0.01 * n_steps  # dt=0.01, n_steps=10
            solve_leapfrog(
                spec,
                grid,
                y0,
                t_span=(0.0, t_end),
                dt=0.01,
                snapshot_interval=t_end,
                order=4,
            )

        # Fused-kick: 1 (pre-loop) + 3 per step = 3N+1
        expected = 3 * n_steps + 1
        assert call_count == expected, (
            f"Force evaluations: {call_count}, expected {expected} (3x{n_steps}+1)"
        )

    def test_verlet_force_eval_count(self) -> None:
        """Verlet should use exactly N+1 force evaluations for N steps."""
        from unittest.mock import patch

        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(16,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        y0 = np.zeros(layout.total_size)
        x = grid.axes_coords(0)
        y0[0:16] = np.sin(x)

        call_count = 0
        original_compute_force = compute_force

        def counting_force(*args, **kwargs) -> np.ndarray:  # noqa: ANN002, ANN003
            nonlocal call_count
            call_count += 1
            return original_compute_force(*args, **kwargs)

        with patch("tidal.solver.leapfrog.compute_force", side_effect=counting_force):
            n_steps = 10
            t_end = 0.01 * n_steps
            solve_leapfrog(
                spec,
                grid,
                y0,
                t_span=(0.0, t_end),
                dt=0.01,
                snapshot_interval=t_end,
                order=2,
            )

        expected = n_steps + 1
        assert call_count == expected, (
            f"Force evaluations: {call_count}, expected {expected} ({n_steps}+1)"
        )

    def test_numerical_equivalence_with_unfused(self) -> None:
        """Fused-kick results should match an explicit unfused implementation.

        Runs both the fused code path (current implementation) and a manually
        unfused Yoshida loop, verifying results match to near machine precision.
        """
        from tidal.solver.fields import FieldSet
        from tidal.solver.leapfrog import _half_kick

        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(32,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        drift_pairs = layout.drift_slot_pairs

        n = grid.num_points
        x = grid.axes_coords(0)
        y0 = np.zeros(layout.total_size)
        y0[0:n] = np.sin(x)
        y0[n : 2 * n] = 0.5 * np.cos(x)

        dt = 0.01
        t_end = 1.0

        # --- Fused version (current implementation) ---
        result_fused = solve_leapfrog(
            spec,
            grid,
            y0.copy(),
            t_span=(0.0, t_end),
            dt=dt,
            snapshot_interval=t_end,
            order=4,
        )
        y_fused = result_fused["y"][-1]

        # --- Unfused version (explicit 6 force evals per step) ---
        import math

        n_steps = max(1, math.ceil(t_end / dt - 1e-10))
        y = y0.copy()
        t = 0.0
        force_buf = np.zeros(layout.total_size)
        fieldset_buf = FieldSet.zeros(layout, grid.shape)

        for _ in range(n_steps):
            t_sub = t
            for w in YOSHIDA_WEIGHTS:
                sub_dt = w * dt
                compute_force(
                    spec,
                    layout,
                    grid,
                    None,
                    y,
                    t_sub,
                    None,
                    out=force_buf,
                    fieldset=fieldset_buf,
                )
                _half_kick(y, force_buf, sub_dt, layout)
                for fs, vs in drift_pairs:
                    y[fs] += sub_dt * y[vs]
                t_sub += sub_dt
                compute_force(
                    spec,
                    layout,
                    grid,
                    None,
                    y,
                    t_sub,
                    None,
                    out=force_buf,
                    fieldset=fieldset_buf,
                )
                _half_kick(y, force_buf, sub_dt, layout)
            t += dt

        # Should match to near machine precision (only FP summation order differs)
        max_diff = float(np.max(np.abs(y_fused - y)))
        assert max_diff < 1e-12, (
            f"Fused vs unfused max difference: {max_diff:.2e} exceeds 1e-12"
        )
