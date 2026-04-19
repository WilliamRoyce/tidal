"""Unit tests for :class:`tidal.solver.perturbative_driver.PerturbativeSolver`.

Validates the Stage 5 orchestration layer: Pass 0 + Pass 1 composition,
validity-monitor thresholds, and robustness to baseline (no-correction)
theories.
"""

from __future__ import annotations

import copy
from typing import TYPE_CHECKING, Any

import numpy as np
import pytest

from tidal.solver.grid import GridInfo
from tidal.solver.perturbative_driver import (
    PerturbativeResult,
    PerturbativeSolver,
)
from tidal.solver.state import StateLayout
from tidal.symbolic.json_loader import EquationSystem

if TYPE_CHECKING:
    from pathlib import Path

_KG_WITH_EPS: dict[str, object] = {
    "metadata": {
        "source": "inline-test",
        "parameters": {"m2": 1.0, "eps": 0.05},
        "perturbation": {"small_parameters": ["eps"], "order": 1},
    },
    "spacetime": {"dimension": 2, "signature": [-1, 1], "coordinates": ["t", "x"]},
    "fields": [{"name": "phi_0", "index": 0, "is_dynamical": True}],
    "equations": [
        {
            "field": "phi_0",
            "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2, "space": 0}},
            "rhs": {
                "type": "linear_combination",
                "terms": [
                    {
                        "coefficient": -1.0,
                        "operator": "identity",
                        "field": "phi_0",
                        "coefficient_symbolic": "-m2",
                    },
                    {"coefficient": 1.0, "operator": "laplacian_x", "field": "phi_0"},
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
    "coupling": {},
}


_KG_BASELINE_NO_PERT: dict[str, object] = {
    "metadata": {"source": "inline-test", "parameters": {"m2": 1.0}},
    "spacetime": {"dimension": 2, "signature": [-1, 1], "coordinates": ["t", "x"]},
    "fields": [{"name": "phi_0", "index": 0, "is_dynamical": True}],
    "equations": [
        {
            "field": "phi_0",
            "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2, "space": 0}},
            "rhs": {
                "type": "linear_combination",
                "terms": [
                    {
                        "coefficient": -1.0,
                        "operator": "identity",
                        "field": "phi_0",
                        "coefficient_symbolic": "-m2",
                    },
                    {"coefficient": 1.0, "operator": "laplacian_x", "field": "phi_0"},
                ],
            },
        }
    ],
    "coupling": {},
}


def _make_spec(data: dict[str, Any]) -> EquationSystem:
    return EquationSystem.from_dict(copy.deepcopy(data))


def _make_ic(spec: EquationSystem, grid: GridInfo) -> np.ndarray:
    layout = StateLayout.from_spec(spec, grid.num_points)
    n = grid.num_points
    x = np.linspace(0.0, 2 * np.pi, n, endpoint=False)
    y0 = np.zeros(layout.num_slots * n)
    y0[:n] = np.sin(x)
    return y0


class TestPerturbativeSolverAPI:
    def test_reports_max_order(self) -> None:
        spec = _make_spec(_KG_WITH_EPS)
        solver = PerturbativeSolver(spec)
        assert solver.max_order == 1
        assert solver.has_corrections()

    def test_baseline_theory_has_no_corrections(self) -> None:
        spec = _make_spec(_KG_BASELINE_NO_PERT)
        solver = PerturbativeSolver(spec)
        assert solver.max_order == 0
        assert solver.has_corrections() is False

    def test_rejects_order_beyond_spec(self) -> None:
        spec = _make_spec(_KG_WITH_EPS)
        solver = PerturbativeSolver(spec)
        grid = GridInfo(shape=(32,), bounds=((0.0, 2 * np.pi),), periodic=(True,))
        y0 = _make_ic(spec, grid)
        with pytest.raises(ValueError, match="order=2"):
            solver.solve(
                y0, grid, (0.0, 1.0), order=2, parameters={"m2": 1.0, "eps": 0.05}
            )


class TestPerturbativeSolverSolve:
    @pytest.fixture
    def solver_setup(self) -> dict[str, Any]:
        spec = _make_spec(_KG_WITH_EPS)
        grid = GridInfo(shape=(64,), bounds=((0.0, 2 * np.pi),), periodic=(True,))
        y0 = _make_ic(spec, grid)
        return {
            "spec": spec,
            "grid": grid,
            "y0": y0,
            "solver": PerturbativeSolver(spec),
        }

    def test_order_zero_returns_base_only(self, solver_setup: dict[str, Any]) -> None:
        res = solver_setup["solver"].solve(
            solver_setup["y0"],
            solver_setup["grid"],
            (0.0, 1.0),
            order=0,
            parameters={"m2": 1.0, "eps": 0.05},
            num_snapshots=11,
        )
        assert isinstance(res, PerturbativeResult)
        assert len(res.orders) == 1  # only Pass 0
        assert res.validity["warn_level"] == "ok"
        # Combined = Pass 0 (no corrections applied)
        np.testing.assert_allclose(res.total["y"], res.orders[0]["y"], atol=1e-14)

    def test_order_one_adds_correction(self, solver_setup: dict[str, Any]) -> None:
        res = solver_setup["solver"].solve(
            solver_setup["y0"],
            solver_setup["grid"],
            (0.0, 1.5),
            order=1,
            parameters={"m2": 1.0, "eps": 0.05},
            num_snapshots=16,
        )
        assert len(res.orders) == 2  # Pass 0 + Pass 1
        assert res.total["y"].shape == res.orders[0]["y"].shape
        # Correction is non-zero
        correction_norm = float(np.max(np.abs(res.orders[1]["y"])))
        assert correction_norm > 1e-6

        # Combined ≠ Pass 0
        diff = float(np.max(np.abs(res.total["y"] - res.orders[0]["y"])))
        assert diff > 1e-6

    def test_validity_monitor_flags_strong_correction(
        self, solver_setup: dict[str, Any]
    ) -> None:
        """Large ε·ω²·t triggers a warn/error band."""
        # ω_max ≈ √(m² + k_max²) ≈ √(1 + 31²) ≈ 31 (Nyquist for N=64, L=2π)
        # With eps = 0.5 and t_end = 10, validity_param ≈ 0.5 × 31² × 10 ≈ 4800
        res = solver_setup["solver"].solve(
            solver_setup["y0"],
            solver_setup["grid"],
            (0.0, 10.0),
            order=1,
            parameters={"m2": 1.0, "eps": 0.5},
            num_snapshots=3,
        )
        assert res.validity["warn_level"] == "error"
        assert res.validity["validity_param"] > 1.0
        assert res.validity["dominant_parameter"] == "eps"


class TestPerturbativeCLIFlag:
    """Verify `tidal simulate --perturbative-order N` routes through the driver."""

    def _write_spec(self, tmp_path: Path) -> Path:
        import json

        path = tmp_path / "kg_with_eps.json"
        path.write_text(json.dumps(_KG_WITH_EPS))
        return path

    def test_cli_explicit_order_one(self, tmp_path: Path) -> None:
        from tidal.cli import main

        spec_path = self._write_spec(tmp_path)
        ret = main(
            [
                "simulate",
                str(spec_path),
                "--param",
                "m2=1.0",
                "--param",
                "eps=0.01",
                "--t-end",
                "0.5",
                "--grid-shape",
                "16",
                "--snapshots",
                "3",
                "--bounds",
                "0:6.283185",
                "--periodic",
                "--perturbative-order",
                "1",
                "--no-plot",
            ]
        )
        assert ret == 0

    def test_cli_default_when_perturbation_metadata_present(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """With [perturbation] metadata in JSON, default is order=1."""
        from tidal.cli import main

        spec_path = self._write_spec(tmp_path)
        ret = main(
            [
                "simulate",
                str(spec_path),
                "--param",
                "m2=1.0",
                "--param",
                "eps=0.01",
                "--t-end",
                "0.5",
                "--grid-shape",
                "16",
                "--snapshots",
                "3",
                "--bounds",
                "0:6.283185",
                "--periodic",
                "--no-plot",
            ]
        )
        assert ret == 0
        # Log goes to stderr; check both streams for the banner.
        captured = capsys.readouterr()
        combined = (captured.out + "\n" + captured.err).lower()
        assert "perturbative modal solver" in combined or "order=1" in combined

    def test_cli_order_zero_skips_driver(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """--perturbative-order 0 routes through plain solve_modal."""
        from tidal.cli import main

        spec_path = self._write_spec(tmp_path)
        ret = main(
            [
                "simulate",
                str(spec_path),
                "--param",
                "m2=1.0",
                "--param",
                "eps=0.01",
                "--t-end",
                "0.5",
                "--grid-shape",
                "16",
                "--snapshots",
                "3",
                "--bounds",
                "0:6.283185",
                "--periodic",
                "--perturbative-order",
                "0",
                "--no-plot",
            ]
        )
        assert ret == 0
        captured = capsys.readouterr()
        combined = (captured.out + "\n" + captured.err).lower()
        assert "perturbative modal" not in combined


class TestPerturbativeSolverValidity:
    def test_validity_ok_for_small_eps(self) -> None:
        spec = _make_spec(_KG_WITH_EPS)
        grid = GridInfo(shape=(16,), bounds=((0.0, 2 * np.pi),), periodic=(True,))
        y0 = _make_ic(spec, grid)
        solver = PerturbativeSolver(spec)
        res = solver.solve(
            y0,
            grid,
            (0.0, 0.5),
            order=1,
            parameters={"m2": 1.0, "eps": 1e-4},
            num_snapshots=5,
        )
        # ω_max ≈ √(1 + 7²) ≈ 7.1 (Nyquist=8 for N=16, L=2π)
        # validity ≈ 1e-4 × 51 × 0.5 ≈ 2.5e-3 → safely < 0.1
        assert res.validity["warn_level"] == "ok"

    def test_validity_records_small_parameter_values(self) -> None:
        spec = _make_spec(_KG_WITH_EPS)
        grid = GridInfo(shape=(16,), bounds=((0.0, 2 * np.pi),), periodic=(True,))
        y0 = _make_ic(spec, grid)
        solver = PerturbativeSolver(spec)
        res = solver.solve(
            y0,
            grid,
            (0.0, 0.1),
            order=1,
            parameters={"m2": 1.0, "eps": 0.02},
            num_snapshots=3,
        )
        assert res.validity["eps_values"] == {"eps": 0.02}
        assert res.validity["omega_max"] > 0


# --- Gap C: Pass 1 constraint recovery (v6 Phase 6A.2) -----------------

_PROCA_WITH_EPS_SPEC: dict[str, object] = {
    "metadata": {
        "source": "inline-test",
        "parameters": {"m2": 1.0, "eps": 0.01},
        "perturbation": {"small_parameters": ["eps"], "order": 1},
    },
    "spacetime": {"dimension": 2, "signature": [-1, 1], "coordinates": ["t", "x"]},
    "fields": [
        {"name": "A_0", "index": 0, "is_dynamical": True},
        {"name": "A_1", "index": 1, "is_dynamical": True},
    ],
    "equations": [
        {
            "field": "A_0",
            "lhs": {"expression": "A_0", "order": {"time": 0, "space": 0}},
            "rhs": {
                "type": "linear_combination",
                "terms": [
                    {
                        "coefficient": 1.0,
                        "operator": "identity",
                        "field": "A_0",
                        "coefficient_symbolic": "m2",
                    },
                    {"coefficient": -1.0, "operator": "laplacian_x", "field": "A_0"},
                    {"coefficient": 1.0, "operator": "gradient_x", "field": "v_A_1"},
                ],
            },
        },
        {
            "field": "A_1",
            "lhs": {"expression": "d2_t(A_1)", "order": {"time": 2, "space": 0}},
            "rhs": {
                "type": "linear_combination",
                "terms": [
                    {
                        "coefficient": -1.0,
                        "operator": "identity",
                        "field": "A_1",
                        "coefficient_symbolic": "-m2",
                    },
                    {"coefficient": 1.0, "operator": "laplacian_x", "field": "A_1"},
                    {"coefficient": 1.0, "operator": "gradient_x", "field": "v_A_0"},
                    # Order-1 correction: eps·A_1 mass shift.  Drives a
                    # non-zero Pass 1 on A_1; through Schur this propagates
                    # back to A_0¹ via the recovery matrix.
                    {
                        "coefficient": -1.0,
                        "operator": "identity",
                        "field": "A_1",
                        "coefficient_symbolic": "-eps",
                        "order_in_eps": 1,
                    },
                ],
            },
        },
    ],
    "coupling": {},
}


class TestPassOneConstraintRecovery:
    """Gap C: Pass 1 constraint-field recovery via Schur (v6 Phase 6A.2)."""

    def test_pass1_constraint_recovery_via_schur(self) -> None:
        """A_0¹_hat should equal recovery_matrix @ A_1¹_hat (the Schur-base
        O(ε¹) recovery) after Pass 1 runs on a Proca-style constraint spec.
        """
        spec = EquationSystem.from_dict(copy.deepcopy(_PROCA_WITH_EPS_SPEC))
        n_grid = 16
        length = 2 * np.pi
        grid = GridInfo(shape=(n_grid,), bounds=((0.0, length),), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        # Non-trivial IC on A_1 so the correction source is non-zero.
        x = np.linspace(0.0, length, n_grid, endpoint=False)
        y0 = np.zeros(layout.num_slots * grid.num_points)
        a1_slot = layout.field_slot_map["A_1"]
        y0[a1_slot * n_grid : (a1_slot + 1) * n_grid] = 0.2 * np.sin(x)

        solver = PerturbativeSolver(spec)
        res = solver.solve(
            y0,
            grid,
            (0.0, 0.5),
            order=1,
            parameters={"m2": 1.0, "eps": 0.01},
            num_snapshots=5,
        )

        assert len(res.orders) == 2
        pass1 = res.orders[1]

        # The recovery check: at each snapshot, extract the Pass 1
        # A_0 component from the full-layout output, FFT it in space,
        # and verify it matches recovery_matrix @ y_hat_dyn[A_1].
        eigendata = res.orders[0]["eigendata"]
        schur_ops = eigendata.get("schur_ops")
        assert schur_ops is not None, "Proca spec must have Schur operators"

        recovery_matrix = schur_ops["recovery_matrix"]
        c_names = schur_ops["constraint_field_names"]
        assert "A_0" in c_names
        c_idx = list(c_names).index("A_0")

        y_hat_dyn = pass1.get("y_hat_dyn")
        assert y_hat_dyn is not None, (
            "Pass 1 result must carry y_hat_dyn for constraint recovery"
        )

        a0_slot = layout.field_slot_map["A_0"]
        rfft_shape = (n_grid // 2 + 1,)

        # Compare recovery-based prediction against the driver's output
        # for one mid-simulation snapshot.
        ti = 2
        predicted_hat = np.einsum("mcj,jm->cm", recovery_matrix, y_hat_dyn[ti])[c_idx]
        predicted_phys = np.real(
            np.fft.irfftn(
                predicted_hat.reshape(rfft_shape), s=grid.shape, axes=[0]
            ).ravel()
        )
        actual_phys = pass1["y"][ti, a0_slot * n_grid : (a0_slot + 1) * n_grid]
        np.testing.assert_allclose(actual_phys, predicted_phys, atol=1e-12)

    def test_pass1_constraint_is_nontrivial(self) -> None:
        """The correction must actually produce a non-zero A_0¹."""
        spec = EquationSystem.from_dict(copy.deepcopy(_PROCA_WITH_EPS_SPEC))
        n_grid = 16
        length = 2 * np.pi
        grid = GridInfo(shape=(n_grid,), bounds=((0.0, length),), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        x = np.linspace(0.0, length, n_grid, endpoint=False)
        y0 = np.zeros(layout.num_slots * grid.num_points)
        a1_slot = layout.field_slot_map["A_1"]
        y0[a1_slot * n_grid : (a1_slot + 1) * n_grid] = 0.2 * np.sin(x)

        solver = PerturbativeSolver(spec)
        res = solver.solve(
            y0,
            grid,
            (0.0, 0.5),
            order=1,
            parameters={"m2": 1.0, "eps": 0.05},
            num_snapshots=5,
        )
        pass1 = res.orders[1]
        a0_slot = layout.field_slot_map["A_0"]
        a0_p1 = pass1["y"][-1, a0_slot * n_grid : (a0_slot + 1) * n_grid]
        assert np.max(np.abs(a0_p1)) > 1e-10, (
            "Gap C recovery produced identically-zero A_0¹; expected "
            "non-trivial recovery from the A_1 correction."
        )
