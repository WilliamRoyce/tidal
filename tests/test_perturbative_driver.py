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
        },
    ],
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
        },
    ],
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
                y0,
                grid,
                (0.0, 1.0),
                order=2,
                parameters={"m2": 1.0, "eps": 0.05},
            )

    def test_gate_order2_raises_not_implemented(self) -> None:
        """Even when the spec carries order_in_eps=2 terms, order>=2 is not
        yet implemented (Pass 2 Duhamel against q¹ missing — #273). We
        must raise NotImplementedError loudly rather than silently reuse
        Pass 0 eigendata at every iteration.
        """
        data = copy.deepcopy(_KG_WITH_EPS)
        # Add an order_in_eps=2 term so max_order = 2.
        terms = data["equations"][0]["rhs"]["terms"]  # type: ignore[index]
        terms.append(
            {
                "coefficient": -1.0,
                "operator": "identity",
                "field": "phi_0",
                "coefficient_symbolic": "-eps2",
                "order_in_eps": 2,
            },
        )
        data["metadata"]["perturbation"] = {  # type: ignore[index]
            "small_parameters": ["eps", "eps2"],
            "order": 2,
        }
        spec = _make_spec(data)
        solver = PerturbativeSolver(spec)
        assert solver.max_order == 2
        grid = GridInfo(shape=(32,), bounds=((0.0, 2 * np.pi),), periodic=(True,))
        y0 = _make_ic(spec, grid)
        with pytest.raises(NotImplementedError, match="#273"):
            solver.solve(
                y0,
                grid,
                (0.0, 1.0),
                order=2,
                parameters={"m2": 1.0, "eps": 0.01, "eps2": 0.001},
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
        self,
        solver_setup: dict[str, Any],
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
            ],
        )
        assert ret == 0

    def test_cli_default_when_perturbation_metadata_present(
        self,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
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
            ],
        )
        assert ret == 0
        # Log goes to stderr; check both streams for the banner.
        captured = capsys.readouterr()
        combined = (captured.out + "\n" + captured.err).lower()
        assert "perturbative modal solver" in combined or "order=1" in combined

    def test_cli_order_zero_skips_driver(
        self,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
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
            ],
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

    def test_validity_formula_is_eps_times_omega_times_t(self) -> None:
        """R4.1 / #284: validity_param = ε·ω·t_end, not ε·ω²·t_end.

        The prior formula over-counted by a factor of ω: safe for ω>1
        but dangerous at low frequencies. Validate the current formula
        directly against omega_max and eps.
        """
        spec = _make_spec(_KG_WITH_EPS)
        grid = GridInfo(shape=(16,), bounds=((0.0, 2 * np.pi),), periodic=(True,))
        y0 = _make_ic(spec, grid)
        solver = PerturbativeSolver(spec)
        t_end = 0.5
        eps = 0.05
        res = solver.solve(
            y0,
            grid,
            (0.0, t_end),
            order=1,
            parameters={"m2": 1.0, "eps": eps},
            num_snapshots=3,
        )
        omega_max = float(res.validity["omega_max"])
        expected = eps * omega_max * t_end  # ε·ω·t (R4.1 formula)
        assert res.validity["validity_param"] == pytest.approx(expected, rel=1e-12), (
            f"validity_param={res.validity['validity_param']!r} "
            f"must equal ε·ω·t_end={expected!r}. See #284."
        )


class TestPerturbativeBaseStability:
    """R4.2 / #285: validity monitor detects base-theory tachyons."""

    def test_base_theory_tachyon_detection(self) -> None:
        """A wrong-sign mass (+m²·φ instead of −m²·φ) produces a
        base-theory tachyon at k=0: Re(λ) = +m > 0. The validity
        monitor must flag this at warn-level even when ε is tiny
        (correction_level = ok).

        We use a mild tachyon (m² = 0.1 → Re(λ) = √0.1 ≈ 0.316) and a
        short t_end so the Pass 0 pre-evolution divergence guard does
        not preempt the run. base_stability_param = 0.316 * 1.0 ≈
        0.316 → warn-level (> 1e-10), not error (< 30).
        """
        data = copy.deepcopy(_KG_WITH_EPS)
        # Flip the sign of the -m²·φ term → +m²·φ (wrong-sign mass).
        terms: list[dict[str, Any]] = data["equations"][0]["rhs"]["terms"]  # type: ignore[index, reportUnknownVariableType]
        for t in terms:  # type: ignore[reportUnknownVariableType]
            if t.get("coefficient_symbolic") == "-m2":
                t["coefficient"] = 1.0  # was -1.0
                t["coefficient_symbolic"] = "m2"
        spec = _make_spec(data)
        grid = GridInfo(shape=(16,), bounds=((0.0, 2 * np.pi),), periodic=(True,))
        # Build IC with a DC component so the k=0 tachyonic mode is
        # excited (modal solver prunes all-zero-IC blocks, which would
        # otherwise hide the λ = +m eigenvalue at k=0).
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points
        x = np.linspace(0.0, 2 * np.pi, n, endpoint=False)
        y0 = np.zeros(layout.num_slots * n)
        y0[:n] = 0.01 + np.sin(x)  # DC offset + sin
        solver = PerturbativeSolver(spec)
        res = solver.solve(
            y0,
            grid,
            (0.0, 1.0),  # short t to stay under the divergence guard
            order=1,
            parameters={"m2": 0.1, "eps": 1e-6},
            num_snapshots=3,
        )
        assert res.validity["max_real_lambda"] > 1e-3, (
            "Wrong-sign mass should produce a positive Re(λ); "
            f"got {res.validity['max_real_lambda']}."
        )
        assert res.validity["base_level"] == "warn", (
            f"Tachyonic base with max Re(λ) = "
            f"{res.validity['max_real_lambda']:.3e} should produce "
            f"base_level='warn' (got {res.validity['base_level']!r}). See #285."
        )
        assert res.validity["warn_level"] == "warn", (
            "Base-stability warn must escalate overall warn_level. See #285."
        )
        assert res.validity["dominant_tachyon_field"] is not None, (
            "dominant_tachyon_field must name the offending slot. See #285."
        )


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
        eigendata: dict[str, Any] = res.orders[0]["eigendata"]
        schur_ops: dict[str, Any] | None = eigendata.get("schur_ops")
        assert schur_ops is not None, "Proca spec must have Schur operators"

        recovery_matrix: Any = schur_ops["recovery_matrix"]
        c_names: Any = schur_ops["constraint_field_names"]
        assert "A_0" in c_names
        c_idx = list(c_names).index("A_0")

        y_hat_dyn: Any = pass1.get("y_hat_dyn")
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
                predicted_hat.reshape(rfft_shape),
                s=grid.shape,
                axes=[0],
            ).ravel(),
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


class TestPerturbativeResultLayoutConsistency:
    """v6 R2.2 / #276: PerturbativeResult carries its own spec so the
    CLI / measurement layer doesn't have to rebind to solver.base_spec.
    """

    def test_result_spec_is_base_spec_not_full(self) -> None:
        spec = _make_spec(_KG_WITH_EPS)
        solver = PerturbativeSolver(spec)
        grid = GridInfo(shape=(32,), bounds=((0.0, 2 * np.pi),), periodic=(True,))
        y0 = _make_ic(spec, grid)
        res = solver.solve(
            y0,
            grid,
            (0.0, 0.5),
            order=1,
            parameters={"m2": 1.0, "eps": 0.05},
        )
        # The result must carry both the base_spec (matches y layout) and
        # the full_spec (for diagnostics).
        assert res.spec is not None
        assert res.full_spec is not None
        assert res.spec is solver.base_spec
        assert res.full_spec is solver.full_spec

    def test_result_layout_matches_total_state(self) -> None:
        """The state vector width equals spec's slot count × grid.num_points.

        Pre-R2.2 the CLI could silently pair total["y"] with the full spec
        (larger slot count) and downstream measurements would read garbage
        slot-by-slot. This confirms the layout invariant.
        """
        from tidal.solver.state import StateLayout

        spec = _make_spec(_KG_WITH_EPS)
        solver = PerturbativeSolver(spec)
        grid = GridInfo(shape=(32,), bounds=((0.0, 2 * np.pi),), periodic=(True,))
        y0 = _make_ic(spec, grid)
        res = solver.solve(
            y0,
            grid,
            (0.0, 0.5),
            order=1,
            parameters={"m2": 1.0, "eps": 0.05},
        )
        assert res.spec is not None
        base_layout = StateLayout.from_spec(res.spec, grid.num_points)
        expected_width = base_layout.num_slots * grid.num_points
        assert res.total["y"].shape[1] == expected_width, (
            f"State width {res.total['y'].shape[1]} != "
            f"spec.num_slots * grid.num_points = {expected_width}. "
            f"See #276."
        )


class TestPerturbativeDefensiveHardening:
    """Defensive diagnostics for unsupported pipeline inputs (#273 resolution)."""

    def test_warns_when_spec_has_higher_order_terms_than_requested(
        self,
        caplog: pytest.LogCaptureFixture,
    ) -> None:
        """When a spec carries order_in_eps=2 terms but solve(order=1) is
        requested, the solver must emit a warning listing the dropped terms
        rather than silently ignoring them.  Under TIDAL's linearized-theory
        constraint this path should never be reached in practice, but the
        diagnostic prevents silent wrong physics if it ever is.
        """
        import logging

        data = copy.deepcopy(_KG_WITH_EPS)
        terms = data["equations"][0]["rhs"]["terms"]  # type: ignore[index]
        terms.append(
            {
                "coefficient": -0.5,
                "operator": "identity",
                "field": "phi_0",
                "coefficient_symbolic": "-eps2",
                "order_in_eps": 2,
            },
        )
        data["metadata"]["perturbation"] = {  # type: ignore[index]
            "small_parameters": ["eps", "eps2"],
            "order": 2,
        }
        spec = _make_spec(data)
        solver = PerturbativeSolver(spec)
        assert solver.max_order == 2

        grid = GridInfo(shape=(16,), bounds=((0.0, 2 * np.pi),), periodic=(True,))
        y0 = _make_ic(spec, grid)

        with caplog.at_level(
            logging.WARNING,
            logger="tidal.solver.perturbative_driver",
        ):
            solver.solve(
                y0,
                grid,
                (0.0, 0.5),
                order=1,
                parameters={"m2": 1.0, "eps": 0.01, "eps2": 0.001},
                num_snapshots=3,
            )

        assert any("order_in_eps > 1" in record.message for record in caplog.records), (
            "Expected a warning about dropped higher-order terms"
        )

    def test_validity_warns_when_small_parameter_missing_from_params(
        self,
        caplog: pytest.LogCaptureFixture,
    ) -> None:
        """If a small parameter declared in [perturbation] is absent from the
        runtime parameters dict, the validity monitor must emit a warning
        rather than silently skipping the check.  A typo in the parameter
        name could otherwise mask a non-small value passing unchecked.

        The correction term uses a purely numeric coefficient so the solve can
        complete even without 'eps' in the parameters dict — the warning fires
        inside _compute_validity when it scans the declared small_parameters.
        """
        import logging

        # Spec where the order-1 correction has a numeric-only coefficient so
        # the modal solve can complete without 'eps' in the parameter dict.
        data: dict[str, Any] = {
            "metadata": {
                "source": "inline-test",
                "parameters": {"m2": 1.0},
                "perturbation": {"small_parameters": ["eps"], "order": 1},
            },
            "spacetime": {
                "dimension": 2,
                "signature": [-1, 1],
                "coordinates": ["t", "x"],
            },
            "fields": [{"name": "phi_0", "index": 0, "is_dynamical": True}],
            "equations": [
                {
                    "field": "phi_0",
                    "lhs": {
                        "expression": "d2_t(phi_0)",
                        "order": {"time": 2, "space": 0},
                    },
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": -1.0,
                                "operator": "identity",
                                "field": "phi_0",
                                "coefficient_symbolic": "-m2",
                            },
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian_x",
                                "field": "phi_0",
                            },
                            {
                                # Purely numeric — no symbolic form. The solver
                                # can evaluate this without 'eps' being defined,
                                # so the solve completes and _compute_validity runs.
                                "coefficient": -0.01,
                                "operator": "identity",
                                "field": "phi_0",
                                "order_in_eps": 1,
                            },
                        ],
                    },
                },
            ],
        }
        spec = _make_spec(data)
        solver = PerturbativeSolver(spec)
        grid = GridInfo(shape=(16,), bounds=((0.0, 2 * np.pi),), periodic=(True,))
        y0 = _make_ic(spec, grid)

        with caplog.at_level(
            logging.WARNING,
            logger="tidal.solver.perturbative_driver",
        ):
            solver.solve(
                y0,
                grid,
                (0.0, 0.5),
                order=1,
                parameters={"m2": 1.0},  # 'eps' intentionally omitted
                num_snapshots=3,
            )

        assert any(
            "small parameter" in record.message and "eps" in record.message
            for record in caplog.records
        ), "Expected a warning about missing small parameter 'eps'"

    def test_pass_n_failure_raises_runtime_error(
        self,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Pass-n solver returning success=False must raise RuntimeError loudly.

        Returning a degraded result (Pass 0 only) would silently mask the
        missing perturbative correction and produce physically wrong physics.
        The caller must see the failure.
        """
        spec = _make_spec(_KG_WITH_EPS)
        solver = PerturbativeSolver(spec)
        grid = GridInfo(shape=(16,), bounds=((0.0, 2 * np.pi),), periodic=(True,))
        y0 = _make_ic(spec, grid)

        # Force solve_modal_pass1 to return a failure result
        from tidal.solver import perturbative_driver as pd_module

        def fake_pass1(eigendata, correction_spec, grid_, t_arr, parameters=None):
            return {
                "t": t_arr,
                "y": np.zeros((len(t_arr), 1)),
                "success": False,
                "message": "synthetic failure for test",
                "correction_drops": [],
            }

        monkeypatch.setattr(pd_module, "solve_modal_pass1", fake_pass1)

        with pytest.raises(RuntimeError, match="Perturbative Pass 1 failed"):
            solver.solve(
                y0,
                grid,
                (0.0, 0.5),
                order=1,
                parameters={"m2": 1.0, "eps": 0.01},
                num_snapshots=3,
            )


# ---------------------------------------------------------------------------
# #301 Phase 3 / #303: Pass 0/Pass 1 hierarchy preservation when a small
# parameter enters kinetic_coefficient_symbolic.
# ---------------------------------------------------------------------------


_KG_WITH_KINETIC_EPS: dict[str, object] = {
    "metadata": {
        "source": "inline-test",
        "parameters": {"m2": 1.0, "alpha": 0.05},
        "perturbation": {"small_parameters": ["alpha"], "order": 1},
    },
    "spacetime": {"dimension": 2, "signature": [-1, 1], "coordinates": ["t", "x"]},
    "fields": [{"name": "phi_0", "index": 0, "is_dynamical": True}],
    "equations": [
        {
            "field": "phi_0",
            "lhs": {
                "expression": "d2_t(phi_0)",
                "order": {"time": 2, "space": 0},
                "kinetic_coefficient_symbolic": "1 + alpha",
            },
            "rhs": {
                "type": "linear_combination",
                "terms": [
                    {"coefficient": 1.0, "operator": "laplacian_x", "field": "phi_0"},
                ],
            },
        },
    ],
}


class TestPhase3HierarchyPreservation:
    """Canonicalizing a spec with a small parameter in the kinetic must leave
    Pass 0 independent of that parameter and let Pass 1 recover the exact
    dispersion. Litmus test for #303.
    """

    def _make_spec(self) -> EquationSystem:
        return _make_spec(_KG_WITH_KINETIC_EPS)

    def _make_ic_wave(self, grid: GridInfo) -> np.ndarray:
        from tidal.solver.state import StateLayout

        spec = self._make_spec()
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points
        x = np.linspace(0.0, 2 * np.pi, n, endpoint=False)
        y0 = np.zeros(layout.num_slots * n)
        y0[:n] = np.cos(x)
        return y0

    def test_canonicalization_strips_small_params(self) -> None:
        """After canonicalize_kinetic_for_perturbation, the kinetic must no
        longer reference any small parameter.
        """
        spec = self._make_spec()
        canon = spec.canonicalize_kinetic_for_perturbation(["alpha"])
        kin = canon.equations[0].kinetic_coefficient_symbolic
        assert kin is not None
        assert "alpha" not in kin
        assert kin.strip() in {"1", "1.0"}

    def test_canonicalization_synthesizes_order1_correction(self) -> None:
        """A synthesized order-1 RHS term must appear covering the
        M₁·M₀⁻¹·K₀ contribution (here: ``-alpha·laplacian_x``).
        """
        spec = self._make_spec()
        canon = spec.canonicalize_kinetic_for_perturbation(["alpha"])
        order1_terms = [t for t in canon.equations[0].rhs_terms if t.order_in_eps == 1]
        assert order1_terms, "expected at least one order-1 synthesized term"
        # Synthesized term has coefficient -1 at unit params (−alpha·1·1/1)
        assert any(
            t.operator == "laplacian_x"
            and t.field == "phi_0"
            and abs(t.coefficient + 1.0) < 1e-12
            for t in order1_terms
        ), f"expected -1·laplacian_x(phi_0) synthesis; got {order1_terms!r}"

    def test_canonicalization_is_idempotent(self) -> None:
        """Running the transform twice must produce the same result — after
        the first pass, the kinetic has no small-parameter dependence left,
        so the second pass is a no-op.
        """
        spec = self._make_spec()
        once = spec.canonicalize_kinetic_for_perturbation(["alpha"])
        twice = once.canonicalize_kinetic_for_perturbation(["alpha"])
        assert once.equations[0].kinetic_coefficient_symbolic == (
            twice.equations[0].kinetic_coefficient_symbolic
        )
        assert len(once.equations[0].rhs_terms) == len(twice.equations[0].rhs_terms)

    def test_pass0_independent_of_small_parameter(self) -> None:
        """The heart of #303: ``solve(order=0, α=anything)`` must be identical
        to ``solve(order=0, α=0)``. If α leaks into Pass 0 the two differ.
        """
        spec = self._make_spec()
        solver = PerturbativeSolver(spec)
        grid = GridInfo(
            shape=(64,),
            bounds=((0.0, 2 * np.pi),),
            periodic=(True,),
        )
        y0 = self._make_ic_wave(grid)

        r_alpha = solver.solve(
            y0,
            grid,
            (0.0, 0.3),
            order=0,
            parameters={"m2": 1.0, "alpha": 0.05},
            num_snapshots=4,
        )
        r_zero = solver.solve(
            y0,
            grid,
            (0.0, 0.3),
            order=0,
            parameters={"m2": 1.0, "alpha": 0.0},
            num_snapshots=4,
        )
        np.testing.assert_allclose(
            r_alpha.orders[0]["y"],
            r_zero.orders[0]["y"],
            rtol=1e-10,
            atol=1e-12,
            err_msg=(
                "Pass 0 differs at α=0.05 vs α=0 — the small parameter "
                "has leaked into the baseline eigendecomposition (#303)."
            ),
        )

    def test_pass1_recovers_exact_modal_dispersion(self) -> None:
        """Pass 0 + Pass 1 at finite α should reproduce the exact dispersion
        given by a direct modal solve on the un-canonicalized spec. Error
        bound is O(α²): the Pass-1 expansion drops terms ≳ α²·(base scale).
        """
        from tidal.solver.modal import solve_modal

        spec = self._make_spec()
        solver = PerturbativeSolver(spec)
        grid = GridInfo(
            shape=(64,),
            bounds=((0.0, 2 * np.pi),),
            periodic=(True,),
        )
        y0 = self._make_ic_wave(grid)
        t_span = (0.0, 0.3)
        alpha = 0.05
        params = {"m2": 1.0, "alpha": alpha}

        r_pert = solver.solve(
            y0,
            grid,
            t_span,
            order=1,
            parameters=params,
            num_snapshots=4,
        )
        r_exact = solve_modal(
            spec,
            grid,
            y0,
            t_span,
            parameters=params,
            num_snapshots=4,
        )

        # Expected accuracy: leading omitted term is O(α²·k²·t_end²) per mode
        # amplitude. Roughly: (0.05)² ≈ 2.5e-3. Use 3% to be safe.
        final_pert = r_pert.total["y"][-1][: grid.num_points]
        final_exact = r_exact["y"][-1][: grid.num_points]
        norm = float(np.linalg.norm(final_exact))
        assert norm > 1e-6
        rel_err = float(np.linalg.norm(final_pert - final_exact) / norm)
        assert rel_err < 3e-2, (
            f"Pass-1 recovery error {rel_err:.2e} exceeds α² bound; "
            "check synthesis formula or canonicalization sign convention."
        )


# ---------------------------------------------------------------------------
# GH #380: Pass 0 must be constant-coefficient even when the kinetic
# coefficient is BOTH coordinate-dependent AND carries a small parameter
# (e.g. Euler-Heisenberg in a localized B background).
# ---------------------------------------------------------------------------

_KG_WITH_POS_DEP_KINETIC_EPS: dict[str, object] = {
    "metadata": {
        "source": "inline-test",
        "parameters": {"m2": 1.0, "alpha": 0.05},
        "perturbation": {"small_parameters": ["alpha"], "order": 1},
    },
    "spacetime": {"dimension": 2, "signature": [-1, 1], "coordinates": ["t", "x"]},
    "fields": [{"name": "phi_0", "index": 0, "is_dynamical": True}],
    "equations": [
        {
            "field": "phi_0",
            "lhs": {
                "expression": "d2_t(phi_0)",
                "order": {"time": 2, "space": 0},
                # Coord-dep + small-param: this is the EH-shape that pre-#380
                # got divided through by Wolfram, hiding the perturbative
                # structure inside 1/(1-alpha*x[]^2).
                "kinetic_coefficient_symbolic": "1 - alpha*x[]**2",
            },
            "rhs": {
                "type": "linear_combination",
                "terms": [
                    {"coefficient": 1.0, "operator": "laplacian_x", "field": "phi_0"},
                ],
            },
        },
    ],
}


class TestGH380CoordDepSmallParamKinetic:
    """Coord-dep + small-param kinetic must canonicalize the same way as
    pure-parameter does (M = M₀ + ε·M₁ split, Pass-1 RHS synthesized). The
    Wolfram side must NOT divide the RHS through by such a kinetic
    coefficient (else everything collapses into order_in_eps=0 RHS terms
    with 1/(1-εX) shapes — see #380).
    """

    def test_canonicalization_strips_small_param_from_kinetic(self) -> None:
        spec = _make_spec(_KG_WITH_POS_DEP_KINETIC_EPS)
        canon = spec.canonicalize_kinetic_for_perturbation(["alpha"])
        kin = canon.equations[0].kinetic_coefficient_symbolic
        assert kin is not None
        assert "alpha" not in kin
        assert kin.strip() in {"1", "1.0"}

    def test_canonicalization_synthesizes_pos_dep_correction(self) -> None:
        spec = _make_spec(_KG_WITH_POS_DEP_KINETIC_EPS)
        canon = spec.canonicalize_kinetic_for_perturbation(["alpha"])
        order1 = [t for t in canon.equations[0].rhs_terms if t.order_in_eps == 1]
        assert order1, "expected synthesized order-1 term from coord-dep kinetic"
        # Synthesized coeff: -((-x[]^2) * alpha * 1) / 1 → +alpha*x[]^2,
        # carrying the coordinate x[] in coefficient_symbolic.
        synth = [
            t for t in order1 if t.operator == "laplacian_x" and t.field == "phi_0"
        ]
        assert synth, f"expected laplacian_x(phi_0) synthesis; got {order1!r}"
        joined = " ".join((t.coefficient_symbolic or "") for t in synth)
        # split_small_parameter_kinetic translates `x[]` → `x` before
        # ast.parse (#380); the renamed form propagates into the synthesized
        # coefficient string. mathematica_to_python at eval time is
        # idempotent on bare coord names, so this is safe.
        import re as _re

        assert _re.search(r"\bx\b", joined), (
            f"synthesized term must carry x coordinate; got {joined!r}"
        )
        assert "alpha" in joined, (
            f"synthesized term must carry alpha small param; got {joined!r}"
        )

    def test_base_spec_is_constant_coefficient(self) -> None:
        """After canonicalization, base_spec must have no coord dependence
        in any order_in_eps=0 RHS coefficient — the heart of the #380 fix.
        If this fails, modal+perturbative will hit NotImplementedError at
        Pass 0.
        """
        import re as _re

        spec = _make_spec(_KG_WITH_POS_DEP_KINETIC_EPS)
        base = spec.canonicalize_kinetic_for_perturbation(["alpha"]).base_spec(
            ["alpha"]
        )
        for eq in base.equations:
            for t in eq.rhs_terms:
                cs = t.coefficient_symbolic or ""
                assert "x[]" not in cs, (
                    f"base spec must be constant-coefficient; "
                    f"term {t!r} has Wolfram x[] in coefficient_symbolic={cs!r}"
                )
                assert not _re.search(r"\bx\b", cs), (
                    f"base spec must be constant-coefficient; "
                    f"term {t!r} has bare x coord in coefficient_symbolic={cs!r}"
                )


# ---------------------------------------------------------------------------
# GH #380 follow-up: PerturbativeSolver must raise an actionable error when
# the base spec has legitimately-position-dependent terms (i.e. not just the
# 1/(1-εX) shape that #380 canonicalizes away).
# ---------------------------------------------------------------------------

_KG_WITH_POS_DEP_BASE: dict[str, object] = {
    "metadata": {
        "source": "inline-test",
        "parameters": {"m2": 1.0, "alpha": 0.05},
        "perturbation": {"small_parameters": ["alpha"], "order": 1},
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
                    # Genuinely position-dependent at order_in_eps=0 (not via a
                    # 1/(1-εX) kinetic shape). canonicalize_kinetic_for_perturbation
                    # will pass this through unchanged.
                    {
                        "coefficient": 1.0,
                        "operator": "laplacian_x",
                        "field": "phi_0",
                        "coefficient_symbolic": "1.0 + 0.1 * x[]**2",
                        "coordinate_dependent": ["x"],
                        "order_in_eps": 0,
                    },
                ],
            },
        },
    ],
}


class TestGH380PerturbativeRejectsPosDepBase:
    """PerturbativeSolver should fail loudly at construction when the base
    spec is position-dependent (after canonicalization), with a message that
    points the user at non-perturbative alternatives instead of letting them
    hit the deep modal.py NotImplementedError after a slow simulate setup.
    """

    def test_init_raises_with_actionable_message(self) -> None:
        spec = _make_spec(_KG_WITH_POS_DEP_BASE)
        with pytest.raises(NotImplementedError) as exc_info:
            PerturbativeSolver(spec)
        msg = str(exc_info.value)
        # Must mention the alternative paths so the user knows what to do.
        assert "modal" in msg.lower()
        assert "cvode" in msg.lower() or "ida" in msg.lower()
        assert "position" in msg.lower()
