"""Tests for the Fourier modal solver (tidal/solver/modal.py).

Tests cover:
1. Eligibility checks (can_use_modal) — flat metric, periodic BCs, operators, etc.
2. Correctness — scalar wave, coupled wave, diffusion, machine precision
3. Resume continuity — split run matches full run
4. Cross-validation — modal vs CVODE agreement
"""

from __future__ import annotations

import copy

import numpy as np
import pytest

from tidal.solver.grid import GridInfo
from tidal.solver.modal import can_use_modal, solve_modal
from tidal.solver.state import StateLayout
from tidal.symbolic.json_loader import EquationSystem

# ---------------------------------------------------------------------------
# Inline specs for testing
# ---------------------------------------------------------------------------

# 1D Klein-Gordon: d2t(phi) = laplacian_x(phi) - m2*phi
_KG_1D_SPEC: dict[str, object] = {
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
    "coupling": {"mass_matrix_symbolic": [["-m2"]]},
}

# Coupled scalars: phi + chi with mass and coupling
_COUPLED_SPEC: dict[str, object] = {
    "metadata": {
        "source": "inline-test",
        "parameters": {"mPhi2": 1.0, "mChi2": 4.0, "gCpl": 0.5},
    },
    "spacetime": {"dimension": 2, "signature": [-1, 1], "coordinates": ["t", "x"]},
    "fields": [
        {"name": "phi_0", "index": 0, "is_dynamical": True},
        {"name": "chi_0", "index": 1, "is_dynamical": True},
    ],
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
                        "field": "chi_0",
                        "coefficient_symbolic": "-gCpl",
                    },
                    {
                        "coefficient": -1.0,
                        "operator": "identity",
                        "field": "phi_0",
                        "coefficient_symbolic": "-mPhi2",
                    },
                    {"coefficient": 1.0, "operator": "laplacian_x", "field": "phi_0"},
                ],
            },
        },
        {
            "field": "chi_0",
            "lhs": {"expression": "d2_t(chi_0)", "order": {"time": 2, "space": 0}},
            "rhs": {
                "type": "linear_combination",
                "terms": [
                    {
                        "coefficient": -1.0,
                        "operator": "identity",
                        "field": "chi_0",
                        "coefficient_symbolic": "-mChi2",
                    },
                    {
                        "coefficient": -1.0,
                        "operator": "identity",
                        "field": "phi_0",
                        "coefficient_symbolic": "-gCpl",
                    },
                    {"coefficient": 1.0, "operator": "laplacian_x", "field": "chi_0"},
                ],
            },
        },
    ],
    "coupling": {
        "mass_matrix_symbolic": [["-mPhi2", None], [None, "-mChi2"]],
        "coupling_matrix_symbolic": [[None, "-gCpl"], ["-gCpl", None]],
    },
}

# Constraint spec (time_order=0) — modal should reject
_CONSTRAINT_SPEC: dict[str, object] = {
    "metadata": {"source": "inline-test"},
    "spacetime": {"dimension": 3, "signature": [-1, 1, 1], "coordinates": ["t", "x", "y"]},
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
                    {"coefficient": -1.0, "operator": "laplacian_x", "field": "A_0"},
                ],
            },
        },
        {
            "field": "A_1",
            "lhs": {"expression": "d2_t(A_1)", "order": {"time": 2, "space": 0}},
            "rhs": {
                "type": "linear_combination",
                "terms": [
                    {"coefficient": 1.0, "operator": "laplacian_x", "field": "A_1"},
                ],
            },
        },
    ],
}

# First-order diffusion: dt(u) = D * laplacian(u)
_DIFFUSION_SPEC: dict[str, object] = {
    "metadata": {"source": "inline-test", "parameters": {"D": 0.1}},
    "spacetime": {"dimension": 2, "signature": [-1, 1], "coordinates": ["t", "x"]},
    "fields": [{"name": "u_0", "index": 0, "is_dynamical": True}],
    "equations": [
        {
            "field": "u_0",
            "lhs": {"expression": "d_t(u_0)", "order": {"time": 1, "space": 0}},
            "rhs": {
                "type": "linear_combination",
                "terms": [
                    {
                        "coefficient": 0.1,
                        "operator": "laplacian_x",
                        "field": "u_0",
                        "coefficient_symbolic": "D",
                    },
                ],
            },
        }
    ],
}


def _make_spec(data: dict) -> EquationSystem:
    return EquationSystem.from_dict(data)


def _make_gaussian_ic(
    spec: EquationSystem,
    grid: GridInfo,
    amplitude: float = 0.1,
    width: float = 1.5,
) -> np.ndarray:
    """Gaussian IC on the first field."""
    layout = StateLayout.from_spec(spec, grid.num_points)
    y0 = np.zeros(layout.num_slots * grid.num_points)
    x = np.linspace(grid.bounds[0][0], grid.bounds[0][1], grid.shape[0], endpoint=False)
    center = (grid.bounds[0][0] + grid.bounds[0][1]) / 2
    y0[: grid.num_points] = amplitude * np.exp(-((x - center) ** 2) / (2 * width**2))
    return y0


# =========================================================================
# Eligibility tests (can_use_modal)
# =========================================================================


class TestModalEligibility:
    """Test can_use_modal correctly identifies eligible/ineligible systems."""

    def test_flat_periodic_eligible(self) -> None:
        """Flat metric + periodic BCs → eligible."""
        spec = _make_spec(_KG_1D_SPEC)
        grid = GridInfo(shape=(64,), bounds=((0.0, 10.0),), periodic=(True,))
        assert can_use_modal(spec, grid, None) is True

    def test_coupled_scalars_eligible(self) -> None:
        """Coupled scalars with periodic BCs → eligible."""
        spec = _make_spec(_COUPLED_SPEC)
        grid = GridInfo(shape=(64,), bounds=((0.0, 10.0),), periodic=(True,))
        assert can_use_modal(spec, grid, None) is True

    def test_non_periodic_rejected(self) -> None:
        """Non-periodic BCs → not eligible."""
        spec = _make_spec(_KG_1D_SPEC)
        grid = GridInfo(
            shape=(64,),
            bounds=((0.0, 10.0),),
            periodic=(False,),
            bc=("neumann",),
        )
        assert can_use_modal(spec, grid, None) is False

    def test_constraints_rejected(self) -> None:
        """Constraint equations (time_order=0) → not eligible."""
        spec = _make_spec(_CONSTRAINT_SPEC)
        grid = GridInfo(
            shape=(32, 32),
            bounds=((0.0, 10.0), (0.0, 10.0)),
            periodic=(True, True),
        )
        assert can_use_modal(spec, grid, None) is False

    def test_dissipation_rejected(self) -> None:
        """first_derivative_t operator → not eligible."""
        spec_data = copy.deepcopy(_KG_1D_SPEC)
        spec_data["equations"][0]["rhs"]["terms"].append(  # type: ignore[index]
            {"coefficient": -0.1, "operator": "first_derivative_t", "field": "phi_0"}
        )
        spec = _make_spec(spec_data)
        grid = GridInfo(shape=(64,), bounds=((0.0, 10.0),), periodic=(True,))
        assert can_use_modal(spec, grid, None) is False

    def test_unsupported_operator_rejected(self) -> None:
        """Operator not in modal multiplier registry → not eligible.

        Uses 'derivative_3_x' (valid in json_loader as dynamic pattern,
        but not in _EXACT_MULTIPLIERS).
        """
        spec_data = copy.deepcopy(_KG_1D_SPEC)
        spec_data["equations"][0]["rhs"]["terms"].append(  # type: ignore[index]
            {"coefficient": 1.0, "operator": "derivative_3_x", "field": "phi_0"}
        )
        spec = _make_spec(spec_data)
        grid = GridInfo(shape=(64,), bounds=((0.0, 10.0),), periodic=(True,))
        assert can_use_modal(spec, grid, None) is False

    def test_time_dependent_rejected(self) -> None:
        """Time-dependent coefficient → not eligible."""
        spec_data = copy.deepcopy(_KG_1D_SPEC)
        spec_data["equations"][0]["rhs"]["terms"][0][  # type: ignore[index]
            "coefficient_symbolic"
        ] = "Sin[t[]]"
        spec_data["equations"][0]["rhs"]["terms"][0]["time_dependent"] = True  # type: ignore[index]
        spec = _make_spec(spec_data)
        grid = GridInfo(shape=(64,), bounds=((0.0, 10.0),), periodic=(True,))
        assert can_use_modal(spec, grid, None) is False

    def test_position_dependent_still_eligible(self) -> None:
        """Position-dependent (but time-independent) coefficient → still eligible.

        Handled via convolution in k-space.
        """
        spec_data = copy.deepcopy(_KG_1D_SPEC)
        spec_data["equations"][0]["rhs"]["terms"][0][  # type: ignore[index]
            "coefficient_symbolic"
        ] = "Sin[x[]]"
        spec = _make_spec(spec_data)
        grid = GridInfo(shape=(64,), bounds=((0.0, 10.0),), periodic=(True,))
        assert can_use_modal(spec, grid, None) is True

    def test_first_order_diffusion_eligible(self) -> None:
        """First-order diffusion with periodic BCs → eligible."""
        spec = _make_spec(_DIFFUSION_SPEC)
        grid = GridInfo(shape=(64,), bounds=((0.0, 10.0),), periodic=(True,))
        assert can_use_modal(spec, grid, None) is True

    def test_curved_metric_rejected(self) -> None:
        """Curved metric (volume_element non-None) → not eligible."""
        spec_data = copy.deepcopy(_KG_1D_SPEC)
        # Add canonical with non-None volume_element (signals curved metric)
        spec_data["canonical"] = {
            "hamiltonian_terms": [],
            "volume_element": "Abs[x[]]",
        }
        spec = _make_spec(spec_data)
        grid = GridInfo(shape=(64,), bounds=((0.0, 10.0),), periodic=(True,))
        assert can_use_modal(spec, grid, None) is False

    def test_explicit_modal_on_ineligible_raises(self) -> None:
        """--scheme modal on non-periodic system → RuntimeError."""
        from tidal.cli._simulate import _resolve_scheme

        spec = _make_spec(_KG_1D_SPEC)
        grid = GridInfo(
            shape=(64,),
            bounds=((0.0, 10.0),),
            periodic=(False,),
            bc=("neumann",),
        )
        with pytest.raises(RuntimeError, match="not eligible"):
            _resolve_scheme("modal", spec, grid, None)


# =========================================================================
# Correctness tests
# =========================================================================


class TestModalCorrectness:
    """Test modal solver produces correct results."""

    def test_scalar_wave_dispersion(self) -> None:
        """Single-mode KG wave: omega^2 = k^2 + m^2.

        Initialize with a single Fourier mode and verify the frequency
        matches the exact dispersion relation.
        """
        spec = _make_spec(_KG_1D_SPEC)
        N = 64
        L = 10.0
        grid = GridInfo(shape=(N,), bounds=((0.0, L),), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        # Single-mode IC: phi = cos(2*pi*n*x/L), velocity = 0
        n_mode = 3
        k = 2 * np.pi * n_mode / L
        x = np.linspace(0, L, N, endpoint=False)
        y0 = np.zeros(layout.num_slots * grid.num_points)
        y0[:N] = np.cos(k * x)  # field

        m2 = 1.0
        omega = np.sqrt(k**2 + m2)

        result = solve_modal(
            spec, grid, y0, t_span=(0, 2 * np.pi / omega),
            parameters={"m2": m2}, num_snapshots=51,
        )
        assert result["success"]

        # After one full period, phi should return to IC
        final_phi = result["y"][-1, :N]
        np.testing.assert_allclose(final_phi, y0[:N], atol=1e-12)

    def test_coupled_wave_modal_vs_cvode(self) -> None:
        """Coupled-scalar system: modal agrees with CVODE to solver tolerance."""
        from tidal.solver.cvode import solve_cvode

        spec = _make_spec(_COUPLED_SPEC)
        grid = GridInfo(shape=(64,), bounds=((0.0, 10.0),), periodic=(True,))
        y0 = _make_gaussian_ic(spec, grid)
        params = {"mPhi2": 1.0, "mChi2": 4.0, "gCpl": 0.5}

        result_modal = solve_modal(
            spec, grid, y0, t_span=(0, 3), parameters=params, num_snapshots=31,
        )
        result_cvode = solve_cvode(
            spec, grid, y0, t_span=(0, 3), parameters=params,
            num_snapshots=31, rtol=1e-10, atol=1e-12,
        )

        assert result_modal["success"]
        assert result_cvode["success"]

        # Modal is exact (machine-precision); CVODE has O(rtol) truncation
        # error, so agreement is limited by CVODE's tolerance, not modal's.
        max_diff = np.max(np.abs(result_modal["y"][-1] - result_cvode["y"][-1]))
        assert max_diff < 1e-3, f"Modal vs CVODE max diff: {max_diff:.2e}"

    def test_diffusion_exponential_decay(self) -> None:
        """Diffusion: single mode decays as exp(-D*k^2*t).

        u(x,0) = cos(k*x) → u(x,t) = exp(-D*k^2*t)*cos(k*x).
        """
        spec = _make_spec(_DIFFUSION_SPEC)
        N = 64
        L = 10.0
        D = 0.1
        grid = GridInfo(shape=(N,), bounds=((0.0, L),), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        n_mode = 2
        k = 2 * np.pi * n_mode / L
        x = np.linspace(0, L, N, endpoint=False)
        y0 = np.zeros(layout.num_slots * grid.num_points)
        y0[:N] = np.cos(k * x)

        t_end = 5.0
        result = solve_modal(
            spec, grid, y0, t_span=(0, t_end),
            parameters={"D": D}, num_snapshots=11,
        )
        assert result["success"]

        # Exact solution at t_end
        exact = np.exp(-D * k**2 * t_end) * np.cos(k * x)
        np.testing.assert_allclose(result["y"][-1, :N], exact, atol=1e-12)

    def test_machine_precision(self) -> None:
        """Modal solver achieves machine-precision (~1e-14) error.

        Compare against exact analytical solution for a single Fourier mode.
        """
        spec = _make_spec(_KG_1D_SPEC)
        N = 32
        L = 2 * np.pi
        grid = GridInfo(shape=(N,), bounds=((0.0, L),), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        # IC: phi = sin(x), v = 0
        k = 1.0
        m2 = 1.0
        omega = np.sqrt(k**2 + m2)
        x = np.linspace(0, L, N, endpoint=False)
        y0 = np.zeros(layout.num_slots * grid.num_points)
        y0[:N] = np.sin(k * x)

        t_end = 1.0
        result = solve_modal(
            spec, grid, y0, t_span=(0, t_end),
            parameters={"m2": m2}, num_snapshots=2,
        )

        # Exact solution: phi(t) = cos(omega*t) * sin(x)
        exact = np.cos(omega * t_end) * np.sin(k * x)
        error = np.max(np.abs(result["y"][-1, :N] - exact))
        assert error < 1e-13, f"Machine-precision test failed: error={error:.2e}"

    def test_resume_continuity(self) -> None:
        """Split run (0→5, 5→10) matches full run (0→10) to machine precision."""
        spec = _make_spec(_KG_1D_SPEC)
        grid = GridInfo(shape=(32,), bounds=((0.0, 10.0),), periodic=(True,))
        y0 = _make_gaussian_ic(spec, grid)
        params = {"m2": 1.0}

        # Full run
        result_full = solve_modal(
            spec, grid, y0, t_span=(0, 10), parameters=params, num_snapshots=21,
        )

        # Split run
        result_a = solve_modal(
            spec, grid, y0, t_span=(0, 5), parameters=params, num_snapshots=11,
        )
        y5 = result_a["y"][-1]
        result_b = solve_modal(
            spec, grid, y5, t_span=(5, 10), parameters=params, num_snapshots=11,
        )

        diff = np.max(np.abs(result_full["y"][-1] - result_b["y"][-1]))
        assert diff < 1e-13, f"Resume continuity diff: {diff:.2e}"

    def test_snapshot_callback(self) -> None:
        """Verify snapshot callback is called at each output time."""
        spec = _make_spec(_KG_1D_SPEC)
        grid = GridInfo(shape=(32,), bounds=((0.0, 10.0),), periodic=(True,))
        y0 = _make_gaussian_ic(spec, grid)

        callback_times: list[float] = []

        def cb(t: float, _y: np.ndarray) -> None:
            callback_times.append(t)

        result = solve_modal(
            spec, grid, y0, t_span=(0, 1), parameters={"m2": 1.0},
            num_snapshots=11, snapshot_callback=cb,
        )
        assert result["success"]
        assert len(callback_times) == 11
        np.testing.assert_allclose(callback_times, result["t"], atol=1e-14)

    def test_output_shape(self) -> None:
        """Verify output shape matches (n_snapshots, state_size)."""
        spec = _make_spec(_KG_1D_SPEC)
        N = 32
        grid = GridInfo(shape=(N,), bounds=((0.0, 10.0),), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        y0 = _make_gaussian_ic(spec, grid)

        n_snap = 21
        result = solve_modal(
            spec, grid, y0, t_span=(0, 1), parameters={"m2": 1.0},
            num_snapshots=n_snap,
        )
        assert result["y"].shape == (n_snap, layout.num_slots * N)
        assert result["t"].shape == (n_snap,)

    def test_zero_ic_stays_zero(self) -> None:
        """Zero initial conditions produce zero output (no spurious modes)."""
        spec = _make_spec(_KG_1D_SPEC)
        grid = GridInfo(shape=(32,), bounds=((0.0, 10.0),), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        y0 = np.zeros(layout.num_slots * grid.num_points)

        result = solve_modal(
            spec, grid, y0, t_span=(0, 5), parameters={"m2": 1.0},
            num_snapshots=11,
        )
        assert result["success"]
        assert np.max(np.abs(result["y"])) < 1e-15


# =========================================================================
# Auto-selection integration tests
# =========================================================================


class TestModalAutoSelection:
    """Test modal solver auto-selection in _resolve_scheme."""

    def test_auto_selects_modal_for_eligible(self) -> None:
        """Auto-selection picks modal for flat+periodic+constant system."""
        from tidal.cli._simulate import _resolve_scheme

        spec = _make_spec(_KG_1D_SPEC)
        grid = GridInfo(shape=(64,), bounds=((0.0, 10.0),), periodic=(True,))
        scheme = _resolve_scheme("auto", spec, grid, None)
        assert scheme == "modal"

    def test_auto_selects_ida_for_constraints(self) -> None:
        """Auto-selection picks IDA for constraint equations."""
        from tidal.cli._simulate import _resolve_scheme

        spec = _make_spec(_CONSTRAINT_SPEC)
        grid = GridInfo(
            shape=(32, 32),
            bounds=((0.0, 10.0), (0.0, 10.0)),
            periodic=(True, True),
        )
        scheme = _resolve_scheme("auto", spec, grid, None)
        assert scheme == "ida"

    def test_auto_selects_cvode_for_non_periodic(self) -> None:
        """Auto-selection falls through to CVODE for non-periodic systems."""
        from tidal.cli._simulate import _resolve_scheme

        spec_data = copy.deepcopy(_KG_1D_SPEC)
        # Add canonical to avoid warning
        spec_data["canonical"] = {"hamiltonian_terms": []}
        spec = _make_spec(spec_data)
        grid = GridInfo(
            shape=(64,),
            bounds=((0.0, 10.0),),
            periodic=(False,),
            bc=("neumann",),
        )
        scheme = _resolve_scheme("auto", spec, grid, None)
        assert scheme == "cvode"

    def test_explicit_modal_allowed_for_eligible(self) -> None:
        """Explicit --scheme modal passes validation for eligible system."""
        from tidal.cli._simulate import _resolve_scheme

        spec = _make_spec(_KG_1D_SPEC)
        grid = GridInfo(shape=(64,), bounds=((0.0, 10.0),), periodic=(True,))
        scheme = _resolve_scheme("modal", spec, grid, None)
        assert scheme == "modal"
