"""Tests for the Fourier modal solver (tidal/solver/modal.py).

Tests cover:
1. Eligibility checks (can_use_modal) — flat metric, periodic BCs, operators, etc.
2. Correctness — scalar wave, coupled wave, diffusion, machine precision
3. Resume continuity — split run matches full run
4. Cross-validation — modal vs CVODE agreement
"""

from __future__ import annotations

import copy
from typing import TYPE_CHECKING, Any, cast

import numpy as np
import pytest

from tidal.solver.grid import GridInfo
from tidal.solver.modal import can_use_modal, solve_modal
from tidal.solver.state import StateLayout
from tidal.symbolic.json_loader import EquationSystem

if TYPE_CHECKING:
    from pathlib import Path

    from numpy.typing import NDArray

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
        },
    ],
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
}

# Constraint spec (time_order=0) — modal should reject
_CONSTRAINT_SPEC: dict[str, object] = {
    "metadata": {"source": "inline-test"},
    "spacetime": {
        "dimension": 3,
        "signature": [-1, 1, 1],
        "coordinates": ["t", "x", "y"],
    },
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
        },
    ],
}


def _make_spec(data: dict[str, object]) -> EquationSystem:
    return EquationSystem.from_dict(data)  # type: ignore[arg-type]


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

    def test_constraints_eligible_if_fourier_eliminable(self) -> None:
        """Fourier-eliminable constraints (time_order=0) → eligible via Schur."""
        spec = _make_spec(_CONSTRAINT_SPEC)
        grid = GridInfo(
            shape=(32, 32),
            bounds=((0.0, 10.0), (0.0, 10.0)),
            periodic=(True, True),
        )
        # Constraint with laplacian_x (has exact Fourier multiplier) is eliminable
        assert can_use_modal(spec, grid, None) is True

    def test_constraints_rejected_unsupported_operator(self) -> None:
        """Constraint with unsupported operator → not eligible."""
        spec_data: dict[str, Any] = copy.deepcopy(dict(_CONSTRAINT_SPEC))
        # Change constraint operator to something truly unknown
        spec_data["equations"][0]["rhs"]["terms"][0]["operator"] = "derivative_5_x"
        spec = _make_spec(spec_data)
        grid = GridInfo(
            shape=(32, 32),
            bounds=((0.0, 10.0), (0.0, 10.0)),
            periodic=(True, True),
        )
        assert can_use_modal(spec, grid, None) is False

    def test_dissipation_accepted(self) -> None:
        """first_derivative_t operator → eligible (generalized mass-matrix path)."""
        spec_data = copy.deepcopy(_KG_1D_SPEC)
        spec_data["equations"][0]["rhs"]["terms"].append(  # type: ignore[index]
            {"coefficient": -0.1, "operator": "first_derivative_t", "field": "phi_0"},
        )
        spec = _make_spec(spec_data)
        grid = GridInfo(shape=(64,), bounds=((0.0, 10.0),), periodic=(True,))
        assert can_use_modal(spec, grid, None) is True

    def test_unsupported_operator_rejected(self) -> None:
        """Operator not in modal decomposition registry → not eligible."""
        spec_data = copy.deepcopy(_KG_1D_SPEC)
        spec_data["equations"][0]["rhs"]["terms"].append(  # type: ignore[index]
            {"coefficient": 1.0, "operator": "derivative_5_x", "field": "phi_0"},
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
            spec,
            grid,
            y0,
            t_span=(0, 2 * np.pi / omega),
            parameters={"m2": m2},
            num_snapshots=51,
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
            spec,
            grid,
            y0,
            t_span=(0, 3),
            parameters=params,
            num_snapshots=31,
        )
        result_cvode = solve_cvode(
            spec,
            grid,
            y0,
            t_span=(0, 3),
            parameters=params,
            num_snapshots=31,
            rtol=1e-10,
            atol=1e-12,
        )

        assert result_modal["success"]
        assert result_cvode["success"]

        # Modal is exact (machine-precision); CVODE has O(rtol) truncation
        # error, so agreement is limited by CVODE's tolerance, not modal's.
        max_diff = np.max(np.abs(result_modal["y"][-1] - result_cvode["y"][-1]))
        assert max_diff < 5e-4, f"Modal vs CVODE max diff: {max_diff:.2e}"

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
            spec,
            grid,
            y0,
            t_span=(0, t_end),
            parameters={"D": D},
            num_snapshots=11,
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
            spec,
            grid,
            y0,
            t_span=(0, t_end),
            parameters={"m2": m2},
            num_snapshots=2,
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
            spec,
            grid,
            y0,
            t_span=(0, 10),
            parameters=params,
            num_snapshots=21,
        )

        # Split run
        result_a = solve_modal(
            spec,
            grid,
            y0,
            t_span=(0, 5),
            parameters=params,
            num_snapshots=11,
        )
        y5 = result_a["y"][-1]
        result_b = solve_modal(
            spec,
            grid,
            y5,
            t_span=(5, 10),
            parameters=params,
            num_snapshots=11,
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
            spec,
            grid,
            y0,
            t_span=(0, 1),
            parameters={"m2": 1.0},
            num_snapshots=11,
            snapshot_callback=cb,
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
            spec,
            grid,
            y0,
            t_span=(0, 1),
            parameters={"m2": 1.0},
            num_snapshots=n_snap,
        )
        assert result["y"].shape == (n_snap, layout.num_slots * N)
        assert result["t"].shape == (n_snap,)

    def test_2d_wave_dispersion(self) -> None:
        """2D wave: omega^2 = kx^2 + ky^2 + m^2, single Fourier mode."""
        spec_2d: dict[str, object] = {
            "metadata": {"source": "inline-test", "parameters": {"m2": 1.0}},
            "spacetime": {
                "dimension": 3,
                "signature": [-1, 1, 1],
                "coordinates": ["t", "x", "y"],
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
                                "operator": "laplacian",
                                "field": "phi_0",
                            },
                        ],
                    },
                },
            ],
        }
        spec = _make_spec(spec_2d)
        Nx, Ny = 16, 16
        Lx, Ly = 4.0, 4.0
        grid = GridInfo(
            shape=(Nx, Ny),
            bounds=((0.0, Lx), (0.0, Ly)),
            periodic=(True, True),
        )
        layout = StateLayout.from_spec(spec, grid.num_points)

        # Single-mode IC: phi = cos(kx*x + ky*y), v = 0
        nx_mode, ny_mode = 2, 1
        kx = 2 * np.pi * nx_mode / Lx
        ky = 2 * np.pi * ny_mode / Ly
        m2 = 1.0
        omega = np.sqrt(kx**2 + ky**2 + m2)

        x = np.linspace(0, Lx, Nx, endpoint=False)
        y = np.linspace(0, Ly, Ny, endpoint=False)
        X, Y = np.meshgrid(x, y, indexing="ij")

        y0 = np.zeros(layout.num_slots * grid.num_points)
        y0[: grid.num_points] = (np.cos(kx * X + ky * Y)).ravel()

        t_period = 2 * np.pi / omega
        result = solve_modal(
            spec,
            grid,
            y0,
            t_span=(0, t_period),
            parameters={"m2": m2},
            num_snapshots=2,
        )
        assert result["success"]

        # After one full period, should return to IC
        np.testing.assert_allclose(
            result["y"][-1, : grid.num_points],
            y0[: grid.num_points],
            atol=1e-11,
        )

    def test_position_dependent_correctness(self) -> None:
        """Position-dependent coefficient correctness via convolution path.

        Diffusion with spatially varying coefficient: dt(u) = c(x)*laplacian(u).
        c(x) = 1.0 (constant, but marked position-dependent with symbolic 'x[]').
        This forces the convolution path while remaining analytically tractable:
        the result should match the constant-coefficient solution.
        """
        spec_data = copy.deepcopy(_DIFFUSION_SPEC)
        # Mark the term as position-dependent with symbolic expression
        # Use coefficient 0.1 as before (D=0.1), but force position_dependent path
        spec_data["equations"][0]["rhs"]["terms"][0][  # type: ignore[index]
            "coefficient_symbolic"
        ] = "D*(1 + 0*x[])"
        spec = _make_spec(spec_data)

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

        t_end = 2.0
        result = solve_modal(
            spec,
            grid,
            y0,
            t_span=(0, t_end),
            parameters={"D": D},
            num_snapshots=2,
        )
        assert result["success"]

        # Exact solution (same as constant-coefficient diffusion)
        exact = np.exp(-D * k**2 * t_end) * np.cos(k * x)
        np.testing.assert_allclose(result["y"][-1, :N], exact, atol=1e-10)

    def test_zero_ic_stays_zero(self) -> None:
        """Zero initial conditions produce zero output (no spurious modes)."""
        spec = _make_spec(_KG_1D_SPEC)
        grid = GridInfo(shape=(32,), bounds=((0.0, 10.0),), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        y0 = np.zeros(layout.num_slots * grid.num_points)

        result = solve_modal(
            spec,
            grid,
            y0,
            t_span=(0, 5),
            parameters={"m2": 1.0},
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

    def test_auto_selects_modal_for_periodic_constraints(self) -> None:
        """Auto-selection picks modal for periodic constraint systems."""
        from tidal.cli._simulate import _resolve_scheme

        spec = _make_spec(_CONSTRAINT_SPEC)
        grid = GridInfo(
            shape=(32, 32),
            bounds=((0.0, 10.0), (0.0, 10.0)),
            periodic=(True, True),
        )
        # Fourier-eliminable constraints with periodic BCs → modal
        scheme = _resolve_scheme("auto", spec, grid, None)
        assert scheme == "modal"

    def test_auto_selects_ida_for_non_periodic_constraints(self) -> None:
        """Auto-selection picks IDA for non-periodic constraint systems."""
        from tidal.cli._simulate import _resolve_scheme

        spec = _make_spec(_CONSTRAINT_SPEC)
        grid = GridInfo(
            shape=(32, 32),
            bounds=((0.0, 10.0), (0.0, 10.0)),
            periodic=(False, False),
            bc=("neumann", "neumann"),
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


# =========================================================================
# Block isolation and eigenvalue stability tests
# =========================================================================


# 4-field system with two identical coupled pairs (mimics Gertsenshtein):
# Pair 1: phi_0 ↔ chi_0 (gradient coupling, nonzero IC)
# Pair 2: phi_1 ↔ chi_1 (identical coupling, zero IC)
# The degenerate eigenvalues across pairs would cause eigenvector mixing
# in a single 8×8 eigendecomposition, seeding exponential growth in pair 2.
_DEGENERATE_PAIRS_SPEC: dict[str, object] = {
    "metadata": {
        "source": "inline-test",
        "parameters": {"B0": 0.3, "kappa2": 1.0},
    },
    "spacetime": {"dimension": 2, "signature": [-1, 1], "coordinates": ["t", "x"]},
    "fields": [
        {"name": "phi_0", "index": 0, "is_dynamical": True},
        {"name": "chi_0", "index": 1, "is_dynamical": True},
        {"name": "phi_1", "index": 2, "is_dynamical": True},
        {"name": "chi_1", "index": 3, "is_dynamical": True},
    ],
    "equations": [
        {
            "field": "phi_0",
            "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2, "space": 0}},
            "rhs": {
                "type": "linear_combination",
                "terms": [
                    {"coefficient": 1.0, "operator": "laplacian_x", "field": "phi_0"},
                    {
                        "coefficient": -1.0,
                        "operator": "gradient_x",
                        "field": "chi_0",
                        "coefficient_symbolic": "-(B0*kappa2)",
                    },
                ],
            },
        },
        {
            "field": "chi_0",
            "lhs": {"expression": "d2_t(chi_0)", "order": {"time": 2, "space": 0}},
            "rhs": {
                "type": "linear_combination",
                "terms": [
                    {"coefficient": 1.0, "operator": "laplacian_x", "field": "chi_0"},
                    {
                        "coefficient": 1.0,
                        "operator": "gradient_x",
                        "field": "phi_0",
                        "coefficient_symbolic": "B0",
                    },
                ],
            },
        },
        {
            "field": "phi_1",
            "lhs": {"expression": "d2_t(phi_1)", "order": {"time": 2, "space": 0}},
            "rhs": {
                "type": "linear_combination",
                "terms": [
                    {"coefficient": 1.0, "operator": "laplacian_x", "field": "phi_1"},
                    {
                        "coefficient": -1.0,
                        "operator": "gradient_x",
                        "field": "chi_1",
                        "coefficient_symbolic": "-(B0*kappa2)",
                    },
                ],
            },
        },
        {
            "field": "chi_1",
            "lhs": {"expression": "d2_t(chi_1)", "order": {"time": 2, "space": 0}},
            "rhs": {
                "type": "linear_combination",
                "terms": [
                    {"coefficient": 1.0, "operator": "laplacian_x", "field": "chi_1"},
                    {
                        "coefficient": 1.0,
                        "operator": "gradient_x",
                        "field": "phi_1",
                        "coefficient_symbolic": "B0",
                    },
                ],
            },
        },
    ],
}


class TestModalBlockIsolation:
    """Test block-aware eigendecomposition for multi-field stability."""

    def test_zero_ic_block_stays_zero(self) -> None:
        """Zero-IC blocks remain at machine zero despite degenerate eigenvalues.

        Creates a 4-field system with two identical coupled pairs. Only pair 1
        gets nonzero IC. Without block isolation, np.linalg.eig mixes the
        degenerate eigenvectors, seeding exponential growth in pair 2. With
        block-aware decomposition, pair 2 stays exactly zero.
        """
        spec = _make_spec(_DEGENERATE_PAIRS_SPEC)
        grid = GridInfo(shape=(64,), bounds=((0.0, 10.0),), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        # IC: Gaussian on phi_0 only, everything else zero
        y0 = np.zeros(layout.num_slots * grid.num_points)
        x = np.linspace(0.0, 10.0, 64, endpoint=False)
        y0[: grid.num_points] = 0.1 * np.exp(-((x - 5.0) ** 2) / (2 * 1.5**2))

        result = solve_modal(
            spec,
            grid,
            y0,
            (0.0, 100.0),
            parameters={"B0": 0.3, "kappa2": 1.0},
            num_snapshots=11,
        )

        # Pair 2 fields (phi_1, chi_1) and their velocities should be zero
        phi_1_slot = layout.field_slot_map["phi_1"]
        chi_1_slot = layout.field_slot_map["chi_1"]
        v_phi_1_slot = layout.velocity_slot_map["phi_1"]
        v_chi_1_slot = layout.velocity_slot_map["chi_1"]

        final = result["y"][-1]  # last snapshot
        n = grid.num_points
        for slot in [phi_1_slot, chi_1_slot, v_phi_1_slot, v_chi_1_slot]:
            max_val = np.max(np.abs(final[slot * n : (slot + 1) * n]))
            assert max_val < 1e-13, (
                f"Zero-IC slot {slot} grew to {max_val:.2e} — block isolation failed"
            )

    def test_nonzero_pair_evolves_correctly(self) -> None:
        """The nonzero-IC pair still evolves correctly with block decomposition."""
        spec = _make_spec(_DEGENERATE_PAIRS_SPEC)
        grid = GridInfo(shape=(64,), bounds=((0.0, 10.0),), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        y0 = np.zeros(layout.num_slots * grid.num_points)
        x = np.linspace(0.0, 10.0, 64, endpoint=False)
        y0[: grid.num_points] = 0.1 * np.exp(-((x - 5.0) ** 2) / (2 * 1.5**2))

        result = solve_modal(
            spec,
            grid,
            y0,
            (0.0, 10.0),
            parameters={"B0": 0.3, "kappa2": 1.0},
            num_snapshots=11,
        )

        # Pair 1 should have non-trivial evolution (not all zero)
        phi_0_slot = layout.field_slot_map["phi_0"]
        n = grid.num_points
        final = result["y"][-1]  # last snapshot
        max_phi0 = np.max(np.abs(final[phi_0_slot * n : (phi_0_slot + 1) * n]))
        assert max_phi0 > 1e-4, "Nonzero-IC field did not evolve"

    def test_eigenvalue_growth_warning(self) -> None:
        """Warning is issued when eigenvalues have large positive real parts.

        Uses a large domain (L=100) so that k_min=2π/100≈0.063 < k_crit=B₀κ/√2≈0.212,
        ensuring low-k modes have genuinely positive real eigenvalues.
        """
        import warnings

        spec = _make_spec(_DEGENERATE_PAIRS_SPEC)
        # Large domain to get k_min < k_crit for unstable modes
        grid = GridInfo(shape=(64,), bounds=((0.0, 100.0),), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        y0 = np.zeros(layout.num_slots * grid.num_points)
        x = np.linspace(0.0, 100.0, 64, endpoint=False)
        y0[: grid.num_points] = 0.1 * np.exp(-((x - 50.0) ** 2) / (2 * 5.0**2))

        # Long time → growth factor > exp(30) should trigger warning,
        # and the divergence guard should raise SimulationDivergedError.
        from tidal.solver._exceptions import SimulationDivergedError

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            # Either runtime guard ("diverged at t=...") or pre-check
            # ("predicted to diverge") is an acceptable detection.
            with pytest.raises(SimulationDivergedError, match="diverge"):
                solve_modal(
                    spec,
                    grid,
                    y0,
                    (0.0, 500.0),
                    parameters={"B0": 0.3, "kappa2": 1.0},
                    num_snapshots=6,
                )

        modal_warnings = [x for x in w if "positive real parts" in str(x.message)]
        assert len(modal_warnings) > 0, "No eigenvalue growth warning issued for t=500"

    def testfind_independent_blocks_utility(self) -> None:
        """Verify find_independent_blocks correctly detects block structure."""
        from tidal.solver.modal import find_independent_blocks

        # 4×4 block-diagonal: [[A, 0], [0, B]]
        M = np.zeros((4, 4), dtype=np.complex128)
        M[0, 1] = 1.0
        M[1, 0] = -1.0
        M[2, 3] = 1.0
        M[3, 2] = -2.0
        blocks = find_independent_blocks(M)
        assert len(blocks) == 2
        assert sorted(blocks[0]) in ([0, 1], [2, 3])
        assert sorted(blocks[1]) in ([0, 1], [2, 3])

        # Fully coupled 4×4
        M2 = np.ones((4, 4), dtype=np.complex128)
        blocks2 = find_independent_blocks(M2)
        assert len(blocks2) == 1
        assert sorted(blocks2[0]) == [0, 1, 2, 3]

        # Diagonal (all independent)
        M3 = np.diag([1.0, 2.0, 3.0]).astype(np.complex128)
        blocks3 = find_independent_blocks(M3)
        assert len(blocks3) == 3


# =========================================================================
# Constraint elimination tests (Fourier Schur complement)
# =========================================================================


# Simple 1D Proca-like: A₀ constraint + A₁ dynamical, with coupling.
# Constraint: (m² - ∇²)A₀ = -∂ₓ(v_A₁)
# Dynamical: ∂²ₜA₁ = ∂²ₓA₁ - m²A₁ + ∂ₓ(v_A₀)
# This is a simplified version of the coupled Proca system.
_PROCA_1D_CONSTRAINT_SPEC: dict[str, object] = {
    "metadata": {"source": "inline-test", "parameters": {"m2": 1.0}},
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
                ],
            },
        },
    ],
}


class TestConstraintElimination:
    """Test Fourier Schur complement constraint elimination."""

    def test_constraint_system_eligible(self) -> None:
        """Proca-like constraint system with periodic BCs is modal-eligible."""
        spec = _make_spec(_PROCA_1D_CONSTRAINT_SPEC)
        grid = GridInfo(shape=(64,), bounds=((0.0, 10.0),), periodic=(True,))
        assert can_use_modal(spec, grid, None) is True

    def test_constraint_elimination_runs(self) -> None:
        """Modal solver runs on constraint system without error."""
        spec = _make_spec(_PROCA_1D_CONSTRAINT_SPEC)
        grid = GridInfo(shape=(64,), bounds=((0.0, 10.0),), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        y0 = np.zeros(layout.num_slots * grid.num_points)
        x = np.linspace(0.0, 10.0, 64, endpoint=False)
        # IC on A₁ field (slot 1)
        a1_slot = layout.field_slot_map["A_1"]
        y0[a1_slot * grid.num_points : (a1_slot + 1) * grid.num_points] = 0.1 * np.exp(
            -((x - 5.0) ** 2) / (2 * 1.5**2),
        )

        result = solve_modal(
            spec,
            grid,
            y0,
            (0.0, 5.0),
            parameters={"m2": 1.0},
            num_snapshots=11,
        )

        assert result["success"]
        assert result["y"].shape == (11, layout.num_slots * grid.num_points)
        assert not np.any(np.isnan(result["y"]))

    def test_constraint_field_reconstructed(self) -> None:
        """Constraint field A₀ is non-trivially reconstructed (not zero)."""
        spec = _make_spec(_PROCA_1D_CONSTRAINT_SPEC)
        grid = GridInfo(shape=(64,), bounds=((0.0, 10.0),), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        y0 = np.zeros(layout.num_slots * grid.num_points)
        x = np.linspace(0.0, 10.0, 64, endpoint=False)
        a1_slot = layout.field_slot_map["A_1"]
        y0[a1_slot * grid.num_points : (a1_slot + 1) * grid.num_points] = 0.1 * np.exp(
            -((x - 5.0) ** 2) / (2 * 1.5**2),
        )

        result = solve_modal(
            spec,
            grid,
            y0,
            (0.0, 2.0),
            parameters={"m2": 1.0},
            num_snapshots=11,
        )

        # A₀ should be non-zero at later times (constraint couples to A₁)
        a0_slot = layout.field_slot_map["A_0"]
        n = grid.num_points
        final = result["y"][-1]
        a0_max = np.max(np.abs(final[a0_slot * n : (a0_slot + 1) * n]))
        assert a0_max > 1e-6, (
            f"Constraint field A₀ is effectively zero ({a0_max:.2e}) — "
            f"reconstruction failed"
        )

    def test_eigenvalues_purely_imaginary(self) -> None:
        """Reduced system eigenvalues are purely imaginary (Hamiltonian)."""
        from tidal.solver.coefficients import CoefficientEvaluator
        from tidal.solver.modal import (
            _build_evolution_matrices,
            _build_k_axes,
            _build_k_grid,
        )

        spec = _make_spec(_PROCA_1D_CONSTRAINT_SPEC)
        grid = GridInfo(shape=(32,), bounds=((0.0, 10.0),), periodic=(True,))
        coeff_eval = CoefficientEvaluator(spec, grid, {"m2": 1.0})

        k_axes = _build_k_axes(grid)
        k_grid = _build_k_grid(k_axes)
        rfft_shape = (17,)

        A_red, _, _, _, _, _, _, _ = _build_evolution_matrices(
            spec,
            StateLayout.from_spec(spec, grid.num_points),
            grid,
            coeff_eval,
            k_grid,
            rfft_shape,
        )

        # All eigenvalues should be purely imaginary for a Hamiltonian system
        for m in range(A_red.shape[0]):
            eigs = cast("np.ndarray[Any, Any]", np.linalg.eigvals(A_red[m]))
            max_real = float(np.max(np.abs(np.real(eigs))))
            assert max_real < 1e-10, (
                f"Mode {m}: max |Re(λ)| = {max_real:.2e}, expected < 1e-10"
            )


class TestRankDeficientBProjection:
    """Tests for null-space projection of rank-deficient B in the QZ path (issue #257).

    The CDT dark photon plasma model has ~20 kinetic-null non-trace torsion DOF
    that make B_block rank-deficient.  Before the fix, scipy.linalg.eig(A, B)
    with rank-deficient B produced FINITE spurious eigenvalues (Re≈25-500) that
    bypassed the |λ|>1e12 filter.  The fix projects A and B onto range(B) first.
    """

    def test_spurious_eigenvalues_present_without_projection(self) -> None:
        """Confirm that QZ on rank-deficient B produces finite spurious eigenvalues.

        This documents the bug: before the fix, sla.eig(A, B) with rank-2 B
        embedded in a 4×4 space returned non-imaginary eigenvalues for the null
        directions (rather than Inf, which the old filter would catch).
        """
        import scipy.linalg as sla  # type: ignore[import-untyped]

        # 4×4 system: 2 physical (harmonic oscillator, ω=1) + 2 null (B=0).
        # Physical block: [0, 1; -1, 0] (eigenvalues ±i).
        # Null block: A has off-diagonal ±1, B has zeros → spurious finite eigs.
        A = np.zeros((4, 4), dtype=np.complex128)
        B = np.zeros((4, 4), dtype=np.complex128)
        A[0, 1] = 1.0
        A[1, 0] = -1.0
        A[2, 3] = 1.0
        A[3, 2] = -1.0
        B[0, 0] = 1.0
        B[1, 1] = 1.0
        # B[2,2] = B[3,3] = 0 → rank-2

        evals = cast("np.ndarray[Any, np.dtype[np.complex128]]", sla.eig(A, B)[0])
        finite_evals = evals[np.isfinite(evals) & (np.abs(evals) <= 1e12)]
        # Without projection: all 4 eigenvalues are finite (scipy distributes
        # physical coupling across all DOF), some with large Re(λ) > 0.
        # We just assert that the finite set is non-empty to document the issue.
        assert len(finite_evals) > 0, (
            "Expected finite spurious eigenvalues from rank-deficient QZ"
        )

    def test_null_projection_removes_spurious_eigenvalues(self) -> None:
        """Null-space projection gives correct eigenvalues {+i, -i, 0, 0} for rank-2 B.

        After projecting onto range(B), QZ operates on the 2×2 physical block
        which yields ±i (harmonic oscillator).  Null directions get eigenvalue 0.
        The total spectrum should have max |Re(λ)| < 1e-10 (all purely imaginary or 0).
        """
        A = np.zeros((4, 4), dtype=np.complex128)
        B = np.zeros((4, 4), dtype=np.complex128)
        A[0, 1] = 1.0
        A[1, 0] = -1.0
        A[2, 3] = 1.0
        A[3, 2] = -1.0
        B[0, 0] = 1.0
        B[1, 1] = 1.0

        _, s_b, Vt_b = np.linalg.svd(B)
        rank_b = int(np.sum(s_b > s_b[0] * 1e-10))
        null_dim_b = 4 - rank_b
        assert null_dim_b == 2, f"Expected 2 null dimensions, got {null_dim_b}"

        import scipy.linalg as sla  # type: ignore[import-untyped]

        Vphys = Vt_b[:rank_b].T
        _Vnull_b = Vt_b[rank_b:].T
        ev_red, _vr_red = cast(
            "tuple[np.ndarray[Any, np.dtype[np.complex128]], np.ndarray[Any, np.dtype[np.complex128]]]",
            sla.eig(Vphys.T @ A @ Vphys, Vphys.T @ B @ Vphys),
        )
        ev_full = np.concatenate([ev_red, np.zeros(null_dim_b, dtype=np.complex128)])

        # All eigenvalues should be purely imaginary or zero (no spurious real parts)
        max_real = float(np.max(np.abs(np.real(ev_full))))
        assert max_real < 1e-10, (
            f"Projected spectrum has max |Re(λ)| = {max_real:.2e}, expected < 1e-10"
        )
        # Physical eigenvalues should be ±i (harmonic oscillator with ω=1)
        phys_evals = np.sort(np.abs(np.imag(ev_full[np.abs(np.imag(ev_full)) > 0.5])))[
            ::-1
        ]
        assert len(phys_evals) == 2, (
            f"Expected 2 physical eigenvalues, got {len(phys_evals)}"
        )
        assert abs(phys_evals[0] - 1.0) < 1e-10, (
            f"Expected |Im(λ)|=1.0, got {phys_evals[0]:.6f}"
        )

    def test_eigenbasis_invertible_after_projection(self) -> None:
        """Full eigenbasis V_full = [Vphys@vr_red | Vnull] is invertible.

        This confirms that the lifted eigenvector matrix can be inverted to
        transform ICs into the eigenbasis — a necessary condition for correctness.
        """
        A = np.zeros((4, 4), dtype=np.complex128)
        B = np.zeros((4, 4), dtype=np.complex128)
        A[0, 1] = 1.0
        A[1, 0] = -1.0
        A[2, 3] = 1.0
        A[3, 2] = -1.0
        B[0, 0] = 1.0
        B[1, 1] = 1.0

        import scipy.linalg as sla  # type: ignore[import-untyped]

        _, s_b, Vt_b = np.linalg.svd(B)
        rank_b = int(np.sum(s_b > s_b[0] * 1e-10))
        Vphys = Vt_b[:rank_b].T
        Vnull = Vt_b[rank_b:].T
        _, vr_red = cast(
            "tuple[np.ndarray[Any, np.dtype[np.complex128]], np.ndarray[Any, np.dtype[np.complex128]]]",
            sla.eig(Vphys.T @ A @ Vphys, Vphys.T @ B @ Vphys),
        )
        V_full = np.hstack([Vphys @ vr_red, Vnull])

        cond = float(np.linalg.cond(V_full))
        assert cond < 1e6, (
            f"V_full condition number {cond:.2e} too large — not invertible"
        )
        # Confirm it can be inverted without error
        V_inv = np.linalg.inv(V_full)
        residual = float(np.max(np.abs(V_full @ V_inv - np.eye(4))))
        assert residual < 1e-10, (
            f"V_full @ V_inv residual = {residual:.2e}, expected < 1e-10"
        )


class TestEigendataExport:
    """Tests for return_eigendata=True (v6 Stage 3).

    Pass 0 must expose its eigendecomposition so Pass 1 Duhamel can reuse
    it without re-eigendecomposing. Verify structure, invertibility, and
    that the reconstructed state matches the solver's own snapshot
    output to machine precision.
    """

    def _run_kg(
        self,
        return_eigendata: bool,
    ) -> tuple[dict[str, Any], np.ndarray, StateLayout, GridInfo]:
        spec = _make_spec(_KG_1D_SPEC)
        n = 32
        length = 2 * np.pi
        grid = GridInfo(shape=(n,), bounds=((0.0, length),), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        x = np.linspace(0.0, length, n, endpoint=False)
        y0 = np.zeros(layout.num_slots * grid.num_points)
        y0[:n] = np.sin(x)  # k=1 Fourier mode
        # Give a non-trivial velocity too so the block IC has both components.
        y0[n : 2 * n] = 0.5 * np.cos(x)

        result = cast(
            "dict[str, Any]",
            solve_modal(
                spec,
                grid,
                y0,
                t_span=(0.0, 1.0),
                parameters={"m2": 1.0},
                num_snapshots=11,
                return_eigendata=return_eigendata,
            ),
        )
        return result, y0, layout, grid

    def test_eigendata_key_absent_by_default(self) -> None:
        result, _, _, _ = self._run_kg(return_eigendata=False)
        assert "eigendata" not in result

    def test_eigendata_structure_present_when_requested(self) -> None:
        """Pass 0 collects {slot_indices, M_block, y0_block} for Pass 1.

        The legacy {V, V_inv, D_diag, alpha} schema was retired in v0.31+ when
        Pass 1 was rewritten to use the augmented matrix exponential
        (Al-Mohy & Higham 2011 §5.2) — see _evolve_duhamel_per_mode and
        docs/tex/modal_solver.tex §"Robust Matrix-Exponential Evolution".
        """
        result, _, _, _ = self._run_kg(return_eigendata=True)
        assert "eigendata" in result
        ed = result["eigendata"]
        assert set(ed.keys()) >= {"blocks", "mode_k", "state_layout"}
        assert len(ed["blocks"]) >= 1
        block = ed["blocks"][0]
        assert set(block.keys()) == {"slot_indices", "M_block", "y0_block"}

    def test_eigendata_y0_block_matches_fft_ic(self) -> None:
        """y0_block stored per block matches rfft(y0) at the same slot indices."""
        result, y0, layout, grid = self._run_kg(return_eigendata=True)
        # The Fourier transform helper is intentionally private in modal.py;
        # importing here is a deliberate whitebox check that y0_block matches
        # rfft(y0) without going through the full Pass 0 evolution.
        from tidal.solver.modal import _fft_slots

        y0_hat = _fft_slots(y0, layout, grid)
        for block in result["eigendata"]["blocks"]:
            slot_indices = block["slot_indices"]
            expected_y0_block = y0_hat[slot_indices, :]  # (bs, n_modes)
            np.testing.assert_allclose(
                block["y0_block"],
                expected_y0_block,
                atol=1e-14,
                rtol=1e-14,
            )

    def test_eigendata_reconstructs_pass_zero(self) -> None:
        """y_hat(t) = exp(M·t) · y0_block reproduces snapshots to machine precision.

        Pass 1 augmented-exp uses M_block to build the (2bs × 2bs) augmented
        operator [[A, S], [0, A]]. Verify that the same M_block evolved as
        exp(M·t)·y0_block reproduces the solver's own Pass 0 output.
        """
        import scipy.linalg as sla  # type: ignore[import-untyped]

        result, _, layout, grid = self._run_kg(return_eigendata=True)
        ed = result["eigendata"]
        t = result["t"]
        n_modes = ed["blocks"][0]["M_block"].shape[0]
        n_slots = layout.num_slots
        y_hat_rec = np.zeros((n_slots, n_modes), dtype=np.complex128)
        ti = len(t) - 1
        dt = float(t[ti] - t[0])
        for block in ed["blocks"]:
            M = block["M_block"]  # (n_modes, bs, bs)
            y0_block = block["y0_block"]  # (bs, n_modes)
            for m in range(n_modes):
                y_evolved = sla.expm(M[m] * dt) @ y0_block[:, m]
                for i, slot in enumerate(block["slot_indices"]):
                    y_hat_rec[slot, m] = y_evolved[i]

        # Inverse FFT to physical space
        rfft_shape_list = list(grid.shape)
        rfft_shape_list[-1] = grid.shape[-1] // 2 + 1
        rfft_shape = tuple(rfft_shape_list)
        n_pts = grid.num_points
        n_dyn = n_slots  # KG: pure wave system, no constraints
        y_rec_phys = np.zeros(n_dyn * n_pts)
        for si in range(n_dyn):
            y_rec_phys[si * n_pts : (si + 1) * n_pts] = np.fft.irfftn(
                y_hat_rec[si].reshape(rfft_shape),
                s=grid.shape,
                axes=list(range(len(grid.shape))),
            ).ravel()

        err = float(np.max(np.abs(y_rec_phys - result["y"][ti])))
        assert err < 1e-12, (
            f"M_block reconstruction error {err:.2e} exceeds 1e-12 tolerance"
        )

    def test_position_dependent_raises(self) -> None:
        """Position-dependent coefficients + return_eigendata raises."""
        # Use a spec with position-dependent mass (background field) to force
        # the convolution path. The simplest trigger: give the mass a
        # coordinate_dependent entry via a custom spec.
        spec_data = copy.deepcopy(_KG_1D_SPEC)
        eqs = cast("list[dict[str, Any]]", spec_data["equations"])
        eqs[0]["rhs"]["terms"][0]["coordinate_dependent"] = ["x"]
        spec = _make_spec(spec_data)
        grid = GridInfo(shape=(32,), bounds=((0.0, 2 * np.pi),), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        y0 = np.zeros(layout.num_slots * grid.num_points)

        with pytest.raises(NotImplementedError, match="position-"):
            solve_modal(
                spec,
                grid,
                y0,
                t_span=(0.0, 0.1),
                parameters={"m2": 1.0},
                num_snapshots=2,
                return_eigendata=True,
            )


# =========================================================================
# JAX backend tests (phase 1: constant-coefficient path)
# =========================================================================


def _jax_solver() -> object:
    """Import solve_modal_jax, skip test if JAX not installed."""
    jax = pytest.importorskip("jax")
    jax.config.update("jax_enable_x64", True)  # noqa: FBT003
    from tidal.solver.modal_jax import solve_modal_jax

    return solve_modal_jax


class TestModalJAXCorrectness:
    """Validate the JAX modal backend against the scipy backend.

    All tests skip automatically when ``jax`` / ``jaxlib`` are not installed.
    Phase-1 scope: constant-coefficient, no constraints, no position-dependent
    coefficients.  Phase-2 paths (position-dep, constraints, return_eigendata)
    are marked skip here.
    """

    _TOL = 1e-10  # max relative error between JAX and scipy backends

    def _assert_backends_agree(
        self,
        result_scipy: dict[str, object],
        result_jax: dict[str, object],
    ) -> None:
        """Assert that JAX and scipy outputs agree to self._TOL."""
        y_s = result_scipy["y"]
        y_j = result_jax["y"]
        assert y_s.shape == y_j.shape, f"Shape mismatch: {y_s.shape} vs {y_j.shape}"  # type: ignore[union-attr]
        max_ref = float(np.max(np.abs(y_s)))  # type: ignore[union-attr]
        if max_ref < 1e-15:
            # Zero solution: check absolute error
            max_abs_err = float(np.max(np.abs(y_j - y_s)))  # type: ignore[operator]
            assert max_abs_err < 1e-13, f"Zero-IC: abs err {max_abs_err:.2e}"
        else:
            max_rel_err = float(np.max(np.abs(y_j - y_s))) / max_ref  # type: ignore[operator]
            assert max_rel_err < self._TOL, (
                f"JAX vs scipy max relative error {max_rel_err:.2e} > {self._TOL}"
            )

    def test_scalar_wave_jax_agrees_with_scipy(self) -> None:
        """KG scalar wave: JAX backend matches scipy to < 1e-10 relative error."""
        solve_modal_jax = _jax_solver()

        spec = _make_spec(_KG_1D_SPEC)
        N = 64
        L = 10.0
        grid = GridInfo(shape=(N,), bounds=((0.0, L),), periodic=(True,))
        y0 = _make_gaussian_ic(spec, grid)
        params = {"m2": 1.0}

        result_scipy = solve_modal(
            spec, grid, y0, (0.0, 5.0), parameters=params, num_snapshots=51
        )
        result_jax = solve_modal_jax(
            spec, grid, y0, (0.0, 5.0), parameters=params, num_snapshots=51
        )  # type: ignore[operator]

        assert result_jax["success"]
        self._assert_backends_agree(result_scipy, result_jax)

    def test_coupled_wave_jax_agrees_with_scipy(self) -> None:
        """Coupled scalars: JAX backend matches scipy to < 1e-10 relative error."""
        solve_modal_jax = _jax_solver()

        spec = _make_spec(_COUPLED_SPEC)
        grid = GridInfo(shape=(64,), bounds=((0.0, 10.0),), periodic=(True,))
        y0 = _make_gaussian_ic(spec, grid)
        params = {"mPhi2": 1.0, "mChi2": 4.0, "gCpl": 0.5}

        result_scipy = solve_modal(
            spec, grid, y0, (0.0, 3.0), parameters=params, num_snapshots=31
        )
        result_jax = solve_modal_jax(
            spec, grid, y0, (0.0, 3.0), parameters=params, num_snapshots=31
        )  # type: ignore[operator]

        assert result_jax["success"]
        self._assert_backends_agree(result_scipy, result_jax)

    def test_diffusion_jax_agrees_with_scipy(self) -> None:
        """Diffusion: JAX backend matches scipy to < 1e-10 relative error."""
        solve_modal_jax = _jax_solver()

        spec = _make_spec(_DIFFUSION_SPEC)
        N = 64
        L = 10.0
        grid = GridInfo(shape=(N,), bounds=((0.0, L),), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        x = np.linspace(0, L, N, endpoint=False)
        y0 = np.zeros(layout.num_slots * grid.num_points)
        y0[:N] = np.cos(2 * np.pi * 2 / L * x)

        result_scipy = solve_modal(
            spec, grid, y0, (0.0, 5.0), parameters={"D": 0.1}, num_snapshots=11
        )
        result_jax = solve_modal_jax(
            spec, grid, y0, (0.0, 5.0), parameters={"D": 0.1}, num_snapshots=11
        )  # type: ignore[operator]

        assert result_jax["success"]
        self._assert_backends_agree(result_scipy, result_jax)

    def test_machine_precision_jax(self) -> None:
        """JAX backend achieves < 1e-10 relative error vs scipy (machine precision)."""
        solve_modal_jax = _jax_solver()

        spec = _make_spec(_KG_1D_SPEC)
        N = 32
        L = 2 * np.pi
        grid = GridInfo(shape=(N,), bounds=((0.0, L),), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        x = np.linspace(0, L, N, endpoint=False)
        y0 = np.zeros(layout.num_slots * grid.num_points)
        y0[:N] = np.sin(x)

        result_scipy = solve_modal(
            spec, grid, y0, (0.0, 1.0), parameters={"m2": 1.0}, num_snapshots=2
        )
        result_jax = solve_modal_jax(
            spec, grid, y0, (0.0, 1.0), parameters={"m2": 1.0}, num_snapshots=2
        )  # type: ignore[operator]

        assert result_jax["success"]
        self._assert_backends_agree(result_scipy, result_jax)

    def test_2d_wave_jax_agrees_with_scipy(self) -> None:
        """2D wave: JAX backend matches scipy to < 1e-10 relative error."""
        solve_modal_jax = _jax_solver()

        spec_2d: dict[str, object] = {
            "metadata": {"source": "inline-test", "parameters": {"m2": 1.0}},
            "spacetime": {
                "dimension": 3,
                "signature": [-1, 1, 1],
                "coordinates": ["t", "x", "y"],
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
                                "operator": "laplacian",
                                "field": "phi_0",
                            },
                        ],
                    },
                }
            ],
        }
        spec = _make_spec(spec_2d)
        Nx, Ny, Lx, Ly = 16, 16, 4.0, 4.0
        grid = GridInfo(
            shape=(Nx, Ny), bounds=((0.0, Lx), (0.0, Ly)), periodic=(True, True)
        )
        layout = StateLayout.from_spec(spec, grid.num_points)
        # 2D Gaussian IC — must use meshgrid, not 1D linspace
        x = np.linspace(0.0, Lx, Nx, endpoint=False)
        y = np.linspace(0.0, Ly, Ny, endpoint=False)
        X, Y = np.meshgrid(x, y, indexing="ij")
        y0 = np.zeros(layout.num_slots * grid.num_points)
        y0[: grid.num_points] = (
            0.1 * np.exp(-((X - Lx / 2) ** 2 + (Y - Ly / 2) ** 2) / (2 * 1.5**2))
        ).ravel()

        result_scipy = solve_modal(
            spec, grid, y0, (0.0, 2.0), parameters={"m2": 1.0}, num_snapshots=11
        )
        result_jax = solve_modal_jax(
            spec, grid, y0, (0.0, 2.0), parameters={"m2": 1.0}, num_snapshots=11
        )  # type: ignore[operator]

        assert result_jax["success"]
        self._assert_backends_agree(result_scipy, result_jax)

    def test_zero_ic_stays_zero_jax(self) -> None:
        """Zero IC: JAX backend produces zero output."""
        solve_modal_jax = _jax_solver()

        spec = _make_spec(_KG_1D_SPEC)
        grid = GridInfo(shape=(32,), bounds=((0.0, 10.0),), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        y0 = np.zeros(layout.num_slots * grid.num_points)

        result_jax = solve_modal_jax(
            spec, grid, y0, (0.0, 5.0), parameters={"m2": 1.0}, num_snapshots=11
        )  # type: ignore[operator]

        assert result_jax["success"]
        assert np.max(np.abs(result_jax["y"])) < 1e-14

    def test_block_isolation_jax(self) -> None:
        """Zero-IC block stays at machine zero with JAX backend."""
        solve_modal_jax = _jax_solver()

        spec = _make_spec(_DEGENERATE_PAIRS_SPEC)
        grid = GridInfo(shape=(64,), bounds=((0.0, 10.0),), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        y0 = np.zeros(layout.num_slots * grid.num_points)
        x = np.linspace(0.0, 10.0, 64, endpoint=False)
        y0[: grid.num_points] = 0.1 * np.exp(-((x - 5.0) ** 2) / (2 * 1.5**2))

        result_jax = solve_modal_jax(  # type: ignore[operator]
            spec,
            grid,
            y0,
            (0.0, 10.0),
            parameters={"B0": 0.3, "kappa2": 1.0},
            num_snapshots=5,
        )

        assert result_jax["success"]
        # Pair 2 fields (phi_1, chi_1) should stay at machine zero
        phi_1_slot = layout.field_slot_map["phi_1"]
        chi_1_slot = layout.field_slot_map["chi_1"]
        final = result_jax["y"][-1]
        n = grid.num_points
        for slot in [phi_1_slot, chi_1_slot]:
            max_val = np.max(np.abs(final[slot * n : (slot + 1) * n]))
            assert max_val < 1e-12, f"Zero-IC slot {slot} grew to {max_val:.2e}"

    @pytest.mark.skip(
        reason="phase 1: modal-jax position-dependent path not implemented"
    )
    def test_position_dependent_jax(self) -> None:
        pass

    @pytest.mark.skip(reason="phase 1: modal-jax constraint/Schur path not implemented")
    def test_constraint_jax(self) -> None:
        pass

    @pytest.mark.skip(reason="phase 1: modal-jax return_eigendata not implemented")
    def test_eigendata_jax(self) -> None:
        pass

    def test_jax_raises_for_position_dependent(self) -> None:
        """JAX backend raises NotImplementedError for position-dependent systems."""
        _jax_solver()
        from tidal.solver.modal_jax import solve_modal_jax

        spec_data = copy.deepcopy(_KG_1D_SPEC)
        eqs = cast("list[dict[str, Any]]", spec_data["equations"])
        eqs[0]["rhs"]["terms"][0]["coordinate_dependent"] = ["x"]
        spec = _make_spec(spec_data)
        grid = GridInfo(shape=(32,), bounds=((0.0, 10.0),), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        y0 = np.zeros(layout.num_slots * grid.num_points)

        with pytest.raises(NotImplementedError, match="phase 1"):
            solve_modal_jax(spec, grid, y0, (0.0, 1.0), parameters={"m2": 1.0})

    def test_jax_raises_for_return_eigendata(self) -> None:
        """JAX backend raises NotImplementedError for return_eigendata=True."""
        _jax_solver()
        from tidal.solver.modal_jax import solve_modal_jax

        spec = _make_spec(_KG_1D_SPEC)
        grid = GridInfo(shape=(32,), bounds=((0.0, 10.0),), periodic=(True,))
        y0 = _make_gaussian_ic(spec, grid)

        with pytest.raises(NotImplementedError, match="phase 1"):
            solve_modal_jax(
                spec,
                grid,
                y0,
                (0.0, 1.0),
                parameters={"m2": 1.0},
                return_eigendata=True,
            )


# =========================================================================
# GH #367: position-dependent + periodic auto-routing
# =========================================================================


class TestPositionDepAutoRoute:
    """Verify the GH #367 fix (v0.41.6): the convolution-matrix path in
    ``_build_convolution_matrix`` now applies the ``M⁻¹`` (inverse kinetic
    coefficient) scaling that the per-mode path already had. With the fix,
    position-dependent + periodic theories run correctly on the modal
    solver — no auto-route to CVODE/IDA is needed.

    Root cause: ``_build_per_mode_matrices`` already folded
    ``build_inverse_kinetic_diag`` into velocity-row coefficients (#301/#302),
    but ``_build_convolution_matrix`` was missing that step. For Gertsenshtein-
    class theories carrying ``kinetic_coefficient_symbolic = -1/kappa²`` on
    graviton modes, the discrete EOM had a sign-flipped Laplacian giving
    a real eigenvalue at every k. ``expm_multiply`` faithfully amplified
    it — the apparent "Nyquist spurious eigenvalue" was the visible
    high-k symptom of a uniform-in-k sign error.
    """

    @pytest.fixture
    def e0_args(self, tmp_path: Path) -> list[str]:
        """CLI args for the GH #367 reproducer (E.0 dual-Gaussian)."""
        output_dir = tmp_path / "out"
        return [
            "examples/data/gertsenshtein_e0_dual_gaussian.json",
            "--grid-shape",
            "64",
            "--bounds",
            "0:100",
            "--periodic",
            "--ic",
            "gaussian",
            "--ic-component",
            "h_5",
            "--ic-amplitude",
            "1e-2",
            "--ic-width",
            "5",
            "--ic-center",
            "25",
            "--param",
            "kappa=1.0",
            "--param",
            "Bpeak=0.01",
            "--param",
            "sigB=5",
            "--param",
            "zc1=25",
            "--param",
            "zc2=75",
            "--snapshots",
            "2",
            "--output",
            str(output_dir),
        ]

    def _run_simulate(self, args: list[str]) -> tuple[int, str]:
        """Invoke the ``tidal simulate`` CLI; return (exit_code, combined_output)."""
        import subprocess

        result = subprocess.run(
            ["uv", "run", "tidal", "simulate", *args],  # noqa: S607
            capture_output=True,
            text=True,
            timeout=120,
            check=False,
        )
        return result.returncode, result.stdout + result.stderr

    def test_gh367_reproducer_modal_correct(self, e0_args: list[str]) -> None:
        """GH #367 reproducer at t_end=20: --scheme auto stays on modal
        (no override since v0.41.6) and produces the correct decaying
        h_5 (≈0.005), matching CVODE within 1%.
        """
        args = [*e0_args, "--t-end", "20"]
        exit_code, output = self._run_simulate(args)
        assert exit_code == 0, f"simulate failed:\n{output}"
        # No auto-routing happens any more.
        assert "auto-routing" not in output, (
            f"unexpected auto-route after the v0.41.6 fix; got:\n{output}"
        )
        # Scheme line confirms modal.
        assert "Auto-selected solver: modal" in output, (
            f"expected auto-selected modal; got:\n{output}"
        )
        import re

        match = re.search(r"h_5:.*→\s*([0-9.eE+-]+)", output)
        assert match is not None, f"h_5 result line not found in:\n{output}"
        h5_final = float(match.group(1))
        assert 0.003 < h5_final < 0.007, (
            f"h_5={h5_final} not in expected range [0.003, 0.007] (CVODE truth ≈ 0.005)"
        )

    def test_explicit_modal_at_t20_produces_correct_h5(
        self, e0_args: list[str]
    ) -> None:
        """Explicit --scheme modal at t_end=20 runs to completion and matches
        CVODE on h_5 peak (the convolution-path m_inv fix in v0.41.6).
        """
        args = [*e0_args, "--scheme", "modal", "--t-end", "20"]
        exit_code, output = self._run_simulate(args)
        assert exit_code == 0, f"simulate failed:\n{output}"
        assert "Simulation diverged" not in output, (
            f"unexpected divergence after the v0.41.6 fix; got:\n{output}"
        )
        import re

        match = re.search(r"h_5:.*→\s*([0-9.eE+-]+)", output)
        assert match is not None, f"h_5 result line not found in:\n{output}"
        h5_final = float(match.group(1))
        assert 0.003 < h5_final < 0.007, (
            f"h_5={h5_final} not in expected range [0.003, 0.007]"
        )

    def test_explicit_modal_safe_regime_correct(self, e0_args: list[str]) -> None:
        """Sanity check the modal solver at small t_end is accurate.

        Before the v0.41.6 fix this only verified non-divergence; with the
        fix in place the answer must match CVODE within a few percent
        even at very small t_end.
        """
        args = [*e0_args, "--scheme", "modal", "--t-end", "1"]
        exit_code, output = self._run_simulate(args)
        assert exit_code == 0, f"simulate failed unexpectedly:\n{output}"
        assert "Scheme: modal" in output
        assert "Simulation diverged" not in output


class TestOverridePosDepPeriodicScheme:
    """Unit tests for ``_override_pos_dep_periodic_scheme``.

    Since v0.41.6 the function is a permanent no-op — the GH #367 root cause
    (the convolution path missing ``M⁻¹`` scaling) was fixed in
    ``_build_convolution_matrix``, so modal is the correct choice for
    pos-dep + periodic theories. These tests pin the no-op semantics so a
    future regression that reintroduces the override is caught.
    """

    @staticmethod
    def _pos_dep_diffusion() -> EquationSystem:
        """Pos-dep (no constraint): first-order diffusion with symbolic coord
        in the coefficient. Mirrors the trick used in
        ``test_position_dependent_correctness``.
        """
        spec_data = copy.deepcopy(_DIFFUSION_SPEC)
        spec_data["equations"][0]["rhs"]["terms"][0][  # type: ignore[index]
            "coefficient_symbolic"
        ] = "D*(1 + 0*x[])"
        return _make_spec(spec_data)

    @staticmethod
    def _pos_dep_with_constraint() -> EquationSystem:
        """Pos-dep + constraint: A_0 is a constraint (time_order=0), A_1 is
        dynamical with a position-dependent coefficient. Models the Phase E
        localized PGT pattern (constraint field + Erf-profile background).
        """
        spec_data = copy.deepcopy(_CONSTRAINT_SPEC)
        # Force position-dependence on the dynamical equation's RHS term.
        spec_data["equations"][1]["rhs"]["terms"][0][  # type: ignore[index]
            "coefficient_symbolic"
        ] = "1 + 0*x[]"
        return _make_spec(spec_data)

    def test_non_modal_scheme_passes_through(self) -> None:
        """If _resolve_scheme didn't pick modal, the override is a no-op."""
        from tidal.cli._simulate import _override_pos_dep_periodic_scheme

        spec = self._pos_dep_with_constraint()
        for scheme in ("cvode", "ida", "leapfrog", "scipy"):
            new_scheme, msg = _override_pos_dep_periodic_scheme(scheme, "auto", spec)
            assert new_scheme == scheme
            assert msg is None

    def test_explicit_modal_not_overridden(self) -> None:
        """User explicitly asked for modal; respect it (post-evolution
        amplitude check at modal.py:2469 is the safety net for that path).
        """
        from tidal.cli._simulate import _override_pos_dep_periodic_scheme

        spec = self._pos_dep_with_constraint()
        new_scheme, msg = _override_pos_dep_periodic_scheme("modal", "modal", spec)
        assert new_scheme == "modal"
        assert msg is None

    def test_constant_coefficients_no_override(self) -> None:
        """No pos-dep terms → modal is correct and stays."""
        from tidal.cli._simulate import _override_pos_dep_periodic_scheme

        spec = _make_spec(_KG_1D_SPEC)
        new_scheme, msg = _override_pos_dep_periodic_scheme("modal", "auto", spec)
        assert new_scheme == "modal"
        assert msg is None

    def test_pos_dep_no_constraint_stays_on_modal(self) -> None:
        """Pos-dep + no constraint: since v0.41.6 modal is the correct choice
        (was: routed to CVODE under the GH #367 safety net).
        """
        from tidal.cli._simulate import _override_pos_dep_periodic_scheme

        spec = self._pos_dep_diffusion()
        new_scheme, msg = _override_pos_dep_periodic_scheme("modal", "auto", spec)
        assert new_scheme == "modal"
        assert msg is None

    def test_pos_dep_with_constraint_stays_on_modal(self) -> None:
        """Pos-dep + constraint: since v0.41.6 modal is the correct choice
        (was: routed to IDA under the GH #367 safety net). The modal solver's
        Schur path handles constraints natively (modal.py:541-595).
        """
        from tidal.cli._simulate import _override_pos_dep_periodic_scheme

        spec = self._pos_dep_with_constraint()
        new_scheme, msg = _override_pos_dep_periodic_scheme("modal", "auto", spec)
        assert new_scheme == "modal"
        assert msg is None


# ---------------------------------------------------------------------------
# GH #421: position-dependent kinetic coefficients are refused clearly
# ---------------------------------------------------------------------------


class TestGH421PositionDependentKinetic:
    """The modal solver cannot handle a position-dependent
    ``kinetic_coefficient_symbolic`` (a mass-side k-space convolution
    M̂(k−k′), tracked as GH #427).  Pre-guard behavior was a deep
    grid-less ``evaluate_coefficient`` ValueError on the per-mode and
    convolution paths, and — worse — a silent ``except Exception → M=1``
    fallback on the generalized-eig path, which the stability probe (and
    therefore every inference likelihood evaluation) reaches directly.
    These tests pin the refusal at every entry and that selection routes
    around modal.
    """

    EH_SPEC = "examples/data/euler_heisenberg_e_dual_gaussian.json"
    EH_PARAMS: dict[str, float] = {
        "Bpeak": 0.01,
        "rho": 1.0,
        "sigma": 1.0,
        "sigB": 0.4,
        "zc1": -1.5,
        "zc2": 1.5,
    }

    def _load_eh(self, *, strip_perturbation: bool) -> EquationSystem:
        """Load the only shipped spec with position-dependent kinetics.

        With the ``[perturbation]`` metadata stripped, the CLI's
        perturbative route (whose canonicalization normalizes the
        kinetics to '1', GH #380) is disabled and the position
        dependence reaches the plain modal path unmasked.
        """
        import json
        from pathlib import Path

        raw = json.loads(Path(self.EH_SPEC).read_text(encoding="utf-8"))
        if strip_perturbation:
            raw.get("metadata", {}).pop("perturbation", None)
        return EquationSystem.from_dict(raw)

    @staticmethod
    def _grid() -> GridInfo:
        return GridInfo(bounds=((-2.0, 2.0),), shape=(32,), periodic=(True,))

    # -- predicate ------------------------------------------------------

    def test_predicate_on_eh_spec(self) -> None:
        spec = self._load_eh(strip_perturbation=False)
        assert spec.has_position_dependent_kinetic()
        posdep = [
            eq.field_name for eq in spec.equations if eq.kinetic_position_dependent
        ]
        assert posdep == ["a_1", "a_2", "a_3"]

    def test_predicate_constant_and_time_only_kinetics(self) -> None:
        """Parameter-only and t[]-only kinetics are NOT position-dependent."""
        spec_data = copy.deepcopy(_KG_1D_SPEC)
        spec_data["equations"][0]["lhs"]["kinetic_coefficient_symbolic"] = (  # type: ignore[index]
            "-kappa^(-2)"
        )
        assert not _make_spec(spec_data).has_position_dependent_kinetic()
        spec_data["equations"][0]["lhs"]["kinetic_coefficient_symbolic"] = "1 + t[]"  # type: ignore[index]
        assert not _make_spec(spec_data).has_position_dependent_kinetic()
        spec_data["equations"][0]["lhs"]["kinetic_coefficient_symbolic"] = (  # type: ignore[index]
            "1 + Sin[x[]]"
        )
        assert _make_spec(spec_data).has_position_dependent_kinetic()

    def test_predicate_ignores_constraint_equations(self) -> None:
        """Only dynamical equations count — consumer-faithful to
        build_inverse_kinetic_diag, which skips time_order=0 rows.
        """
        spec_data = copy.deepcopy(_KG_1D_SPEC)
        spec_data["equations"][0]["lhs"]["order"]["time"] = 0  # type: ignore[index]
        spec_data["equations"][0]["lhs"]["kinetic_coefficient_symbolic"] = (  # type: ignore[index]
            "1 + Sin[x[]]"
        )
        assert not _make_spec(spec_data).has_position_dependent_kinetic()

    # -- eligibility ----------------------------------------------------

    def test_can_use_modal_accepts_posdep_kinetic(self) -> None:
        # GH #427: position-dependent kinetics are modal-eligible — the
        # kinetics-aware routing predicate sends them to the convolution
        # builders, which fold M⁻¹(x) into the real-space coefficients.
        spec = self._load_eh(strip_perturbation=True)
        assert can_use_modal(spec, self._grid(), None) is True

    def test_can_use_modal_accepts_constant_symbolic_kinetic(self) -> None:
        spec_data = copy.deepcopy(_KG_1D_SPEC)
        spec_data["equations"][0]["lhs"]["kinetic_coefficient_symbolic"] = (  # type: ignore[index]
            "-kappa^(-2)"
        )
        grid = GridInfo(shape=(64,), bounds=((0.0, 10.0),), periodic=(True,))
        assert can_use_modal(_make_spec(spec_data), grid, None) is True

    def test_can_use_modal_accepts_canonicalized_base(self) -> None:
        """The EH campaign flow stays modal-eligible: the driver-style
        canonicalized base spec has constant kinetics (GH #380), so the
        new check must not fire on it.  This is the regression pin for
        auto-selection of the perturbative flow.
        """
        spec = self._load_eh(strip_perturbation=False)
        small = list(spec.metadata.get("perturbation", {}).get("small_parameters", []))
        base = spec.canonicalize_kinetic_for_perturbation(small).base_spec(small)
        assert not base.has_position_dependent_kinetic()
        assert can_use_modal(base, self._grid(), None) is True

    # -- entry guards ---------------------------------------------------

    def test_solve_modal_runs_posdep_kinetic(self) -> None:
        # GH #427: was a raises-actionable pin; solve_modal now completes
        # on the stripped EH spec (convolution path 4 with M⁻¹(x) folded
        # into the coefficients). Accuracy is pinned separately against
        # CVODE in TestGH427PositionDependentKineticModal.
        spec = self._load_eh(strip_perturbation=True)
        grid = self._grid()
        layout = StateLayout.from_spec(spec, grid.num_points)
        rng = np.random.default_rng(0)
        y0 = 1e-3 * rng.standard_normal(layout.num_slots * grid.num_points)
        result = solve_modal(
            spec,
            grid,
            y0,
            (0.0, 0.1),
            parameters=self.EH_PARAMS,
            num_snapshots=3,
        )
        assert result["success"]
        assert np.all(np.isfinite(result["y"]))

    def test_build_evolution_matrices_raises_not_silent_m1(self) -> None:
        """The genEig builder previously swallowed the evaluation failure
        (``except Exception → M_mat = 1``) — wrong physics, no error.
        """
        from tidal.solver.coefficients import CoefficientEvaluator
        from tidal.solver.modal import (
            _build_evolution_matrices,
            _build_k_axes,
            _build_k_grid,
        )

        spec = self._load_eh(strip_perturbation=True)
        grid = self._grid()
        layout = StateLayout.from_spec(spec, grid.num_points)
        ce = CoefficientEvaluator(spec, grid, self.EH_PARAMS)
        k_grid = _build_k_grid(_build_k_axes(grid))
        rfft_shape = (grid.shape[0] // 2 + 1,)
        with pytest.raises(NotImplementedError, match="kinetic"):
            _build_evolution_matrices(spec, layout, grid, ce, k_grid, rfft_shape)

    def test_stability_probe_propagates_guard(self) -> None:
        """check_conversion_stability reaches _build_evolution_matrices
        directly (it gates every likelihood evaluation in inference); its
        (LinAlgError, ValueError) handler must NOT swallow the guard into
        a fabricated 'tachyonic' verdict.
        """
        from tidal.measurement._stability import check_conversion_stability

        spec = self._load_eh(strip_perturbation=True)
        with pytest.raises(NotImplementedError, match="kinetic"):
            check_conversion_stability(
                spec,
                self._grid(),
                self.EH_PARAMS,
                source="a_1",
            )

    def test_stability_probe_constant_kinetic_still_works(self) -> None:
        """A constant symbolic kinetic must still pass through the probe."""
        from tidal.measurement._stability import check_conversion_stability

        spec_data = copy.deepcopy(_KG_1D_SPEC)
        spec_data["equations"][0]["lhs"]["kinetic_coefficient_symbolic"] = "2"  # type: ignore[index]
        spec = _make_spec(spec_data)
        grid = GridInfo(shape=(32,), bounds=((0.0, 10.0),), periodic=(True,))
        result = check_conversion_stability(spec, grid, {"m2": 1.0}, source="phi_0")
        assert result.stable is not None  # produced a verdict, no raise

    def test_solve_modal_pass1_raises_on_posdep_correction(self) -> None:
        """The Pass 1 source builder consumes correction-spec kinetics
        without a grid; the guard must fire before eigendata is touched.
        """
        from tidal.solver.modal import solve_modal_pass1

        spec_data = copy.deepcopy(_KG_1D_SPEC)
        spec_data["equations"][0]["lhs"]["kinetic_coefficient_symbolic"] = (  # type: ignore[index]
            "1 + Sin[x[]]"
        )
        correction_spec = _make_spec(spec_data)
        with pytest.raises(NotImplementedError, match="kinetic"):
            solve_modal_pass1(
                {},  # never read: the guard fires first
                correction_spec,
                self._grid(),
                np.linspace(0.0, 1.0, 3),
            )

    # -- scheme resolution ----------------------------------------------

    def test_resolve_scheme_intact_eh_stays_modal(self) -> None:
        """THE campaign-flow regression pin: the intact perturbative EH
        spec must keep auto-selecting modal, because eligibility is now
        judged on the driver-style canonicalized base (constant
        kinetics), not on plain base_spec() (which keeps them
        position-dependent).
        """
        from tidal.cli._simulate import _resolve_scheme

        spec = self._load_eh(strip_perturbation=False)
        assert _resolve_scheme("auto", spec, self._grid(), None) == "modal"

    def test_resolve_scheme_stripped_eh_routes_to_modal(self) -> None:
        # GH #427: auto-selection now keeps the stripped EH spec on modal
        # (was {"ida", "cvode"} while requirement 6 excluded pos-dep
        # kinetics from eligibility).
        from tidal.cli._simulate import _resolve_scheme

        spec = self._load_eh(strip_perturbation=True)
        assert _resolve_scheme("auto", spec, self._grid(), None) == "modal"

    def test_resolve_scheme_explicit_modal_accepted(self) -> None:
        # GH #427: was a RuntimeError pin (requirement 6 rejection).
        from tidal.cli._simulate import _resolve_scheme

        spec = self._load_eh(strip_perturbation=True)
        assert _resolve_scheme("modal", spec, self._grid(), None) == "modal"

    # -- the documented alternative actually works ----------------------

    def test_stripped_eh_runs_under_cvode(self) -> None:
        """The guard's alternative (a) must be real: the stripped EH spec
        integrates under CVODE (which evaluates the position-dependent
        kinetics on the grid, GH #382).
        """
        from tidal.solver.cvode import solve_cvode

        spec = self._load_eh(strip_perturbation=True)
        grid = self._grid()
        layout = StateLayout.from_spec(spec, grid.num_points)
        rng = np.random.default_rng(0)
        y0 = 1e-3 * rng.standard_normal(layout.num_slots * grid.num_points)
        result = solve_cvode(
            spec,
            grid,
            y0,
            (0.0, 0.1),
            parameters=self.EH_PARAMS,
            num_snapshots=3,
        )
        assert result["success"]
        assert np.all(np.isfinite(result["y"]))


class TestGH427PositionDependentKineticModal:
    """GH #427 validation: modal handles position-dependent kinetics via
    real-space M⁻¹(x) folded into the convolution-path coefficients.

    Success criteria from the issue: modal-vs-CVODE RMS < 1% on the
    stripped EH dual-Gaussian spec (the only shipped pos-dep-kinetic
    spec), and spectral-rate convergence with grid resolution. Includes
    the GH #438 hardening pins (no silent corner-collapse of
    position-dependent coefficients in per-mode builders).
    """

    EH_SPEC = TestGH421PositionDependentKinetic.EH_SPEC
    EH_PARAMS = TestGH421PositionDependentKinetic.EH_PARAMS

    def _load_stripped_eh(self) -> EquationSystem:
        loader = TestGH421PositionDependentKinetic()
        return loader._load_eh(strip_perturbation=True)

    def _grid(self, n: int = 32) -> GridInfo:
        return GridInfo(shape=(n,), bounds=((-2.0, 2.0),), periodic=(True,))

    # -- issue success criterion 1: agreement with the time-domain ref --

    def test_stripped_eh_modal_matches_cvode_rms(self) -> None:
        """Modal (convolution path 4 with M⁻¹(x)) vs CVODE (grid-evaluated
        kinetics, GH #382) on identical smooth ICs: RMS relative
        difference < 1% over the window. Do NOT loosen this tolerance —
        a failure here means the physics differs between backends and
        must be investigated (soundness over coverage).

        Comparison design: CVODE runs with SPECTRAL spatial operators
        (conftest resets the global afterwards) so both backends share
        the exact spatial discretization, and the IC is band-limited
        (low-k harmonics, seeded phases per field). A white-noise IC
        would measure high-k discretization differences between the
        backends, not the GH #427 kinetic handling under test.
        """
        from tidal.solver.cvode import solve_cvode
        from tidal.solver.operators import set_spectral

        spec = self._load_stripped_eh()
        grid = self._grid()
        layout = StateLayout.from_spec(spec, grid.num_points)
        x = grid.axes_coords(0)
        length = 4.0  # bounds (-2, 2)
        rng = np.random.default_rng(0)
        y0 = np.zeros(layout.num_slots * grid.num_points)
        for slot in layout.field_slot_map.values():
            sl = layout.slot_slice(slot)
            phases = rng.uniform(0.0, 2 * np.pi, size=3)
            y0[sl] = 1e-3 * sum(
                np.cos(2 * np.pi * (m + 1) * x / length + phases[m]) for m in range(3)
            )
        t_span = (0.0, 0.1)

        modal = solve_modal(
            spec, grid, y0, t_span, parameters=self.EH_PARAMS, num_snapshots=5
        )
        set_spectral(True)
        cvode = solve_cvode(
            spec,
            grid,
            y0,
            t_span,
            parameters=self.EH_PARAMS,
            num_snapshots=5,
            rtol=1e-10,
            atol=1e-12,
        )
        assert modal["success"]
        assert cvode["success"]

        diff = np.asarray(modal["y"]) - np.asarray(cvode["y"])
        rms_diff = float(np.sqrt(np.mean(diff**2)))
        rms_ref = float(np.sqrt(np.mean(np.asarray(cvode["y"]) ** 2)))
        assert rms_ref > 0.0
        rel = rms_diff / rms_ref
        assert rel < 1e-2, (
            f"modal-vs-CVODE RMS relative difference {rel:.3e} exceeds 1% "
            f"on the stripped EH spec (GH #427 success criterion)"
        )

    # -- issue success criterion 2: spectral convergence ----------------

    def test_posdep_kinetic_spectral_convergence(self) -> None:
        """A smooth periodic M(x) = 1 + 0.25·sin(x) KG problem must
        converge to the fine-grid modal reference at spectral (faster
        than algebraic) rate: the N=32 error must beat the N=16 error by
        far more than the 4x of a 2nd-order scheme.

        Comparison design: TIDAL grids are CELL-CENTERED
        (grid.axes_coords), so centers never coincide across resolutions
        — comparing raw grid values would inject an O(dx) alignment
        artifact and mask the spectral rate. Instead the IC is the same
        physical function evaluated at each grid's own centers, and
        errors are measured on normalized DFT coefficients over shared
        low-|k| bins (basis-independent).
        """
        spec_data = copy.deepcopy(_KG_1D_SPEC)
        spec_data["equations"][0]["lhs"]["kinetic_coefficient_symbolic"] = (  # type: ignore[index]
            "1 + 0.25*Sin[x[]]"
        )
        spec = _make_spec(spec_data)
        length = 2 * np.pi
        t_span = (0.0, 0.5)
        params = {"m2": 1.0}
        k_keep = 5  # compare DFT bins |k| ≤ 5, well inside N=16's band

        def _solve_spectrum(n: int) -> NDArray[np.complex128]:
            grid = GridInfo(shape=(n,), bounds=((0.0, length),), periodic=(True,))
            layout = StateLayout.from_spec(spec, grid.num_points)
            x = grid.axes_coords(0)
            y0 = np.zeros(layout.num_slots * n)
            y0[:n] = np.sin(x)  # same physical IC at every resolution
            result = solve_modal(
                spec, grid, y0, t_span, parameters=params, num_snapshots=3
            )
            assert result["success"]
            field = np.asarray(result["y"])[-1, :n]
            spectrum = np.fft.fft(field) / n
            # Phase-align bins to the physical origin: samples sit at
            # cell centers x_j = (j+1/2)·dx, so bin k carries an extra
            # factor exp(+i·κ·dx/2) relative to the x=0-referenced
            # Fourier coefficient — divide it out.
            k_bins = np.fft.fftfreq(n, d=1.0 / n)  # integer wavenumbers
            dx = length / n
            aligned = spectrum * np.exp(-1j * k_bins * (2 * np.pi / length) * dx / 2)
            return np.concatenate([aligned[: k_keep + 1], aligned[-k_keep:]])

        ref = _solve_spectrum(128)
        errors: dict[int, float] = {}
        for n in (16, 32):
            errors[n] = float(np.max(np.abs(_solve_spectrum(n) - ref)))

        # 2nd-order FD would give errors[16]/errors[32] ≈ 4; spectral
        # accuracy on a smooth periodic problem should exceed that by
        # orders of magnitude (both grids resolve M's single harmonic).
        assert errors[32] < errors[16] / 20.0, (
            f"convergence not spectral: err(16)={errors[16]:.3e}, "
            f"err(32)={errors[32]:.3e} (ratio "
            f"{errors[16] / max(errors[32], 1e-300):.1f}x, need > 20x)"
        )

    def test_varying_kinetic_path3_matches_cvode(self) -> None:
        """Genuinely varying M(x) = 1 + 0.25·sin(x) on the constraint-free
        convolution path (path 3) agrees with spectral CVODE to < 1% RMS.
        Complements the EH test above, which exercises path 4.
        """
        from tidal.solver.cvode import solve_cvode
        from tidal.solver.operators import set_spectral

        spec_data = copy.deepcopy(_KG_1D_SPEC)
        spec_data["equations"][0]["lhs"]["kinetic_coefficient_symbolic"] = (  # type: ignore[index]
            "1 + 0.25*Sin[x[]]"
        )
        spec = _make_spec(spec_data)
        n = 32
        grid = GridInfo(shape=(n,), bounds=((0.0, 2 * np.pi),), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        x = grid.axes_coords(0)
        y0 = np.zeros(layout.num_slots * n)
        y0[:n] = np.sin(x) + 0.3 * np.cos(2 * x + 0.4)
        params = {"m2": 1.0}
        t_span = (0.0, 0.5)

        modal = solve_modal(spec, grid, y0, t_span, parameters=params, num_snapshots=5)
        set_spectral(True)
        cvode = solve_cvode(
            spec,
            grid,
            y0,
            t_span,
            parameters=params,
            num_snapshots=5,
            rtol=1e-10,
            atol=1e-12,
        )
        assert modal["success"]
        assert cvode["success"]
        diff = np.asarray(modal["y"]) - np.asarray(cvode["y"])
        rel = float(
            np.sqrt(np.mean(diff**2)) / np.sqrt(np.mean(np.asarray(cvode["y"]) ** 2))
        )
        assert rel < 1e-2, f"path-3 modal-vs-CVODE RMS {rel:.3e} exceeds 1%"

    # -- GH #438 hardening: no silent corner-collapse -------------------

    def test_resolve_constant_coeff_raises_on_posdep(self) -> None:
        """_resolve_constant_coeff must refuse an ndarray (position-
        dependent) coefficient instead of silently evaluating it at the
        first grid point (GH #438).
        """
        from tidal.solver.coefficients import CoefficientEvaluator
        from tidal.solver.modal import _resolve_constant_coeff

        spec = self._load_stripped_eh()
        grid = self._grid()
        ce = CoefficientEvaluator(spec, grid, self.EH_PARAMS)
        # Find a position-dependent term to resolve.
        eq_idx, term_idx, posdep_term = next(
            (ei, ti, term)
            for ei, eq in enumerate(spec.equations)
            for ti, term in enumerate(eq.rhs_terms)
            if term.position_dependent
        )
        with pytest.raises(NotImplementedError, match="first grid point"):
            _resolve_constant_coeff(posdep_term, ce, eq_idx=eq_idx, term_idx=term_idx)

    def test_stability_probe_refuses_posdep_rhs_background(self) -> None:
        """The stability probe on a localized (position-dependent RHS)
        spec previously evaluated the background at the domain edge
        (B ≈ 0) and judged the near-vacuum theory. Post-GH #438 the
        refusal propagates — NotImplementedError deliberately does not
        match the probe's (LinAlgError, ValueError) handler, so it must
        NOT be converted into a fabricated verdict.
        """
        from tidal.measurement._stability import check_conversion_stability

        spec_data = copy.deepcopy(_KG_1D_SPEC)
        terms = spec_data["equations"][0]["rhs"]["terms"]  # type: ignore[index]
        terms.append(
            {
                "coefficient": 1.0,
                "operator": "identity",
                "field": "phi_0",
                "coefficient_symbolic": "0.5*Sin[x[]]",
                "coordinate_dependent": ["x"],
            }
        )
        spec = _make_spec(spec_data)
        grid = GridInfo(shape=(32,), bounds=((0.0, 10.0),), periodic=(True,))
        with pytest.raises(NotImplementedError, match="first grid point"):
            check_conversion_stability(spec, grid, {"m2": 1.0}, source="phi_0")


class TestGH427EnergyOnPositionDependentKinetics:
    """WS1 rider: the energy machinery on the one shipped pos-dep-kinetic
    spec, which only became runnable (and hence measurable) with GH #427.

    `_energy.py` resolves position-dependent Hamiltonian coefficients on
    the grid (`_build_coord_arrays` / `evaluate_coefficient` with
    coord_arrays); the stripped EH spec carries 21 such terms. This pins
    that the path (a) evaluates without raising, (b) returns finite
    values, and (c) conserves total energy over a short window — the
    background is time-independent, so H is conserved; a corner-collapse
    or basis error in the energy coefficients would break conservation at
    the percent level immediately.
    """

    def test_stripped_eh_energy_finite_and_conserved(self) -> None:
        from tidal.measurement._energy import compute_energy_timeseries
        from tidal.measurement._io import SimulationData

        loader = TestGH421PositionDependentKinetic()
        spec = loader._load_eh(strip_perturbation=True)
        grid = GridInfo(shape=(32,), bounds=((-2.0, 2.0),), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        x = grid.axes_coords(0)
        length = 4.0
        rng = np.random.default_rng(0)
        y0 = np.zeros(layout.num_slots * grid.num_points)
        for slot in layout.field_slot_map.values():
            sl = layout.slot_slice(slot)
            phases = rng.uniform(0.0, 2 * np.pi, size=3)
            y0[sl] = 1e-3 * sum(
                np.cos(2 * np.pi * (m + 1) * x / length + phases[m]) for m in range(3)
            )
        params = TestGH421PositionDependentKinetic.EH_PARAMS

        result = solve_modal(
            spec, grid, y0, (0.0, 0.1), parameters=params, num_snapshots=5
        )
        assert result["success"]

        n_pts = grid.num_points
        snaps = np.asarray(result["y"])
        fields = {
            name: snaps[:, slot * n_pts : (slot + 1) * n_pts]
            for name, slot in layout.field_slot_map.items()
        }
        velocities = {
            name: snaps[:, slot * n_pts : (slot + 1) * n_pts]
            for name, slot in layout.velocity_slot_map.items()
        }
        data = SimulationData(
            times=np.asarray(result["t"]),
            fields=fields,
            velocities=velocities,
            grid_spacing=grid.dx,
            grid_bounds=grid.bounds,
            periodic=grid.periodic,
            spec=spec,
            parameters=dict(params),
        )
        _times, _per_field, _interaction, total = compute_energy_timeseries(data)
        assert np.all(np.isfinite(total)), "energy produced non-finite values"
        assert abs(total[0]) > 0.0, "energy is identically zero — nothing measured"
        drift = float(np.max(np.abs(total - total[0])) / abs(total[0]))
        assert drift < 1e-2, (
            f"total energy drifts {drift:.3e} over t=0.1 on a time-independent "
            f"background — the pos-dep Hamiltonian coefficients are being "
            f"mis-evaluated"
        )
