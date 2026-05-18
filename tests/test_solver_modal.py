"""Tests for the Fourier modal solver (tidal/solver/modal.py).

Tests cover:
1. Eligibility checks (can_use_modal) — flat metric, periodic BCs, operators, etc.
2. Correctness — scalar wave, coupled wave, diffusion, machine precision
3. Resume continuity — split run matches full run
4. Cross-validation — modal vs CVODE agreement
"""

from __future__ import annotations

import copy
from typing import Any, cast

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
        },
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
    """Verify the GH #367 fix: position-dependent + periodic theories auto-route to
    CVODE under ``--scheme auto`` because the modal solver's Fourier-convolution
    path has discretization artifacts that give wrong answers (spurious eigenvalue
    tracks Nyquist wavenumber k_max = π/dx; verified empirically across N).

    The CLI auto-routing lives in ``tidal/cli/_simulate.py`` after
    ``_resolve_scheme``; these tests exercise that path via the ``tidal`` CLI.
    """

    @pytest.fixture
    def e0_args(self, tmp_path: Path) -> list[str]:  # noqa: F821
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

    def test_gh367_reproducer_auto_routes_to_cvode(self, e0_args: list[str]) -> None:
        """GH #367 reproducer at t_end=20: --scheme auto silently routes to CVODE
        and produces the correct decaying h_5 (≈0.005), not 25 million.
        """
        args = [*e0_args, "--t-end", "20"]
        exit_code, output = self._run_simulate(args)
        assert exit_code == 0, f"simulate failed:\n{output}"
        # Auto-routing notice must appear
        assert "auto-routing to CVODE" in output, (
            f"expected auto-route notice; got:\n{output}"
        )
        # Scheme line confirms CVODE
        assert "Auto-selected solver: cvode" in output, (
            f"expected auto-selected cvode; got:\n{output}"
        )
        # The physical result: h_5 should be ≈0.005 (CVODE truth)
        # Format from CLI: "h_5: peak 0.0099 → 0.0050 (ratio: 0.5035)"
        import re

        match = re.search(r"h_5:.*→\s*([0-9.eE+-]+)", output)
        assert match is not None, f"h_5 result line not found in:\n{output}"
        h5_final = float(match.group(1))
        assert 0.003 < h5_final < 0.007, (
            f"h_5={h5_final} not in expected CVODE range [0.003, 0.007]"
        )

    def test_explicit_modal_at_t20_raises_diverged(self, e0_args: list[str]) -> None:
        """--scheme modal at the broken regime hits the post-evolution divergence
        check (commit a2cdaa5) and surfaces the cvode workaround message.
        """
        args = [*e0_args, "--scheme", "modal", "--t-end", "20"]
        exit_code, output = self._run_simulate(args)
        # CLI returns non-zero on simulation failure, but the error message must
        # mention the workaround (we don't assert on exit_code because the CLI
        # may return 0 with [ERROR] printed — check the message itself).
        del exit_code  # exit code semantics vary; message check is authoritative
        assert "Simulation diverged" in output, (
            f"expected diverged message; got:\n{output}"
        )
        assert "--scheme cvode" in output, (
            f"expected workaround message; got:\n{output}"
        )

    def test_explicit_modal_safe_regime_does_not_raise(
        self, e0_args: list[str]
    ) -> None:
        """--scheme modal at small t_end (within the modal-solver's
        amplitude-growth window) runs to completion without raising
        SimulationDivergedError. The result will not match CVODE for
        position-dependent backgrounds (the discretization artifact exists at
        all t_end > 0), but the user has explicitly opted in to --scheme modal,
        and the post-evolution amplitude-growth check (10⁶× threshold) does
        not fire at small t_end. This documents the override-still-runs
        behaviour.
        """
        args = [*e0_args, "--scheme", "modal", "--t-end", "1"]
        exit_code, output = self._run_simulate(args)
        assert exit_code == 0, f"simulate failed unexpectedly:\n{output}"
        assert "Scheme: modal" in output
        assert "Simulation diverged" not in output


class TestOverridePosDepPeriodicScheme:
    """Unit tests for ``_override_pos_dep_periodic_scheme`` (GH #367 auto-route).

    Position-dependent + periodic theories must NOT run through the modal
    solver's Fourier-convolution path under ``--scheme auto``. The override
    sends them to IDA when algebraic constraints (``time_order==0``) are
    present, or CVODE otherwise.
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

    def test_pos_dep_no_constraint_routes_to_cvode(self) -> None:
        """Pos-dep + no constraint: CVODE is the correct fallback (faster
        than IDA and correct in physical space).
        """
        from tidal.cli._simulate import _override_pos_dep_periodic_scheme

        spec = self._pos_dep_diffusion()
        new_scheme, msg = _override_pos_dep_periodic_scheme("modal", "auto", spec)
        assert new_scheme == "cvode"
        assert msg is not None
        assert "auto-routing to CVODE" in msg
        assert "GH #367" in msg

    def test_pos_dep_with_constraint_routes_to_ida(self) -> None:
        """Pos-dep + constraint: IDA is required (CVODE would freeze the
        algebraic constraint at IC and silently produce wrong physics).
        """
        from tidal.cli._simulate import _override_pos_dep_periodic_scheme

        spec = self._pos_dep_with_constraint()
        new_scheme, msg = _override_pos_dep_periodic_scheme("modal", "auto", spec)
        assert new_scheme == "ida"
        assert msg is not None
        assert "auto-routing to IDA" in msg
        assert "constraints" in msg
        assert "GH #367" in msg
