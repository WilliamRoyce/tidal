"""Tests for FFT-based spectral operators (Phase 3).

Verifies that spectral operators achieve machine-precision accuracy for
band-limited functions and exponential convergence for smooth problems.

Reference: Burns et al. (2020), "Dedalus: A flexible framework for numerical
simulations with spectral methods", Physical Review Research 2:023068.
"""

from __future__ import annotations

import numpy as np
import pytest

from tidal.solver.grid import GridInfo
from tidal.solver.operators import (
    biharmonic,
    cross_derivative,
    directional_laplacian,
    get_spectral,
    gradient,
    laplacian,
    set_spectral,
)
from tidal.symbolic.json_loader import EquationSystem

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _periodic_grid_1d(n: int, length: float = 2 * np.pi) -> GridInfo:
    """Create a 1D periodic grid on [0, length) with n points."""
    return GridInfo(
        bounds=((0.0, length),),
        shape=(n,),
        periodic=(True,),
    )


def _periodic_grid_2d(
    nx: int,
    ny: int,
    lx: float = 2 * np.pi,
    ly: float = 2 * np.pi,
) -> GridInfo:
    """Create a 2D periodic grid."""
    return GridInfo(
        bounds=((0.0, lx), (0.0, ly)),
        shape=(nx, ny),
        periodic=(True, True),
    )


# ---------------------------------------------------------------------------
# Test: spectral state management
# ---------------------------------------------------------------------------


class TestSpectralState:
    """Tests for set_spectral / get_spectral module state."""

    def test_default_off(self) -> None:
        """Spectral mode is off by default."""
        assert not get_spectral()

    def test_toggle_on_off(self) -> None:
        """Can enable and disable spectral mode."""
        set_spectral(True)
        assert get_spectral()
        set_spectral(False)
        assert not get_spectral()


# ---------------------------------------------------------------------------
# Test: spectral gradient (1D)
# ---------------------------------------------------------------------------


class TestSpectralGradient1D:
    """Spectral gradient d/dx on 1D periodic grids."""

    def test_sin_exact(self) -> None:
        """d/dx sin(x) = cos(x) — exact for single Fourier mode."""
        set_spectral(True)
        grid = _periodic_grid_1d(16)
        x = grid.axes_coords(0)
        f = np.sin(x)
        result = gradient(f, 0, grid)
        expected = np.cos(x)
        np.testing.assert_allclose(result, expected, atol=1e-13)

    def test_cos_exact(self) -> None:
        """d/dx cos(x) = -sin(x) — exact."""
        set_spectral(True)
        grid = _periodic_grid_1d(16)
        x = grid.axes_coords(0)
        f = np.cos(x)
        result = gradient(f, 0, grid)
        expected = -np.sin(x)
        np.testing.assert_allclose(result, expected, atol=1e-13)

    def test_sin_3x_exact(self) -> None:
        """d/dx sin(3x) = 3*cos(3x) — exact with N=8 (Nyquist = 4)."""
        set_spectral(True)
        grid = _periodic_grid_1d(8)
        x = grid.axes_coords(0)
        f = np.sin(3 * x)
        result = gradient(f, 0, grid)
        expected = 3.0 * np.cos(3 * x)
        np.testing.assert_allclose(result, expected, atol=1e-13)

    def test_constant_zero(self) -> None:
        """d/dx (constant) = 0."""
        set_spectral(True)
        grid = _periodic_grid_1d(16)
        f = np.ones(grid.shape) * 7.5
        result = gradient(f, 0, grid)
        np.testing.assert_allclose(result, 0.0, atol=1e-13)

    def test_sum_of_modes(self) -> None:
        """d/dx [sin(x) + 0.5*sin(2x)] = cos(x) + cos(2x)."""
        set_spectral(True)
        grid = _periodic_grid_1d(32)
        x = grid.axes_coords(0)
        f = np.sin(x) + 0.5 * np.sin(2 * x)
        result = gradient(f, 0, grid)
        expected = np.cos(x) + np.cos(2 * x)
        np.testing.assert_allclose(result, expected, atol=1e-13)

    def test_non_unit_domain(self) -> None:
        """Gradient on [0, L) with L != 2*pi."""
        set_spectral(True)
        L = 10.0
        grid = _periodic_grid_1d(64, length=L)
        x = grid.axes_coords(0)
        k = 2 * np.pi / L
        f = np.sin(k * x)
        result = gradient(f, 0, grid)
        expected = k * np.cos(k * x)
        np.testing.assert_allclose(result, expected, atol=1e-12)


# ---------------------------------------------------------------------------
# Test: spectral Laplacian (1D)
# ---------------------------------------------------------------------------


class TestSpectralLaplacian1D:
    """Spectral Laplacian d^2/dx^2 on 1D periodic grids."""

    def test_sin_exact(self) -> None:
        """d^2/dx^2 sin(x) = -sin(x) — exact."""
        set_spectral(True)
        grid = _periodic_grid_1d(16)
        x = grid.axes_coords(0)
        f = np.sin(x)
        result = directional_laplacian(f, 0, grid)
        expected = -np.sin(x)
        np.testing.assert_allclose(result, expected, atol=1e-13)

    def test_sin_3x_exact(self) -> None:
        """d^2/dx^2 sin(3x) = -9*sin(3x)."""
        set_spectral(True)
        grid = _periodic_grid_1d(16)
        x = grid.axes_coords(0)
        f = np.sin(3 * x)
        result = directional_laplacian(f, 0, grid)
        expected = -9.0 * np.sin(3 * x)
        np.testing.assert_allclose(result, expected, atol=1e-12)

    def test_full_laplacian_1d(self) -> None:
        """Full Laplacian = directional in 1D."""
        set_spectral(True)
        grid = _periodic_grid_1d(32)
        x = grid.axes_coords(0)
        f = np.sin(2 * x) + np.cos(3 * x)
        result_full = laplacian(f, grid)
        result_dir = directional_laplacian(f, 0, grid)
        np.testing.assert_allclose(result_full, result_dir, atol=1e-14)


# ---------------------------------------------------------------------------
# Test: spectral operators in 2D
# ---------------------------------------------------------------------------


class TestSpectral2D:
    """Spectral operators on 2D periodic grids."""

    def test_gradient_x(self) -> None:
        """d/dx sin(x)*cos(y) = cos(x)*cos(y)."""
        set_spectral(True)
        grid = _periodic_grid_2d(16, 16)
        x = grid.axes_coords(0)
        y = grid.axes_coords(1)
        xx, yy = np.meshgrid(x, y, indexing="ij")
        f = np.sin(xx) * np.cos(yy)
        result = gradient(f, 0, grid)
        expected = np.cos(xx) * np.cos(yy)
        np.testing.assert_allclose(result, expected, atol=1e-13)

    def test_gradient_y(self) -> None:
        """d/dy sin(x)*cos(y) = -sin(x)*sin(y)."""
        set_spectral(True)
        grid = _periodic_grid_2d(16, 16)
        x = grid.axes_coords(0)
        y = grid.axes_coords(1)
        xx, yy = np.meshgrid(x, y, indexing="ij")
        f = np.sin(xx) * np.cos(yy)
        result = gradient(f, 1, grid)
        expected = -np.sin(xx) * np.sin(yy)
        np.testing.assert_allclose(result, expected, atol=1e-13)

    def test_laplacian_2d(self) -> None:
        """nabla^2 sin(x)*sin(y) = -2*sin(x)*sin(y)."""
        set_spectral(True)
        grid = _periodic_grid_2d(16, 16)
        x = grid.axes_coords(0)
        y = grid.axes_coords(1)
        xx, yy = np.meshgrid(x, y, indexing="ij")
        f = np.sin(xx) * np.sin(yy)
        result = laplacian(f, grid)
        expected = -2.0 * np.sin(xx) * np.sin(yy)
        np.testing.assert_allclose(result, expected, atol=1e-12)

    def test_cross_derivative_xy(self) -> None:
        """d^2/(dxdy) sin(x)*sin(y) = cos(x)*cos(y)."""
        set_spectral(True)
        grid = _periodic_grid_2d(16, 16)
        x = grid.axes_coords(0)
        y = grid.axes_coords(1)
        xx, yy = np.meshgrid(x, y, indexing="ij")
        f = np.sin(xx) * np.sin(yy)
        result = cross_derivative(f, 0, 1, grid)
        expected = np.cos(xx) * np.cos(yy)
        np.testing.assert_allclose(result, expected, atol=1e-12)

    def test_biharmonic_2d(self) -> None:
        """nabla^4 sin(x)*sin(y) = 4*sin(x)*sin(y)."""
        set_spectral(True)
        grid = _periodic_grid_2d(16, 16)
        x = grid.axes_coords(0)
        y = grid.axes_coords(1)
        xx, yy = np.meshgrid(x, y, indexing="ij")
        f = np.sin(xx) * np.sin(yy)
        result = biharmonic(f, grid)
        # nabla^4 sin(x)sin(y) = nabla^2(-2 sin(x)sin(y)) = 4 sin(x)sin(y)
        expected = 4.0 * np.sin(xx) * np.sin(yy)
        np.testing.assert_allclose(result, expected, atol=1e-11)


# ---------------------------------------------------------------------------
# Test: exponential convergence
# ---------------------------------------------------------------------------


class TestSpectralConvergence:
    """Verify exponential convergence of spectral operators.

    For smooth (C^inf, periodic) functions, spectral error decreases
    exponentially with N, unlike algebraic O(dx^p) convergence for FD.
    """

    def test_exponential_convergence_gradient(self) -> None:
        """Gradient error decreases exponentially with N."""
        set_spectral(True)
        errors = []
        ns = [8, 16, 32, 64]
        for n in ns:
            grid = _periodic_grid_1d(n)
            x = grid.axes_coords(0)
            # Use a smooth function with multiple modes
            f = np.sin(x) + 0.3 * np.sin(2 * x) + 0.1 * np.sin(3 * x)
            result = gradient(f, 0, grid)
            exact = np.cos(x) + 0.6 * np.cos(2 * x) + 0.3 * np.cos(3 * x)
            errors.append(np.max(np.abs(result - exact)))

        # With N=8 and modes up to k=3 (Nyquist = 4), already near-exact.
        # N=16+: machine precision.
        assert errors[-1] < 1e-12, (
            f"Expected machine precision at N=64, got {errors[-1]}"
        )

        # For a band-limited function with max k=3, N>=8 should be exact
        # (Nyquist = N/2 >= 4 > 3). The error should be near machine epsilon.
        assert errors[0] < 1e-12, (
            f"Expected near-exact at N=8 for k_max=3, got {errors[0]}"
        )

    def test_exponential_convergence_laplacian(self) -> None:
        """Laplacian error decreases exponentially with N for Gaussian."""
        set_spectral(True)
        errors = []
        ns = [16, 32, 64, 128]

        for n in ns:
            L = 2 * np.pi
            grid = _periodic_grid_1d(n, length=L)
            x = grid.axes_coords(0)
            # Periodic Gaussian-like: exp(-sigma * sin^2(x/2))
            sigma = 4.0
            f = np.exp(-sigma * np.sin(x / 2) ** 2)
            result = laplacian(f, grid)

            # Exact Laplacian via chain rule for f = exp(-sigma*sin^2(x/2)):
            # f'' = -(sigma/2)*f * [cos(x) - (sigma/2)*sin^2(x)]
            exact = -(sigma / 2.0) * f * (np.cos(x) - (sigma / 2.0) * np.sin(x) ** 2)
            errors.append(np.max(np.abs(result - exact)))

        # Exponential convergence: each doubling of N should dramatically
        # reduce the error (much faster than algebraic 4x for O(dx^2))
        assert errors[-1] < 1e-10, f"Expected high accuracy at N=128, got {errors[-1]}"

        # Check that the convergence ratio is much better than algebraic
        if errors[0] > 1e-10:
            ratio = errors[0] / errors[1]
            assert ratio > 10, (
                f"Expected exponential convergence ratio > 10, got {ratio:.1f} "
                f"(errors: {errors})"
            )


# ---------------------------------------------------------------------------
# Test: spectral vs FD comparison
# ---------------------------------------------------------------------------


class TestSpectralVsFD:
    """Compare spectral and FD operators — spectral should be more accurate."""

    @pytest.mark.parametrize("fd_order", [2, 4, 6])
    def test_gradient_accuracy_comparison(self, fd_order: int) -> None:
        """Spectral gradient more accurate than any FD order at same N."""
        from tidal.solver.operators import set_fd_order

        n = 32
        grid = _periodic_grid_1d(n)
        x = grid.axes_coords(0)
        f = np.sin(x)
        exact = np.cos(x)

        # FD
        set_fd_order(fd_order)
        set_spectral(False)
        fd_result = gradient(f, 0, grid)
        fd_error = np.max(np.abs(fd_result - exact))

        # Spectral
        set_spectral(True)
        sp_result = gradient(f, 0, grid)
        sp_error = np.max(np.abs(sp_result - exact))

        assert sp_error < fd_error, (
            f"Spectral error {sp_error:.2e} should be less than "
            f"FD order {fd_order} error {fd_error:.2e}"
        )

    def test_laplacian_accuracy_comparison(self) -> None:
        """Spectral Laplacian more accurate than FD-2 at same N."""
        from tidal.solver.operators import set_fd_order

        n = 32
        grid = _periodic_grid_1d(n)
        x = grid.axes_coords(0)
        f = np.sin(2 * x)
        exact = -4.0 * np.sin(2 * x)

        # FD order 2
        set_fd_order(2)
        set_spectral(False)
        fd_result = laplacian(f, grid)
        fd_error = np.max(np.abs(fd_result - exact))

        # Spectral
        set_spectral(True)
        sp_result = laplacian(f, grid)
        sp_error = np.max(np.abs(sp_result - exact))

        assert sp_error < 1e-12, f"Expected machine precision, got {sp_error}"
        assert sp_error < fd_error


# ---------------------------------------------------------------------------
# Test: CLI integration
# ---------------------------------------------------------------------------


class TestSpectralCLI:
    """Test CLI --spectral flag wiring."""

    def test_simulate_spectral_1d(
        self,
        inline_kg_1d_json: pytest.fixture,
        tmp_path: pytest.fixture,
    ) -> None:
        """1D KG simulation with --spectral runs without error."""
        from tidal.cli import main

        output = tmp_path / "spectral_out"
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "m2=1.0",
                "--t-end",
                "1.0",
                "--grid-shape",
                "64",
                "--periodic",
                "--scheme",
                "leapfrog",
                "--spectral",
                "--output",
                str(output),
                "--quiet",
            ]
        )
        assert ret == 0, "Spectral simulation failed"
        assert output.is_dir()

    def test_spectral_rejects_nonperiodic(
        self,
        inline_kg_1d_json: pytest.fixture,
        tmp_path: pytest.fixture,
    ) -> None:
        """--spectral with non-periodic BC should fail."""
        from tidal.cli import main

        output = tmp_path / "spectral_fail"
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "m2=1.0",
                "--t-end",
                "1.0",
                "--grid-shape",
                "64",
                "--no-periodic",
                "--spectral",
                "--output",
                str(output),
                "--quiet",
            ]
        )
        assert ret == 1, "Should fail with non-periodic + spectral"

    def test_spectral_ida_auto_switch(
        self,
        inline_kg_1d_json: pytest.fixture,
        tmp_path: pytest.fixture,
    ) -> None:
        """--spectral with auto scheme should avoid IDA."""
        from tidal.cli import main

        output = tmp_path / "spectral_auto"
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "m2=1.0",
                "--t-end",
                "1.0",
                "--grid-shape",
                "32",
                "--periodic",
                "--spectral",
                "--output",
                str(output),
                "--quiet",
            ]
        )
        assert ret == 0

    def test_spectral_explicit_ida_error(
        self,
        inline_kg_1d_json: pytest.fixture,
        tmp_path: pytest.fixture,
    ) -> None:
        """--spectral with explicit --scheme ida should fail."""
        from tidal.cli import main

        output = tmp_path / "spectral_ida"
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "m2=1.0",
                "--t-end",
                "1.0",
                "--grid-shape",
                "32",
                "--periodic",
                "--scheme",
                "ida",
                "--spectral",
                "--output",
                str(output),
                "--quiet",
            ]
        )
        assert ret == 1


# ---------------------------------------------------------------------------
# Test: spectral metadata persistence
# ---------------------------------------------------------------------------


class TestSpectralMetadata:
    """Test that spectral flag is saved and restored from metadata."""

    def test_spectral_saved_in_metadata(
        self,
        inline_kg_1d_json: pytest.fixture,
        tmp_path: pytest.fixture,
    ) -> None:
        """Running with --spectral saves spectral=True in metadata.json."""
        import json

        from tidal.cli import main

        output = tmp_path / "spectral_meta"
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "m2=1.0",
                "--t-end",
                "0.5",
                "--grid-shape",
                "32",
                "--periodic",
                "--scheme",
                "leapfrog",
                "--spectral",
                "--output",
                str(output),
                "--quiet",
            ]
        )
        assert ret == 0

        metadata_path = output / "metadata.json"
        assert metadata_path.exists()
        metadata = json.loads(metadata_path.read_text())
        assert metadata.get("spectral") is True

    def test_no_spectral_explicit(
        self,
        inline_kg_1d_json: pytest.fixture,
        tmp_path: pytest.fixture,
    ) -> None:
        """--no-spectral forces FD operators, no spectral in metadata."""
        import json

        from tidal.cli import main

        output = tmp_path / "fd_meta"
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "m2=1.0",
                "--t-end",
                "0.5",
                "--grid-shape",
                "32",
                "--periodic",
                "--scheme",
                "leapfrog",
                "--no-spectral",
                "--output",
                str(output),
                "--quiet",
            ]
        )
        assert ret == 0

        metadata_path = output / "metadata.json"
        metadata = json.loads(metadata_path.read_text())
        assert "spectral" not in metadata  # --no-spectral forces FD


# ---------------------------------------------------------------------------
# Test: spectral integration accuracy
# ---------------------------------------------------------------------------


class TestSpectralIntegration:
    """End-to-end: spectral at low N matches FD at high N."""

    def test_spectral_n64_vs_fd2_n512(self) -> None:
        """Spectral N=64 accuracy comparable to or better than FD-2 N=512.

        Both runs use the SAME dt to isolate spatial error differences.
        The analytic solution sin(x)*cos(omega*t) is a single Fourier mode,
        so spectral operators differentiate it exactly — only temporal error
        remains.  FD-2 at N=512 also has small spatial error but measurable.
        """
        from tidal.solver.grid import GridInfo
        from tidal.solver.leapfrog import solve_leapfrog
        from tidal.solver.operators import set_fd_order
        from tidal.symbolic.json_loader import EquationSystem

        spec_dict = {
            "metadata": {"source": "test"},
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
                            },
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian_x",
                                "field": "phi_0",
                            },
                        ],
                    },
                }
            ],
            "coupling": {"mass_matrix_symbolic": [["-1"]]},
        }
        spec = EquationSystem.from_dict(spec_dict)

        L = 2 * np.pi
        t_end = 2.0
        # Use same small dt for both to isolate spatial accuracy
        dt = 0.005

        # FD order 2, N=512
        set_fd_order(2)
        set_spectral(False)
        grid_fd = GridInfo(bounds=((0.0, L),), shape=(512,), periodic=(True,))
        x_fd = grid_fd.axes_coords(0)
        y0_fd = np.zeros(2 * 512)
        y0_fd[:512] = np.sin(x_fd)
        res_fd = solve_leapfrog(
            spec,
            grid_fd,
            y0_fd,
            (0.0, t_end),
            dt,
            parameters={"m2": 1.0},
            snapshot_interval=t_end,
        )

        # Spectral N=64
        set_spectral(True)
        grid_sp = GridInfo(bounds=((0.0, L),), shape=(64,), periodic=(True,))
        x_sp = grid_sp.axes_coords(0)
        y0_sp = np.zeros(2 * 64)
        y0_sp[:64] = np.sin(x_sp)
        res_sp = solve_leapfrog(
            spec,
            grid_sp,
            y0_sp,
            (0.0, t_end),
            dt,
            parameters={"m2": 1.0},
            snapshot_interval=t_end,
        )

        # Analytic: sin(x)*cos(omega*t) where omega = sqrt(k^2 + m^2) = sqrt(2)
        omega = np.sqrt(2.0)

        phi_fd_final = res_fd["y"][-1][:512]
        error_fd = np.max(np.abs(phi_fd_final - np.sin(x_fd) * np.cos(omega * t_end)))

        phi_sp_final = res_sp["y"][-1][:64]
        error_sp = np.max(np.abs(phi_sp_final - np.sin(x_sp) * np.cos(omega * t_end)))

        # Both should have similar temporal error (same dt, same integrator).
        # Spectral has zero spatial error; FD-2 has very small spatial error
        # at N=512.  Temporal error constant depends on omega, so there may
        # be a small constant-factor difference — within 5x means comparable.
        # The key result: spectral N=64 (8x fewer DOFs) matches FD-2 N=512.
        assert error_sp < error_fd * 5, (
            f"Spectral N=64 error {error_sp:.2e} should be within 5x of "
            f"FD-2 N=512 error {error_fd:.2e}"
        )
        # Both should be small (sub 1e-4)
        assert error_sp < 1e-4, f"Spectral N=64 error {error_sp:.2e} too large"
        assert error_fd < 1e-4, f"FD-2 N=512 error {error_fd:.2e} too large"


# ---------------------------------------------------------------------------
# Test: auto-detection of solver settings
# ---------------------------------------------------------------------------


class TestAutoDetection:
    """Verify auto-detection of spectral mode and leapfrog order."""

    def test_spectral_auto_enabled_periodic(
        self,
        inline_kg_1d_json: pytest.fixture,
        tmp_path: pytest.fixture,
    ) -> None:
        """Spectral auto-enables when all BCs are periodic."""
        import json

        from tidal.cli import main

        output = tmp_path / "auto_spectral"
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "m2=1.0",
                "--t-end",
                "0.5",
                "--grid-shape",
                "32",
                "--periodic",
                "--scheme",
                "leapfrog",
                "--output",
                str(output),
                "--quiet",
            ]
        )
        assert ret == 0

        metadata = json.loads((output / "metadata.json").read_text())
        assert metadata.get("spectral") is True, (
            "Spectral should auto-enable for all-periodic BCs"
        )

    def test_spectral_auto_disabled_nonperiodic(
        self,
        inline_kg_1d_json: pytest.fixture,
        tmp_path: pytest.fixture,
    ) -> None:
        """Spectral stays off when BCs are non-periodic (auto mode)."""
        import json

        from tidal.cli import main

        output = tmp_path / "auto_fd"
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "m2=1.0",
                "--t-end",
                "0.5",
                "--grid-shape",
                "32",
                "--no-periodic",
                "--output",
                str(output),
                "--quiet",
            ]
        )
        assert ret == 0

        metadata = json.loads((output / "metadata.json").read_text())
        assert "spectral" not in metadata, (
            "Spectral should not auto-enable with non-periodic BCs"
        )

    def test_spectral_auto_disabled_for_ida(
        self,
        inline_kg_1d_json: pytest.fixture,
        tmp_path: pytest.fixture,
    ) -> None:
        """Spectral auto-disabled when IDA is needed (auto mode)."""
        from tidal.cli import main

        output = tmp_path / "auto_ida"
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "m2=1.0",
                "--t-end",
                "0.5",
                "--grid-shape",
                "32",
                "--periodic",
                "--scheme",
                "ida",
                "--output",
                str(output),
                "--quiet",
            ]
        )
        assert ret == 0  # should not error — spectral silently disabled

    def test_no_spectral_overrides_auto(
        self,
        inline_kg_1d_json: pytest.fixture,
        tmp_path: pytest.fixture,
    ) -> None:
        """--no-spectral disables auto-detection even for periodic BCs."""
        import json

        from tidal.cli import main

        output = tmp_path / "no_spectral"
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "m2=1.0",
                "--t-end",
                "0.5",
                "--grid-shape",
                "32",
                "--periodic",
                "--scheme",
                "leapfrog",
                "--no-spectral",
                "--output",
                str(output),
                "--quiet",
            ]
        )
        assert ret == 0

        metadata = json.loads((output / "metadata.json").read_text())
        assert "spectral" not in metadata  # forced off by --no-spectral


# ---------------------------------------------------------------------------
# Test: system property detection helpers
# ---------------------------------------------------------------------------


class TestSystemPropertyDetection:
    """Test _has_dissipation() and _has_time_dependent_coeffs() helpers."""

    def test_standard_kg_no_dissipation(self) -> None:
        """Standard KG system has no dissipation."""
        from tidal.cli._simulate import _has_dissipation

        spec = EquationSystem.from_dict(
            {
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
                                {
                                    "coefficient": -1.0,
                                    "operator": "identity",
                                    "field": "phi_0",
                                },
                            ],
                        },
                    }
                ],
            }
        )
        assert not _has_dissipation(spec)

    def test_dissipative_system_detected(self) -> None:
        """System with first_derivative_t is detected as dissipative."""
        from tidal.cli._simulate import _has_dissipation

        spec = EquationSystem.from_dict(
            {
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
                                {
                                    "coefficient": -0.1,
                                    "operator": "first_derivative_t",
                                    "field": "phi_0",
                                },
                            ],
                        },
                    }
                ],
            }
        )
        assert _has_dissipation(spec)

    def test_time_independent_system(self) -> None:
        """Standard KG has no time-dependent coefficients."""
        from tidal.cli._simulate import _has_time_dependent_coeffs

        spec = EquationSystem.from_dict(
            {
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
            }
        )
        assert not _has_time_dependent_coeffs(spec)

    def test_time_dependent_detected(self) -> None:
        """System with time_dependent=True term is detected."""
        from tidal.cli._simulate import _has_time_dependent_coeffs

        spec = EquationSystem.from_dict(
            {
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
                                    "time_dependent": True,
                                },
                            ],
                        },
                    }
                ],
            }
        )
        assert _has_time_dependent_coeffs(spec)
