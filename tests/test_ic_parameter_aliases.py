"""Tests for parameter name aliases in initial conditions."""

from __future__ import annotations

import numpy as np
import pytest

from torsion_gertsenshtein.kgsim import (
    GaussianPulse,
    GridConfig,
    RingPulse2D,
    make_grid,
)


def test_gaussian_pulse_width_sigma_aliases() -> None:
    """Test that width and sigma parameters are equivalent for GaussianPulse."""
    grid = make_grid(GridConfig(dim=1, shape=(64,), bounds=((0.0, 50.0),)))

    # Create with width parameter
    ic_width = GaussianPulse(amplitude=1.0, width=2.0, center=[25.0])
    state_width = ic_width.build(grid)

    # Create with sigma parameter (alias)
    ic_sigma = GaussianPulse(amplitude=1.0, sigma=2.0, center=[25.0])
    state_sigma = ic_sigma.build(grid)

    # Should produce identical results
    assert np.allclose(state_width[0].data, state_sigma[0].data)
    assert np.allclose(state_width[1].data, state_sigma[1].data)

    # Both should store as width internally
    assert ic_width.width == 2.0  # noqa: PLR2004
    assert ic_sigma.width == 2.0  # noqa: PLR2004


def test_gaussian_pulse_both_params_raises_error() -> None:
    """Test that specifying both width and sigma raises ValueError."""
    with pytest.raises(ValueError, match="Cannot specify both 'width' and 'sigma'"):
        GaussianPulse(amplitude=1.0, width=2.0, sigma=3.0, center=[25.0])


def test_ring_pulse_2d_width_sigma_aliases() -> None:
    """Test that width and sigma parameters are equivalent for RingPulse2D."""
    grid = make_grid(
        GridConfig(dim=2, shape=(32, 32), bounds=((-10.0, 10.0), (-10.0, 10.0)))
    )

    # Create with sigma parameter (default name)
    ic_sigma = RingPulse2D(amplitude=1.0, initial_radius=5.0, sigma=1.5)
    state_sigma = ic_sigma.build(grid)

    # Create with width parameter (alias)
    ic_width = RingPulse2D(amplitude=1.0, initial_radius=5.0, width=1.5)
    state_width = ic_width.build(grid)

    # Should produce identical results
    assert np.allclose(state_sigma[0].data, state_width[0].data)
    assert np.allclose(state_sigma[1].data, state_width[1].data)

    # Both should store as width internally (consistent with GaussianPulse)
    assert ic_sigma.width == 1.5  # noqa: PLR2004
    assert ic_width.width == 1.5  # noqa: PLR2004


def test_ring_pulse_2d_both_params_raises_error() -> None:
    """Test that specifying both sigma and width raises ValueError."""
    with pytest.raises(ValueError, match="Cannot specify both 'sigma' and 'width'"):
        RingPulse2D(amplitude=1.0, initial_radius=5.0, sigma=1.5, width=2.5)


@pytest.mark.parametrize("use_sigma", [True, False])
def test_consistent_naming_across_ics(use_sigma: bool) -> None:  # noqa: FBT001
    """Test that both ICs can use either width or sigma consistently."""
    grid_1d = make_grid(GridConfig(dim=1, shape=(64,), bounds=((0.0, 50.0),)))
    grid_2d = make_grid(
        GridConfig(dim=2, shape=(32, 32), bounds=((-10.0, 10.0), (-10.0, 10.0)))
    )

    if use_sigma:
        # Both use sigma
        ic_gaussian = GaussianPulse(amplitude=1.0, sigma=2.0)
        ic_ring = RingPulse2D(amplitude=1.0, initial_radius=5.0, sigma=2.0)
    else:
        # Both use width
        ic_gaussian = GaussianPulse(amplitude=1.0, width=2.0)
        ic_ring = RingPulse2D(amplitude=1.0, initial_radius=5.0, width=2.0)

    # Both should build successfully
    state_gaussian = ic_gaussian.build(grid_1d)
    state_ring = ic_ring.build(grid_2d)

    # Verify they have expected structure
    assert len(state_gaussian) == 2  # noqa: PLR2004
    assert len(state_ring) == 2  # noqa: PLR2004

    # Both should store as .width for consistency
    assert ic_gaussian.width == 2.0  # noqa: PLR2004
    assert ic_ring.width == 2.0  # noqa: PLR2004


class TestNaNInfValidation:
    """Tests for NaN and infinity validation in initial conditions."""

    def test_gaussian_nan_amplitude_raises(self) -> None:
        """Test that NaN amplitude raises ValueError."""
        with pytest.raises(ValueError, match="amplitude must be finite"):
            GaussianPulse(amplitude=float("nan"), width=2.0)

    def test_gaussian_inf_amplitude_raises(self) -> None:
        """Test that infinite amplitude raises ValueError."""
        with pytest.raises(ValueError, match="amplitude must be finite"):
            GaussianPulse(amplitude=float("inf"), width=2.0)

    def test_gaussian_nan_width_raises(self) -> None:
        """Test that NaN width raises ValueError."""
        with pytest.raises(ValueError, match="width must be finite"):
            GaussianPulse(amplitude=1.0, width=float("nan"))

    def test_gaussian_inf_width_raises(self) -> None:
        """Test that infinite width raises ValueError."""
        with pytest.raises(ValueError, match="width must be finite"):
            GaussianPulse(amplitude=1.0, width=float("inf"))

    def test_gaussian_nan_initial_velocity_raises(self) -> None:
        """Test that NaN initial_velocity raises ValueError."""
        with pytest.raises(ValueError, match="initial_velocity must be finite"):
            GaussianPulse(amplitude=1.0, width=2.0, initial_velocity=float("nan"))

    def test_ring_nan_amplitude_raises(self) -> None:
        """Test that NaN amplitude raises ValueError for RingPulse2D."""
        with pytest.raises(ValueError, match="amplitude must be finite"):
            RingPulse2D(amplitude=float("nan"), initial_radius=5.0, sigma=1.0)

    def test_ring_inf_sigma_raises(self) -> None:
        """Test that infinite sigma raises ValueError for RingPulse2D."""
        with pytest.raises(ValueError, match="sigma must be finite"):
            RingPulse2D(amplitude=1.0, initial_radius=5.0, sigma=float("inf"))

    def test_ring_nan_initial_radius_raises(self) -> None:
        """Test that NaN initial_radius raises ValueError for RingPulse2D."""
        with pytest.raises(ValueError, match="initial_radius must be finite"):
            RingPulse2D(amplitude=1.0, initial_radius=float("nan"), sigma=1.0)

    def test_negative_inf_also_rejected(self) -> None:
        """Test that negative infinity is also rejected."""
        with pytest.raises(ValueError, match="amplitude must be finite"):
            GaussianPulse(amplitude=float("-inf"), width=2.0)
