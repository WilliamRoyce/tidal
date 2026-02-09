"""Tests for validation edge cases added to improve robustness."""

from __future__ import annotations

import pytest
from pde import CartesianGrid

from tidal.kgsim.animation_builder import (
    AnimationBuilder,
    AnimationConfig,
)
from tidal.kgsim.config import GridConfig
from tidal.kgsim.grids import make_grid


class TestEmptyGridValidation:
    """Test that empty grids are properly rejected."""

    def test_multi_gaussian_rejects_empty_grid(self) -> None:
        """multi_gaussian should reject grids with zero points (via py-pde)."""
        # Create grid with zero shape - py-pde will reject this
        grid_cfg = GridConfig(dim=1, shape=(0,), bounds=((0.0, 1.0),), periodic=True)

        # py-pde rejects empty grids at construction time
        with pytest.raises(ValueError, match="not a valid number of support points"):
            make_grid(grid_cfg)

    def test_multi_gaussian_2d_rejects_empty_grid(self) -> None:
        """multi_gaussian_2d should reject grids with zero points in any dimension."""
        # Create grid with zero in one dimension - py-pde will reject this
        grid_cfg = GridConfig(
            dim=2, shape=(10, 0), bounds=((0.0, 1.0), (0.0, 1.0)), periodic=True
        )

        with pytest.raises(ValueError, match="not a valid number of support points"):
            make_grid(grid_cfg)


class TestBoundsValidation:
    """Test that invalid bounds are properly rejected."""

    def test_grid_config_rejects_invalid_bounds(self) -> None:
        """GridConfig should reject bounds where lower >= upper."""
        with pytest.raises(ValueError, match=r"lower bound.*>= upper bound"):
            GridConfig(dim=1, shape=(100,), bounds=((10.0, 5.0),), periodic=True)

    def test_grid_config_rejects_equal_bounds(self) -> None:
        """GridConfig should reject bounds where lower == upper."""
        with pytest.raises(ValueError, match=r"lower bound.*>= upper bound"):
            GridConfig(dim=1, shape=(100,), bounds=((5.0, 5.0),), periodic=True)

    def test_grid_config_accepts_valid_bounds(self) -> None:
        """GridConfig should accept valid bounds."""
        # This should not raise
        cfg = GridConfig(dim=1, shape=(100,), bounds=((0.0, 10.0),), periodic=True)
        assert cfg.bounds[0][0] < cfg.bounds[0][1]


class TestPathTraversalValidation:
    """Test that path traversal is prevented in AnimationConfig."""

    def test_animation_config_rejects_path_traversal(self) -> None:
        """AnimationConfig should reject paths with .. traversal."""
        with pytest.raises(ValueError, match="path traversal pattern"):
            AnimationConfig(output_path="../etc/passwd")

    def test_animation_config_rejects_nested_traversal(self) -> None:
        """AnimationConfig should reject nested .. patterns."""
        with pytest.raises(ValueError, match="path traversal pattern"):
            AnimationConfig(output_path="output/../../sensitive/file.mp4")

    def test_animation_config_accepts_valid_path(self) -> None:
        """AnimationConfig should accept valid paths without traversal."""
        # This should not raise
        cfg = AnimationConfig(output_path="output/animation.mp4")
        assert str(cfg.output_path) == "output/animation.mp4"


class TestDivisionByZeroFix:
    """Test that division by zero is handled in FPS calculation."""

    def test_choose_writer_handles_zero_t_end(self) -> None:
        """_choose_writer should handle t_end=0 without division by zero."""
        # This should not raise ZeroDivisionError
        writer, ext = AnimationBuilder._choose_writer(snap_count=10, t_end=0.0)  # pyright: ignore[reportPrivateUsage]
        assert writer.fps == 1
        assert ext in {".mp4", ".gif"}

    def test_choose_writer_normal_case(self) -> None:
        """_choose_writer should work normally for positive t_end."""
        writer, ext = AnimationBuilder._choose_writer(snap_count=100, t_end=20.0)  # pyright: ignore[reportPrivateUsage]
        # fps should be calculated based on snap_count and t_end
        assert writer.fps >= 1
        assert ext in {".mp4", ".gif"}


class TestStorageBoundsCheck:
    """Test that storage bounds are validated before access."""

    def test_empty_storage_handled(self) -> None:
        """Animation methods should validate storage is non-empty."""
        from pde import MemoryStorage

        grid = CartesianGrid([(0, 10), (0, 10)], 10)
        storage = MemoryStorage()  # Empty storage

        builder = AnimationBuilder(storage, grid)
        config = AnimationConfig(output_path="test.mp4")

        with pytest.raises(ValueError, match="storage is empty"):
            builder.create_2d_heatmap_animation(config)
