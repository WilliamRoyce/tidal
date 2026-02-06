"""Tests for center parameter mutation prevention.

This module verifies that GaussianPulse properly handles the center
parameter, including defensive copying to prevent external mutation
and type validation.
"""

from __future__ import annotations

import pytest

from torsion_gertsenshtein.kgsim.initial_conditions import GaussianPulse

# ==================== Center Parameter Tests ====================


class TestGaussianPulseCenter:
    """Tests for GaussianPulse center parameter handling."""

    def test_not_mutated_by_external_change(self) -> None:
        """Verify modifying the original center list doesn't affect the IC.

        The GaussianPulse should make a defensive copy of the center parameter
        to prevent external mutations from affecting internal state.
        """
        center: list[float] = [10.0, 15.0]
        ic = GaussianPulse(amplitude=1.0, width=2.0, center=center)

        # Modify the original list
        center[0] = 999.0
        center[1] = 888.0

        # IC should have its own copy
        assert ic.center == [10.0, 15.0], "IC center was mutated by external change"

    def test_tuple_conversion(self) -> None:
        """Verify tuples are converted to lists as center parameters."""
        center: tuple[float, float] = (10.0, 15.0)
        ic = GaussianPulse(amplitude=1.0, width=2.0, center=center)

        # Should be converted to list internally
        assert ic.center == [10.0, 15.0]
        assert isinstance(ic.center, list)

    def test_validates_type(self) -> None:
        """Verify non-iterable types are rejected with TypeError."""
        with pytest.raises(TypeError, match="not iterable"):
            GaussianPulse(
                amplitude=1.0,
                width=2.0,
                center=42.0,  # type: ignore[arg-type]
            )
