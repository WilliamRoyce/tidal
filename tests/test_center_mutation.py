"""Tests for center parameter mutation prevention."""

import pytest

from torsion_gertsenshtein.kgsim.initial_conditions import GaussianPulse


def test_center_parameter_not_mutated() -> None:
    """Verify that modifying the original center list doesn't affect the IC."""
    center = [10.0, 15.0]
    ic = GaussianPulse(amplitude=1.0, width=2.0, center=center)

    # Modify the original list
    center[0] = 999.0
    center[1] = 888.0

    # IC should have its own copy
    assert ic.center == [10.0, 15.0], "IC center was mutated by external change"


def test_center_parameter_tuple_works() -> None:
    """Verify that tuples work as center parameters."""
    center = (10.0, 15.0)
    ic = GaussianPulse(amplitude=1.0, width=2.0, center=center)

    # Should be converted to list internally
    assert ic.center == [10.0, 15.0]
    assert isinstance(ic.center, list)


def test_center_parameter_validates_type() -> None:
    """Verify that non-list/tuple types are rejected with TypeError."""
    with pytest.raises(TypeError, match="not iterable"):
        GaussianPulse(amplitude=1.0, width=2.0, center=42.0)  # type: ignore[arg-type]
