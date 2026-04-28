"""Tests for the versioned stability-probe profile registry.

Pinned regression fixtures for the campaign's pre-flight tachyonic
guard. Every change to the probe IC construction, threshold semantics,
or borderline detection must keep these tests passing or explicitly
update them with reasoning recorded in #323.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pytest

from tidal.measurement._stability import check_conversion_stability
from tidal.measurement._stability_profile import (
    DEFAULT_PROFILE_NAME,
    PROFILES,
    StabilityProfile,
    adhoc_profile,
    get_profile,
)

# ---------------------------------------------------------------------
# Registry sanity
# ---------------------------------------------------------------------


class TestProfileRegistry:
    def test_default_profile_is_legacy(self) -> None:
        """Default must remain v1 until Stage C validation completes (#323)."""
        assert DEFAULT_PROFILE_NAME == "v1-unit-0.3"
        assert DEFAULT_PROFILE_NAME in PROFILES

    def test_registry_contains_v1_and_v2(self) -> None:
        """Both probe modes available so they can run in parallel during Stage C."""
        assert "v1-unit-0.3" in PROFILES
        assert "v2-consistent-0.1" in PROFILES

    def test_v1_legacy_thresholds(self) -> None:
        p = get_profile("v1-unit-0.3")
        assert p.ic_mode == "unit"
        assert p.threshold == pytest.approx(0.3)
        assert p.t_test == pytest.approx(10.0)
        assert p.perturbativity_p_max == pytest.approx(0.5)

    def test_v2_consistent_thresholds(self) -> None:
        p = get_profile("v2-consistent-0.1")
        assert p.ic_mode == "consistent"
        assert p.threshold == pytest.approx(0.1)

    def test_unknown_profile_raises(self) -> None:
        with pytest.raises(KeyError, match="Unknown stability profile"):
            get_profile("not-a-real-profile")

    def test_borderline_strip_is_top_30_percent(self) -> None:
        """borderline_low = 0.7·threshold by definition; pin the contract."""
        p = get_profile("v1-unit-0.3")
        assert p.borderline_low == pytest.approx(0.21)

    def test_adhoc_profile_synthesizes_unit_default(self) -> None:
        p = adhoc_profile(threshold=0.5, t_test=5.0)
        assert p.ic_mode == "unit"
        assert p.threshold == pytest.approx(0.5)
        assert p.t_test == pytest.approx(5.0)
        # Ad-hoc names mark the profile as non-canonical.
        assert p.name.startswith("adhoc-")

    def test_profile_is_immutable(self) -> None:
        p = get_profile("v1-unit-0.3")
        with pytest.raises((AttributeError, TypeError)):
            p.threshold = 0.5  # type: ignore[misc]


# ---------------------------------------------------------------------
# Probe behaviour parity & divergence
# ---------------------------------------------------------------------


def _load_t1_spec_and_grid() -> tuple[Any, Any]:
    from tidal.solver.grid import GridInfo
    from tidal.symbolic.json_loader import load_equation_system

    repo = Path(__file__).resolve().parents[1]
    spec = load_equation_system(repo / "examples/data/dark_photon_plasma.json")
    grid = GridInfo(shape=(32,), bounds=((0.0, 50.0),), periodic=(True,))
    return spec, grid


class TestProbeBehaviour:
    """Pin the v1 vs v2 verdicts on canonical fixtures.

    Updates here must be paired with an updated comment recording the
    physics reason and the simulation-truth result that justifies the
    change.
    """

    def test_v1_unit_ic_decoupling_stable(self) -> None:
        """T1 deltam=0 (decoupling limit) → both profiles stable."""
        from tidal.cli._simulate import (
            _parse_params,  # pyright: ignore[reportPrivateUsage]
        )

        spec, grid = _load_t1_spec_and_grid()
        params = _parse_params(["kappa=1.0", "B0=0.01"], spec)
        params.update({"mA2": 0.5, "deltam": 0.0, "xi": 0.5, "alpha3": 0.1})

        for name in ("v1-unit-0.3", "v2-consistent-0.1"):
            profile = get_profile(name)
            result = check_conversion_stability(
                spec,
                grid,
                params,
                source="h_5",
                target="a_1",
                profile=profile,
            )
            assert result.stable, f"{name}: decoupling must be stable"
            assert result.profile_name == name

    def test_v1_catches_known_ghost(self) -> None:
        """v4 amplify MAP point (ghost-contaminated, A~1454) → still rejected."""
        from tidal.cli._simulate import (
            _parse_params,  # pyright: ignore[reportPrivateUsage]
        )

        spec, grid = _load_t1_spec_and_grid()
        params = _parse_params(["kappa=1.0", "B0=0.01"], spec)
        params.update(
            {
                "mA2": 0.001018,
                "deltam": 0.283592,
                "xi": 0.320470,
                "alpha3": 0.002783,
            },
        )

        result = check_conversion_stability(
            spec,
            grid,
            params,
            source="h_5",
            target="a_1",
            profile=get_profile("v1-unit-0.3"),
        )
        assert not result.stable
        assert result.profile_name == "v1-unit-0.3"

    def test_back_compat_no_profile_uses_unit_mode(self) -> None:
        """Legacy callers (no profile arg) get unit-IC + adhoc profile."""
        from tidal.cli._simulate import (
            _parse_params,  # pyright: ignore[reportPrivateUsage]
        )

        spec, grid = _load_t1_spec_and_grid()
        params = _parse_params(["kappa=1.0", "B0=0.01"], spec)
        params.update({"mA2": 0.5, "deltam": 0.0, "xi": 0.5, "alpha3": 0.1})

        result = check_conversion_stability(
            spec,
            grid,
            params,
            source="h_5",
            target="a_1",
            t_test=10.0,
        )
        assert result.profile_name.startswith("adhoc-unit-")

    def test_borderline_flag_clean_in_decoupling(self) -> None:
        """Decoupling limit has γ_eff well below borderline strip → flag False."""
        from tidal.cli._simulate import (
            _parse_params,  # pyright: ignore[reportPrivateUsage]
        )

        spec, grid = _load_t1_spec_and_grid()
        params = _parse_params(["kappa=1.0", "B0=0.01"], spec)
        params.update({"mA2": 0.5, "deltam": 0.0, "xi": 0.5, "alpha3": 0.1})

        result = check_conversion_stability(
            spec,
            grid,
            params,
            source="h_5",
            target="a_1",
            profile=get_profile("v1-unit-0.3"),
        )
        assert result.borderline is False

    def test_consistent_mode_skips_high_k_no_ic_modes(self) -> None:
        """v2 consistent IC must not flag high-k modes the IC doesn't excite.

        Probe at v2 with IC at fundamental k → the rfft IC weight at high k
        is ≈ 0, so those modes are skipped. This is the architectural fix
        that the legacy v1 all-k unit-IC probe lacks (#323).
        """
        from tidal.cli._simulate import (
            _parse_params,  # pyright: ignore[reportPrivateUsage]
        )

        spec, grid = _load_t1_spec_and_grid()
        params = _parse_params(["kappa=1.0", "B0=0.01"], spec)
        params.update({"mA2": 0.5, "deltam": 0.0, "xi": 0.5, "alpha3": 0.1})

        result = check_conversion_stability(
            spec,
            grid,
            params,
            source="h_5",
            target="a_1",
            ic_wavevector=2.0 * np.pi / 50.0,
            profile=get_profile("v2-consistent-0.1"),
        )
        # The decoupling-limit run is stable under both profiles; under v2
        # the per-k IC weights from rfft(cos) decay to ≈ 0 above the IC
        # wavenumber, so the probe must not over-reject.
        assert result.stable


# ---------------------------------------------------------------------
# Result schema
# ---------------------------------------------------------------------


class TestResultSchema:
    def test_result_has_borderline_and_profile_fields(self) -> None:
        """Schema additions for the post-hoc audit pipeline (#323)."""
        from tidal.cli._simulate import (
            _parse_params,  # pyright: ignore[reportPrivateUsage]
        )

        spec, grid = _load_t1_spec_and_grid()
        params = _parse_params(["kappa=1.0", "B0=0.01"], spec)
        params.update({"mA2": 0.5, "deltam": 0.0, "xi": 0.5, "alpha3": 0.1})

        result = check_conversion_stability(
            spec,
            grid,
            params,
            source="h_5",
            target="a_1",
            profile=get_profile("v1-unit-0.3"),
        )
        assert hasattr(result, "borderline")
        assert hasattr(result, "profile_name")
        assert isinstance(result.borderline, bool)
        assert isinstance(result.profile_name, str)


# ---------------------------------------------------------------------
# Type-safety guards
# ---------------------------------------------------------------------


class TestProfileTypeSafety:
    def test_non_profile_arg_rejected(self) -> None:
        from tidal.cli._simulate import (
            _parse_params,  # pyright: ignore[reportPrivateUsage]
        )

        spec, grid = _load_t1_spec_and_grid()
        params = _parse_params(["kappa=1.0", "B0=0.01"], spec)
        params.update({"mA2": 0.5})

        with pytest.raises(TypeError, match="must be a StabilityProfile"):
            check_conversion_stability(
                spec,
                grid,
                params,
                source="h_5",
                target="a_1",
                profile="v1-unit-0.3",  # type: ignore[arg-type]
            )

    def test_profile_dataclass_field_count_pinned(self) -> None:
        """Adding a field is a behaviour change — bump version, update tests."""
        import dataclasses

        fields = {f.name for f in dataclasses.fields(StabilityProfile)}
        assert fields == {
            "name",
            "ic_mode",
            "threshold",
            "t_test",
            "perturbativity_p_max",
            "notes",
        }
