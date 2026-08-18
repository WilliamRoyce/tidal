"""Curved-metric derivation coverage and derivation-integrity guards.

Regression cover for two defects that were silent at derivation time and only
surfaced later as wrong physics:

* **GH #394** -- ``ComponentEulerLagrange`` varied the scalar component
  Lagrangian instead of the density ``sqrt|g| * L``, dropping the Christoffel
  first-derivative terms for curvilinear metrics.  ``spherical_kg_1d`` lost its
  ``(2/r) d_r`` term and stopped conserving energy (1.6e-01), while the energy
  integral still used the correct ``r^2`` volume element.
* **GH #397/#381** -- pre-fix exports left a vector field's temporal component
  un-normalized, giving ``d2_t a_0 = -laplacian a_0``: a temporal-only tachyon.

The guard tests below are as much a regression suite for the *detection* logic
as for the specs: every exclusion corresponds to a case that a naive check
misreports.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pytest

from tidal.cli._validate import (
    _check_temporal_component_sign,
    _check_volume_element_consistency,
)
from tidal.symbolic.json_loader import load_equation_system

if TYPE_CHECKING:
    from tidal.symbolic.json_loader import EquationSystem

_EXAMPLES = Path(__file__).parent.parent / "examples" / "data"


def _load(name: str) -> EquationSystem:
    path = _EXAMPLES / f"{name}.json"
    if not path.exists():
        pytest.skip(f"{name}.json not found")
    return load_equation_system(path)


def _self_terms(spec: EquationSystem, field: str, operator_prefix: str) -> list[Any]:
    for eq in spec.equations:
        if eq.field_name != field:
            continue
        return [
            t
            for t in eq.rhs_terms
            if t.field == field and t.operator.startswith(operator_prefix)
        ]
    return []


# --------------------------------------------------------------------------
# GH #394 -- the Christoffel terms supplied by the sqrt|g| measure
# --------------------------------------------------------------------------


class TestCurvilinearChristoffelTerms:
    """The first-derivative terms that only the volume measure produces."""

    def test_spherical_kg_1d_has_radial_christoffel_term(self) -> None:
        """Radial KG must carry (2/r) d_r, else it is the flat wave operator."""
        spec = _load("spherical_kg_1d")
        terms = _self_terms(spec, "phi_0", "gradient")
        assert terms, "spherical_kg_1d lost its gradient term (GH #394)"
        assert any(
            t.coefficient_symbolic
            and "2/x[]" in t.coefficient_symbolic.replace(" ", "")
            for t in terms
        ), f"expected 2/x[], got {[t.coefficient_symbolic for t in terms]}"

    def test_polar_kg_has_radial_christoffel_term(self) -> None:
        """Polar KG must carry (1/r) d_r, matching the historical reference."""
        spec = _load("polar_kg")
        terms = _self_terms(spec, "phi_0", "gradient")
        assert terms, "polar_kg lost its gradient term (GH #394)"
        assert any(
            t.coefficient_symbolic and "x[]" in t.coefficient_symbolic for t in terms
        ), (
            f"expected an x[]-dependent coefficient, got {[t.coefficient_symbolic for t in terms]}"
        )

    @pytest.mark.parametrize("name", ["spherical_kg_1d", "polar_kg"])
    def test_curvilinear_specs_have_non_constant_volume_element(
        self, name: str
    ) -> None:
        spec = _load(name)
        assert spec.canonical is not None
        ve = spec.canonical.volume_element
        assert ve is not None, f"{name} lost its volume_element"
        assert any(f"{c}[]" in str(ve) for c in spec.coordinates), (
            f"{name} volume_element {ve!r} is constant; expected coordinate dependence"
        )

    @pytest.mark.parametrize("name", ["spherical_kg_1d", "polar_kg"])
    def test_measure_and_equations_agree(self, name: str) -> None:
        """The #393 failure mode: curved energy weight, flat evolution operator."""
        assert _check_volume_element_consistency(_load(name)) == []


# --------------------------------------------------------------------------
# GH #394 guard -- non-constant measure implies Christoffel gradient terms
# --------------------------------------------------------------------------


class TestVolumeElementGuard:
    def test_fires_when_measure_present_but_gradients_absent(self) -> None:
        """A curved measure with a purely flat operator must be rejected."""
        spec = _load("spherical_kg_1d")
        stripped = json.loads((_EXAMPLES / "spherical_kg_1d.json").read_text())
        for eq in stripped["equations"]:
            eq["rhs"]["terms"] = [
                t
                for t in eq["rhs"]["terms"]
                if not t["operator"].startswith("gradient")
            ]
        import tempfile

        with tempfile.NamedTemporaryFile(
            "w", encoding="utf-8", suffix=".json", delete=False
        ) as fh:
            json.dump(stripped, fh)
            broken = Path(fh.name)
        try:
            errors = _check_volume_element_consistency(load_equation_system(broken))
            assert errors, "guard failed to fire on a spec with the measure dropped"
            assert "#394" in errors[0]
        finally:
            broken.unlink()
        assert spec is not None


# --------------------------------------------------------------------------
# GH #397 guard -- temporal vs spatial component signs
# --------------------------------------------------------------------------

# Specs re-derived after the GH #381 fix, or never affected.  Each entry after
# the first block is a case that a naive check misreports -- see the docstrings.
KNOWN_GOOD = [
    # re-derived this session
    "dark_photon_plasma",
    "torsion_dark_photon",
    "torsion_dark_photon_fv",
    "gertsenshtein_localized",
    "gertsenshtein_ungauged",
    "gertsenshtein_ungauged_e_dual_gaussian",
    "torsion_gertsenshtein_b5_zero",
    # re-derived by the original GH #381 rollout
    "gertsenshtein",
    "gertsenshtein_proca",
    "gertsenshtein_e0_dual_gaussian",
    # Euler-Heisenberg: a_1..a_3 carry lap=-1 with kin=-1+2*B0^2*rho, i.e.
    # orientation +1.  Comparing raw coefficients marks these inconsistent.
    "euler_heisenberg",
    "euler_heisenberg_e_dual_gaussian",
    "gertsenshtein_eh",
    "gertsenshtein_eh_top",
    # constraint LHS (a_0 is algebraic, not evolved): an overall sign is
    # conventional there, so the guard must not look at them.
    "chern_simons_3d",
    "proca_background",
    "scalar_vector_coupling",
]


class TestTemporalSignGuard:
    @pytest.mark.parametrize("name", KNOWN_GOOD)
    def test_no_false_positive(self, name: str) -> None:
        """Every one of these misreports under a naive sign comparison."""
        assert _check_temporal_component_sign(_load(name)) == [], (
            f"{name} incorrectly flagged"
        )

    def test_detects_an_inverted_temporal_component(self) -> None:
        """Invert a_0 in a known-good spec; the guard must catch it."""
        data = json.loads((_EXAMPLES / "dark_photon_plasma.json").read_text())
        for eq in data["equations"]:
            if eq["field"] != "a_0":
                continue
            for term in eq["rhs"]["terms"]:
                if isinstance(term.get("coefficient"), (int, float)):
                    term["coefficient"] = -term["coefficient"]
                if term.get("coefficient_symbolic"):
                    term["coefficient_symbolic"] = f"-({term['coefficient_symbolic']})"
        import tempfile

        with tempfile.NamedTemporaryFile(
            "w", encoding="utf-8", suffix=".json", delete=False
        ) as fh:
            json.dump(data, fh)
            broken = Path(fh.name)
        try:
            errors = _check_temporal_component_sign(load_equation_system(broken))
            assert errors, "guard failed to detect an inverted temporal component"
            assert "a_0" in errors[0]
            assert "#397" in errors[0]
        finally:
            broken.unlink()

    def test_torsion_family_non_uniformity_is_not_flagged(self) -> None:
        """Rank-3 torsion components legitimately differ in normalization.

        ``torsion_dark_photon`` came through re-derivation with its torsion
        sector unchanged, so that structure is what correct code produces.
        Flagging it would make the guard unusable on every PGT theory.
        """
        errors = _check_temporal_component_sign(_load("torsion_dark_photon"))
        assert not any("'t'" in e for e in errors), (
            f"torsion family incorrectly flagged: {errors}"
        )
