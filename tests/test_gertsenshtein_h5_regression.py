"""Regression tests for the h_5 graviton kinetic term fix (issue #250).

These tests load the committed JSONs and verify that the structural
properties expected from a correct derivation are preserved. They run fast
(no wolframscript) and catch regressions in the canonical pipeline that
would silently misclassify graviton or photon fields.

Background: commit 9f496ef (#250) fixed an issue where h_5 = h_{xy}
(graviton cross-polarization) was being classified as a constraint
(time_order=0) instead of a dynamical wave field (time_order=2). The fix
introduced ``LagrangianFullExpand`` — a fully-expanded copy of L^(2)
processed by ToCanonical, which preserves the non-trace graviton kinetic
structure needed for off-diagonal graviton components.

These tests protect against future changes to:
- ``_wls_deferred_field_expand`` (which controls F → d(a) substitution)
- The ``LagrangianFullExpand`` Module block in ``_wls_linearize_from_lagrangian``
- ``_wls_canonical_phase_b`` (which uses LagrangianFullExpand for lagForCanon)
- The xPert/ToCanonical interaction with abstract deferred fields
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest


def _load_spec(name: str) -> dict:
    """Load a committed example JSON spec."""
    path = Path(__file__).parent.parent / "examples" / "data" / f"{name}.json"
    if not path.exists():
        pytest.skip(f"Example JSON {name}.json not present")
    return json.loads(path.read_text())


def _equation_for(spec: dict, field_name: str) -> dict:
    """Return the equation dict for a given field name."""
    for eq in spec["equations"]:
        if eq["field"] == field_name:
            return eq
    msg = f"Field {field_name} not in spec equations"
    raise AssertionError(msg)


class TestGertsenshteinH5Regression:
    """Regression tests for #250: h_5 kinetic terms in base Gertsenshtein."""

    def test_h5_is_dynamical(self) -> None:
        """h_5 must have a wave equation (d2_t LHS), not a constraint."""
        spec = _load_spec("gertsenshtein")
        eq = _equation_for(spec, "h_5")
        lhs_expr = eq["lhs"]["expression"]
        assert lhs_expr.startswith("d2_t"), (
            f"h_5 should be dynamical (d2_t) not a constraint: lhs={lhs_expr}. "
            "If this fails, the LagrangianFullExpand fix for #250 has regressed."
        )

    def test_h7_is_dynamical(self) -> None:
        """h_7 must have a wave equation (d2_t LHS)."""
        spec = _load_spec("gertsenshtein")
        eq = _equation_for(spec, "h_7")
        lhs_expr = eq["lhs"]["expression"]
        assert lhs_expr.startswith("d2_t"), (
            f"h_7 should be dynamical (d2_t): lhs={lhs_expr}"
        )

    def test_h5_has_kinetic_term(self) -> None:
        """h_5 EOM must contain a laplacian_x(h_5) kinetic term."""
        spec = _load_spec("gertsenshtein")
        eq = _equation_for(spec, "h_5")
        rhs_terms = eq.get("rhs", {}).get("terms", [])
        has_laplacian = any(
            t.get("operator") == "laplacian_x" and t.get("field") == "h_5"
            for t in rhs_terms
        )
        assert has_laplacian, (
            "h_5 EOM missing laplacian_x(h_5) kinetic term. "
            "This indicates the canonical pipeline lost the non-trace "
            "graviton kinetic structure (#250 regression)."
        )

    def test_h5_has_mass_term(self) -> None:
        """h_5 EOM must contain a B0^2/2 identity mass term (from background EM)."""
        spec = _load_spec("gertsenshtein")
        eq = _equation_for(spec, "h_5")
        rhs_terms = eq.get("rhs", {}).get("terms", [])
        has_mass = any(
            t.get("operator") == "identity"
            and t.get("field") == "h_5"
            and "B0" in (t.get("coefficient_symbolic") or "")
            for t in rhs_terms
        )
        assert has_mass, "h_5 EOM missing B0^2 mass term (background EM coupling)"

    def test_h5_couples_to_a1(self) -> None:
        """h_5 EOM must contain a B0 * gradient_x(a_1) coupling term."""
        spec = _load_spec("gertsenshtein")
        eq = _equation_for(spec, "h_5")
        rhs_terms = eq.get("rhs", {}).get("terms", [])
        has_coupling = any(
            t.get("operator") == "gradient_x"
            and t.get("field") == "a_1"
            and "B0" in (t.get("coefficient_symbolic") or "")
            for t in rhs_terms
        )
        assert has_coupling, (
            "h_5 EOM missing B0 * gradient_x(a_1) Gertsenshtein coupling"
        )

    def test_photons_all_dynamical(self) -> None:
        """All photon fields a_0..a_3 must have time_order=2 (#218 protection).

        If a photon becomes a constraint, it indicates ToCanonical is merging
        Maxwell contractions incorrectly — the original #218 failure mode.
        """
        spec = _load_spec("gertsenshtein")
        for i in range(4):
            eq = _equation_for(spec, f"a_{i}")
            lhs_expr = eq["lhs"]["expression"]
            assert lhs_expr.startswith("d2_t"), (
                f"a_{i} should be dynamical (d2_t), got {lhs_expr}. "
                "This is the #218 regression: photon EOM has been "
                "constrained by ToCanonical merging Maxwell contractions."
            )


class TestTorsionDarkPhotonRegression:
    """Regression tests for torsion_dark_photon: tests #218, #250, #255 jointly.

    This theory has BOTH matter-deferred (F = dA) AND torsion-deferred
    (Ftorsion = dT) fields. It exercises:
      - #218: photon EOM correctness (a_* dynamical)
      - #250: h_5 kinetic structure (h_5 dynamical)
      - #255: no coefficient corruption from mixed abstract tensors
    """

    def test_photons_dynamical(self) -> None:
        """All photon a_0..a_3 dynamical (#218)."""
        spec = _load_spec("torsion_dark_photon")
        for i in range(4):
            eq = _equation_for(spec, f"a_{i}")
            assert eq["lhs"]["expression"].startswith("d2_t"), (
                f"a_{i} not dynamical: {eq['lhs']['expression']}"
            )

    def test_gravitons_dynamical(self) -> None:
        """h_5 and h_7 dynamical (#250)."""
        spec = _load_spec("torsion_dark_photon")
        for fname in ("h_5", "h_7"):
            eq = _equation_for(spec, fname)
            assert eq["lhs"]["expression"].startswith("d2_t"), (
                f"{fname} not dynamical: {eq['lhs']['expression']}"
            )

    def test_torsion_propagating_modes_dynamical(self) -> None:
        """The propagating torsion modes must be dynamical.

        The remaining torsion components (t_6, t_13, t_20) are
        gauge-eliminated constraints — this is a physical property of the
        propagating dark photon model. If a propagating mode becomes a
        constraint, it indicates that the canonical pipeline lost the
        torsion kinetic structure (#255 / torsion-deferred regression in
        the LagrangianFullExpand handling).
        """
        spec = _load_spec("torsion_dark_photon")
        # Modes expected to be dynamical (verified against committed JSON)
        propagating_modes = (
            "t_0",
            "t_1",
            "t_2",
            "t_9",
            "t_10",
            "t_15",
            "t_17",
            "t_22",
            "t_23",
        )
        for fname in propagating_modes:
            eq = _equation_for(spec, fname)
            assert eq["lhs"]["expression"].startswith("d2_t"), (
                f"{fname} expected dynamical but got: {eq['lhs']['expression']}"
            )
