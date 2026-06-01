"""Tests for tidal.symbolic.kinetic_matrix — kinetic-matrix assembly (issue #372)."""

from __future__ import annotations

from pathlib import Path

import pytest

from tidal.symbolic.json_loader import load_equation_system
from tidal.symbolic.kinetic_matrix import (
    _negate_coefficient,
    build_kinetic_matrix,
)

_EXAMPLES = Path(__file__).resolve().parent.parent / "examples" / "data"
_BASELINE = _EXAMPLES / "gertsenshtein.json"
_needs_baseline = pytest.mark.skipif(
    not _BASELINE.exists(),
    reason="examples/data/gertsenshtein.json not derived",
)


class TestNegateCoefficient:
    def test_positive_to_negative(self) -> None:
        assert _negate_coefficient("B0/2") == "-B0/2"

    def test_negative_to_positive(self) -> None:
        assert _negate_coefficient("-kappa^(-2)") == "kappa^(-2)"

    def test_one_to_minus_one(self) -> None:
        assert _negate_coefficient("1") == "-1"

    def test_minus_one_to_one(self) -> None:
        assert _negate_coefficient("-1") == "1"

    def test_sum_wrapped_in_parens(self) -> None:
        # x+y -> -(x+y) so the leading sign applies to the whole sum.
        assert _negate_coefficient("a+b") == "-(a+b)"


@_needs_baseline
class TestBaselineKineticMatrix:
    """Acceptance: Einstein-Maxwell baseline produces the canonical structure."""

    def setup_method(self) -> None:
        self.spec = load_equation_system(_BASELINE)
        self.km = build_kinetic_matrix(self.spec)

    def test_shape(self) -> None:
        # 6 dynamical fields in the constraint-eliminated baseline.
        assert self.km.n == len(self.km.fields)
        assert all(len(row) == self.km.n for row in self.km.cells)

    def test_gertsenshtein_off_diagonal_h_xy_A_x(self) -> None:
        """K[h_xy, A_x] and its transpose carry the -B_0 ∂_z Gertsenshtein coupling."""
        i = self.km.field_index("h_5")  # h_xy
        j = self.km.field_index("a_1")  # A_x
        assert i is not None
        assert j is not None

        # K[h_xy, A_x]: should have a single (gradient_x, -B0) entry
        # (gradient_x is the propagation-axis gradient before axis remap).
        cell = self.km.cells[i][j]
        assert len(cell.entries) == 1
        op, coeff = cell.entries[0]
        assert op == "gradient_x"
        assert coeff == "-B0"

        # K[A_x, h_xy]: symmetric coupling, same entry.
        cell_T = self.km.cells[j][i]
        assert cell_T.entries == ((op, coeff),)

    def test_orthogonal_polarization_decoupled(self) -> None:
        """K[h_xy, h_yy] and K[A_x, A_y] are empty — the two polarizations
        evolve in independent channels at the kinetic level.
        """
        i_hxy = self.km.field_index("h_5")
        i_hyy = self.km.field_index("h_7")
        i_ax = self.km.field_index("a_1")
        i_ay = self.km.field_index("a_2")
        assert i_hxy is not None
        assert i_hyy is not None
        assert i_ax is not None
        assert i_ay is not None
        assert self.km.cells[i_hxy][i_hyy].is_zero()
        assert self.km.cells[i_ax][i_ay].is_zero()

    def test_diagonal_h_xy_has_expected_structure(self) -> None:
        """K[h_xy, h_xy] carries the LHS kinetic prefactor (-1/κ²)∂_t²,
        the propagation-axis Laplacian (+1/κ²)∂_z², and the
        magnetic mass (-B₀²/2)·1.
        """
        i = self.km.field_index("h_5")
        assert i is not None
        entries = self.km.cells[i][i].entries
        ops_present = {op for op, _ in entries}
        assert "d2_t" in ops_present
        assert (
            "laplacian_x" in ops_present
        )  # propagation-axis (renamed to z at render time)
        assert "identity" in ops_present


@_needs_baseline
class TestKineticMatrixToLatex:
    """The rendered LaTeX uses a labelled array and includes the expected entries."""

    def test_baseline_render(self) -> None:
        import tidal.symbolic.latex as latex_mod
        from tidal.symbolic.latex import kinetic_matrix_to_latex, system_to_latex

        spec = load_equation_system(_BASELINE)
        # system_to_latex mutates the module-level metric symbol as a
        # side effect; save and restore so this test does not leak
        # state into subsequent tests that assert flat-metric \eta.
        saved_metric = latex_mod._metric_symbol
        try:
            _ = system_to_latex(spec, output_format="raw")
            km = build_kinetic_matrix(spec)
            out = kinetic_matrix_to_latex(km, spec)
        finally:
            latex_mod._metric_symbol = saved_metric

        assert r"\mathcal{K}" in out
        # Labelled array wrapper (replaces the earlier bare bmatrix).
        assert r"\begin{array}" in out
        assert r"\end{array}" in out
        # Gertsenshtein off-diagonal: -B_0 ∂_z (after axis remap to z).
        assert r"-B_{0} \, \partial_z" in out
        # Diagonal photon kinetic: ∂_t^2 - ∂_z^2 entries present.
        assert r"\partial_t^2 - \partial_z^2" in out
        # Field labels: \mathcal{A}_{t} appears as both row and column
        # label for the baseline (calligraphic at render time; the
        # appendix-side sed re-letters to italic A in the published
        # listing).
        assert r"\mathcal{A}_{t}" in out


_NONMINIMAL = _EXAMPLES / "torsion_gertsenshtein_nonminimal.json"
_needs_nonminimal = pytest.mark.skipif(
    not _NONMINIMAL.exists(),
    reason="examples/data/torsion_gertsenshtein_nonminimal.json not derived",
)


@_needs_nonminimal
class TestVelocityPairColumns:
    """Acceptance for the rectangular extension: v-fields show up as columns
    so first-time-derivative couplings (gradient_z(v_X), etc.) are not
    silently dropped from the displayed matrix.
    """

    def setup_method(self) -> None:
        self.spec = load_equation_system(_NONMINIMAL)
        self.km = build_kinetic_matrix(self.spec)

    def test_matrix_is_rectangular(self) -> None:
        """The nonminimal theory uses velocity-pair fields on the RHS;
        the matrix grows in columns (not rows).
        """
        assert self.km.n_cols > self.km.n_rows
        assert not self.km.is_square

    def test_all_rhs_v_fields_have_columns(self) -> None:
        """Every v_<X> referenced on any RHS now has a dedicated
        column — no skipped terms.
        """
        v_refs = {
            term.field
            for eq in self.spec.equations
            for term in eq.rhs_terms
            if term.field.startswith("v_")
        }
        column_set = set(self.km.column_fields)
        assert v_refs.issubset(column_set), f"missing v-columns: {v_refs - column_set}"

    def test_v_columns_appear_after_q_columns(self) -> None:
        """Convention: column ordering is q-fields first (one per
        equation), then v-fields ordered by first RHS appearance.
        """
        n_q = self.km.n_rows
        head = self.km.column_fields[:n_q]
        tail = self.km.column_fields[n_q:]
        # q head matches row_fields exactly.
        assert head == self.km.row_fields
        # All trailing columns are v-fields.
        assert all(name.startswith("v_") for name in tail)
