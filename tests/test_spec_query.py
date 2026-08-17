"""Tests for :mod:`tidal.symbolic.spec_query`.

Every test goes through the public accessors rather than raw JSON, so the
suite doubles as usage documentation — and an awkward API shows up here first.

The corpus-level tests matter most.  All six #401 misreadings were claims
*about the whole corpus* ("28 JSONs show the signature"), which no single-file
unit test can catch.  #394 is the cautionary precedent: its test asserted
against a hand-written dict that already contained the term under test, so it
passed for months while the pipeline emitted nothing.  The tests here read the
committed specs through the loader instead.
"""

from __future__ import annotations

import json
import subprocess
import unittest
import warnings
from pathlib import Path

import pytest

from tidal.symbolic.json_loader import EquationSystem, load_equation_system
from tidal.symbolic.sign_algebra import Sign
from tidal.symbolic.spec_query import (
    coefficient_provenance,
    compare_equations,
    diff_systems,
    effective_coefficient,
    field_families,
    matrix_encoding_agrees,
    self_terms,
    sibling_sign_conflicts,
    terms_for,
)

REPO = Path(__file__).resolve().parent.parent
EXAMPLES = REPO / "examples" / "data"

# The 19 committed specs whose photon components carry the #397 opposite-sign
# signature. Pinned as a set so a re-derivation that fixes one shows up here.
EXPECTED_CONFLICT_FILES = 19


def _load(name: str) -> EquationSystem:
    """Load a committed example spec.

    ``strict_v6=False`` matches what ``tidal inspect`` already does for
    read-only paths: analysis never evolves the system, so the higher-time-order
    guard would only block inspection of specs that are perfectly readable.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return load_equation_system(EXAMPLES / f"{name}.json", strict_v6=False)


class TestEffectiveCoefficient(unittest.TestCase):
    """Summing all matching terms and dividing out the kinetic coefficient."""

    def test_all_matching_terms_are_summed(self) -> None:
        """`gertsenshtein_eh` a_0 has two self-laplacians, not one (#401 row 3)."""
        spec = _load("gertsenshtein_eh")
        equation = spec.equations[spec.equation_map["a_0"]]
        eff = effective_coefficient(equation, "a_0", "laplacian_x")
        self.assertEqual(eff.term_count, 2)
        self.assertIn("-2*B0^2*rho", eff.expression)

    def test_kinetic_coefficient_is_divided_out(self) -> None:
        """The LHS divisor is part of the coefficient, not decoration."""
        spec = _load("gertsenshtein_ungauged")
        equation = spec.equations[spec.equation_map["h_5"]]
        eff = effective_coefficient(equation, "h_5", "laplacian_x")
        self.assertEqual(eff.kinetic, "-kappa^(-2)")
        self.assertEqual(eff.expression, "(-kappa^(-2))/(-kappa^(-2))")
        self.assertIs(eff.sign().sign, Sign.POSITIVE)

    def test_absent_kinetic_coefficient_means_one(self) -> None:
        """A bare LHS leaves the summed numerator unchanged."""
        spec = _load("gertsenshtein_ungauged")
        equation = spec.equations[spec.equation_map["a_0"]]
        eff = effective_coefficient(equation, "a_0", "laplacian_x")
        self.assertIsNone(eff.kinetic)
        self.assertEqual(eff.expression, eff.numerator)

    def test_missing_term_reports_absent(self) -> None:
        """A key with no matching term is empty rather than zero-by-accident."""
        spec = _load("gertsenshtein_ungauged")
        equation = spec.equations[spec.equation_map["a_0"]]
        eff = effective_coefficient(equation, "a_0", "biharmonic")
        self.assertFalse(eff.exists)
        self.assertEqual(eff.term_count, 0)

    def test_terms_for_and_self_terms_agree_on_self_keys(self) -> None:
        """``self_terms`` is ``terms_for`` specialised to the equation's own field."""
        spec = _load("gertsenshtein_eh")
        equation = spec.equations[spec.equation_map["a_0"]]
        self.assertEqual(
            terms_for(equation, "a_0", "laplacian_x"),
            self_terms(equation, "laplacian_x"),
        )

    def test_effective_mass_divides_by_the_kinetic_coefficient(self) -> None:
        """The defect made five times over, pinned at the accessor level.

        Ignoring ``kinetic_coefficient_symbolic`` has been independently
        re-implemented wrong by the modal solver (#258), the non-modal solvers
        (#302), the ``--require-stable`` mass check (#237, closed with the
        ``--no-require-stable`` workaround rather than a fix), and twice in the
        #401 analysis scans.

        #237's case: propagating torsion has a kinetic coefficient of ``-xi``,
        so the effective mass² is ``pot[i][j] / xi``, not ``pot[i][j]``.  The
        accessor must carry that divisor, so #395's validators inherit the
        right quantity instead of repeating the mistake a sixth time.
        """
        spec = _load("torsion_gertsenshtein_propagating")
        equation = next(
            eq
            for eq in spec.equations
            if eq.field_name.startswith("t_") and eq.time_derivative_order >= 2
        )
        eff = effective_coefficient(equation, equation.field_name, "identity")
        self.assertEqual(eff.kinetic, "-xi")
        self.assertIn("/(-xi)", eff.expression)
        # Several contributions, so a single-term read would also be wrong.
        self.assertGreater(eff.term_count, 1)

    def test_symbolic_coefficients_are_not_skipped(self) -> None:
        """`torsion_dark_photon_fv` t_1 is symbolic throughout (#401 row 2)."""
        spec = _load("torsion_dark_photon_fv")
        equation = spec.equations[spec.equation_map["t_1"]]
        eff = effective_coefficient(equation, "t_1", "laplacian_x")
        self.assertEqual(eff.kinetic, "-xi")
        # -xi / -xi is exactly +1 without any value for xi.
        self.assertIs(eff.sign().sign, Sign.POSITIVE)


class TestFieldFamilies(unittest.TestCase):
    """Grouping by tensor metadata rather than by parsing the name suffix."""

    def test_families_use_tensor_head(self) -> None:
        """Components group by head, and the grouping is marked exact."""
        spec = _load("torsion_dark_photon")
        families = {f.head: f for f in field_families(spec)}
        self.assertEqual(set(families), {"a", "h", "t"})
        self.assertTrue(all(f.exact for f in families.values()))
        self.assertEqual(families["a"].rank, 1)
        self.assertEqual(families["t"].rank, 3)

    def test_rank3_temporal_slots_replace_suffix_reading(self) -> None:
        """`t_13` is ``[2, 0, 2]`` — 'component 0 is temporal' is meaningless (#401 row 4).

        The flat ``t_0..t_23`` numbering interleaves the classes arbitrarily,
        so only the index tuple can tell them apart.
        """
        spec = _load("torsion_dark_photon")
        family = next(f for f in field_families(spec) if f.head == "t")
        self.assertEqual(family.temporal_slots("t_0"), 2)  # [0, 0, 1]
        self.assertEqual(family.temporal_slots("t_13"), 1)  # [2, 0, 2]
        self.assertEqual(family.temporal_slots("t_15"), 0)  # [2, 1, 2]
        groups = family.group_by_temporal_slots()
        self.assertEqual(set(groups), {0, 1, 2})

    def test_rank1_temporal_slots_recover_the_familiar_case(self) -> None:
        """For a rank-1 field the classification reduces to 'a_0 is temporal'."""
        spec = _load("torsion_dark_photon")
        family = next(f for f in field_families(spec) if f.head == "a")
        self.assertEqual(family.group_by_temporal_slots()[1], ("a_0",))
        self.assertEqual(family.group_by_temporal_slots()[0], ("a_1", "a_2", "a_3"))

    def test_specs_without_metadata_are_marked_inexact(self) -> None:
        """Older specs fall back to name splitting, and say so.

        12 committed specs predate ``78374c1``; a caller must be able to tell a
        guessed grouping from a derived one.
        """
        spec = _load("torsion_gertsenshtein")
        families = field_families(spec)
        self.assertTrue(families)
        self.assertFalse(any(f.exact for f in families))


class TestSiblingSignConflicts(unittest.TestCase):
    """The #397 invariant, and the partition choice it depends on."""

    def test_detects_the_known_defect_across_the_corpus(self) -> None:
        """Exactly the known-affected specs are flagged, corpus-wide."""
        flagged = []
        for path in sorted(EXAMPLES.glob("*.json")):
            spec = _load(path.stem)
            if sibling_sign_conflicts(spec):
                flagged.append(path.name)
        self.assertEqual(len(flagged), EXPECTED_CONFLICT_FILES)

    def test_euler_heisenberg_is_not_flagged(self) -> None:
        """The #401 row-1 false positive must not reappear.

        `gertsenshtein_eh` looks inconsistent only if the kinetic coefficient
        is ignored; with it divided out the components agree.
        """
        for name in ("gertsenshtein_eh", "gertsenshtein_eh_top"):
            with self.subTest(spec=name):
                self.assertEqual(sibling_sign_conflicts(_load(name)), ())

    def test_constraints_are_excluded(self) -> None:
        """A constraint's overall sign is conventional, so it is never compared.

        In `torsion_dark_photon` the one-temporal-slot torsion components are
        all ``time_order=0`` constraints; comparing them against the dynamical
        components would manufacture conflicts.
        """
        spec = _load("torsion_dark_photon")
        constrained = [
            eq.field_name
            for eq in spec.equations
            if eq.time_derivative_order == 0 and eq.field_name.startswith("t_")
        ]
        self.assertTrue(constrained, "fixture should contain torsion constraints")
        conflicts = sibling_sign_conflicts(spec)
        reported = {name for _, left, right in conflicts for name in (left, right)}
        self.assertEqual(reported & set(constrained), set())

    def test_comparison_spans_the_whole_family(self) -> None:
        """Temporal-slot grouping must not partition the comparison.

        The defect is a temporal component disagreeing with its *spatial*
        siblings, so restricting comparisons to like-indexed components would
        never look at that pair. Measured on this corpus, such a partition
        misses 13 of the 19 affected files — this test pins that the photon
        family is compared as a whole.
        """
        spec = _load("torsion_gertsenshtein_complete_even")
        conflicts = sibling_sign_conflicts(spec)
        pairs = {(left, right) for _, left, right in conflicts}
        self.assertIn(
            ("a_0", "a_1"),
            pairs,
            "a_0 (1 temporal slot) must be compared against a_1 (0 slots)",
        )


class TestCoefficientProvenance(unittest.TestCase):
    """Where a coefficient is recorded, and how those places relate."""

    def test_reports_parts_duplicates_and_related_quantities(self) -> None:
        """One `h_5` mass term appears in three places with three conventions."""
        spec = _load("gertsenshtein_ungauged")
        prov = coefficient_provenance(spec, "h_5", "h_5", "identity")

        # (a) parts: the RHS term plus the LHS kinetic coefficient
        self.assertEqual(prov.effective.numerator, "B0^2/2")
        self.assertEqual(prov.effective.kinetic, "-kappa^(-2)")

        # (b) duplicate encoding: the mass matrix, un-normalised, negated
        self.assertIsNotNone(prov.matrix_entry)

        # (c) related but distinct: the Hamiltonian carries B0^2/4, not B0^2/2
        symbolics = {t.coefficient_symbolic for t in prov.hamiltonian_terms}
        self.assertIn("B0^2/4", symbolics)

    def test_numeric_and_symbolic_encodings_agree(self) -> None:
        """The stored numeric equals the symbolic evaluated at all-ones."""
        spec = _load("gertsenshtein_ungauged")
        prov = coefficient_provenance(spec, "h_5", "h_5", "identity")
        statuses = {c.status for c in prov.checks}
        self.assertNotIn("mismatch", statuses)

    def test_unknown_equation_raises(self) -> None:
        """A mistyped component fails cleanly rather than guessing."""
        spec = _load("gertsenshtein_ungauged")
        with pytest.raises(KeyError):
            coefficient_provenance(spec, "not_a_field", "h_5", "identity")

    def test_matrix_encoding_agrees_corpus_wide(self) -> None:
        """The redundant matrix encoding matches the summed identity terms everywhere.

        This replaces a characterization test that pinned 26 specs as failing
        while GH #403 was open: ``_compute_matrices_from_terms`` accumulated
        into the numeric matrix (``+=``) but overwrote the symbolic one (``=``),
        so a component with several ``identity`` self-terms kept only the last.
        Both now accumulate, and the invariant holds across the corpus.

        Nothing else enforces this redundancy — the loader recomputes both
        matrices and ignores the values stored in the JSON, and #274 records
        them drifting apart once already — so this is the enforcement.
        """
        for path in sorted(EXAMPLES.glob("*.json")):
            with self.subTest(spec=path.name):
                self.assertEqual(matrix_encoding_agrees(_load(path.stem)), ())

    def test_multi_term_identity_is_summed_not_truncated(self) -> None:
        """A component with several identity terms keeps all of them (GH #403).

        ``coupled_scalars`` ``h_0`` carries both ``B0^2`` and ``mg2/kappa^2``;
        before the fix the symbolic matrix held only the second.
        """
        spec = _load("coupled_scalars")
        index = spec.equation_map["h_0"]
        entry = spec.mass_matrix_symbolic[index][index]
        self.assertIsNotNone(entry)
        assert entry is not None
        for part in ("B0^2", "mg2/kappa^2"):
            self.assertIn(part, entry)


class TestCorpusGolden(unittest.TestCase):
    """The committed corpus summary must stay in step with the accessors."""

    def test_committed_report_is_up_to_date(self) -> None:
        """``tests/data/spec_semantics.txt`` matches a fresh regeneration.

        The report is the committed answer to corpus-level questions, so that
        answering one means reading a file rather than writing another
        throwaway scan — which is how all six #401 misreadings happened.  Its
        value depends entirely on being current, so drift is a test failure.

        Regenerate with ``python -m scripts.spec_semantics_report``.
        """
        from scripts.spec_semantics_report import GOLDEN, build_report

        self.assertTrue(GOLDEN.exists(), f"{GOLDEN} is missing; regenerate it")
        self.assertEqual(
            GOLDEN.read_text(),
            build_report(),
            "spec_semantics.txt is stale; run "
            "`python -m scripts.spec_semantics_report` and review the diff",
        )

    def test_report_records_the_known_conflict_count(self) -> None:
        """The summary line agrees with the accessor, guarding against a silent drop."""
        from scripts.spec_semantics_report import GOLDEN

        text = GOLDEN.read_text()
        self.assertIn(
            f"specs with a proven sign conflict: {EXPECTED_CONFLICT_FILES}",
            text,
        )


class TestDiffing(unittest.TestCase):
    """Separating a genuine change from a rewritten one."""

    @staticmethod
    def _spec_at_revision(revision: str, relative: str) -> EquationSystem:
        """Load a spec as it was at a git revision."""
        blob = subprocess.run(
            ["git", "show", f"{revision}:{relative}"],  # noqa: S607
            capture_output=True,
            text=True,
            check=True,
            cwd=REPO,
        ).stdout
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return EquationSystem.from_dict(json.loads(blob))

    def test_rescaling_is_representational_and_rhs_flip_is_real(self) -> None:
        """The `gertsenshtein_ungauged` re-derivation, classified correctly.

        ``h_5``/``h_6``/``h_8`` were multiplied through by ``-kappa^(-2)`` on
        both sides — physically unchanged, but read as three separate "fixes"
        by a naive diff (#401 row 5).  ``a_0`` flipped its RHS with the kinetic
        coefficient untouched, which *is* the real #397 defect.

        Both directions are asserted: a tool that called everything
        representational would pass a test of invariance alone.
        """
        relative = "examples/data/gertsenshtein_ungauged.json"
        old = self._spec_at_revision("7b779288^", relative)
        new = self._spec_at_revision("7b779288", relative)

        diff = diff_systems(old, new)
        verdicts = {c.field: c.verdict for c in diff.comparisons}

        for field in ("h_5", "h_6", "h_8"):
            with self.subTest(field=field):
                self.assertEqual(
                    verdicts[field],
                    "representational",
                    f"{field} was rescaled on both sides, not fixed",
                )
        self.assertEqual(verdicts["a_0"], "real")
        self.assertTrue(diff.has_real_changes)

    def test_identical_systems_show_no_changes(self) -> None:
        """A spec compared with itself yields no real or representational change."""
        spec = _load("gertsenshtein_ungauged")
        diff = diff_systems(spec, spec)
        self.assertFalse(diff.has_real_changes)
        self.assertEqual(diff.representational, ())
        self.assertTrue(all(c.verdict == "identical" for c in diff.comparisons))

    def test_comparison_names_the_changed_keys(self) -> None:
        """A real change reports which term keys moved."""
        relative = "examples/data/gertsenshtein_ungauged.json"
        old = self._spec_at_revision("7b779288^", relative)
        new = self._spec_at_revision("7b779288", relative)
        comparison = compare_equations(
            old.equations[old.equation_map["a_0"]],
            new.equations[new.equation_map["a_0"]],
        )
        self.assertEqual(comparison.verdict, "real")
        self.assertIn("laplacian_x(a_0)", comparison.changed_keys)


if __name__ == "__main__":
    unittest.main()
