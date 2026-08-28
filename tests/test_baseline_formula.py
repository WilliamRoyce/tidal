"""Baseline formulas are refused at launch, not swallowed (GH #407).

`--baseline-formula` referencing spec constants that were never passed via
`--param` used to fail per evaluation, be swallowed to ``None``, become
``-inf``, and land on the v3 soft floor -- so an inference chain ran to
completion and produced a posterior built **entirely from floor values**.

Only 6 of 48 committed specs carry a ``metadata.parameters`` block, and
none is a PGT or dark-photon theory, so for every campaign spec ``--param``
is the only channel that supplies ``kappa``/``B0`` at all. #270's
prescribed remedy -- raise if the formula references an unresolvable name
-- was never implemented until now.
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from tidal.cli._simulate import (  # pyright: ignore[reportPrivateUsage]
    FORMULA_NAMESPACE,
    unresolved_formula_names,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
SPEC = REPO_ROOT / "examples/data/gertsenshtein.json"
CAMPAIGN_FORMULA = "sin(kappa * B0 * t_end / 2)**2"


class TestUnresolvedFormulaNames:
    """The shared check itself."""

    def test_reports_every_missing_name_at_once(self) -> None:
        """Not just the first: fixing one and re-running is a poor loop.

        ``_validate_formula_ast`` raises on the first offender, which is
        right for evaluation but wrong for a launch-time report.
        """
        missing = unresolved_formula_names(
            CAMPAIGN_FORMULA, set(FORMULA_NAMESPACE) | {"t_end"}
        )
        assert missing == {"kappa", "B0"}

    def test_satisfiable_formula_reports_nothing(self) -> None:
        available = set(FORMULA_NAMESPACE) | {"t_end", "kappa", "B0"}
        assert unresolved_formula_names(CAMPAIGN_FORMULA, available) == set()

    def test_math_functions_are_not_missing_names(self) -> None:
        assert unresolved_formula_names(
            "sqrt(pi) * exp(1)", set(FORMULA_NAMESPACE)
        ) == (set())

    def test_syntax_error_is_not_this_check_s_business(self) -> None:
        """A malformed formula is ``_validate_formula_ast``'s error to report."""
        assert unresolved_formula_names("sin(", set(FORMULA_NAMESPACE)) == set()


class TestSampleRefusesAtLaunch:
    """`tidal sample` must not start a chain it cannot score."""

    @staticmethod
    def _args(tmp_path: Path, *extra: str) -> object:
        from tidal.cli import _build_parser

        return _build_parser().parse_args(
            [
                "sample",
                str(SPEC),
                "--prior",
                "B0=uniform:0.001:0.1",
                "--likelihood",
                "P_max:maximize",
                "--method",
                "mc",
                "--n-samples",
                "1",
                "--output",
                str(tmp_path),
                "--baseline-formula",
                CAMPAIGN_FORMULA,
                *extra,
            ]
        )

    def test_missing_constant_is_refused_before_any_simulation(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        from tidal.cli._sample import sample_command

        rc = sample_command(self._args(tmp_path))
        # error_with_hint writes to stderr.
        err = capsys.readouterr().err

        assert rc == 1
        assert "kappa" in err
        assert "#407" in err
        # Nothing ran: no results were written.
        assert not (tmp_path / "results.csv").exists()

    def test_a_sampled_name_is_not_missing(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """The failure mode this check risks is over-rejection.

        ``B0`` is supplied by ``--prior``, so it resolves per evaluation
        and must not be reported. Rejecting it would break the campaign's
        normal usage, where the formula references a swept coupling.
        """
        from tidal.cli._sample import sample_command

        rc = sample_command(self._args(tmp_path, "--param", "kappa=1.0"))
        err = capsys.readouterr().err

        assert "baseline-formula references" not in err, (
            "B0 is supplied by --prior and must not be reported missing"
        )
        assert rc == 0, "a satisfiable formula must not be refused"

    def test_no_formula_no_check(self, tmp_path: Path) -> None:
        from tidal.cli import _build_parser
        from tidal.cli._sample import sample_command

        args = _build_parser().parse_args(
            [
                "sample",
                str(SPEC),
                "--prior",
                "B0=uniform:0.001:0.1",
                "--likelihood",
                "P_max:maximize",
                "--method",
                "mc",
                "--n-samples",
                "1",
                "--output",
                str(tmp_path),
                "--param",
                "kappa=1.0",
            ]
        )
        assert sample_command(args) == 0


class TestEvalBaselineNoLongerConflates:
    """A config error and a zero baseline are different things."""

    def test_unresolvable_name_raises(self) -> None:
        """It used to return None, which the caller read as 'baseline zero'.

        Both then produced ``-inf`` -> soft floor, so a missing ``--param``
        was indistinguishable from physics.
        """
        from tidal.inference._likelihood import (
            LikelihoodConfig,
            compute_log_likelihood,
        )

        cfg = LikelihoodConfig(
            metric="P_max",
            likelihood_type="maximize",
            baseline_formula=CAMPAIGN_FORMULA,
        )
        with pytest.raises(ValueError, match="kappa"):
            compute_log_likelihood(0.5, cfg, {"t_end": 10.0})

    def test_non_positive_baseline_still_soft_floors(self) -> None:
        """The legitimate half must be untouched: it is a physics outcome."""
        from tidal.inference._likelihood import (
            LikelihoodConfig,
            compute_log_likelihood,
        )

        cfg = LikelihoodConfig(
            metric="P_max",
            likelihood_type="maximize",
            baseline_formula="0.0",
        )
        assert compute_log_likelihood(0.5, cfg, {}) == -math.inf

    def test_satisfiable_formula_scores_normally(self) -> None:
        from tidal.inference._likelihood import (
            LikelihoodConfig,
            compute_log_likelihood,
        )

        cfg = LikelihoodConfig(
            metric="P_max",
            likelihood_type="maximize",
            baseline_formula=CAMPAIGN_FORMULA,
        )
        logl = compute_log_likelihood(
            0.5, cfg, {"t_end": 10.0, "kappa": 1.0, "B0": 0.01}
        )
        assert math.isfinite(logl)
        assert logl > 0.0


class TestSpecDefaultsAreInertForCampaignSpecs:
    """Why the launch check is needed at all (#270's fix does not cover us)."""

    def test_no_campaign_spec_carries_parameter_defaults(self) -> None:
        import json

        with_params = [
            p.name
            for p in sorted((REPO_ROOT / "examples/data").glob("*.json"))
            if (json.loads(p.read_text()).get("metadata") or {}).get("parameters")
        ]
        # The six that do are toy/analytic specs; none is PGT or dark-photon,
        # so --param is the only channel for kappa/B0 on every campaign run.
        assert "gertsenshtein.json" not in with_params
        assert "dark_photon_plasma.json" not in with_params
        assert len(with_params) < 10, (
            "if many specs now carry defaults, revisit the #407 premise"
        )
