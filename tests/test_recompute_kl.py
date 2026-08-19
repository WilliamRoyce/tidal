"""Unit tests for scripts/analysis/recompute_parameter_kl.py helpers.

The two defect classes PR #431 found by eye — column-rename prior
misassignment and unapplied/malformed overrides — live in pure functions
that previously had zero pytest coverage (post-merge review, prevention
gap 5).  The script is loaded by file path (scripts/ is not a package).
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pytest

if TYPE_CHECKING:
    from types import ModuleType

_SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "analysis"
    / "recompute_parameter_kl.py"
)


@pytest.fixture(scope="module")
def rk() -> ModuleType:
    spec = importlib.util.spec_from_file_location("recompute_parameter_kl", _SCRIPT)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    # Register before exec: the script's @dataclass needs its module
    # resolvable through sys.modules for annotation lookup.
    sys.modules["recompute_parameter_kl"] = module
    spec.loader.exec_module(module)
    return module


def _scalar(name: str, dist: str = "arctan_uniform") -> dict[str, Any]:
    return {
        "kind": "scalar",
        "name": name,
        "distribution": dist,
        "low": -89.0,
        "high": 89.0,
    }


class TestAlignPriorsToParams:
    def test_renamed_columns_matched_positionally(self, rk: ModuleType) -> None:
        """The ricci_em case: TOML relabels alpha* -> beta*; records must
        follow position, not name.
        """
        priors = [
            _scalar("alpha1"),
            _scalar("alpha2"),
            _scalar("alpha3", "log_uniform"),
            _scalar("delta1"),
        ]
        aligned = rk._align_priors_to_params(
            priors, ["beta1", "beta2", "beta3", "delta1"], "test"
        )
        assert [r["name"] for r in aligned] == ["beta1", "beta2", "beta3", "delta1"]
        # Kinds travel with position, not with the old name.
        assert aligned[2]["distribution"] == "log_uniform"

    def test_identity_alignment_is_noop(self, rk: ModuleType) -> None:
        priors = [_scalar("a"), _scalar("b")]
        aligned = rk._align_priors_to_params(priors, ["a", "b"], "test")
        assert aligned == priors

    def test_radial_angular_consumes_block(self, rk: ModuleType) -> None:
        joint = {
            "kind": "radial_angular",
            "names": ["c1", "c2", "c3"],
            "r_lo": 1e-3,
            "r_hi": 1e3,
        }
        aligned = rk._align_priors_to_params(
            [joint, _scalar("delta1")], ["x1", "x2", "x3", "delta1"], "test"
        )
        assert aligned[0]["names"] == ["x1", "x2", "x3"]
        assert aligned[1]["name"] == "delta1"

    @pytest.mark.parametrize(
        ("priors", "params"),
        [
            # too few columns for the joint record
            (
                [{"kind": "radial_angular", "names": ["c1", "c2", "c3"]}],
                ["x1", "x2"],
            ),
            # more scalar records than columns
            ([_scalar("a"), _scalar("b")], ["a"]),
            # fewer records than columns
            ([_scalar("a")], ["a", "b"]),
        ],
    )
    def test_arity_mismatch_raises(
        self, rk: ModuleType, priors: list[dict[str, Any]], params: list[str]
    ) -> None:
        """Misalignment must refuse loudly — proceeding scores every column
        against the wrong (or no) prior, the silent-misassignment class
        the #420 campaign was about.
        """
        with pytest.raises(ValueError, match="misaligned priors"):
            rk._align_priors_to_params(priors, params, "test")


class TestApplyPriorOverrides:
    def test_override_replaces_record(self, rk: ModuleType) -> None:
        """The T7/T6/Barker case: fabricated arctan xi replaced by the
        sbatch-template log-uniform prior.
        """
        priors = [_scalar("beta1"), _scalar("xi")]
        out = rk._apply_prior_overrides(priors, {"xi": "log_uniform:1e-3:1e3"}, "test")
        assert out[1] == {
            "kind": "scalar",
            "name": "xi",
            "distribution": "log_uniform",
            "low": 1e-3,
            "high": 1e3,
        }
        assert out[0] == priors[0]  # untouched

    def test_malformed_override_raises_with_context(self, rk: ModuleType) -> None:
        priors = [_scalar("xi")]
        with pytest.raises(ValueError, match=r"malformed prior_overrides.*xi"):
            rk._apply_prior_overrides(priors, {"xi": "log_uniform:1e-3"}, "test")
        with pytest.raises(ValueError, match="malformed prior_overrides"):
            rk._apply_prior_overrides(priors, {"xi": "log_uniform:abc:1"}, "test")

    def test_unknown_name_raises(self, rk: ModuleType) -> None:
        """An override that matches nothing means the intended correction
        silently would not happen — refuse instead of warn.
        """
        priors = [_scalar("beta1")]
        with pytest.raises(ValueError, match="no scalar record"):
            rk._apply_prior_overrides(priors, {"xii": "log_uniform:1e-3:1e3"}, "test")


class TestReadPriorsProvenance:
    def test_recorded(self, rk: ModuleType, tmp_path: Path) -> None:
        import json

        (tmp_path / "inference.json").write_text(json.dumps({"priors": [_scalar("a")]}))
        records, provenance = rk._read_priors(tmp_path, ["a"])
        assert provenance == "recorded"
        assert records[0]["name"] == "a"

    def test_missing_priors_key_fabricates_with_warning(
        self, rk: ModuleType, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        """The Barker-amp failure mode (#434): inference.json exists but
        carries no priors key — fabrication must be loud and tagged.
        """
        import json
        import logging

        (tmp_path / "inference.json").write_text(json.dumps({"reconstructed": True}))
        with caplog.at_level(logging.WARNING, logger="recompute_kl"):
            records, provenance = rk._read_priors(tmp_path, ["xi"])
        assert provenance == "fabricated"
        assert records[0]["distribution"] == "arctan_uniform"
        assert "NO priors key" in caplog.text

    def test_no_file_fabricates_with_warning(
        self, rk: ModuleType, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        import logging

        with caplog.at_level(logging.WARNING, logger="recompute_kl"):
            _, provenance = rk._read_priors(tmp_path, ["xi"])
        assert provenance == "fabricated"
        assert "no inference.json" in caplog.text
