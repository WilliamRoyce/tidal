"""Integration tests for Phase C perturbative reduction (v6).

Verifies order-tagging integration:

* ``_WlsContext`` carries the ``perturbative_reduction`` field
* ``_validate_perturbation_config`` rejects malformed TOML sections
* ``_audit_higher_derivative_lagrangian`` warns on R^2/EH patterns
* ``generate_wls`` injects ``small_parameters`` into the Wolfram
  ``metadata`` association when ``[perturbation]`` is configured

The v5 JLM module (``PerturbativeReduction.wl``) was removed in Stage 1
of the v6 plan. Perturbative work now happens at simulate time via the
Python ``PerturbativeSolver``; derive-time is limited to emitting
``order_in_eps`` tags on each OperatorTerm in the JSON output.

End-to-end derivation tests that require running ``wolframscript`` live
under ``tests/integration/`` or are marked ``slow`` — not included here
because they take minutes per case.
"""

from __future__ import annotations

import re
import tomllib
from pathlib import Path
from typing import Any

import pytest

from tidal.cli._derive import (
    _audit_higher_derivative_lagrangian,
    _validate_perturbation_config,
    generate_wls,
)

# === TOML perturbation validator ===


class TestValidatePerturbationConfig:
    def test_none_returns_none(self) -> None:
        assert _validate_perturbation_config(None, []) is None

    def test_disabled_returns_none(self) -> None:
        cfg = {"small_parameters": ["b5"], "enabled": False}
        assert _validate_perturbation_config(cfg, ["b5"]) is None

    def test_missing_small_parameters_raises(self) -> None:
        with pytest.raises(ValueError, match="requires 'small_parameters'"):
            _validate_perturbation_config({"order": 1}, ["b5"])

    def test_empty_small_parameters_raises(self) -> None:
        with pytest.raises(ValueError, match="non-empty list"):
            _validate_perturbation_config({"small_parameters": []}, [])

    def test_non_list_small_parameters_raises(self) -> None:
        with pytest.raises(ValueError, match="non-empty list"):
            _validate_perturbation_config({"small_parameters": "b5"}, ["b5"])

    def test_non_string_entry_raises_typeerror(self) -> None:
        with pytest.raises(TypeError, match="entries must be strings"):
            _validate_perturbation_config({"small_parameters": [5]}, [])

    def test_undeclared_parameter_raises(self) -> None:
        with pytest.raises(ValueError, match="not declared in \\[constants\\]"):
            _validate_perturbation_config({"small_parameters": ["missing"]}, ["kappa"])

    def test_underscore_name_raises(self) -> None:
        with pytest.raises(ValueError, match="underscore"):
            _validate_perturbation_config({"small_parameters": ["b_5"]}, ["b_5"])

    def test_zero_order_raises(self) -> None:
        with pytest.raises(ValueError, match="positive integer"):
            _validate_perturbation_config(
                {"small_parameters": ["b5"], "order": 0}, ["b5"]
            )

    def test_negative_order_raises(self) -> None:
        with pytest.raises(ValueError, match="positive integer"):
            _validate_perturbation_config(
                {"small_parameters": ["b5"], "order": -1}, ["b5"]
            )

    def test_non_int_order_raises(self) -> None:
        with pytest.raises(ValueError, match="positive integer"):
            _validate_perturbation_config(
                {"small_parameters": ["b5"], "order": 1.5}, ["b5"]
            )

    def test_valid_single_param(self) -> None:
        result = _validate_perturbation_config(
            {"small_parameters": ["b5"]}, ["kappa", "b5"]
        )
        assert result is not None
        assert result["small_parameters"] == ["b5"]
        assert result["order"] == 1  # default
        assert result["enabled"] is True  # default
        assert result["validity_warnings"] is True  # default

    def test_valid_multi_param(self) -> None:
        result = _validate_perturbation_config(
            {"small_parameters": ["b1", "b5"], "order": 2},
            ["kappa", "b1", "b5"],
        )
        assert result is not None
        assert result["small_parameters"] == ["b1", "b5"]
        assert result["order"] == 2

    def test_validity_warnings_override(self) -> None:
        result = _validate_perturbation_config(
            {"small_parameters": ["b5"], "validity_warnings": False},
            ["b5"],
        )
        assert result is not None
        assert result["validity_warnings"] is False


# === Derive-time audit warning ===


class TestAuditHigherDerivativeLagrangian:
    @pytest.mark.parametrize(
        ("label", "expr"),
        [
            ("R^2", "(1/kappa^2) * RicciScalar[] + a1 * RicciScalar[]^2"),
            ("R-tilde^2", "(1/kappa^2) * RicciScalarCD[] + b5 * RicciCDT[-a,-b]^2"),
            (
                "Ricci_uv * Ricci^uv",
                "a1 * RicciCD[-a,-b] * RicciCD[a,b] + (1/kappa^2) * RicciScalarCD[]",
            ),
            (
                "Riemann 4-cov",
                "c3 * RiemannCD[-a,-b,-c,-d] * RiemannCD[a,b,c,d]",
            ),
            ("RiemannCDT standalone", "b * RiemannCDT[-a,-b] * RiemannCDT[a,b]"),
            ("EH F^4 power", "(alpha^2) * F[-a,-b]**4"),
            ("EH (F.F)^2", "(alpha^2) * (F[-a,-b] * F[a,b])^2"),
            (
                "EH F^4 as product",
                "(alpha^2/mass^4) * FStrength[-a,-b] * FStrength[a,b] * "
                "FStrength[-c,-d] * FStrength[c,d]",
            ),
            ("FStrength^4", "(alpha^2) * FStrength[-a,-b]^4"),
            ("Torsion cubed", "c * TorsionCD[-a,-b,-c]^3"),
        ],
    )
    def test_flags_higher_derivative_patterns(
        self, capsys: pytest.CaptureFixture[str], label: str, expr: str
    ) -> None:
        """Each pattern should trigger the audit warning."""
        config: dict[str, Any] = {"lagrangian": {"expression": expr}}
        _audit_higher_derivative_lagrangian(config)
        captured = capsys.readouterr()
        # Warning goes to stderr via _cwarn; the hint goes to stdout.
        assert "higher-derivative patterns" in captured.err, (
            f"[{label}] expected warning for expression {expr!r}, got: {captured.err!r}"
        )
        assert "[perturbation]" in captured.out

    def test_clean_lagrangian_no_warning(
        self, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """Pure Einstein-Maxwell (no higher derivatives) should produce no output."""
        config: dict[str, Any] = {
            "lagrangian": {
                "expression": "(1/kappa^2) * RicciScalarCD[] - (1/4) * F[-a,-b] * F[a,b]"
            }
        }
        _audit_higher_derivative_lagrangian(config)
        captured = capsys.readouterr()
        assert "[WARN]" not in captured.err
        assert "[perturbation]" not in captured.out

    def test_perturbation_configured_silences_warning(
        self, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """When [perturbation] is configured, no warning is emitted."""
        config: dict[str, Any] = {
            "lagrangian": {"expression": "b5 * RicciScalar[]^2"},
            "perturbation": {
                "small_parameters": ["b5"],
                "order": 1,
                "enabled": True,
            },
        }
        _audit_higher_derivative_lagrangian(config)
        captured = capsys.readouterr()
        assert "[WARN]" not in captured.err
        assert "[perturbation]" not in captured.out

    def test_empty_lagrangian_no_warning(
        self, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """Missing or empty Lagrangian expression should not crash."""
        _audit_higher_derivative_lagrangian({})
        _audit_higher_derivative_lagrangian({"lagrangian": {}})
        _audit_higher_derivative_lagrangian({"lagrangian": {"expression": ""}})
        captured = capsys.readouterr()
        assert "[WARN]" not in captured.err


# === WLS generation ===


@pytest.fixture
def gertsenshtein_config() -> dict[str, Any]:
    """Load the gertsenshtein theory TOML (a standard small-theory fixture)."""
    path = Path(__file__).parent.parent / "examples" / "gertsenshtein" / "theory.toml"
    with path.open("rb") as f:
        return tomllib.load(f)


class TestWlsGeneration:
    def test_no_perturbation_section_omits_metadata_keys(
        self, gertsenshtein_config: dict[str, Any]
    ) -> None:
        """Theory without [perturbation] omits small_parameters from metadata."""
        wls = generate_wls(gertsenshtein_config)
        assert '"small_parameters"' not in wls
        assert '"perturbation_order"' not in wls
        # The v5 JLM module should no longer be loaded or referenced anywhere.
        assert "PerturbativeReduction.wl" not in wls
        assert "PerturbativeReduce[" not in wls

    def test_with_perturbation_injects_small_parameters_into_metadata(
        self, gertsenshtein_config: dict[str, Any]
    ) -> None:
        """TOML with [perturbation] injects small_parameters into metadata."""
        cfg = dict(gertsenshtein_config)
        cfg["perturbation"] = {"small_parameters": ["kappa"], "order": 1}
        wls = generate_wls(cfg)
        # The small parameter symbol list must appear as a Wolfram list
        # inside the metadata Association.
        match = re.search(r'"small_parameters"\s*->\s*\{kappa\}', wls)
        assert match is not None, (
            "metadata must carry small_parameters -> {kappa} when [perturbation] "
            f"is configured. WLS excerpt: {wls[:500]}"
        )

    def test_multi_parameter_injection(
        self, gertsenshtein_config: dict[str, Any]
    ) -> None:
        """Multi-parameter config passes all names into the metadata list."""
        cfg = dict(gertsenshtein_config)
        cfg["perturbation"] = {"small_parameters": ["kappa", "B0"], "order": 1}
        wls = generate_wls(cfg)
        match = re.search(r'"small_parameters"\s*->\s*\{kappa,\s*B0\}', wls)
        assert match is not None, (
            "metadata must carry small_parameters -> {kappa, B0} for multi-param "
            f"config. WLS excerpt: {wls[:500]}"
        )

    def test_order_propagated_into_metadata(
        self, gertsenshtein_config: dict[str, Any]
    ) -> None:
        """Default truncation order flows into metadata.perturbation_order."""
        cfg = dict(gertsenshtein_config)
        cfg["perturbation"] = {"small_parameters": ["kappa"], "order": 2}
        wls = generate_wls(cfg)
        match = re.search(r'"perturbation_order"\s*->\s*(\d+)', wls)
        assert match is not None, (
            f"could not locate perturbation_order in metadata: {wls[:500]}"
        )
        assert match.group(1) == "2"

    def test_disabled_is_noop(self, gertsenshtein_config: dict[str, Any]) -> None:
        """enabled=false suppresses metadata injection entirely."""
        cfg = dict(gertsenshtein_config)
        cfg["perturbation"] = {
            "small_parameters": ["kappa"],
            "order": 1,
            "enabled": False,
        }
        wls = generate_wls(cfg)
        assert '"small_parameters"' not in wls
        assert '"perturbation_order"' not in wls

    def test_invalid_config_raises_at_generate(
        self, gertsenshtein_config: dict[str, Any]
    ) -> None:
        """Malformed [perturbation] is rejected at generate_wls, before Wolfram."""
        cfg = dict(gertsenshtein_config)
        cfg["perturbation"] = {"small_parameters": ["undeclared"], "order": 1}
        with pytest.raises(ValueError, match="not declared in \\[constants\\]"):
            generate_wls(cfg)
