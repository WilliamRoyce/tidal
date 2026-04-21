"""Tests for tidal.symbolic.latex — LaTeX equation export."""

from __future__ import annotations

from pathlib import Path
from typing import ClassVar

import pytest

from tidal.symbolic.latex import (
    coefficient_to_latex,
    equation_to_latex,
    field_to_latex,
    hamiltonian_to_latex,
    lagrangian_to_latex,
    operator_to_latex,
    system_to_latex,
)

# ---------------------------------------------------------------------------
# Test data paths
# ---------------------------------------------------------------------------

_EXAMPLES = Path(__file__).resolve().parent.parent / "examples" / "data"

# Derived JSON files live in gitignored examples/data/ — only present after
# running `tidal derive`.  Skip integration tests when absent.
_CS_JSON = _EXAMPLES / "coupled_scalars.json"
_needs_cs = pytest.mark.skipif(
    not _CS_JSON.exists(), reason="coupled_scalars.json not derived",
)


# ---------------------------------------------------------------------------
# coefficient_to_latex
# ---------------------------------------------------------------------------


class TestCoefficientToLatex:
    """Unit tests for Mathematica symbolic → LaTeX coefficient conversion."""

    def test_empty(self) -> None:
        assert not coefficient_to_latex("")

    def test_simple_param(self) -> None:
        result = coefficient_to_latex("m2")
        assert "m" in result

    def test_greek(self) -> None:
        result = coefficient_to_latex("kappa")
        assert r"\kappa" in result

    def test_negation(self) -> None:
        result = coefficient_to_latex("-m2")
        assert result.startswith("-")

    def test_product(self) -> None:
        result = coefficient_to_latex("B0^2*kappa^2")
        assert r"\kappa" in result
        assert "B" in result

    def test_simple_fraction(self) -> None:
        result = coefficient_to_latex("1/2")
        assert r"\frac" in result

    def test_negative_fraction(self) -> None:
        result = coefficient_to_latex("-1/2")
        assert r"\frac" in result
        assert "-" in result

    def test_fraction_with_param(self) -> None:
        result = coefficient_to_latex("1/(2*kappa^2)")
        assert r"\frac" in result
        assert r"\kappa" in result

    def test_prefix_fraction_product(self) -> None:
        """Ensure -1/2*(B0^2*kappa^2) renders as fraction times product."""
        result = coefficient_to_latex("-1/2*(B0^2*kappa^2)")
        assert r"\frac{1}{2}" in result
        assert r"\kappa" in result

    def test_coordinate_call_stripped(self) -> None:
        result = coefficient_to_latex("x[]^(-1)")
        assert "[]" not in result

    def test_sqrt(self) -> None:
        result = coefficient_to_latex("Sqrt[2]")
        assert r"\sqrt" in result

    def test_e_power(self) -> None:
        result = coefficient_to_latex("E^(x)")
        assert "e^" in result

    def test_rational(self) -> None:
        result = coefficient_to_latex("Rational[1, 3]")
        assert r"\frac{1}{3}" in result

    def test_pi(self) -> None:
        result = coefficient_to_latex("Pi")
        assert r"\pi" in result

    def test_tanh(self) -> None:
        result = coefficient_to_latex("Tanh[x/W]")
        assert r"\tanh" in result
        assert "x/W" in result

    def test_sin(self) -> None:
        result = coefficient_to_latex("Sin[2*x]")
        assert r"\sin" in result

    def test_abs(self) -> None:
        result = coefficient_to_latex("Abs[x]")
        assert r"\left|" in result
        assert r"\right|" in result

    def test_lam_as_lambda(self) -> None:
        result = coefficient_to_latex("lam")
        assert r"\lambda" in result

    def test_greek_prefix_split(self) -> None:
        """omegaP2 should split Greek prefix from trailing label."""
        result = coefficient_to_latex("omegaP2")
        assert r"\omega" in result


# ---------------------------------------------------------------------------
# field_to_latex
# ---------------------------------------------------------------------------


class TestFieldToLatex:
    """Unit tests for field name → LaTeX conversion."""

    def test_greek_field(self) -> None:
        assert r"\phi" in field_to_latex("phi_0")

    def test_roman_field_calligraphic(self) -> None:
        result = field_to_latex("h_5")
        assert r"\mathcal{H}" in result
        assert "5" in result

    def test_vector_field_calligraphic(self) -> None:
        result = field_to_latex("a_1")
        assert r"\mathcal{A}" in result

    def test_velocity_prefix(self) -> None:
        result = field_to_latex("v_phi_0")
        assert r"\dot" in result
        assert r"\phi" in result

    def test_chi_field(self) -> None:
        assert r"\chi" in field_to_latex("chi_0")

    def test_with_tensor_meta_scalar(self) -> None:
        meta: dict[str, str | int | list[int]] = {
            "tensor_head": "phi",
            "tensor_rank": 0,
            "tensor_indices": [],
        }
        result = field_to_latex("phi_0", tensor_meta=meta)
        assert r"\phi" in result

    def test_with_tensor_meta_vector(self) -> None:
        meta = {"tensor_head": "a", "tensor_rank": 1, "tensor_indices": [1]}
        result = field_to_latex(
            "a_1", tensor_meta=meta, coordinates=("t", "x", "y", "z"),
        )
        assert r"\mathcal{A}" in result
        assert "x" in result

    def test_with_tensor_meta_rank2(self) -> None:
        meta = {"tensor_head": "h", "tensor_rank": 2, "tensor_indices": [2, 2]}
        result = field_to_latex(
            "h_5", tensor_meta=meta, coordinates=("t", "x", "y", "z"),
        )
        assert r"\mathcal{H}" in result
        assert "yy" in result

    def test_with_tensor_meta_rank2_no_coords(self) -> None:
        """Falls back to numeric indices when coordinates not provided."""
        meta = {"tensor_head": "h", "tensor_rank": 2, "tensor_indices": [2, 2]}
        result = field_to_latex("h_5", tensor_meta=meta)
        assert "22" in result

    def test_velocity_with_tensor_meta(self) -> None:
        meta = {"tensor_head": "h", "tensor_rank": 2, "tensor_indices": [2, 2]}
        result = field_to_latex(
            "v_h_5", tensor_meta=meta, coordinates=("t", "x", "y", "z"),
        )
        assert r"\dot" in result


# ---------------------------------------------------------------------------
# operator_to_latex
# ---------------------------------------------------------------------------


class TestOperatorToLatex:
    """Unit tests for operator → LaTeX rendering."""

    def test_identity(self) -> None:
        assert operator_to_latex("identity", r"\phi") == r"\phi"

    def test_laplacian(self) -> None:
        result = operator_to_latex("laplacian", r"\phi")
        assert r"\nabla^2" in result

    def test_laplacian_x(self) -> None:
        result = operator_to_latex("laplacian_x", r"\phi")
        assert r"\partial_x^2" in result

    def test_gradient_x(self) -> None:
        result = operator_to_latex("gradient_x", r"\phi")
        assert r"\partial_x" in result

    def test_biharmonic(self) -> None:
        result = operator_to_latex("biharmonic", r"\phi")
        assert r"\nabla^4" in result

    def test_first_derivative_t(self) -> None:
        result = operator_to_latex("first_derivative_t", r"\phi")
        assert r"\partial_t" in result

    def test_cross_derivative(self) -> None:
        result = operator_to_latex("cross_derivative_xy", r"\phi")
        assert r"\partial_x" in result
        assert r"\partial_y" in result

    def test_time_derivative(self) -> None:
        result = operator_to_latex("time_derivative", r"\phi")
        assert r"\dot" in result

    def test_d2_t(self) -> None:
        result = operator_to_latex("d2_t", r"\phi")
        assert r"\partial_t^2" in result

    def test_d3_t(self) -> None:
        result = operator_to_latex("d3_t", r"\phi")
        assert r"\partial_t^{3}" in result

    def test_mixed_new_format(self) -> None:
        result = operator_to_latex("mixed_T2_S1x", r"\phi")
        assert r"\partial_t" in result
        assert r"\partial_x" in result

    def test_mixed_old_format(self) -> None:
        result = operator_to_latex("mixed_1_0_0_1", r"\phi")
        assert r"\partial_t" in result
        assert r"\partial_z" in result

    def test_mixed_old_time_only(self) -> None:
        """mixed_2_0_0_0 = pure second time derivative."""
        result = operator_to_latex("mixed_2_0_0_0", r"\phi")
        assert r"\partial_t" in result
        assert "2" in result

    def test_dynamic_single_axis(self) -> None:
        result = operator_to_latex("derivative_3_x", r"\phi")
        assert r"\partial_x" in result
        assert "3" in result

    def test_unknown_operator(self) -> None:
        result = operator_to_latex("unknown_op", r"\phi")
        assert r"\mathrm{unknown_op}" in result


# ---------------------------------------------------------------------------
# lagrangian_to_latex
# ---------------------------------------------------------------------------


class TestLagrangianToLatex:
    """Tests for xAct Lagrangian → LaTeX with tensor index notation."""

    def test_empty(self) -> None:
        assert not lagrangian_to_latex("")

    def test_klein_gordon(self) -> None:
        expr = "-1/2 CD[-a][phi[]] eta[a, b] CD[-b][phi[]] - m2/2 phi[]^2"
        result = lagrangian_to_latex(expr)
        assert r"\nabla" in result
        assert r"\phi" in result
        assert r"\tensor" in result
        assert r"\eta" in result

    def test_field_strength(self) -> None:
        expr = "F[-a,-b] eta[a,c] eta[b,d] F[-c,-d]"
        result = lagrangian_to_latex(expr)
        assert r"\tensor{F}{_a_b}" in result
        assert r"\tensor{\eta}" in result

    def test_ricci_scalar(self) -> None:
        result = lagrangian_to_latex("RicciScalarCD[]")
        assert result == r"\mathcal{R}"

    def test_ricci_scalar_torsion(self) -> None:
        result = lagrangian_to_latex("RicciScalarCDT[]")
        assert r"\tilde{\mathcal{R}}" in result

    def test_ricci_scalar_torsion_squared(self) -> None:
        result = lagrangian_to_latex("b5 RicciScalarCDT[]^2")
        assert r"\tilde{\mathcal{R}}" in result
        assert "b_{5}" in result

    def test_parameter_subscript_splitting(self) -> None:
        result = lagrangian_to_latex("alpha1 * phi[]^2")
        assert r"\alpha_{1}" in result
        assert r"\phi" in result

    def test_levi_civita(self) -> None:
        expr = "epsiloneta[a, b, c] * A[-a] * CD[-b][A[-c]]"
        result = lagrangian_to_latex(expr)
        assert r"\tensor{\epsilon}{^a^b^c}" in result

    def test_torsion(self) -> None:
        expr = "TorsionCDT[-a,-b,-c] eta[a,d] eta[b,e] eta[c,f] TorsionCDT[-d,-e,-f]"
        result = lagrangian_to_latex(expr)
        assert r"\tensor{T}{_a_b_c}" in result

    def test_vector_field_indices(self) -> None:
        expr = "A[-a] eta[a,b] A[-b]"
        result = lagrangian_to_latex(expr)
        assert r"\tensor" in result
        assert "_a" in result

    def test_scalar_field_no_indices(self) -> None:
        expr = "phi[]^2"
        result = lagrangian_to_latex(expr)
        assert r"\phi" in result
        assert "[]" not in result

    def test_background_metric(self) -> None:
        expr = "bg[-a, -b]"
        result = lagrangian_to_latex(expr)
        assert r"\bar{g}" in result

    def test_position_dependent_scalar(self) -> None:
        expr = "G[] * phi[]"
        result = lagrangian_to_latex(expr)
        assert "G" in result
        assert r"\phi" in result

    def test_coupled_scalars(self) -> None:
        expr = "B0 * h[] * n[a] * CD[-a][a[]]"
        result = lagrangian_to_latex(expr)
        assert r"\nabla" in result
        assert r"\tensor{n}{^a}" in result

    def test_3form(self) -> None:
        expr = "C[-a, -b, -c] eta[a, e] eta[b, f] eta[c, g] C[-e, -f, -g]"
        result = lagrangian_to_latex(expr)
        assert r"\tensor{C}{_a_b_c}" in result

    def test_fraction_in_parens(self) -> None:
        expr = "(1/kappa^2) RicciScalarCD[]"
        result = lagrangian_to_latex(expr)
        assert r"\frac" in result
        assert "R" in result


# ---------------------------------------------------------------------------
# equation_to_latex (integration tests with real JSONs)
# ---------------------------------------------------------------------------


class TestEquationToLatex:
    """Integration tests for equation rendering with real JSON specs."""

    @_needs_cs
    def test_single_field(self) -> None:
        from tidal.symbolic.json_loader import load_equation_system

        spec = load_equation_system(_EXAMPLES / "coupled_scalars.json")
        result = equation_to_latex(spec.equations[0], spec)
        assert r"\partial_t" in result
        assert r"\mathcal" in result
        assert "&=" in result

    @_needs_cs
    def test_coupled_scalars(self) -> None:
        from tidal.symbolic.json_loader import load_equation_system

        spec = load_equation_system(_EXAMPLES / "coupled_scalars.json")
        # Should have 2 equations
        assert len(spec.equations) == 2
        for eq in spec.equations:
            result = equation_to_latex(eq, spec)
            assert "&=" in result
            assert r"\partial_t" in result


# ---------------------------------------------------------------------------
# hamiltonian_to_latex
# ---------------------------------------------------------------------------


class TestHamiltonianToLatex:
    """Integration tests for Hamiltonian rendering."""

    @_needs_cs
    def test_hamiltonian(self) -> None:
        from tidal.symbolic.json_loader import load_equation_system

        spec = load_equation_system(_EXAMPLES / "coupled_scalars.json")
        assert spec.canonical is not None
        result = hamiltonian_to_latex(list(spec.canonical.hamiltonian_terms), spec)
        assert r"\mathscr{H}" in result
        assert "&=" in result

    @_needs_cs
    def test_coupled_scalars_hamiltonian(self) -> None:
        from tidal.symbolic.json_loader import load_equation_system

        spec = load_equation_system(_EXAMPLES / "coupled_scalars.json")
        assert spec.canonical is not None
        result = hamiltonian_to_latex(list(spec.canonical.hamiltonian_terms), spec)
        assert r"\mathscr{H}" in result


# ---------------------------------------------------------------------------
# system_to_latex
# ---------------------------------------------------------------------------


@_needs_cs
class TestSystemToLatex:
    """Tests for full system rendering and output formats."""

    def test_align_format(self) -> None:
        from tidal.symbolic.json_loader import load_equation_system

        spec = load_equation_system(_EXAMPLES / "coupled_scalars.json")
        result = system_to_latex(spec, output_format="align")
        assert r"\begin{align}" in result
        assert r"\end{align}" in result
        assert r"\mathcal{L}" in result

    def test_document_format(self) -> None:
        from tidal.symbolic.json_loader import load_equation_system

        spec = load_equation_system(_EXAMPLES / "coupled_scalars.json")
        result = system_to_latex(spec, output_format="document")
        assert r"\documentclass{article}" in result
        assert r"\usepackage{tensor}" in result
        assert r"\begin{document}" in result
        assert r"\end{document}" in result
        assert r"\begin{align}" in result

    def test_raw_format(self) -> None:
        from tidal.symbolic.json_loader import load_equation_system

        spec = load_equation_system(_EXAMPLES / "coupled_scalars.json")
        result = system_to_latex(spec, output_format="raw")
        assert r"\begin{align}" not in result
        # Should have Lagrangian + equation + Hamiltonian lines
        lines = result.strip().split("\n")
        assert len(lines) >= 2

    def test_no_hamiltonian(self) -> None:
        from tidal.symbolic.json_loader import load_equation_system

        spec = load_equation_system(_EXAMPLES / "coupled_scalars.json")
        result = system_to_latex(spec, include_hamiltonian=False)
        assert r"\mathscr{H}" not in result

    def test_no_lagrangian(self) -> None:
        from tidal.symbolic.json_loader import load_equation_system

        spec = load_equation_system(_EXAMPLES / "coupled_scalars.json")
        result = system_to_latex(spec, include_lagrangian=False)
        assert r"\mathcal{L}" not in result

    def test_package_comment(self) -> None:
        from tidal.symbolic.json_loader import load_equation_system

        spec = load_equation_system(_EXAMPLES / "coupled_scalars.json")
        result = system_to_latex(spec, output_format="align")
        assert "tensor" in result


# ---------------------------------------------------------------------------
# CLI integration test
# ---------------------------------------------------------------------------


@_needs_cs
class TestCLIIntegration:
    """Tests that the CLI --latex flag works end-to-end."""

    def test_inspect_latex_flag(self) -> None:
        from tidal.cli import main

        exit_code = main(
            ["inspect", str(_EXAMPLES / "coupled_scalars.json"), "--latex"],
        )
        assert exit_code == 0

    def test_inspect_latex_document_format(self) -> None:
        from tidal.cli import main

        exit_code = main(
            [
                "inspect",
                str(_EXAMPLES / "coupled_scalars.json"),
                "--latex",
                "--latex-format",
                "document",
            ],
        )
        assert exit_code == 0


# ---------------------------------------------------------------------------
# Crash test: all example JSONs must render without exceptions
# ---------------------------------------------------------------------------


class TestAllExamplesRender:
    """Ensure every example JSON renders without errors."""

    # v6 re-derivation migration complete for all three shipped higher-
    # derivative theories (Phase 6B.3/6B.4 of #267).  Kept as an empty
    # class-level set for future migrations.
    _PENDING_V6_REDERIVATION: ClassVar[set[str]] = set()

    @pytest.mark.parametrize(
        "json_file",
        sorted(_EXAMPLES.glob("*.json")),
        ids=lambda p: p.stem,
    )
    def test_no_crash(self, json_file: Path, request: pytest.FixtureRequest) -> None:
        if json_file.stem in self._PENDING_V6_REDERIVATION:
            request.applymarker(
                pytest.mark.xfail(
                    reason=(
                        "Pending Phase 6B.3/6B.4 re-derivation with "
                        "[perturbation] section (v6 #267)."
                    ),
                    raises=ValueError,
                    strict=True,
                ),
            )
        from tidal.symbolic.json_loader import load_equation_system

        spec = load_equation_system(json_file)
        result = system_to_latex(spec, output_format="align")
        assert r"\begin{align}" in result
        assert r"\end{align}" in result
