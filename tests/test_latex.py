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
    load_symbol_overrides,
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
    not _CS_JSON.exists(),
    reason="coupled_scalars.json not derived",
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
        assert r"\tfrac" in result

    def test_negative_fraction(self) -> None:
        result = coefficient_to_latex("-1/2")
        assert r"\tfrac" in result
        assert "-" in result

    def test_fraction_with_param(self) -> None:
        result = coefficient_to_latex("1/(2*kappa^2)")
        assert r"\tfrac" in result
        assert r"\kappa" in result

    def test_prefix_fraction_product(self) -> None:
        """Ensure -1/2*(B0^2*kappa^2) renders as fraction times product."""
        result = coefficient_to_latex("-1/2*(B0^2*kappa^2)")
        assert r"\tfrac{1}{2}" in result
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
        assert r"\tfrac{1}{3}" in result

    def test_pi(self) -> None:
        result = coefficient_to_latex("Pi")
        assert r"\pi" in result

    def test_tanh(self) -> None:
        result = coefficient_to_latex("Tanh[x/W]")
        assert r"\tanh" in result
        # The bare-denominator regex in _coefficient_inner now converts the
        # ``x/W`` slash form to a proper fraction inside the \tanh argument.
        assert r"\tfrac{x}{W}" in result

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
            "a_1",
            tensor_meta=meta,
            coordinates=("t", "x", "y", "z"),
        )
        assert r"\mathcal{A}" in result
        assert "x" in result

    def test_with_tensor_meta_rank2(self) -> None:
        meta = {"tensor_head": "h", "tensor_rank": 2, "tensor_indices": [2, 2]}
        result = field_to_latex(
            "h_5",
            tensor_meta=meta,
            coordinates=("t", "x", "y", "z"),
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
            "v_h_5",
            tensor_meta=meta,
            coordinates=("t", "x", "y", "z"),
        )
        assert r"\dot" in result

    def test_velocity_subscript_outside_dot(self) -> None:
        # Issue #371: \dot{\mathcal{X}_{N}} triggers stix-mathcal missing-digit
        # warnings under \usepackage{accents}; subscript must sit outside.
        meta = {"tensor_head": "a", "tensor_rank": 1, "tensor_indices": [1]}
        result = field_to_latex(
            "v_a_1",
            tensor_meta=meta,
            coordinates=("t", "x", "y", "z"),
        )
        assert result == r"\dot{\mathcal{A}}_{x}"

    def test_velocity_subscript_outside_dot_no_coords(self) -> None:
        meta = {"tensor_head": "h", "tensor_rank": 2, "tensor_indices": [2, 2]}
        result = field_to_latex("v_h_5", tensor_meta=meta)
        assert result == r"\dot{\mathcal{H}}_{22}"

    def test_velocity_no_subscript_unchanged(self) -> None:
        assert field_to_latex("v_phi") == r"\dot{\phi}"


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

    def test_time_derivative_subscript_outside_dot(self) -> None:
        # Issue #371: lift trailing _{N} / ^{N} out of \dot{...}.
        assert (
            operator_to_latex("time_derivative", r"\mathcal{A}_{2}")
            == r"\dot{\mathcal{A}}_{2}"
        )
        assert (
            operator_to_latex("time_derivative", r"\mathcal{A}^{7}")
            == r"\dot{\mathcal{A}}^{7}"
        )

    def test_time_derivative_preserves_braced_index(self) -> None:
        # Indices with nested braces (e.g. \mu\nu inside another group)
        # are NOT lifted — keep accent placement correct for tensor bases.
        result = operator_to_latex("time_derivative", r"\mathcal{T}_{\mu\nu}")
        assert result == r"\dot{\mathcal{T}}_{\mu\nu}"

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
        assert r"\tfrac" in result
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

    def test_lhs_kinetic_coefficient_renders(self) -> None:
        """LHS gains the symbolic kinetic coefficient when one is set.

        The Wolfram pipeline strips a non-unity kinetic coefficient
        from the LHS (to avoid divide-by-zero on the constraint side)
        and stores it under ``equation.lhs.kinetic_coefficient_symbolic``;
        the renderer must restore it so the equation is mathematically
        complete.
        """
        from tidal.symbolic.json_loader import load_equation_system

        # gertsenshtein.json's h_5 / h_7 EOMs have kinetic_coefficient_symbolic
        # = "-kappa^(-2)". The bare-slash converter turns that into a fraction
        # which must appear before \partial_t in the rendered LHS.
        spec = load_equation_system(_EXAMPLES / "gertsenshtein.json")
        eq = next(
            e for e in spec.equations if e.kinetic_coefficient_symbolic == "-kappa^(-2)"
        )
        result = equation_to_latex(eq, spec)
        assert r"\tfrac{1}{\kappa^{2}}" in result, result
        assert r"\partial_t" in result
        # Ordering: the kinetic prefactor sits to the left of \partial_t.
        idx_frac = result.index(r"\tfrac{1}{\kappa^{2}}")
        idx_partial = result.index(r"\partial_t")
        assert idx_frac < idx_partial, (
            f"kinetic coefficient must precede \\partial_t: {result!r}"
        )


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
        # The rendered Hamiltonian is the self-GW + self-EM sector
        # restriction relevant to the conversion measurement; cross-sector
        # interaction terms and torsion contributions are dropped, so the
        # LHS uses \supset (contains) rather than = (equality).
        assert r"&\supset" in result

    @_needs_cs
    def test_coupled_scalars_hamiltonian(self) -> None:
        from tidal.symbolic.json_loader import load_equation_system

        spec = load_equation_system(_EXAMPLES / "coupled_scalars.json")
        assert spec.canonical is not None
        result = hamiltonian_to_latex(list(spec.canonical.hamiltonian_terms), spec)
        assert r"\mathscr{H}" in result


class TestSectorFilter:
    """Tests for the self-GW + self-EM sector filter applied at render time.

    ``_sector`` and ``_is_self_sector_term`` are module-private helpers; we
    import them inside each test to avoid auto-formatter rules that strip
    underscore-prefixed names from the module-level import list.
    """

    def test_sector_gw_prefixes(self) -> None:
        from tidal.symbolic.latex import _sector

        assert _sector("h_xy") == "gw"
        assert _sector("h_yy") == "gw"
        assert _sector("H_tz") == "gw"

    def test_sector_em_prefixes(self) -> None:
        from tidal.symbolic.latex import _sector

        assert _sector("A_t") == "em"
        assert _sector("A_x") == "em"
        assert _sector("a_phi") == "em"

    def test_sector_torsion_prefixes(self) -> None:
        from tidal.symbolic.latex import _sector

        assert _sector("t_xyz") == "torsion"
        assert _sector("T_xtx") == "torsion"

    def test_sector_other(self) -> None:
        from tidal.symbolic.latex import _sector

        # Constraint multipliers, residuals, or unrecognized field names
        # fall through to "other".
        assert _sector("lambda_0") == "other"
        assert _sector("phi_1") == "other"

    def test_self_sector_term_same_sector(self) -> None:
        """GW-GW and EM-EM terms pass the filter."""
        from types import SimpleNamespace

        from tidal.symbolic.latex import _is_self_sector_term

        gw_gw = SimpleNamespace(
            factor_a=SimpleNamespace(field="h_xy"),
            factor_b=SimpleNamespace(field="h_yy"),
        )
        em_em = SimpleNamespace(
            factor_a=SimpleNamespace(field="A_x"),
            factor_b=SimpleNamespace(field="A_y"),
        )
        assert _is_self_sector_term(gw_gw)
        assert _is_self_sector_term(em_em)

    def test_self_sector_term_cross_sector_dropped(self) -> None:
        """GW-EM cross terms (the conversion-driving interactions) are dropped."""
        from types import SimpleNamespace

        from tidal.symbolic.latex import _is_self_sector_term

        gw_em = SimpleNamespace(
            factor_a=SimpleNamespace(field="h_xy"),
            factor_b=SimpleNamespace(field="A_x"),
        )
        assert not _is_self_sector_term(gw_em)

    def test_self_sector_term_other_kept(self) -> None:
        """Terms involving an "other"-sector field are kept unconditionally."""
        from types import SimpleNamespace

        from tidal.symbolic.latex import _is_self_sector_term

        other_gw = SimpleNamespace(
            factor_a=SimpleNamespace(field="lambda_0"),
            factor_b=SimpleNamespace(field="h_xy"),
        )
        assert _is_self_sector_term(other_gw)


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
        # Lagrangian density now renders as \mathscr{L} (corpus convention,
        # matching \lag in manuscript/macros.tex and the \mathscr{H} that
        # hamiltonian_to_latex() already emits).
        assert r"\mathscr{L}" in result

    def test_gather_format(self) -> None:
        """`gather` wraps each equation in its own `aligned` block.

        Used by the Appendix-E driver so each equation centers on its own
        natural width — no global `&=` column drag across the listing.
        """
        from tidal.symbolic.json_loader import load_equation_system

        spec = load_equation_system(_EXAMPLES / "coupled_scalars.json")
        result = system_to_latex(spec, output_format="gather")
        assert r"\begin{gather*}" in result
        assert r"\end{gather*}" in result
        # Each section gets its own \begin{aligned}...\end{aligned} block.
        assert result.count(r"\begin{aligned}") >= 2
        assert result.count(r"\begin{aligned}") == result.count(r"\end{aligned}")
        # Equations are separated by `\\[1ex]` to give visible inter-equation
        # spacing distinct from the tight intra-equation `\jot`.
        assert r"\\[1ex]" in result

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

        # Read-only render test — relax the v6 time_order > 2 guard so
        # theories whose [perturbation] block is intentionally disabled
        # (to extract the exact non-LPS EOMs) can still be rendered.
        spec = load_equation_system(json_file, strict_v6=False)
        result = system_to_latex(spec, output_format="align")
        assert r"\begin{align}" in result
        assert r"\end{align}" in result


# ---------------------------------------------------------------------------
# load_symbol_overrides — tensor head & parameter rebinding
# ---------------------------------------------------------------------------


class TestSymbolOverrides:
    """Project-scoped overrides for tensor heads and scalar parameters.

    Loaded via ``load_symbol_overrides(path)`` from a TOML file (typically
    ``manuscript/latex_symbols.toml``). State is module-level; each test
    re-loads from a tmp_path so state is bounded.
    """

    @staticmethod
    def _write(tmp_path: Path, content: str) -> Path:
        f = tmp_path / "symbols.toml"
        f.write_text(content, encoding="utf-8")
        return f

    def test_clear_when_path_missing(self, tmp_path: Path) -> None:
        # Loading a file that does not exist clears existing overrides.
        present = self._write(
            tmp_path,
            """
[parameters]
deltam = "\\\\delta_m"
""",
        )
        load_symbol_overrides(present)
        # Sanity: override active
        assert "\\delta_m" in coefficient_to_latex("deltam")
        # Now clear by pointing at a path that doesn't exist
        load_symbol_overrides(tmp_path / "does_not_exist.toml")
        # Override gone; raw "deltam" passes through (built-ins know no such name)
        assert "\\delta_m" not in coefficient_to_latex("deltam")

    def test_tensor_head_override(self, tmp_path: Path) -> None:
        cfg = self._write(
            tmp_path,
            """
[tensor_heads]
Ftorsion = "\\\\tilde{F}"
""",
        )
        load_symbol_overrides(cfg)
        # Lagrangian-side rendering of a tensor with overridden head.
        out = lagrangian_to_latex("Ftorsion[-a,-b]")
        assert "\\tilde{F}" in out
        assert "Ftorsion" not in out
        load_symbol_overrides(tmp_path / "clear.toml")  # cleanup

    def test_parameter_override_bare(self, tmp_path: Path) -> None:
        cfg = self._write(
            tmp_path,
            """
[parameters]
deltam = "\\\\delta_m"
""",
        )
        load_symbol_overrides(cfg)
        # Bare parameter token (no trailing digits) gets the override verbatim.
        out = coefficient_to_latex("deltam")
        assert "\\delta_m" in out
        assert "deltam" not in out
        load_symbol_overrides(tmp_path / "clear.toml")

    def test_parameter_override_with_subscript_auto_braces(
        self,
        tmp_path: Path,
    ) -> None:
        cfg = self._write(
            tmp_path,
            """
[parameters]
zt = "\\\\zeta_T"
""",
        )
        load_symbol_overrides(cfg)
        # ``zt1`` is split into base ``zt`` + digits ``1`` by the existing
        # subscript handler; the override target is auto-braced so the
        # subscript binds to the whole ``\zeta_T`` rather than producing the
        # ambiguous ``\zeta_T_{1}``.
        out = coefficient_to_latex("zt1")
        assert "{\\zeta_T}_{1}" in out
        load_symbol_overrides(tmp_path / "clear.toml")

    def test_override_does_not_break_greek(self, tmp_path: Path) -> None:
        # Loading an override TOML must not interfere with Greek substitution
        # of unrelated identifiers.
        cfg = self._write(
            tmp_path,
            """
[parameters]
deltam = "\\\\delta_m"
""",
        )
        load_symbol_overrides(cfg)
        assert "\\kappa" in coefficient_to_latex("kappa")
        assert "\\alpha_{1}" in coefficient_to_latex("alpha1")
        load_symbol_overrides(tmp_path / "clear.toml")

    def test_repeat_load_replaces_state(self, tmp_path: Path) -> None:
        first = self._write(
            tmp_path,
            """
[parameters]
deltam = "\\\\delta_m"
""",
        )
        load_symbol_overrides(first)
        assert "\\delta_m" in coefficient_to_latex("deltam")

        second = tmp_path / "second.toml"
        second.write_text(
            """
[parameters]
otherparam = "\\\\xi"
""",
            encoding="utf-8",
        )
        load_symbol_overrides(second)
        # The first override is gone (replaced, not merged).
        assert "\\delta_m" not in coefficient_to_latex("deltam")
        assert "\\xi" in coefficient_to_latex("otherparam")
        load_symbol_overrides(tmp_path / "clear.toml")
