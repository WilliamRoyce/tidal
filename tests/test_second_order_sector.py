"""GH #457 Stage 1: the second-order-sector classification (json_loader).

``EquationSystem.second_order_sector`` is the ONE definition of the
promoted/residual split of order-0 rows (see :class:`SecondOrderSector`
for the one-definition rule). These tests are its regression contract:

* unit tests pin the predicate (edge seeding, closure, deliberate
  non-edges, demotion interaction);
* the corpus scan pins the promoted set of every shipped spec — it FAILS
  if a second classification path or a predicate change silently moves
  the boundary;
* the index-2 pre-scan records, ahead of the Stage-3 solver build, which
  promoted fields carry velocity references with no mass support
  anywhere (the D_cc-only class the velocity-aware Schur must handle
  rather than refuse — h_3 on the E-family class).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from tidal.symbolic.json_loader import EquationSystem, operator_time_order

REPO_ROOT = Path(__file__).resolve().parent.parent


# ---------------------------------------------------------------------------
# operator_time_order — the symbolic-side source of truth
# ---------------------------------------------------------------------------


class TestOperatorTimeOrder:
    @pytest.mark.parametrize(
        ("name", "order"),
        [
            ("identity", 0),
            ("laplacian_x", 0),
            ("gradient_x", 0),
            ("biharmonic", 0),
            ("derivative_3_x", 0),
            ("cross_derivative_xy", 0),
            ("first_derivative_t", 1),
            ("d2_t", 2),
            ("d3_t", 3),
            ("d4_t", 4),
            ("mixed_T_S1x", 1),
            ("mixed_T1_S1x", 1),
            ("mixed_T2_S1x", 2),
            ("mixed_T2_S2x", 2),
            ("mixed_1_0_0_1", 1),
            ("mixed_2_1_0", 2),
        ],
    )
    def test_orders(self, name: str, order: int) -> None:
        assert operator_time_order(name) == order


# ---------------------------------------------------------------------------
# Synthetic classification tests (predicate + closure)
# ---------------------------------------------------------------------------


def _mini_spec(rows: dict[str, tuple[int, list[tuple[str, str]]]]) -> EquationSystem:
    """Minimal spec factory: field -> (lhs time order, [(operator, target)]).

    All coefficients 1.0; every row also gets an identity self-term
    (mirroring real constraint rows and the demotion-injected self-term,
    which must never create a classification edge).
    """
    fields = [{"name": f, "index": i} for i, f in enumerate(rows)]
    equations: list[dict[str, Any]] = []
    for f, (t_order, terms) in rows.items():
        expr = f"d2_t({f})" if t_order >= 2 else f
        rhs = [{"coefficient": 1.0, "operator": "identity", "field": f}]
        rhs += [{"coefficient": 1.0, "operator": op, "field": tgt} for op, tgt in terms]
        equations.append(
            {
                "field": f,
                "lhs": {"expression": expr, "order": {"time": t_order, "space": 0}},
                "rhs": {"type": "linear_combination", "terms": rhs},
            }
        )
    return EquationSystem.from_dict(
        {
            "metadata": {"name": "sector_synthetic"},
            "spacetime": {
                "dimension": 2,
                "signature": [-1, 1],
                "coordinates": ["t", "x"],
            },
            "fields": fields,
            "equations": equations,
        }
    )


class TestClassificationPredicate:
    def test_m_cc_seed_promotes_both_endpoints(self) -> None:
        spec = _mini_spec(
            {
                "c_a": (0, [("d2_t", "c_b")]),
                "c_b": (0, []),
                "phi": (2, [("laplacian_x", "phi")]),
            }
        )
        s = spec.second_order_sector
        assert s.promoted == frozenset({"c_a", "c_b"})
        assert s.reasons["c_a"] == "carries M_cc"
        assert s.reasons["c_b"] == "mass-targeted"

    def test_d_cc_seed_via_velocity_prefix(self) -> None:
        spec = _mini_spec(
            {
                "c_a": (0, [("gradient_x", "v_c_b")]),
                "c_b": (0, []),
                "phi": (2, []),
            }
        )
        s = spec.second_order_sector
        assert s.promoted == frozenset({"c_a", "c_b"})
        assert s.reasons["c_a"] == "carries D_cc"
        assert s.reasons["c_b"] == "velocity-targeted"

    def test_d_cc_seed_via_first_derivative_t(self) -> None:
        spec = _mini_spec(
            {
                "c_a": (0, [("first_derivative_t", "c_b")]),
                "c_b": (0, []),
                "phi": (2, []),
            }
        )
        assert spec.second_order_sector.promoted == frozenset({"c_a", "c_b"})

    def test_velocity_of_velocity_counts_as_mass(self) -> None:
        # first_derivative_t(v_C) = ẍ_C: total time order 2 → M_cc
        spec = _mini_spec(
            {
                "c_a": (0, [("first_derivative_t", "v_c_b")]),
                "c_b": (0, []),
                "phi": (2, []),
            }
        )
        s = spec.second_order_sector
        assert s.reasons["c_a"] == "carries M_cc"
        assert s.reasons["c_b"] == "mass-targeted"

    def test_closure_is_transitive(self) -> None:
        # a→b time edge, b→c time edge, d isolated: {a,b,c} promoted, d not.
        spec = _mini_spec(
            {
                "c_a": (0, [("d2_t", "c_b")]),
                "c_b": (0, [("gradient_x", "v_c_c")]),
                "c_c": (0, []),
                "c_d": (0, []),
                "phi": (2, []),
            }
        )
        assert spec.second_order_sector.promoted == frozenset({"c_a", "c_b", "c_c"})

    def test_spatial_only_coupling_never_promotes(self) -> None:
        # Inter-constraint SPATIAL coupling is S_cc's legitimate job.
        spec = _mini_spec(
            {
                "c_a": (0, [("laplacian_x", "c_b")]),
                "c_b": (0, [("identity", "c_a")]),
                "phi": (2, []),
            }
        )
        assert spec.second_order_sector.promoted == frozenset()

    def test_constraint_time_ref_of_dynamical_never_promotes(self) -> None:
        # Handled by S_cd velocity slots / deferred substitution — a
        # deliberate non-edge (measured-healthy machinery, WS2 oracle).
        spec = _mini_spec(
            {
                "c_a": (0, [("d2_t", "phi"), ("gradient_x", "v_phi")]),
                "phi": (2, []),
            }
        )
        assert spec.second_order_sector.promoted == frozenset()

    def test_dynamical_row_refs_never_seed(self) -> None:
        # Dynamical-row references to constraints route through A_dc —
        # their defects are GH #458, not classification.
        spec = _mini_spec(
            {
                "c_a": (0, []),
                "phi": (2, [("first_derivative_t", "c_a"), ("gradient_x", "v_c_a")]),
            }
        )
        assert spec.second_order_sector.promoted == frozenset()

    def test_identity_self_term_creates_no_edge(self) -> None:
        # Every row of the factory carries identity(self) — including the
        # demotion-injected self-term shape. Must never promote alone.
        spec = _mini_spec({"c_a": (0, []), "c_b": (0, []), "phi": (2, [])})
        assert spec.second_order_sector.promoted == frozenset()

    def test_time_self_edge_promotes_single_field(self) -> None:
        # A constraint row carrying d2_t of ITSELF is second-order content.
        spec = _mini_spec({"c_a": (0, [("d2_t", "c_a")]), "phi": (2, [])})
        s = spec.second_order_sector
        assert s.promoted == frozenset({"c_a"})
        assert s.reasons["c_a"] == "carries M_cc+mass-targeted"


# ---------------------------------------------------------------------------
# Corpus scan — the regression contract over shipped specs
# ---------------------------------------------------------------------------

#: Expected promoted sets (measured 2026-08-25, commit of Stage 1).
#: The E-family class {h_3, h_4, h_7, h_9} is the metric trace/longitudinal
#: sector: h_9 M_cc-only, h_4/h_7 M_cc+D_cc, h_3 D_cc-only. A change here
#: is a classification-boundary change and must be deliberate.
E_CLASS = frozenset({"h_3", "h_4", "h_7", "h_9"})

CORPUS_EXPECTED: dict[str, frozenset[str]] = {
    # E-family class
    "gertsenshtein_ungauged": E_CLASS,
    "gertsenshtein_ungauged_e_dual_gaussian": E_CLASS,
    "dark_photon_plasma_e_dual_gaussian": E_CLASS,
    "torsion_gertsenshtein_b5_zero": E_CLASS,
    "torsion_gertsenshtein_complete_even": E_CLASS,
    "torsion_gertsenshtein_complete_even_e_dual_gaussian": E_CLASS,
    "torsion_gertsenshtein_complete_even_full_xi": E_CLASS,
    "torsion_gertsenshtein_einstein_cartan_e_dual_gaussian": E_CLASS,
    "torsion_gertsenshtein_minimal_propagating": E_CLASS,
    "torsion_gertsenshtein_minimal_propagating_e0_dual_gaussian": E_CLASS,
    # E-class + torsion D_cc extension
    "torsion_gertsenshtein_general_nonminimal": E_CLASS
    | frozenset({"t_5", "t_14", "t_19"}),
    "torsion_gertsenshtein_general_nonminimal_e_dual_gaussian": E_CLASS
    | frozenset({"t_5", "t_14", "t_19"}),
    "torsion_gertsenshtein_propagating": E_CLASS | frozenset({"t_5", "t_14", "t_19"}),
    # larger closures
    "torsion_gertsenshtein_nonminimal": E_CLASS
    | frozenset(
        {"h_2", "t_1", "t_2", "t_5", "t_9", "t_10", "t_14", "t_17", "t_19", "t_23"}
    ),
    "torsion_gertsenshtein_nonminimal_e_dual_gaussian": E_CLASS
    | frozenset(
        {"h_2", "t_1", "t_2", "t_5", "t_9", "t_10", "t_14", "t_17", "t_19", "t_23"}
    ),
    "torsion_gertsenshtein_parity_odd": E_CLASS
    | frozenset(
        {"h_2", "t_3", "t_4", "t_5", "t_7", "t_8", "t_12", "t_14", "t_18", "t_19"}
    ),
    "torsion_gertsenshtein_parity_odd_minimal": E_CLASS
    | frozenset(
        {"h_2", "t_0", "t_3", "t_4", "t_7", "t_8", "t_12", "t_15", "t_18", "t_22"}
    ),
    "torsion_gertsenshtein_complete_odd": E_CLASS
    | frozenset(
        {
            "h_0",
            "h_1",
            "h_2",
            "t_3",
            "t_4",
            "t_5",
            "t_6",
            "t_7",
            "t_8",
            "t_11",
            "t_12",
            "t_13",
            "t_14",
            "t_16",
            "t_18",
            "t_19",
            "t_20",
            "t_21",
        }
    ),
    "torsion_gertsenshtein_complete_odd_e_dual_gaussian": E_CLASS
    | frozenset(
        {
            "h_0",
            "h_1",
            "h_2",
            "t_3",
            "t_4",
            "t_5",
            "t_6",
            "t_7",
            "t_8",
            "t_11",
            "t_12",
            "t_13",
            "t_14",
            "t_16",
            "t_18",
            "t_19",
            "t_20",
            "t_21",
        }
    ),
    # FULL specs with time_order-4 rows (h_4/h_7/h_9 are order 4 here, not
    # 0, so they are outside C₀; the solver consumes the DEMOTED base spec
    # — see TestDemotionInteraction)
    "torsion_gertsenshtein": frozenset({"h_0", "h_3", "t_2", "t_10", "t_17"}),
    "graviton_torsion": frozenset({"h_0", "h_3", "t_2", "t_10", "t_17"}),
    "torsion_gertsenshtein_combined": frozenset({"h_0", "h_1", "h_2", "h_3"}),
    "massive_gravity_3d": frozenset({"h_1", "h_2", "h_3", "h_5"}),
}

#: Every other loadable shipped spec must classify EMPTY. Notably
#: dark_photon_plasma (uniform): its constraint defect is the rank-1 S_cc
#: degeneracy (#459), NOT #457 — no inter-constraint time references.
CORPUS_EMPTY: frozenset[str] = frozenset(
    {
        "chern_simons_3d",
        "conformal_kg_static",
        "coupled_scalars",
        "coupled_scattering",
        "cylindrical_kg_1d",
        "dark_photon_plasma",
        "euler_heisenberg",
        "euler_heisenberg_e_dual_gaussian",
        "gertsenshtein",
        "gertsenshtein_e0_dual_gaussian",
        "gertsenshtein_eh",
        "gertsenshtein_eh_top",
        "gertsenshtein_localized",
        "gertsenshtein_proca",
        "gw_plane_wave_1d",
        "linearized_gravity",
        "massive_3form",
        "navier_cauchy_2d",
        "polar_kg",
        "proca_background",
        "scalar_vector_coupling",
        "spherical_kg_1d",
        "torsion_dark_photon",
        "torsion_dark_photon_fv",
    }
)


def _load_spec(name: str) -> EquationSystem:
    p = REPO_ROOT / f"examples/data/{name}.json"
    if not p.exists():
        pytest.skip(f"{name}.json not present")
    return EquationSystem.from_dict(json.loads(p.read_text()))


class TestCorpusScan:
    @pytest.mark.parametrize("name", sorted(CORPUS_EXPECTED))
    def test_promoted_set_pinned(self, name: str) -> None:
        spec = _load_spec(name)
        assert spec.second_order_sector.promoted == CORPUS_EXPECTED[name]

    @pytest.mark.parametrize("name", sorted(CORPUS_EMPTY))
    def test_empty_specs_stay_empty(self, name: str) -> None:
        spec = _load_spec(name)
        assert spec.second_order_sector.promoted == frozenset()

    def test_corpus_is_fully_pinned(self) -> None:
        """Every loadable spec in examples/data/ has a pinned expectation.

        A new spec must be classified here deliberately; this test makes
        forgetting impossible. (torsion_gertsenshtein_exact predates the
        v6 time_order>2 guard and does not load — excluded.)
        """
        known = (
            set(CORPUS_EXPECTED) | set(CORPUS_EMPTY) | {"torsion_gertsenshtein_exact"}
        )
        present = {p.stem for p in (REPO_ROOT / "examples/data").glob("*.json")}
        unpinned = present - known
        assert not unpinned, f"specs without a pinned classification: {unpinned}"


class TestDemotionInteraction:
    """The accessor is provenance-agnostic: ε-demoted base-spec rows
    classify exactly like JSON-native order-0 rows.
    """

    def test_torsion_base_spec_promotes_demoted_mass_rows(self) -> None:
        spec = _load_spec("torsion_gertsenshtein")
        base = spec.base_spec(small_parameters=["b5"])
        s = base.second_order_sector
        # Demotion moves h_4/h_7/h_9 (order 4 → 0, keeping their κ⁻²
        # ε⁰ d2_t cross terms) into C₀; they join the closure. The
        # demoted t_6/t_13/t_20 rows are pure identity — no edges.
        assert s.promoted >= E_CLASS
        assert {"t_6", "t_13", "t_20"}.isdisjoint(s.promoted)

    def test_gap_b_single_constraint_stays_empty(self) -> None:
        """The single-constraint demotion shape (Gap-B h/φ class): no
        inter-constraint edges → P = ∅ → provably zero behavior change,
        and test_perturbative_algebraic_constraint's order==0 assertion
        keeps holding.
        """
        spec = _mini_spec(
            {
                "h": (0, [("d2_t", "phi"), ("laplacian_x", "phi")]),
                "phi": (2, [("laplacian_x", "phi")]),
            }
        )
        assert spec.second_order_sector.promoted == frozenset()


class TestIndexTwoPreScan:
    """Stage-1 diagnostic (recorded for the Stage-3 build): promoted
    fields that are velocity-targeted but never mass-targeted — their ż
    appears in the coupled system while their ḧ never does (the D_cc-only
    shape). The Stage-3 velocity-aware Schur must HANDLE this shape (it
    is present in essentially every affected spec — h_3 on the E-family
    class); the index-2 refusal is only for ż of a column that ends up
    GAUGE after the per-k SVD, which cannot be decided spec-statically.
    """

    @staticmethod
    def _risky(spec: EquationSystem) -> frozenset[str]:
        s = spec.second_order_sector
        return frozenset(
            f
            for f, r in s.reasons.items()
            if "velocity-targeted" in r and "mass-targeted" not in r
        )

    def test_e_class_risk_is_h3(self) -> None:
        spec = _load_spec("gertsenshtein_ungauged")
        assert self._risky(spec) == frozenset({"h_3"})

    def test_risk_set_is_recorded_for_all_promoted_specs(self) -> None:
        # Documentation-of-record: every promoted spec's D_cc-only set is
        # non-empty somewhere in the corpus — the Stage-3 build cannot
        # treat the D_cc-only shape as a corner case.
        risky_specs = 0
        for name in CORPUS_EXPECTED:
            p = REPO_ROOT / f"examples/data/{name}.json"
            if not p.exists():
                continue
            spec = EquationSystem.from_dict(json.loads(p.read_text()))
            if self._risky(spec):
                risky_specs += 1
        assert risky_specs > 0


# ---------------------------------------------------------------------------
# Stage 2a: refusal consistency + display (accessor-only, no layout change)
# ---------------------------------------------------------------------------


class TestStage2Refusals:
    """Promoted specs must be refused EARLY and actionably by backends that
    cannot represent the second-order sector, instead of failing late
    (missing-v_-slot KeyError) or silently building wrong operators.

    The state-layout flip and the modal builders' promotion land together
    in the Stage-3/4 commit — until then the modal paths keep their
    measured (defective) behavior, pinned by the oracle xfail marks.
    """

    def test_rhs_evaluator_refuses_promoted_spec(self) -> None:
        import numpy as np

        from tidal.solver.coefficients import CoefficientEvaluator
        from tidal.solver.grid import GridInfo
        from tidal.solver.rhs import RHSEvaluator

        spec = _load_spec("gertsenshtein_ungauged")
        grid = GridInfo(bounds=((0.0, 100.0),), shape=(16,), periodic=(True,))
        ce = CoefficientEvaluator(spec, grid, {"kappa": 1.0, "B0": 0.01})
        with pytest.raises(NotImplementedError, match=r"h_3.*#457") as exc:
            RHSEvaluator(spec, grid, ce)
        assert "modal" in str(exc.value)
        del np

    def test_rhs_evaluator_accepts_ordinary_constraint_spec(self) -> None:
        from tidal.solver.coefficients import CoefficientEvaluator
        from tidal.solver.grid import GridInfo
        from tidal.solver.rhs import RHSEvaluator

        spec = _mini_spec(
            {
                "c_a": (0, [("laplacian_x", "phi")]),
                "phi": (2, [("laplacian_x", "phi")]),
            }
        )
        grid = GridInfo(bounds=((0.0, 6.0),), shape=(16,), periodic=(True,))
        ce = CoefficientEvaluator(spec, grid, {})
        RHSEvaluator(spec, grid, ce)  # must not raise

    def test_modal_jax_refuses_promoted_spec(self) -> None:
        import numpy as np

        from tidal.solver.coefficients import CoefficientEvaluator
        from tidal.solver.grid import GridInfo
        from tidal.solver.modal_jax import (
            _solve_modal_jax_constrained,  # pyright: ignore[reportPrivateUsage]
        )
        from tidal.solver.state import StateLayout

        spec = _load_spec("gertsenshtein_ungauged")
        grid = GridInfo(bounds=((0.0, 100.0),), shape=(16,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        ce = CoefficientEvaluator(spec, grid, {"kappa": 1.0, "B0": 0.01})
        # The refusal fires before any matrix work or JAX use, so dummy
        # evolution arguments are never touched.
        with pytest.raises(NotImplementedError, match=r"h_3.*#457"):
            _solve_modal_jax_constrained(
                spec,
                layout,
                grid,
                ce,
                [np.zeros(9)],
                (9,),
                np.zeros((layout.num_slots, 9), dtype=np.complex128),
                np.linspace(0.0, 1.0, 2),
                None,
                2,
                None,
                None,
            )

    def test_inspect_annotates_promoted_rows(
        self, capsys: pytest.CaptureFixture[str]
    ) -> None:
        from tidal.cli import main

        path = REPO_ROOT / "examples/data/gertsenshtein_ungauged.json"
        if not path.exists():
            pytest.skip("spec not present")
        assert main(["inspect", str(path), "--equation", "h_3"]) == 0
        out = capsys.readouterr().out
        assert "promoted to second-order sector" in out
        assert "GH #457" in out
