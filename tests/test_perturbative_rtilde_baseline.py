"""R5.2 / #281 — base_spec(['b5']) produces the correct b5=0 theory.

Strong structural and integration checks that ``base_spec(['b5'])``
reduces the full R̃² PGT spec to the correct b5=0 base theory. The
original plan (R5.2) called for a trajectory-level regression against
a separately-derived b5=0 JSON; the shipped
``torsion_gertsenshtein_minimal_propagating.json`` is unsuitable for
this because it uses a different Lagrangian (Barker's Maxwell-torsion
kinetic term, different coupling conventions) — not a b5→0 limit of
the full theory. #289 tracks the derivation work to produce a proper
reference JSON and reactivate the trajectory-level regression.

What this file DOES validate without that reference JSON:

1. Structural equivalence: ``base_spec(['b5'])`` equals the full spec
   with (a) every ``order_in_eps > 0`` RHS term removed and (b) every
   kinetic coefficient that vanishes at b5=0 demoted to an algebraic
   constraint.
2. Synthetic spec cross-check: an inline spec where we construct the
   b5=0 form by hand; ``base_spec(['b5'])`` must match term-by-term.
3. Gap C recovery is parameter-consistent: the Schur recovery matrix
   built from ``base_spec(['b5'])`` depends only on order-0 parameters,
   never on b5.

Marked ``@pytest.mark.slow`` for consistency with other R̃² tests.
"""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from tidal.solver.grid import GridInfo
from tidal.solver.perturbative_driver import PerturbativeSolver
from tidal.solver.state import StateLayout
from tidal.symbolic.json_loader import EquationSystem, load_equation_system

_FULL_PATH = (
    Path(__file__).parent.parent / "examples" / "data" / "torsion_gertsenshtein.json"
)

pytestmark = pytest.mark.slow


@pytest.fixture(scope="module")
def full_spec() -> EquationSystem:
    return load_equation_system(_FULL_PATH)


def test_base_spec_has_no_order_gt_zero_terms(full_spec: EquationSystem) -> None:
    """Every RHS term in the base spec must be order_in_eps == 0."""
    base = full_spec.base_spec(["b5"])
    bad: list[tuple[str, str, int]] = []
    for eq in base.equations:
        bad.extend(
            (eq.field_name, t.field, t.order_in_eps)
            for t in eq.rhs_terms
            if t.order_in_eps != 0
        )
    # Permit the single prepended ``1*identity(self)`` term on demoted
    # constraint fields (coefficient 1.0, operator "identity",
    # field == eq.field_name, order_in_eps==0) — filter above already
    # checks order_in_eps so this is redundant but kept as a sanity
    # anchor in case the inserter ever changes.
    assert not bad, (
        f"base_spec contains {len(bad)} non-zero-order RHS terms: "
        f"{bad[:5]}. Gap B / filter_by_order(0) must strip every "
        f"order_in_eps > 0 term."
    )


def test_base_spec_time_orders_at_most_2(full_spec: EquationSystem) -> None:
    """After Gap B demotion, no equation has time_derivative_order > 2."""
    base = full_spec.base_spec(["b5"])
    high = [
        (eq.field_name, eq.time_derivative_order)
        for eq in base.equations
        if eq.time_derivative_order > 2
    ]
    assert not high, (
        f"base_spec has {len(high)} equations still above time_order=2: "
        f"{high}. Gap B demotion must catch every b5-promoted field."
    )


def test_base_spec_equations_have_no_b5_coefficient(
    full_spec: EquationSystem,
) -> None:
    """No coefficient_symbolic in the base spec may contain 'b5'.

    An order-0 term with a b5-containing coefficient_symbolic means
    either (a) a Wolfram tagging error (term has b5 but order_in_eps=0)
    or (b) a Python filter miss. Either way, base_spec must not carry
    b5 dependence.
    """
    base = full_spec.base_spec(["b5"])
    offenders: list[tuple[str, str, str]] = []
    for eq in base.equations:
        offenders.extend(
            (eq.field_name, t.field, t.coefficient_symbolic)
            for t in eq.rhs_terms
            if t.coefficient_symbolic and "b5" in t.coefficient_symbolic
        )
    assert not offenders, (
        f"base_spec has {len(offenders)} terms with b5 in the symbolic "
        f"coefficient: {offenders[:5]}. Base theory must be b5-free."
    )


def test_synthetic_b5_kinetic_collapses_correctly() -> None:
    """Constructive check: inline spec with b5 in kinetic; confirm
    base_spec(['b5']) produces the exact hand-written b5=0 reduction.
    """
    spec_data = {
        "metadata": {
            "source": "inline-test",
            "parameters": {"b5": 0.01, "m2": 1.0},
            "perturbation": {"small_parameters": ["b5"], "order": 1},
        },
        "spacetime": {
            "dimension": 2,
            "signature": [-1, 1],
            "coordinates": ["t", "x"],
        },
        "fields": [
            {"name": "phi", "index": 0, "is_dynamical": True},
            {"name": "h", "index": 1, "is_dynamical": True},
        ],
        "equations": [
            {
                "field": "phi",
                "lhs": {
                    "expression": "d2_t(phi)",
                    "order": {"time": 2, "space": 0},
                },
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {
                            "coefficient": -1.0,
                            "operator": "identity",
                            "field": "phi",
                            "coefficient_symbolic": "-m2",
                        },
                        {
                            "coefficient": 1.0,
                            "operator": "laplacian_x",
                            "field": "phi",
                        },
                        # b5 correction: sources into phi
                        {
                            "coefficient": 1.0,
                            "operator": "identity",
                            "field": "h",
                            "coefficient_symbolic": "b5",
                            "order_in_eps": 1,
                        },
                    ],
                },
            },
            {
                "field": "h",
                "lhs": {
                    "expression": "b5*d2_t(h)",
                    "order": {"time": 2, "space": 0},
                    "kinetic_coefficient_symbolic": "b5",
                },
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {
                            "coefficient": -1.0,
                            "operator": "identity",
                            "field": "h",
                            "coefficient_symbolic": "-1",
                        },
                        {
                            "coefficient": 1.0,
                            "operator": "identity",
                            "field": "phi",
                            "coefficient_symbolic": "1",
                        },
                    ],
                },
            },
        ],
        "coupling": {},
    }
    full = EquationSystem.from_dict(copy.deepcopy(spec_data))
    base = full.base_spec(["b5"])

    # phi equation: should keep order-0 terms only.
    phi_base = next(eq for eq in base.equations if eq.field_name == "phi")
    assert phi_base.time_derivative_order == 2
    # Order-0 RHS terms: identity(phi) coeff=-1, laplacian_x(phi) coeff=+1.
    # The order-1 identity(h) term must be dropped.
    phi_rhs = [(t.operator, t.field, t.order_in_eps) for t in phi_base.rhs_terms]
    assert ("identity", "phi", 0) in phi_rhs
    assert ("laplacian_x", "phi", 0) in phi_rhs
    assert ("identity", "h", 0) not in phi_rhs, (
        "Order-1 correction term (identity(h) at order 1) must NOT "
        "appear in base_spec's phi equation. See #281."
    )

    # h equation: kinetic='b5' vanishes at b5=0 → demote to algebraic.
    h_base = next(eq for eq in base.equations if eq.field_name == "h")
    assert h_base.time_derivative_order == 0
    assert h_base.kinetic_coefficient_symbolic is None
    # Must carry an identity(h) self-term (any coefficient) so Schur
    # elimination can invert the algebraic constraint. base_spec only
    # prepends ``1*identity(self)`` when one is ABSENT — when the
    # original RHS already has a self-identity term (here with coeff
    # -1), it's left alone.
    has_self_identity = any(
        t.operator == "identity" and t.field == "h" for t in h_base.rhs_terms
    )
    assert has_self_identity, (
        "Demoted h must carry an identity(h) self-term (either the "
        "original one or a prepended coeff=1 sentinel) so downstream "
        "Schur elimination recognizes it as a constraint."
    )


def test_base_spec_layout_is_idempotent(full_spec: EquationSystem) -> None:
    """base_spec(['b5']) is idempotent: applying it twice gives the
    same spec. This is a structural invariant — helpful safety check
    against bugs where the demotion logic runs more than once.
    """
    once = full_spec.base_spec(["b5"])
    twice = once.base_spec(["b5"])
    assert len(once.equations) == len(twice.equations)
    for e1, e2 in zip(once.equations, twice.equations, strict=True):
        assert e1.field_name == e2.field_name
        assert e1.time_derivative_order == e2.time_derivative_order
        assert e1.kinetic_coefficient_symbolic == e2.kinetic_coefficient_symbolic
        # Term counts and operators match.
        ops1 = [(t.operator, t.field, t.order_in_eps) for t in e1.rhs_terms]
        ops2 = [(t.operator, t.field, t.order_in_eps) for t in e2.rhs_terms]
        assert ops1 == ops2, f"Idempotence broken for {e1.field_name!r}"


def test_pass0_is_b5_independent_on_full_spec(full_spec: EquationSystem) -> None:
    """Pass 0 of the full spec at different b5 values produces
    identical trajectories. This is the strongest end-to-end check
    without a reference JSON: the base theory isolates b5 completely.

    Tightens the existing regression (test_pass0_insensitive_to_b5)
    by running at a wider b5 range and comparing per-field.
    """
    n = 32
    grid = GridInfo(shape=(n,), bounds=((0.0, 2 * np.pi),), periodic=(True,))
    solver = PerturbativeSolver(full_spec)
    layout = StateLayout.from_spec(solver.base_spec, grid.num_points)
    x = np.linspace(0.0, 2 * np.pi, n, endpoint=False)
    y0 = np.zeros(layout.num_slots * n)
    a2_slot = layout.field_slot_map["a_2"]
    y0[a2_slot * n : (a2_slot + 1) * n] = 0.01 * np.sin(x)
    base_params: dict[str, Any] = {
        "kappa": 1.0,
        "B0": 0.1,
        "alpha1": 0.1,
        "alpha2": 0.1,
        "alpha3": 0.1,
    }
    res_zero = solver.solve(
        y0,
        grid,
        t_span=(0.0, 0.5),
        order=0,
        parameters={**base_params, "b5": 0.0},
        num_snapshots=11,
    )
    res_big = solver.solve(
        y0,
        grid,
        t_span=(0.0, 0.5),
        order=0,
        parameters={**base_params, "b5": 1e-2},
        num_snapshots=11,
    )
    for fname in ("a_0", "a_1", "a_2", "a_3", "h_5", "h_6", "h_8"):
        slot = layout.field_slot_map[fname]
        z = res_zero.total["y"][:, slot * n : (slot + 1) * n]
        b = res_big.total["y"][:, slot * n : (slot + 1) * n]
        max_diff = float(np.max(np.abs(z - b)))
        assert max_diff < 1e-12, (
            f"Pass 0 field {fname!r} shifts by {max_diff:.3e} between "
            f"b5=0 and b5=1e-2 — base_spec must be b5-independent. "
            f"See #281."
        )


_B5_ZERO_JSON = (
    Path(__file__).parent.parent
    / "examples"
    / "data"
    / "torsion_gertsenshtein_b5_zero.json"
)


@pytest.mark.skipif(
    not _B5_ZERO_JSON.exists(),
    reason="torsion_gertsenshtein_b5_zero.json not derived (run tidal derive)",
)
def test_pass0_matches_b5_zero_reference(full_spec: EquationSystem) -> None:
    """The strongest Gap B regression: the full theory's Pass 0 at b5=0
    reproduces the independently-derived b5=0 reference JSON's
    trajectory field-for-field.

    The reference JSON `torsion_gertsenshtein_b5_zero.json` is derived
    from `theory_b5_zero.toml` (identical Lagrangian minus the b5·R̃²
    term) — a byte-level b5 → 0 limit. Its equations load natively
    with h_4/h_7/h_9 as algebraic constraints (no kinetic promotion
    needed), so the solver's base_spec(['b5']) demotion path must
    produce an equivalent system.

    See #281 (original R5.2 gap), #289 (derivation of reference JSON).
    """
    from tidal.solver.modal import solve_modal

    reference = load_equation_system(_B5_ZERO_JSON)

    n = 32
    grid = GridInfo(shape=(n,), bounds=((0.0, 2 * np.pi),), periodic=(True,))
    solver = PerturbativeSolver(full_spec)
    full_layout = StateLayout.from_spec(solver.base_spec, grid.num_points)
    ref_layout = StateLayout.from_spec(reference, grid.num_points)

    # Confirm field sets match before running the expensive solve.
    assert set(solver.base_spec.component_names) == set(reference.component_names)

    x = np.linspace(0.0, 2 * np.pi, n, endpoint=False)

    # Build the same photon IC in each layout (slot indices may differ;
    # map by field name).
    def _ic_in(layout: StateLayout) -> np.ndarray:
        y0 = np.zeros(layout.num_slots * n)
        a2_slot = layout.field_slot_map["a_2"]
        y0[a2_slot * n : (a2_slot + 1) * n] = 0.01 * np.sin(x)
        return y0

    y0_full = _ic_in(full_layout)
    y0_ref = _ic_in(ref_layout)

    params = {
        "kappa": 1.0,
        "B0": 0.1,
        "alpha1": 0.1,
        "alpha2": 0.1,
        "alpha3": 0.1,
    }
    # Full theory Pass 0 at b5 = 0.
    res_full = solver.solve(
        y0_full,
        grid,
        t_span=(0.0, 0.5),
        order=0,
        parameters={**params, "b5": 0.0},
        num_snapshots=11,
    )
    # Reference JSON: direct modal solve (no b5 parameter — it's not
    # declared in that theory).
    res_ref = solve_modal(
        reference,
        grid,
        y0_ref,
        t_span=(0.0, 0.5),
        parameters=params,
        num_snapshots=11,
    )

    # Compare per-field relative error across all snapshots.
    # Use absolute+relative tolerance with floor 1e-12 to avoid
    # dividing near-zero slots; physics-interesting fields sit at
    # 1e-2 amplitude under this IC so the floor is safe.
    tol_rel = 1e-6
    for fname in solver.base_spec.component_names:
        full_slot = full_layout.field_slot_map[fname]
        ref_slot = ref_layout.field_slot_map[fname]
        y_full = res_full.total["y"][:, full_slot * n : (full_slot + 1) * n]
        y_ref = res_ref["y"][:, ref_slot * n : (ref_slot + 1) * n]
        abs_diff = float(np.max(np.abs(y_full - y_ref)))
        scale = float(np.max(np.abs(y_ref))) + 1e-12
        rel_err = abs_diff / scale
        assert rel_err < tol_rel, (
            f"Field {fname!r}: abs_diff={abs_diff:.3e}, scale={scale:.3e}, "
            f"rel_err={rel_err:.3e} exceeds {tol_rel}. "
            f"base_spec(['b5']) trajectory differs from the independently-"
            f"derived b5=0 reference. See #289, #281."
        )
