"""Ostrogradsky reduction for higher-derivative field equations.

Converts 4th-order-in-time equations to 2nd-order by introducing
auxiliary fields.  For ``d⁴h/dt⁴ = RHS``, defines ``w = d²h/dt²``
and rewrites as two coupled 2nd-order equations:

    d²h/dt² = w           (kinematic)
    d²w/dt² = RHS'        (dynamics, with substituted operators)

This is the standard Ostrogradsky reduction (Ostrogradsky 1850,
"Mémoires sur les équations différentielles").  Modern review:
Woodard (2015, arXiv:1506.02210).

**Physics warning**: Non-degenerate 4th-order Lagrangians generically
have Ostrogradsky ghosts (unbounded Hamiltonian).  Ghost-free parameter
windows exist for PGT (Barker 2024, arXiv:2406.12826).  Full ghost
detection is tracked in issue #164.

The reduction is applied in-memory during JSON loading, transparent
to all solvers.  The JSON on disk is unchanged.
"""

from __future__ import annotations

import logging
import re
import warnings
from dataclasses import dataclass, replace
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from tidal.symbolic.json_loader import (
        ComponentEquation,
        EquationSystem,
        OperatorTerm,
    )

logger = logging.getLogger(__name__)

#: Maximum time-derivative order that solvers support natively.
#: Equations with higher order are reduced via Ostrogradsky.
MAX_NATIVE_TIME_ORDER = 2

#: Time-derivative order constants used in Ostrogradsky reduction logic.
_T_ORDER_3 = 3  # 3rd time derivative: d³/dt³
_T_ORDER_4 = 4  # 4th time derivative: Ostrogradsky threshold
#: Spatial derivative orders for operator name mapping.
_S_ORDER_GRADIENT = 1  # 1st spatial derivative → gradient_*
_S_ORDER_LAPLACIAN = 2  # 2nd spatial derivative → laplacian_*

# --- Operator substitution for d⁴ → d² reduction ---
# Maps (operator_on_original_field) → (new_operator, use_auxiliary_field)
# "auxiliary" means the term references the auxiliary field w instead of h.
_OPERATOR_REDUCTION_MAP: dict[str, tuple[str, bool]] = {
    # Pure time derivatives → auxiliary field
    "d2_t": ("identity", True),  # d²_t(h) = w → identity(w)
    "d3_t": ("first_derivative_t", True),  # d³_t(h) = d_t(w) → v_w
    # d4_t: handled specially (cross-field substitution)
    # Mixed time-space → spatial on auxiliary
    "mixed_T2_S1x": ("gradient_x", True),  # ∂²_t ∂_x(h) = ∂_x(w)
    "mixed_T2_S2x": ("laplacian_x", True),  # ∂²_t ∂²_x(h) = ∂²_x(w)
    "mixed_T2_S1y": ("gradient_y", True),
    "mixed_T2_S2y": ("laplacian_y", True),
    "mixed_T2_S1z": ("gradient_z", True),
    "mixed_T2_S2z": ("laplacian_z", True),
}

# Pattern for mixed_T[n]_S[m][axis] operators
_MIXED_RE = re.compile(r"^mixed_T(\d+)_(.+)$")


@dataclass(frozen=True)
class AuxiliaryField:
    """Metadata for an Ostrogradsky auxiliary field.

    Attributes
    ----------
    original_field : str
        Name of the higher-order field (e.g., "h_3").
    auxiliary_name : str
        Name of the auxiliary field (e.g., "w_h_3").
    original_time_order : int
        Time derivative order of the original equation (4, 6, ...).
    """

    original_field: str
    auxiliary_name: str
    original_time_order: int


def apply_ostrogradsky_reduction(  # noqa: C901, PLR0914
    spec: EquationSystem,
) -> EquationSystem:
    """Apply Ostrogradsky reduction to higher-derivative equations.

    For each equation with ``time_derivative_order > 2``, introduces
    auxiliary fields to reduce to 2nd-order.  The reduction is applied
    in-memory — the JSON on disk is unchanged.

    Parameters
    ----------
    spec : EquationSystem
        Original equation system, possibly with higher-order equations.

    Returns
    -------
    EquationSystem
        Modified system where all equations have ``time_order <= 2``.
        Auxiliary fields are added to the component list.

    Raises
    ------
    ValueError
        If equations have odd time-derivative order (not physically
        meaningful for Euler-Lagrange equations from a Lagrangian).

    References
    ----------
        - Ostrogradsky M.V. (1850), Mem. Acad. St. Petersbourg VI 4, 385.
        - Woodard R.P. (2015), arXiv:1506.02210.
        - Barker W.E.V. et al. (2024), arXiv:2406.12826.
    """
    from tidal.symbolic.json_loader import (  # noqa: PLC0415
        ComponentEquation,
        OperatorTerm,
    )

    # Identify higher-order equations
    fourth_order = [
        eq for eq in spec.equations if eq.time_derivative_order > MAX_NATIVE_TIME_ORDER
    ]
    if not fourth_order:
        return spec  # Nothing to reduce

    # Validate: only even orders (from E-L equations)
    for eq in fourth_order:
        if eq.time_derivative_order % 2 != 0:
            msg = (
                f"Field {eq.field_name}: odd time-derivative order "
                f"{eq.time_derivative_order} not supported by Ostrogradsky "
                f"reduction (E-L equations should have even order)"
            )
            raise ValueError(msg)

    # Build auxiliary field mapping
    aux_map: dict[str, str] = {}  # original_name → auxiliary_name
    aux_fields: list[AuxiliaryField] = []
    for eq in fourth_order:
        aux_name = f"w_{eq.field_name}"
        aux_map[eq.field_name] = aux_name
        aux_fields.append(
            AuxiliaryField(
                original_field=eq.field_name,
                auxiliary_name=aux_name,
                original_time_order=eq.time_derivative_order,
            )
        )

    logger.info(
        "Ostrogradsky reduction: %d fields with order > 2 → %d auxiliary fields",
        len(fourth_order),
        len(aux_fields),
    )
    warnings.warn(
        f"Higher-derivative theory: {len(fourth_order)} equation(s) with "
        f"time order > 2 reduced via Ostrogradsky. Ghost instability may "
        f"be present (see issue #164 for propagator analysis).",
        stacklevel=2,
    )

    # Build time-order lookup for ALL fields (needed for mixed operator resolution)
    {eq.field_name: eq.time_derivative_order for eq in spec.equations}

    # Pre-process: substitute d4_t cross-references
    # If d4_t(h_5) appears in h_3's RHS, and h_5 is also 4th-order,
    # then d4_t(h_5) = RHS_h5. Substitute RHS_h5 into h_3's equation.
    rhs_lookup: dict[str, tuple[OperatorTerm, ...]] = {
        eq.field_name: eq.rhs_terms for eq in fourth_order
    }

    # Build new equation list
    new_equations: list[ComponentEquation] = []
    next_idx = max(eq.field_index for eq in spec.equations) + 1

    for eq in spec.equations:
        if eq.time_derivative_order <= MAX_NATIVE_TIME_ORDER:
            # Non-4th-order: pass through, but substitute any d2_t/d3_t
            # references to 4th-order fields with auxiliary references
            new_terms = _substitute_cross_field_refs(eq.rhs_terms, aux_map, rhs_lookup)
            new_equations.append(replace(eq, rhs_terms=new_terms))
            continue

        # 4th-order equation: split into two 2nd-order equations
        aux_name = aux_map[eq.field_name]

        # Equation 1: d²h/dt² = w (kinematic)
        eq_kinematic = ComponentEquation(
            field_name=eq.field_name,
            field_index=eq.field_index,
            time_derivative_order=2,
            rhs_terms=(
                OperatorTerm(
                    operator="identity",
                    field=aux_name,
                    coefficient=1.0,
                    coefficient_symbolic="1",
                ),
            ),
            constraint_solver=eq.constraint_solver,
        )
        new_equations.append(eq_kinematic)

        # Equation 2: d²w/dt² = RHS' (dynamics)
        # Substitute operators: d2_t(h) → identity(w), etc.
        new_rhs = _reduce_rhs_terms(eq.rhs_terms, eq.field_name, aux_map, rhs_lookup)
        eq_dynamics = ComponentEquation(
            field_name=aux_name,
            field_index=next_idx,
            time_derivative_order=2,
            rhs_terms=new_rhs,
            constraint_solver=eq.constraint_solver,
        )
        new_equations.append(eq_dynamics)
        next_idx += 1

    # Build new component names and field list
    new_names = list(spec.component_names)
    new_names.extend(af.auxiliary_name for af in aux_fields)

    # Expand mass/coupling matrices to include auxiliary fields (zero rows/cols)
    import numpy as np  # noqa: PLC0415

    n_new = len(new_names)
    n_old = spec.n_components
    n_aux = n_new - n_old

    def _expand_matrix(
        mat: tuple | np.ndarray | None,
        n_old: int,
        n_aux: int,
    ) -> tuple | np.ndarray | None:
        if mat is None:
            return None
        arr = np.array(mat, dtype=float)
        expanded = np.zeros((n_old + n_aux, n_old + n_aux))
        expanded[:n_old, :n_old] = arr
        # Return as tuple-of-tuples to match original format
        return tuple(tuple(row) for row in expanded)

    new_mass = _expand_matrix(spec.mass_matrix, n_old, n_aux)
    new_coupling = _expand_matrix(spec.coupling_matrix, n_old, n_aux)

    # Post-reduction: eliminate time-derivative operators on constraint fields.
    # For constraint fields (time_order=0), ∂²_t(field) = 0 in the linearized
    # theory. Any d2_t, d3_t, d4_t, or mixed_T* operator on a constraint
    # field produces zero → drop the term.
    constraint_fields = {
        eq.field_name for eq in spec.equations if eq.time_derivative_order == 0
    }
    time_ops = {"d2_t", "d3_t", "d4_t", "first_derivative_t"}

    cleaned_equations: list[ComponentEquation] = []
    for eq in new_equations:
        cleaned_terms = tuple(
            t
            for t in eq.rhs_terms
            if not (
                t.field in constraint_fields
                and (t.operator in time_ops or _MIXED_RE.match(t.operator))
            )
        )
        if cleaned_terms != eq.rhs_terms:
            n_dropped = len(eq.rhs_terms) - len(cleaned_terms)
            logger.info(
                "Ostrogradsky: dropped %d time-derivative terms on constraint "
                "fields in %s equation",
                n_dropped,
                eq.field_name,
            )
        cleaned_equations.append(replace(eq, rhs_terms=cleaned_terms))

    # Reconstruct EquationSystem with reduced equations
    return replace(
        spec,
        n_components=n_new,
        component_names=tuple(new_names),
        equations=tuple(cleaned_equations),
        mass_matrix=new_mass,
        coupling_matrix=new_coupling,
    )


def _reduce_rhs_terms(
    terms: tuple[OperatorTerm, ...],
    field_name: str,
    aux_map: dict[str, str],
    rhs_lookup: dict[str, tuple[OperatorTerm, ...]],
) -> tuple[OperatorTerm, ...]:
    """Reduce RHS terms for a 4th-order equation's auxiliary dynamics.

    Substitutes time-derivative operators with auxiliary field references.
    """
    result: list[OperatorTerm] = []

    for term in terms:
        new_term = _reduce_single_term(term, field_name, aux_map, rhs_lookup)
        if new_term is not None:
            if isinstance(new_term, list):
                result.extend(new_term)
            else:
                result.append(new_term)

    return tuple(result)


def _reduce_single_term(  # noqa: C901, PLR0911
    term: OperatorTerm,
    own_field: str,
    aux_map: dict[str, str],
    _rhs_lookup: dict[str, tuple[OperatorTerm, ...]],
) -> OperatorTerm | list[OperatorTerm] | None:
    """Reduce a single RHS term via Ostrogradsky substitution.

    Raises
    ------
    ValueError
        If ``d4_t`` appears on a field that is not a 4th-order field.
    """
    from tidal.symbolic.json_loader import OperatorTerm  # noqa: PLC0415

    op = term.operator
    field = term.field

    # Check if the field is a 4th-order field with auxiliary
    field_has_aux = field in aux_map
    aux_name = aux_map.get(field, "")

    # d4_t cross-reference: substitute with RHS of target equation
    if op == "d4_t" and field_has_aux:
        # d4_t(field) = d²_t(w_field) = the RHS of w_field's equation
        # This is circular for own field (d4_t(h) in h's own equation)
        # For cross-field: d4_t(h_5) in h_3's equation = RHS_h5
        if field == own_field:
            msg = (
                f"d4_t({field}) in own equation — circular reference. "
                f"This shouldn't happen in correctly derived E-L equations."
            )
            logger.warning(msg)
            return None
        # Substitute with the target field's RHS (scaled by coefficient)
        # d4_t(h_5) = d²_t(w_h_5), which becomes the RHS of w_h_5
        # But we can just keep it as d2_t(w_h_5) — the dynamics will handle it
        return OperatorTerm(
            operator="d2_t",
            field=aux_name,
            coefficient=term.coefficient,
            coefficient_symbolic=term.coefficient_symbolic,
        )

    # d4_t on non-4th-order field: error
    if op == "d4_t" and not field_has_aux:
        msg = f"d4_t({field}) but {field} is not 4th-order"
        raise ValueError(msg)

    # Check lookup table for direct reduction
    if op in _OPERATOR_REDUCTION_MAP and field_has_aux:
        new_op, use_aux = _OPERATOR_REDUCTION_MAP[op]
        target = aux_name if use_aux else field
        return OperatorTerm(
            operator=new_op,
            field=target,
            coefficient=term.coefficient,
            coefficient_symbolic=term.coefficient_symbolic,
        )

    # Check d2_t/d3_t on non-auxiliary fields (2nd-order fields)
    if op == "d2_t" and not field_has_aux:
        # d2_t on a 2nd-order field = the LHS of that field's equation
        # This is an implicit coupling — keep as-is for now
        # The solver's RHS evaluator handles d2_t via acceleration
        return term

    if op == "d3_t" and not field_has_aux:
        # d3_t on a 2nd-order field — needs velocity of acceleration
        # For now, keep as-is (solver needs special handling)
        return term

    # Handle mixed_T[n]_S[m] operators
    m = _MIXED_RE.match(op)
    if m and field_has_aux:
        t_order = int(m.group(1))
        s_part = m.group(2)

        if t_order == MAX_NATIVE_TIME_ORDER:
            # mixed_T2_Sx(h) = ∂_x(∂²_t h) = ∂_x(w) → spatial_op(w)
            spatial_op = _spatial_part_to_operator(s_part)
            return OperatorTerm(
                operator=spatial_op,
                field=aux_name,
                coefficient=term.coefficient,
                coefficient_symbolic=term.coefficient_symbolic,
            )

        if t_order == _T_ORDER_3:
            # mixed_T3_Sx(h) = ∂_x(∂³_t h) = ∂_x(∂_t w)
            # → mixed_T1_Sx(w) (one less time order)
            new_op = f"mixed_T{t_order - MAX_NATIVE_TIME_ORDER}_{s_part}"
            return OperatorTerm(
                operator=new_op,
                field=aux_name,
                coefficient=term.coefficient,
                coefficient_symbolic=term.coefficient_symbolic,
            )

        if t_order >= _T_ORDER_4:
            # Higher-order mixed: reduce by 2
            new_op = f"mixed_T{t_order - MAX_NATIVE_TIME_ORDER}_{s_part}"
            return OperatorTerm(
                operator=new_op,
                field=aux_name,
                coefficient=term.coefficient,
                coefficient_symbolic=term.coefficient_symbolic,
            )

    # Mixed operators on non-4th-order fields: check if resolvable.
    # For constraint fields (time_order=0), ∂²_t(field) = 0 in the
    # linearized theory → mixed_T2_S*(field) = ∂_x(∂²_t field) = 0.
    # For 2nd-order dynamical fields, mixed_T2_S1x(field) = ∂_x(d²_t field)
    # — this is a genuine RHS term that needs solver support.
    if m and not field_has_aux:
        return term

    # Pure spatial operators: pass through unchanged
    return term


def _substitute_cross_field_refs(
    terms: tuple[OperatorTerm, ...],
    aux_map: dict[str, str],
    rhs_lookup: dict[str, tuple[OperatorTerm, ...]],
) -> tuple[OperatorTerm, ...]:
    """Substitute d2_t/d3_t references to 4th-order fields in non-4th-order equations."""
    from tidal.symbolic.json_loader import OperatorTerm  # noqa: PLC0415

    result: list[OperatorTerm] = []
    for term in terms:
        if term.field in aux_map:
            # This term references a 4th-order field
            if term.operator == "d2_t":
                # d2_t(4th_order_field) → identity(w_field)
                result.append(
                    OperatorTerm(
                        operator="identity",
                        field=aux_map[term.field],
                        coefficient=term.coefficient,
                        coefficient_symbolic=term.coefficient_symbolic,
                    )
                )
            elif term.operator == "d3_t":
                # d3_t(4th_order_field) → first_derivative_t(w_field)
                result.append(
                    OperatorTerm(
                        operator="first_derivative_t",
                        field=aux_map[term.field],
                        coefficient=term.coefficient,
                        coefficient_symbolic=term.coefficient_symbolic,
                    )
                )
            elif term.operator == "d4_t":
                # d4_t(4th_order_field) → d2_t(w_field)
                result.append(
                    OperatorTerm(
                        operator="d2_t",
                        field=aux_map[term.field],
                        coefficient=term.coefficient,
                        coefficient_symbolic=term.coefficient_symbolic,
                    )
                )
            elif _MIXED_RE.match(term.operator):
                # Mixed operator on 4th-order field — reduce time order
                reduced = _reduce_single_term(term, "", aux_map, rhs_lookup)
                if reduced is not None:
                    if isinstance(reduced, list):
                        result.extend(reduced)
                    else:
                        result.append(reduced)
            else:
                # Spatial-only operator on 4th-order field: pass through
                result.append(term)
        else:
            # Non-4th-order field: pass through
            result.append(term)
    return tuple(result)


def _spatial_part_to_operator(s_part: str) -> str:
    """Convert spatial part code to standard operator name.

    Examples: "S1x" → "gradient_x", "S2x" → "laplacian_x"

    Raises
    ------
    ValueError
        If ``s_part`` does not match the expected ``S[order][axis]`` pattern.
    """
    # Parse S[order][axis] pattern
    m = re.match(r"S(\d+)(\w+)", s_part)
    if not m:
        msg = f"Cannot parse spatial part '{s_part}'"
        raise ValueError(msg)

    order = int(m.group(1))
    axis = m.group(2)

    if order == _S_ORDER_GRADIENT:
        return f"gradient_{axis}"
    if order == _S_ORDER_LAPLACIAN:
        return f"laplacian_{axis}"
    # Higher-order spatial: use generic name
    return f"derivative_{order}_{axis}"
