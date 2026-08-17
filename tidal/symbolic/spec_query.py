"""Semantic accessors for equation specifications.

:meth:`~tidal.symbolic.json_loader.EquationSystem.from_dict` gives a *typed*
model of a spec.  This module adds the *semantic* one: the physics-level
questions a reader actually asks, each with a single vetted implementation.

The questions, and the trap in each (GH #401 records one wrong answer per row):

``effective_coefficient``
    What is the coefficient of ``operator(field)`` in an equation?  It is the
    **sum of every matching term**, divided by the LHS kinetic coefficient.
    Components routinely carry more than one matching term — the
    Euler-Heisenberg photons carry a base ``laplacian_x`` *and* a ``B0^2*rho``
    correction — so taking "the" matching term is wrong, and ignoring the
    kinetic coefficient is the single most repeated mistake in this codebase
    (#237, #258, #302, and twice in #401).

``field_families``
    Which components belong together?  Grouped by ``tensor_head`` and
    classified by ``tensor_indices``, never by parsing the numeric suffix.
    Index 0 is the temporal component of a rank-1 field, but rank-3 torsion
    has components like ``t_13 = [2, 0, 2]`` where that reading is meaningless.

``coefficient_provenance``
    Where is this coefficient written down, and do those places agree?  Parts
    of it live across several RHS terms and the LHS; it is *duplicated* in the
    mass/coupling matrices; and a *related but distinct* quantity lives in
    ``canonical.hamiltonian_terms``.  These three relationships are reported
    separately, because presenting the Hamiltonian counterpart as "the same
    coefficient elsewhere" would be its own misreading.

``compare_equations`` / ``diff_systems``
    Did re-derivation change the physics?  Multiplying an equation through by
    a constant — both the kinetic coefficient and every RHS term — changes
    nothing physical.  A naive diff reports such a rescaling as three separate
    "fixes" (the ``gertsenshtein_ungauged`` case).  Flipping only the RHS,
    with the kinetic coefficient unchanged, *is* a real change (the #397
    defect).  The two are distinguished here rather than by eye.

All sign and equality reasoning delegates to :mod:`tidal.symbolic.sign_algebra`,
which answers only when it can prove an answer.  Anything undecided is reported
as such and never guessed.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from tidal.symbolic.sign_algebra import (
    Sign,
    are_equal,
    constant_ratio,
    evaluate_numeric,
    ratio_sign,
    sign_of,
)

if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping

    from tidal.symbolic.json_loader import (
        ComponentEquation,
        EquationSystem,
        HamiltonianTerm,
        OperatorTerm,
    )
    from tidal.symbolic.sign_algebra import SignResult

__all__ = [
    "CoefficientProvenance",
    "ConsistencyCheck",
    "EffectiveCoefficient",
    "EquationComparison",
    "FieldFamily",
    "SystemDiff",
    "coefficient_provenance",
    "compare_equations",
    "diff_systems",
    "effective_coefficient",
    "field_families",
    "matrix_matches_summed_terms",
    "self_terms",
    "terms_for",
]


# --- Effective coefficients ---


def terms_for(
    equation: ComponentEquation,
    field: str,
    operator: str,
) -> tuple[OperatorTerm, ...]:
    """Return **all** RHS terms of *equation* matching *field* and *operator*.

    Plural by design.  A component may carry several terms with the same
    ``(field, operator)`` key — for example a numeric base term and a symbolic
    background correction — and summing them is the caller's whole question.

    Parameters
    ----------
    equation : ComponentEquation
        The equation to search.
    field : str
        Field the operator acts on.
    operator : str
        Operator name, e.g. ``"laplacian_x"``.

    Returns
    -------
    tuple[OperatorTerm, ...]
        Matching terms in their original order; empty when none match.
    """
    return tuple(
        term
        for term in equation.rhs_terms
        if term.field == field and term.operator == operator
    )


def self_terms(
    equation: ComponentEquation,
    operator: str | None = None,
) -> tuple[OperatorTerm, ...]:
    """Return the terms where *equation* acts on its own field.

    Parameters
    ----------
    equation : ComponentEquation
        The equation to search.
    operator : str | None
        Restrict to one operator, or ``None`` for every self-term.

    Returns
    -------
    tuple[OperatorTerm, ...]
        Matching self-terms.
    """
    return tuple(
        term
        for term in equation.rhs_terms
        if term.field == equation.field_name
        and (operator is None or term.operator == operator)
    )


def _term_expression(term: OperatorTerm) -> str:
    """Return a term's coefficient as an expression string."""
    if term.coefficient_symbolic is not None:
        return term.coefficient_symbolic
    return repr(term.coefficient)


@dataclass(frozen=True)
class EffectiveCoefficient:
    """The coefficient of ``operator(field)`` in one equation, fully resolved.

    "Effective" means two things have already been applied: every matching RHS
    term has been summed, and the LHS kinetic coefficient has been divided out.
    The *value* is the primary result; its sign is one derived view of it.

    Attributes
    ----------
    equation_field : str
        Field whose equation this came from.
    field : str
        Field the operator acts on.
    operator : str
        Operator name.
    terms : tuple[OperatorTerm, ...]
        The individual contributing terms, kept so a reader can see the parts.
    numerator : str
        Their sum, as an expression string.
    kinetic : str | None
        The LHS kinetic coefficient, or ``None`` when the LHS is bare (= 1).
    """

    equation_field: str
    field: str
    operator: str
    terms: tuple[OperatorTerm, ...]
    numerator: str
    kinetic: str | None

    @property
    def exists(self) -> bool:
        """Whether any term contributed."""
        return bool(self.terms)

    @property
    def expression(self) -> str:
        """The effective coefficient as a single expression string."""
        if self.kinetic is None:
            return self.numerator
        return f"({self.numerator})/({self.kinetic})"

    @property
    def term_count(self) -> int:
        """How many RHS terms were summed — more than one is common."""
        return len(self.terms)

    def sign(
        self,
        *,
        assume_positive: Iterable[str] | None = None,
        assume_nonzero: Iterable[str] | None = None,
        parameters: Mapping[str, float] | None = None,
    ) -> SignResult:
        """Return the proven sign of this effective coefficient.

        Parameters
        ----------
        assume_positive, assume_nonzero : Iterable[str] | None
            Caller-declared parameter facts; see :func:`sign_algebra.sign_of`.
        parameters : Mapping[str, float] | None
            Optional values, used only for corroboration.

        Returns
        -------
        SignResult
            The verdict and the tactic that decided it.
        """
        return sign_of(
            self.expression,
            assume_positive=assume_positive,
            assume_nonzero=assume_nonzero,
            parameters=parameters,
        )

    def value(self, parameters: Mapping[str, float]) -> float | None:
        """Evaluate numerically at *parameters*, or ``None`` if unresolved."""
        return evaluate_numeric(self.expression, parameters)


def effective_coefficient(
    equation: ComponentEquation,
    field: str,
    operator: str,
) -> EffectiveCoefficient:
    """Return the effective coefficient of ``operator(field)`` in *equation*.

    Sums every matching term and records the kinetic coefficient, so callers
    cannot accidentally use one term of several, or forget the LHS divisor.

    Parameters
    ----------
    equation : ComponentEquation
        Equation to read.
    field : str
        Field the operator acts on.
    operator : str
        Operator name.

    Returns
    -------
    EffectiveCoefficient
        The resolved coefficient; check :attr:`EffectiveCoefficient.exists`.
    """
    matched = terms_for(equation, field, operator)
    if not matched:
        numerator = "0"
    elif len(matched) == 1:
        numerator = _term_expression(matched[0])
    else:
        numerator = "(" + ") + (".join(_term_expression(t) for t in matched) + ")"
    return EffectiveCoefficient(
        equation_field=equation.field_name,
        field=field,
        operator=operator,
        terms=matched,
        numerator=numerator,
        kinetic=equation.kinetic_coefficient_symbolic,
    )


# --- Field families ---


@dataclass(frozen=True)
class FieldFamily:
    """Components sharing a tensor head, classified by their index structure.

    Attributes
    ----------
    head : str
        Tensor head, e.g. ``"a"``, ``"h"``, ``"t"``.
    rank : int
        Tensor rank; 0 when unknown.
    members : tuple[str, ...]
        Component names in spec order.
    indices : dict[str, tuple[int, ...]]
        Each component's index tuple, when known.
    exact : bool
        ``True`` when grouping used exported ``tensor_head`` metadata.  ``False``
        means it fell back to splitting the name on its last underscore, which
        is a guess — 12 of the committed example specs predate the metadata.
    """

    head: str
    rank: int
    members: tuple[str, ...]
    indices: dict[str, tuple[int, ...]]
    exact: bool

    def temporal_slots(self, component: str) -> int | None:
        """Return how many of *component*'s indices are temporal (zero).

        This replaces "index 0 is the temporal component", which holds for a
        rank-1 field but not for rank-3 torsion, where ``t_13`` has indices
        ``(2, 0, 2)``.  Returns ``None`` when index metadata is unavailable.
        """
        idx = self.indices.get(component)
        if idx is None:
            return None
        return sum(1 for i in idx if i == 0)

    def group_by_temporal_slots(self) -> dict[int, tuple[str, ...]]:
        """Group members by their number of temporal indices.

        This is **descriptive metadata**, not a comparison partition.  It is
        the correct replacement for "component 0 is the temporal one" (#401
        row 4): for the photon it recovers ``a_0`` alone against ``a_1..a_3``,
        and for rank-3 torsion it separates the twelve components into three
        classes that the flat ``t_0..t_23`` numbering completely obscures.

        Do **not** restrict sign comparisons to within a group.  Each equation
        is normalised by its own kinetic coefficient, so every evolution
        equation in a family should agree in sign regardless of index
        structure — and the #397 defect is precisely a temporal component
        disagreeing with its spatial siblings, which a within-group comparison
        would never look at.  Measured on the committed corpus, partitioning by
        temporal slots misses 13 of the 19 files carrying that defect.
        """
        grouped: dict[int, list[str]] = {}
        for member in self.members:
            slots = self.temporal_slots(member)
            if slots is None:
                continue
            grouped.setdefault(slots, []).append(member)
        return {k: tuple(v) for k, v in sorted(grouped.items())}


def field_families(spec: EquationSystem) -> tuple[FieldFamily, ...]:
    """Group a spec's components into tensor families.

    Uses the ``tensor_head`` / ``tensor_rank`` / ``tensor_indices`` metadata
    exported since ``78374c1``.  Specs predating it fall back to splitting the
    component name on its final underscore, and the resulting families are
    marked ``exact=False`` so callers can tell a guess from a fact.

    Parameters
    ----------
    spec : EquationSystem
        The loaded equation system.

    Returns
    -------
    tuple[FieldFamily, ...]
        Families in order of first appearance.
    """
    metadata: dict[str, dict[str, Any]] = spec.metadata.get("tensor_metadata", {})

    ordered: dict[str, list[str]] = {}
    ranks: dict[str, int] = {}
    indices: dict[str, dict[str, tuple[int, ...]]] = {}
    exactness: dict[str, bool] = {}

    for name in spec.component_names:
        entry = metadata.get(name)
        if entry is not None:
            head = str(entry["tensor_head"])
            exact = True
        else:
            # Fallback: strip a trailing numeric suffix. A guess, flagged as one.
            head = name.rsplit("_", maxsplit=1)[0] if "_" in name else name
            exact = False
        ordered.setdefault(head, []).append(name)
        exactness[head] = exactness.get(head, True) and exact
        if entry is not None:
            ranks[head] = int(entry.get("tensor_rank", 0))
            indices.setdefault(head, {})[name] = tuple(
                int(i) for i in entry.get("tensor_indices", ())
            )

    return tuple(
        FieldFamily(
            head=head,
            rank=ranks.get(head, 0),
            members=tuple(members),
            indices=indices.get(head, {}),
            exact=exactness.get(head, False),
        )
        for head, members in ordered.items()
    )


# --- Coefficient provenance ---


@dataclass(frozen=True)
class ConsistencyCheck:
    """One cross-representation check, and whether it could be settled.

    Attributes
    ----------
    name : str
        Short identifier, e.g. ``"numeric-vs-symbolic"``.
    status : str
        ``"ok"``, ``"mismatch"``, or ``"undecided"``.  ``"undecided"`` is a
        first-class outcome, not a failure.
    detail : str
        Human-readable explanation.
    """

    name: str
    status: str
    detail: str


@dataclass(frozen=True)
class CoefficientProvenance:
    """Every place one coefficient is recorded, grouped by *how* it relates.

    The three groups are kept apart deliberately.  Only :attr:`effective` and
    its parts are this coefficient; :attr:`matrix_entry` is a redundant
    re-encoding of it; :attr:`hamiltonian_terms` is a *different quantity* in a
    different formalism, related by a derivation and carrying its own factor
    convention.  Flattening them into one list would invite the very confusion
    this module exists to prevent.

    Attributes
    ----------
    effective : EffectiveCoefficient
        The coefficient itself, summed and kinetic-normalised.
    order_spread : dict[int, tuple[str, ...]]
        Contributing term expressions grouped by ``order_in_eps``.
    matrix_entry : str | float | None
        The mass/coupling matrix encoding, for ``identity`` operators.  Note
        the matrix convention is ``matrix[i][j] = -(coefficient)`` and is
        **not** normalised by the kinetic coefficient.
    hamiltonian_terms : tuple[HamiltonianTerm, ...]
        Hamiltonian terms mentioning the same field pair — related, distinct.
    checks : tuple[ConsistencyCheck, ...]
        Only checks that can be settled without unproven factor reasoning.
    """

    effective: EffectiveCoefficient
    order_spread: dict[int, tuple[str, ...]]
    matrix_entry: str | float | None
    hamiltonian_terms: tuple[HamiltonianTerm, ...]
    checks: tuple[ConsistencyCheck, ...]


def _matrix_entry(
    spec: EquationSystem,
    equation: ComponentEquation,
    field: str,
    operator: str,
) -> str | float | None:
    """Return the mass/coupling matrix encoding of an identity coefficient."""
    if operator != "identity":
        return None
    try:
        row = spec.equation_map[equation.field_name]
        col = spec.component_names.index(field)
    except (KeyError, ValueError):
        return None
    if spec.mass_matrix_symbolic:
        symbolic = spec.mass_matrix_symbolic[row][col]
        if symbolic is not None:
            return symbolic
    if row < len(spec.mass_matrix) and col < len(spec.mass_matrix[row]):
        return spec.mass_matrix[row][col]
    return None


def _hamiltonian_terms_for(
    spec: EquationSystem,
    equation_field: str,
    field: str,
) -> tuple[HamiltonianTerm, ...]:
    """Return Hamiltonian terms coupling *equation_field* and *field*."""
    if spec.canonical is None:
        return ()
    wanted = {equation_field, field}
    return tuple(
        term
        for term in spec.canonical.hamiltonian_terms
        if {term.factor_a.field, term.factor_b.field} & wanted == wanted
        or (
            equation_field == field
            and term.factor_a.field == field
            and term.factor_b.field == field
        )
    )


def _numeric_symbolic_checks(
    terms: tuple[OperatorTerm, ...],
) -> list[ConsistencyCheck]:
    """Check each term's numeric coefficient against its symbolic form.

    The exporter's convention is that the numeric value equals the symbolic
    expression evaluated with every free name set to 1.0, so a disagreement is
    detectable without knowing any parameter.
    """
    checks: list[ConsistencyCheck] = []
    for term in terms:
        if term.coefficient_symbolic is None:
            continue
        from tidal.symbolic.sign_algebra import free_names  # noqa: PLC0415

        names = free_names(term.coefficient_symbolic)
        value = evaluate_numeric(term.coefficient_symbolic, dict.fromkeys(names, 1.0))
        if value is None:
            checks.append(
                ConsistencyCheck(
                    name="numeric-vs-symbolic",
                    status="undecided",
                    detail=f"{term.coefficient_symbolic!r} did not evaluate at all-ones",
                ),
            )
        elif not _close(value, term.coefficient):
            checks.append(
                ConsistencyCheck(
                    name="numeric-vs-symbolic",
                    status="mismatch",
                    detail=(
                        f"{term.coefficient_symbolic!r} is {value:g} at all-ones "
                        f"but the stored numeric coefficient is {term.coefficient:g}"
                    ),
                ),
            )
        else:
            checks.append(
                ConsistencyCheck(
                    name="numeric-vs-symbolic",
                    status="ok",
                    detail=f"{term.coefficient_symbolic!r} agrees at all-ones",
                ),
            )
    return checks


def _close(left: float, right: float, tol: float = 1e-9) -> bool:
    """Whether two floats agree to a relative tolerance."""
    return abs(left - right) <= tol * max(1.0, abs(left), abs(right))


def coefficient_provenance(
    spec: EquationSystem,
    equation_field: str,
    field: str,
    operator: str,
) -> CoefficientProvenance:
    """Gather every recorded form of one coefficient, with its relationships.

    Parameters
    ----------
    spec : EquationSystem
        The loaded system.
    equation_field : str
        Which equation to read.
    field : str
        Field the operator acts on.
    operator : str
        Operator name.

    Returns
    -------
    CoefficientProvenance
        The parts, the duplicate encodings, and the related-but-distinct
        Hamiltonian terms, kept separate.

    Raises
    ------
    KeyError
        If *equation_field* is not a component of *spec*.
    """
    index = spec.equation_map.get(equation_field)
    if index is None:
        msg = f"unknown equation field {equation_field!r}"
        raise KeyError(msg)
    equation = spec.equations[index]

    effective = effective_coefficient(equation, field, operator)

    order_spread: dict[int, list[str]] = {}
    for term in effective.terms:
        order_spread.setdefault(term.order_in_eps, []).append(_term_expression(term))

    return CoefficientProvenance(
        effective=effective,
        order_spread={k: tuple(v) for k, v in sorted(order_spread.items())},
        matrix_entry=_matrix_entry(spec, equation, field, operator),
        hamiltonian_terms=_hamiltonian_terms_for(spec, equation_field, field),
        checks=tuple(_numeric_symbolic_checks(effective.terms)),
    )


# --- Diffing two specs ---


@dataclass(frozen=True)
class EquationComparison:
    """How one equation differs between two specs.

    Attributes
    ----------
    field : str
        The component compared.
    verdict : str
        ``"identical"`` — byte-equal effective coefficients.
        ``"representational"`` — every effective coefficient is unchanged
        even though the written form differs, i.e. the equation was rescaled.
        ``"real"`` — at least one effective coefficient provably changed.
        ``"undecided"`` — a difference exists that could not be settled.
    changed_keys : tuple[str, ...]
        ``operator(field)`` keys whose effective coefficient changed.
    undecided_keys : tuple[str, ...]
        Keys that could neither be proven equal nor proven different.
    detail : str
        Human-readable summary.
    """

    field: str
    verdict: str
    changed_keys: tuple[str, ...]
    undecided_keys: tuple[str, ...]
    detail: str


def _term_keys(equation: ComponentEquation) -> set[tuple[str, str]]:
    """Return the distinct ``(operator, field)`` keys on an equation's RHS."""
    return {(term.operator, term.field) for term in equation.rhs_terms}


def compare_equations(
    left: ComponentEquation,
    right: ComponentEquation,
) -> EquationComparison:
    """Compare two versions of one equation, separating real from cosmetic change.

    An equation multiplied through by any non-zero constant — the kinetic
    coefficient *and* every RHS term — is physically unchanged.  Comparing
    effective coefficients makes that rescaling invisible automatically, while
    a flip of the RHS alone still shows up as real.

    Parameters
    ----------
    left, right : ComponentEquation
        The two versions.

    Returns
    -------
    EquationComparison
        The verdict and the keys responsible for it.
    """
    keys = sorted(_term_keys(left) | _term_keys(right))
    changed: list[str] = []
    undecided: list[str] = []

    for operator, field in keys:
        left_eff = effective_coefficient(left, field, operator)
        right_eff = effective_coefficient(right, field, operator)
        equal = are_equal(left_eff.expression, right_eff.expression)
        label = f"{operator}({field})"
        if equal is True:
            continue
        if equal is False:
            changed.append(label)
        else:
            undecided.append(label)

    if changed:
        verdict = "real"
        detail = f"{len(changed)} effective coefficient(s) changed"
    elif undecided:
        verdict = "undecided"
        detail = f"{len(undecided)} coefficient(s) could not be compared"
    elif left.rhs_terms == right.rhs_terms and (
        left.kinetic_coefficient_symbolic == right.kinetic_coefficient_symbolic
    ):
        verdict = "identical"
        detail = "byte-identical"
    else:
        verdict = "representational"
        detail = (
            "written differently but every effective coefficient is unchanged "
            "(an overall rescaling of the equation)"
        )

    return EquationComparison(
        field=left.field_name,
        verdict=verdict,
        changed_keys=tuple(changed),
        undecided_keys=tuple(undecided),
        detail=detail,
    )


@dataclass(frozen=True)
class SystemDiff:
    """Result of comparing two equation systems.

    Attributes
    ----------
    comparisons : tuple[EquationComparison, ...]
        Per-equation results for fields present in both systems.
    only_left, only_right : tuple[str, ...]
        Components present in just one system.
    """

    comparisons: tuple[EquationComparison, ...]
    only_left: tuple[str, ...]
    only_right: tuple[str, ...]

    @property
    def real(self) -> tuple[EquationComparison, ...]:
        """Equations whose physics provably changed."""
        return tuple(c for c in self.comparisons if c.verdict == "real")

    @property
    def representational(self) -> tuple[EquationComparison, ...]:
        """Equations rewritten without changing their physics."""
        return tuple(c for c in self.comparisons if c.verdict == "representational")

    @property
    def undecided(self) -> tuple[EquationComparison, ...]:
        """Equations whose difference could not be settled."""
        return tuple(c for c in self.comparisons if c.verdict == "undecided")

    @property
    def has_real_changes(self) -> bool:
        """Whether any equation provably changed, or a component appeared/vanished."""
        return bool(self.real or self.only_left or self.only_right)


def diff_systems(left: EquationSystem, right: EquationSystem) -> SystemDiff:
    """Compare two equation systems, separating real changes from rescalings.

    Parameters
    ----------
    left, right : EquationSystem
        The two systems, e.g. a committed spec and a re-derived one.

    Returns
    -------
    SystemDiff
        Per-equation verdicts plus components unique to either side.
    """
    left_map = left.equation_map
    right_map = right.equation_map
    shared = [name for name in left.component_names if name in right_map]

    comparisons = tuple(
        compare_equations(
            left.equations[left_map[name]],
            right.equations[right_map[name]],
        )
        for name in shared
    )
    return SystemDiff(
        comparisons=comparisons,
        only_left=tuple(n for n in left.component_names if n not in right_map),
        only_right=tuple(n for n in right.component_names if n not in left_map),
    )


# --- Sibling sign consistency (the #397 invariant) ---

# An equation is a genuine evolution equation from second time order upward;
# below that it is a constraint, whose overall sign is conventional.
_EVOLUTION_TIME_ORDER = 2


def sibling_sign_conflicts(
    spec: EquationSystem,
    operator: str = "laplacian_x",
    *,
    assume_positive: Iterable[str] | None = None,
    assume_nonzero: Iterable[str] | None = None,
) -> tuple[tuple[str, str, str], ...]:
    """Find components whose self-``operator`` sign provably opposes a sibling's.

    Compares every component in a tensor family against the first, using
    *effective* coefficients so each equation's own kinetic coefficient is
    divided out first — which is what makes components with different index
    structure legitimately comparable.

    Two guards matter, and only two:

    * **Evolution equations only.** A constraint (``time_order == 0``) carries
      no kinetic term, so an overall sign on it is conventional rather than
      physical; #397 excludes those by hand and this does it structurally.
      For rank-3 torsion that is exactly what drops ``t_6``, ``t_13`` and
      ``t_20``.
    * **Proven conflicts only.** Pairs the sign algebra cannot decide are
      skipped rather than guessed at.

    Note the comparison deliberately spans the *whole* family rather than
    grouping by index structure: the #397 defect is a temporal component
    disagreeing with its spatial siblings, so restricting to like-indexed
    components would miss it — measurably, 13 of the 19 affected files.

    Parameters
    ----------
    spec : EquationSystem
        The system to check.
    operator : str
        Self-operator to compare, defaulting to the spatial laplacian.
    assume_positive, assume_nonzero : Iterable[str] | None
        Caller-declared parameter facts.

    Returns
    -------
    tuple[tuple[str, str, str], ...]
        ``(family_head, reference_component, conflicting_component)`` triples.
    """
    conflicts: list[tuple[str, str, str]] = []
    equations = {eq.field_name: eq for eq in spec.equations}

    for family in field_families(spec):
        reference: tuple[str, EffectiveCoefficient] | None = None
        for member in family.members:
            equation = equations.get(member)
            # Constraints carry no kinetic term, so their overall sign is
            # conventional; comparing them would manufacture false conflicts.
            if (
                equation is None
                or equation.time_derivative_order < _EVOLUTION_TIME_ORDER
            ):
                continue
            eff = effective_coefficient(equation, member, operator)
            if not eff.exists:
                continue
            if reference is None:
                reference = (member, eff)
                continue
            # Cross-multiplied so no division is performed: the ratio of the
            # two effective coefficients is (num_a * kin_b) / (num_b * kin_a).
            verdict = ratio_sign(
                f"({eff.numerator})*({reference[1].kinetic or 1})",
                f"({reference[1].numerator})*({eff.kinetic or 1})",
                assume_positive=assume_positive,
                assume_nonzero=assume_nonzero,
            )
            if verdict.sign is Sign.NEGATIVE:
                conflicts.append((family.head, reference[0], member))
    return tuple(conflicts)


def matrix_matches_summed_terms(spec: EquationSystem) -> tuple[ConsistencyCheck, ...]:
    """Cross-check the derived mass matrix against the summed identity terms.

    Two Python derivations of the same quantity are compared:
    :meth:`EquationSystem._compute_matrices_from_terms`, which builds
    ``mass_matrix_symbolic`` while scanning the equations, and
    :func:`effective_coefficient`, which sums the matching terms on demand.
    They should never disagree, and when they did it was because the matrix
    builder overwrote multi-term coefficients instead of accumulating them
    (GH #403).

    This reads the **derived** ``spec.mass_matrix_symbolic``, not the JSON: the
    loader has not read a stored ``coupling`` section since c407240d, and the
    exporter no longer writes one, so there is no stored encoding left to
    compare against.

    The symbolic matrix stores ``coefficient_symbolic`` verbatim while the
    numeric one negates it, so the comparison here is against the un-negated
    sum.

    Parameters
    ----------
    spec : EquationSystem
        The system to check.

    Returns
    -------
    tuple[ConsistencyCheck, ...]
        One check per disagreeing entry, empty when all agree.
    """
    if not spec.mass_matrix_symbolic:
        return ()
    problems: list[ConsistencyCheck] = []
    for equation in spec.equations:
        row = spec.equation_map[equation.field_name]
        for col, field in enumerate(spec.component_names):
            eff = effective_coefficient(equation, field, "identity")
            derived = spec.mass_matrix_symbolic[row][col]
            if derived is None or not eff.exists:
                continue
            if constant_ratio(derived, eff.numerator) != 1:
                problems.append(
                    ConsistencyCheck(
                        name="matrix-vs-terms",
                        status="undecided",
                        detail=(
                            f"{equation.field_name}/{field}: derived matrix entry "
                            f"{derived!r} is not provably equal to the summed "
                            f"identity term {eff.numerator!r}"
                        ),
                    ),
                )
    return tuple(problems)
