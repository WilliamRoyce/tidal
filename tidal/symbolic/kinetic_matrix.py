r"""Kinetic-matrix assembly from a derived EquationSystem.

Reorganises the per-equation JSON output into the
$\mathcal{K}(\partial_t, \partial_z)$ wave-operator matrix defined in
`manuscript/sections/theory.tex` (label `KineticMatrix`):

.. math::

    \mathcal{L}^{(2)} = \tfrac{1}{2}\,\xi^\mathsf{T}\,
        \mathcal{K}(\partial_t, \partial_z)\,\xi

i.e., each row $i$ and column $j$ of $\mathcal{K}$ is the
differential-operator polynomial acting on $\phi_j$ in the
linearised equation of motion for $\phi_i$. The diagonal $i = j$
collects the LHS kinetic prefactor and any RHS self-terms; the
off-diagonal entries collect inter-field couplings (Gertsenshtein
$h \leftrightarrow A$, dark-photon $A \leftrightarrow T$, etc.).

This module operates purely on the existing JSON `equations[]`
block (no Wolfram changes, no JSON-schema changes, no
re-derivation) — see GitHub issue #372 for the design discussion.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from tidal.symbolic.json_loader import EquationSystem

_log = logging.getLogger(__name__)


# Each kinetic-matrix cell is a list of (operator_label,
# coefficient_symbolic) pairs. The cell value is the sum
# `Σ coeff_k · op_k(\phi_j)` — i.e., a differential-operator
# polynomial. Empty list → zero entry. Coefficients carry their
# own sign (no separate sign field).
CellEntry = tuple[str, str]


@dataclass(frozen=True)
class KineticMatrixCell:
    r"""One entry $\\mathcal{K}_{ij}$ in the kinetic matrix.

    A *cell* aggregates the contributions to the operator
    polynomial acting on column-field $\\phi_j$ in equation $i$:
    each ``(operator_label, coefficient_symbolic)`` pair encodes
    one term ``coeff · op(\\phi_j)``.

    Operator labels follow the same vocabulary the JSON RHS uses
    (``identity``, ``gradient_z``, ``d2_t``, ...); the LaTeX
    renderer reuses
    :func:`tidal.symbolic.latex.operator_to_latex` to map them.
    """

    entries: tuple[CellEntry, ...] = ()

    def is_zero(self) -> bool:
        return len(self.entries) == 0


@dataclass(frozen=True)
class KineticMatrix:
    """Assembled kinetic matrix for an EquationSystem.

    The matrix may be rectangular for theories that carry
    velocity-pair fields (``v_<X>``) on the right-hand side of
    EOMs. Rows are indexed by the q-fields with equations
    (one per :attr:`EquationSystem.equations` entry); columns
    are indexed by the union of q-fields and any v-fields
    referenced on a RHS, with q-fields first.

    Attributes
    ----------
    row_fields
        Field name for each row (length ``n_rows``). Matches
        the q-field of the corresponding equation.
    column_fields
        Field name for each column (length ``n_cols``). The
        first ``n_rows`` entries match ``row_fields`` (q-fields
        in equation order); any trailing entries are v-fields
        ordered by first appearance on a RHS.
    cells
        ``n_rows x n_cols`` grid of :class:`KineticMatrixCell`.
    """

    row_fields: tuple[str, ...]
    column_fields: tuple[str, ...]
    cells: tuple[tuple[KineticMatrixCell, ...], ...]

    def row_index(self, name: str) -> int | None:
        """Return the row index of ``name``, or ``None`` if absent."""
        try:
            return self.row_fields.index(name)
        except ValueError:
            return None

    def column_index(self, name: str) -> int | None:
        """Return the column index of ``name``, or ``None`` if absent."""
        try:
            return self.column_fields.index(name)
        except ValueError:
            return None

    def get(self, row: int, col: int) -> KineticMatrixCell:
        return self.cells[row][col]

    @property
    def n_rows(self) -> int:
        return len(self.row_fields)

    @property
    def n_cols(self) -> int:
        return len(self.column_fields)

    @property
    def is_square(self) -> bool:
        return self.n_rows == self.n_cols

    # Backwards-compatible aliases for the pre-rectangular API.
    # ``fields`` and ``field_index`` previously did double duty
    # for both row and column lookup; map them to the column
    # variants (the wider list).
    @property
    def fields(self) -> tuple[str, ...]:
        return self.column_fields

    def field_index(self, name: str) -> int | None:
        return self.column_index(name)

    @property
    def n(self) -> int:
        return self.n_cols


# Map LHS time-order → operator label that the existing
# `operator_to_latex` knows how to render.
_LHS_OPERATOR_BY_TIME_ORDER: dict[int, str] = {
    0: "identity",
    1: "first_derivative_t",
    2: "d2_t",
    3: "d3_t",
    4: "d4_t",
}


def _negate_coefficient(coeff_symbolic: str) -> str:
    r"""Flip the leading sign of a Mathematica-style coefficient string.

    Used to move RHS terms across the equality when assembling the
    kinetic matrix $\\mathcal{K}\\xi = 0$ from the original
    ``LHS = RHS`` emission. Examples:

    - ``"B0/2"`` → ``"-B0/2"``
    - ``"-kappa^(-2)"`` → ``"kappa^(-2)"``
    - ``"1"`` → ``"-1"``
    - ``"-1"`` → ``"1"``
    """
    s = coeff_symbolic.strip()
    if s.startswith("-"):
        return s[1:].lstrip()
    # Wrap in parens if the body contains '+' / '-' so that '-x+y'
    # negates correctly as '-(x+y)' rather than the lexical
    # '--x+y'.
    if any(ch in s[1:] for ch in "+-"):
        return f"-({s})"
    return f"-{s}"


def _coefficient_string(coefficient: float, coefficient_symbolic: str | None) -> str:
    """Return the coefficient to thread into a kinetic-matrix cell.

    Prefer the symbolic form (carries the parameter names the
    appendix is meant to expose); fall back to the numeric form
    when no symbolic is attached. Numeric ``1.0`` / ``-1.0`` are
    stringified bare (so the leading ``op`` is clean).
    """
    if coefficient_symbolic:
        return coefficient_symbolic
    if abs(coefficient - 1.0) < 1e-12:
        return "1"
    if abs(coefficient + 1.0) < 1e-12:
        return "-1"
    return repr(coefficient)


def build_kinetic_matrix(spec: EquationSystem) -> KineticMatrix:
    r"""Assemble $\mathcal{K}$ from the equation system.

    Rows are the dynamical EOMs in :attr:`spec.equations`
    (one per equation, indexed by ``eq.field_name``). Columns are
    the union of:

    1. the same q-fields as the rows (in equation order), and
    2. any *velocity-pair* fields ``v_<X>`` referenced on a RHS,
       ordered by first appearance.

    The matrix is rectangular when (2) is non-empty: theories
    with first-time-derivative couplings (e.g. the
    $\delta_1$-weighted ``gradient_z(v_X)``
    $= \partial_z\partial_t X$ terms in the R²/nonminimal
    sectors) acquire extra columns labelled $\dot{X}$ in the
    rendered matrix.

    Cell semantics:

    - **LHS contribution** lands on the q-diagonal cell
      $\mathcal{K}_{ii_q}$ where $i_q$ is the column index of
      ``eq.field_name``. The LHS reads
      ``kinetic_coeff * d^{time_order}_t phi_i``; encoded as
      ``(d^k_t, +kinetic_coeff)``.
    - **RHS contributions** are negated when moved to the LHS to
      put the EOM in the form $\mathcal{K}\xi = 0$. Each RHS
      term ``(coeff, op, field=j)`` contributes
      ``(op, -coeff)`` to $\mathcal{K}_{ij}$ where $j$ is the
      column index of ``term.field`` (possibly a v-field).
    """
    row_fields = tuple(eq.field_name for eq in spec.equations)
    q_set = set(row_fields)

    # Pass 1: discover any v-fields referenced on RHS but not
    # already in the row/q list. Order by first appearance.
    extra_columns: list[str] = []
    seen_extras: set[str] = set()
    for eq in spec.equations:
        for term in eq.rhs_terms:
            if term.field in q_set or term.field in seen_extras:
                continue
            extra_columns.append(term.field)
            seen_extras.add(term.field)

    column_fields = row_fields + tuple(extra_columns)
    col_to_idx: dict[str, int] = {name: i for i, name in enumerate(column_fields)}

    # Pass 2: assemble. Accumulator is n_rows × n_cols.
    n_rows = len(row_fields)
    n_cols = len(column_fields)
    accum: list[list[list[CellEntry]]] = [
        [[] for _ in range(n_cols)] for _ in range(n_rows)
    ]

    for i, eq in enumerate(spec.equations):
        # LHS → q-diagonal cell K[i, col_of(eq.field_name)].
        diag_col = col_to_idx[eq.field_name]
        lhs_op = _LHS_OPERATOR_BY_TIME_ORDER.get(
            eq.time_derivative_order,
            f"d{eq.time_derivative_order}_t",
        )
        lhs_coeff = eq.kinetic_coefficient_symbolic or "1"
        accum[i][diag_col].append((lhs_op, lhs_coeff))

        # RHS → entries K[i, j] with negated coefficient.
        for term in eq.rhs_terms:
            j = col_to_idx.get(term.field)
            if j is None:
                # Defensive: should be unreachable after Pass 1.
                _log.warning(
                    "kinetic_matrix: equation %r RHS references field %r "
                    "which is missing from column list; skipping term.",
                    eq.field_name,
                    term.field,
                )
                continue
            coeff_str = _coefficient_string(term.coefficient, term.coefficient_symbolic)
            negated = _negate_coefficient(coeff_str)
            accum[i][j].append((term.operator, negated))

    cells = tuple(
        tuple(KineticMatrixCell(entries=tuple(col)) for col in row) for row in accum
    )
    return KineticMatrix(
        row_fields=row_fields, column_fields=column_fields, cells=cells
    )
