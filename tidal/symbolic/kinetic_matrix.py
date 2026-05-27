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

    ``fields`` is the row/column label list (one entry per
    dynamical equation in the spec). ``cells`` is an
    ``len(fields) x len(fields)`` 2D tuple of
    :class:`KineticMatrixCell` instances.
    """

    fields: tuple[str, ...]
    cells: tuple[tuple[KineticMatrixCell, ...], ...]

    def field_index(self, name: str) -> int | None:
        """Return the matrix index of ``name``, or ``None`` if absent."""
        try:
            return self.fields.index(name)
        except ValueError:
            return None

    def get(self, row: int, col: int) -> KineticMatrixCell:
        return self.cells[row][col]

    @property
    def n(self) -> int:
        return len(self.fields)


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

    The row/column labels are taken from the field of each
    component equation in :attr:`spec.equations`. Each cell
    aggregates the operator-polynomial contributions to that
    (row-equation, column-field) pair:

    - **LHS contribution** lands on the diagonal $\mathcal{K}_{ii}$.
      The LHS reads `kinetic_coeff * d^{time_order}_t phi_i`; we
      encode it as a single ``(d^k_t, +kinetic_coeff)`` entry.
    - **RHS contributions** are negated when moved to the LHS to
      put the EOM in the form $\mathcal{K}\xi = 0$. So each
      RHS term ``(coeff, op, field=j)`` contributes
      ``(op, -coeff)`` to $\mathcal{K}_{ij}$.

    Sign convention follows the LHS=RHS form emitted by the
    Wolfram pipeline; rearranging to $\mathcal{K}\xi = 0$ flips
    the RHS sign.
    """
    field_names = tuple(eq.field_name for eq in spec.equations)
    n = len(field_names)
    field_to_idx: dict[str, int] = {name: i for i, name in enumerate(field_names)}

    # Accumulator: mutable list of entries per (i, j) cell.
    accum: list[list[list[CellEntry]]] = [[[] for _ in range(n)] for _ in range(n)]

    for i, eq in enumerate(spec.equations):
        # LHS → diagonal cell K_ii.
        lhs_op = _LHS_OPERATOR_BY_TIME_ORDER.get(
            eq.time_derivative_order,
            f"d{eq.time_derivative_order}_t",
        )
        lhs_coeff = eq.kinetic_coefficient_symbolic or "1"
        accum[i][i].append((lhs_op, lhs_coeff))

        # RHS → entries K_ij with negated coefficient.
        for term in eq.rhs_terms:
            j = field_to_idx.get(term.field)
            if j is None:
                _log.warning(
                    "kinetic_matrix: equation %r RHS references field %r "
                    "which is not in spec.equations (eliminated or "
                    "constraint-only); skipping term.",
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
    return KineticMatrix(fields=field_names, cells=cells)
