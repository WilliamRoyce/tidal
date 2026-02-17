"""Build py-pde PDEBase subclasses from equation specifications.

This module provides the core functionality for converting symbolically-derived
field equations (loaded from JSON) into executable py-pde PDE classes.

The key principle is that NO physics is hardcoded here - all equation structure
comes from the specification that was derived from the Lagrangian.
"""

from __future__ import annotations

import logging
import math
import operator
import re
import warnings
from dataclasses import dataclass
from functools import partial
from typing import TYPE_CHECKING, Any, ClassVar, SupportsFloat, cast

import numpy as np
from pde import FieldCollection, PDEBase, ScalarField
from scipy import (  # type: ignore[reportMissingTypeStubs]
    sparse,
    special,
)
from scipy.sparse.linalg import spsolve  # type: ignore[reportUnknownVariableType]
from typing_extensions import override

from tidal.symbolic.json_loader import (
    _CUSTOM_OPERATORS,  # type: ignore[reportPrivateUsage]
    AXIS_LETTERS,
    BoundaryCondition,
    ComponentEquation,
    ConstraintSolverConfig,
    EquationSystem,
    OperatorTerm,
    load_equation_system,
)
from tidal.utils import infer_bc_from_grid

if TYPE_CHECKING:
    from collections.abc import Callable
    from pathlib import Path

    from numpy.typing import NDArray
    from pde.grids.base import GridBase
    from pde.pdes.base import TState

    from tidal.utils import BCDescriptor

    NumericArray = NDArray[np.float64]

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Pre-compiled regex patterns for field name parsing (avoid re-compilation
# on every call to ParsedFieldName.parse / parse_momentum_field_name)
# ---------------------------------------------------------------------------
_STANDARD_FORMAT_RE = re.compile(r"^([a-zA-Z]+)_([0-9]+)$")
_TENSOR_FORMAT_RE = re.compile(r"^(.+)_([0-9]+)$")
_COMPACT_FORMAT_RE = re.compile(r"^([a-zA-Z]+)([0-9]+)$")
_SIMPLE_FORMAT_RE = re.compile(r"^[a-zA-Z]+$")
_MOMENTUM_STANDARD_RE = re.compile(r"^pi_([0-9]+)$")
_MOMENTUM_COMPACT_RE = re.compile(r"^pi([0-9]+)$")


# ---------------------------------------------------------------------------
# Operator registry: maps operator name -> (handler, min_grid_dimension)
# Each handler takes (field: ScalarField, bc: BCDescriptor) -> ScalarField
# ---------------------------------------------------------------------------


def _op_laplacian(field: ScalarField, bc: BCDescriptor) -> ScalarField:
    return field.laplace(bc=bc)


def _op_identity(field: ScalarField, _bc: BCDescriptor) -> ScalarField:
    return field.copy()


def _op_gradient(axis: int) -> Callable[[ScalarField, BCDescriptor], ScalarField]:
    """Create a gradient handler for a specific axis."""

    def _handler(field: ScalarField, bc: BCDescriptor) -> ScalarField:
        grad = field.gradient(bc=bc)
        component = grad[axis]
        if not isinstance(component, ScalarField):  # pyright: ignore[reportUnnecessaryIsInstance]
            msg = f"Expected ScalarField from gradient, got {type(component).__name__}"
            raise TypeError(msg)
        return component

    return _handler


def _op_directional_laplacian(
    axis: int,
) -> Callable[[ScalarField, BCDescriptor], ScalarField]:
    """Create a directional Laplacian handler (∂²/∂x_i²)."""

    def _handler(field: ScalarField, bc: BCDescriptor) -> ScalarField:
        grad = field.gradient(bc=bc)[axis]
        if not isinstance(grad, ScalarField):  # pyright: ignore[reportUnnecessaryIsInstance]
            msg = f"Expected ScalarField from gradient, got {type(grad).__name__}"
            raise TypeError(msg)
        d2 = grad.gradient(bc=bc)[axis]
        if not isinstance(d2, ScalarField):  # pyright: ignore[reportUnnecessaryIsInstance]
            msg = f"Expected ScalarField from gradient, got {type(d2).__name__}"
            raise TypeError(msg)
        return d2

    return _handler


def _op_cross_derivative(
    axis1: int, axis2: int
) -> Callable[[ScalarField, BCDescriptor], ScalarField]:
    """Create a cross derivative handler (∂²/∂x_i ∂x_j)."""

    def _handler(field: ScalarField, bc: BCDescriptor) -> ScalarField:
        grad_j = field.gradient(bc=bc)[axis2]
        if not isinstance(grad_j, ScalarField):  # pyright: ignore[reportUnnecessaryIsInstance]
            msg = f"Expected ScalarField from gradient, got {type(grad_j).__name__}"
            raise TypeError(msg)
        grad_ij = grad_j.gradient(bc=bc)[axis1]
        if not isinstance(grad_ij, ScalarField):  # pyright: ignore[reportUnnecessaryIsInstance]
            msg = f"Expected ScalarField from gradient, got {type(grad_ij).__name__}"
            raise TypeError(msg)
        return grad_ij

    return _handler


def _op_biharmonic(field: ScalarField, bc: BCDescriptor) -> ScalarField:
    """Biharmonic operator: ∇⁴f = ∇²(∇²f).

    Raises
    ------
    TypeError
        If intermediate results are not ``ScalarField``.
    """
    lap = field.laplace(bc=bc)
    if not isinstance(lap, ScalarField):  # pyright: ignore[reportUnnecessaryIsInstance]
        msg = f"Expected ScalarField from laplace, got {type(lap).__name__}"
        raise TypeError(msg)
    bilap = lap.laplace(bc=bc)
    if not isinstance(bilap, ScalarField):  # pyright: ignore[reportUnnecessaryIsInstance]
        msg = f"Expected ScalarField from laplace, got {type(bilap).__name__}"
        raise TypeError(msg)
    return bilap


def _op_nth_derivative(
    axis: int, order: int
) -> Callable[[ScalarField, BCDescriptor], ScalarField]:
    """Create a handler for Nth-order derivative along a single axis.

    Applies gradient in the given axis direction ``order`` times.
    E.g., order=3, axis=0 computes ∂³f/∂x³.
    """

    def _handler(field: ScalarField, bc: BCDescriptor) -> ScalarField:
        result: ScalarField = field
        for _ in range(order):
            grad = result.gradient(bc=bc)
            component = grad[axis]
            if not isinstance(component, ScalarField):  # pyright: ignore[reportUnnecessaryIsInstance]
                msg = f"Expected ScalarField from gradient, got {type(component).__name__}"
                raise TypeError(msg)
            result = component
        return result

    return _handler


def _parse_multi_axis_spec(spec: str) -> list[tuple[int, int]]:
    """Parse a multi-axis derivative spec into (axis_index, order) pairs.

    Parameters
    ----------
    spec : str
        Multi-axis spec like ``"2x_1y"``, ``"3x_2z"``, ``"1x_1y_1z"``.

    Returns
    -------
    list[tuple[int, int]]
        List of (axis_index, order) pairs, sorted by axis index.

    Raises
    ------
    ValueError
        If the spec cannot be parsed or contains invalid axis names.
    """
    axis_pattern = re.compile(r"^(\d+)(" + _AXIS_RE_CLASS + r")$")
    parts = spec.split("_")
    result: list[tuple[int, int]] = []
    for part in parts:
        match = axis_pattern.match(part)
        if not match:
            msg = (
                f"Invalid multi-axis derivative part: '{part}'. "
                f"Expected format like '2x', '1y', '3z'."
            )
            raise ValueError(msg)
        order = int(match.group(1))
        axis_letter = match.group(2)
        axis = _AXIS_INDEX[axis_letter]
        result.append((axis, order))

    # Detect duplicate axes — e.g. "2x_1x_1y" has x twice
    seen_axes: set[int] = set()
    for axis, _order in result:
        if axis in seen_axes:
            axis_name = _AXIS_LETTER[axis]
            msg = (
                f"Duplicate axis '{axis_name}' in multi-axis derivative spec: '{spec}'. "
                f"Combine orders for the same axis into a single entry."
            )
            raise ValueError(msg)
        seen_axes.add(axis)

    return sorted(result, key=operator.itemgetter(0))


def _op_multi_axis_derivative(
    axes_and_orders: list[tuple[int, int]],
) -> Callable[[ScalarField, BCDescriptor], ScalarField]:
    """Create a handler for multi-axis mixed spatial derivatives.

    Applies gradients sequentially for each (axis, order) pair.
    E.g., for ``derivative_2x_1y``: applies gradient in y once, then
    gradient in x twice, giving d^2/dx^2 d/dy.
    """

    def _handler(field: ScalarField, bc: BCDescriptor) -> ScalarField:
        result: ScalarField = field
        for axis, order in axes_and_orders:
            for _ in range(order):
                grad = result.gradient(bc=bc)
                component = grad[axis]
                if not isinstance(component, ScalarField):  # pyright: ignore[reportUnnecessaryIsInstance]
                    msg = f"Expected ScalarField from gradient, got {type(component).__name__}"
                    raise TypeError(msg)
                result = component
        return result

    return _handler


#: Map axis letter to axis index (derived from AXIS_LETTERS in json_loader).
_AXIS_INDEX: dict[str, int] = {letter: i for i, letter in enumerate(AXIS_LETTERS)}

#: Minimum grid dimension required for each axis.
_AXIS_MIN_DIM: dict[str, int] = {letter: i + 1 for i, letter in enumerate(AXIS_LETTERS)}

#: Map axis index back to letter.
_AXIS_LETTER: dict[int, str] = dict(enumerate(AXIS_LETTERS))

#: Character class matching all known axis letters.
_AXIS_RE_CLASS = "[" + "".join(AXIS_LETTERS) + "]"

#: Regex for parsing generic single-axis derivative names.
_GENERIC_SINGLE_RE = re.compile(r"^derivative_(\d+)_(" + _AXIS_RE_CLASS + r")$")

#: Regex for parsing generic multi-axis derivative names.
_GENERIC_MULTI_RE = re.compile(
    r"^derivative_(\d+" + _AXIS_RE_CLASS + r"(?:_\d+" + _AXIS_RE_CLASS + r")*)$"
)


#: Registry mapping operator names to (handler, min_dimension) pairs.
#: To add a new operator, simply add an entry here.
_OPERATOR_REGISTRY: dict[str, tuple[Any, int]] = {
    "laplacian": (_op_laplacian, 1),
    "identity": (_op_identity, 1),
    "gradient_x": (_op_gradient(0), 1),
    "gradient_y": (_op_gradient(1), 2),
    "gradient_z": (_op_gradient(2), 3),
    "laplacian_x": (_op_directional_laplacian(0), 1),
    "laplacian_y": (_op_directional_laplacian(1), 2),
    "laplacian_z": (_op_directional_laplacian(2), 3),
    "cross_derivative_xy": (_op_cross_derivative(0, 1), 2),
    "cross_derivative_xz": (_op_cross_derivative(0, 2), 3),
    "cross_derivative_yz": (_op_cross_derivative(1, 2), 3),
    "biharmonic": (_op_biharmonic, 1),
    # first_derivative_t is handled specially in _compute_rhs_for_component.
    # This sentinel entry ensures _get_operator gives a clear error if called directly.
    "first_derivative_t": (None, 1),
}


def register_operator(
    name: str,
    handler: Callable[[ScalarField, BCDescriptor], ScalarField],
    min_dim: int = 1,
) -> None:
    """Register a custom spatial operator for use in the PDE pipeline.

    The operator becomes available both for JSON validation (``is_known_operator``)
    and for runtime evaluation (``PDEFromSpec._get_operator``).

    Parameters
    ----------
    name : str
        Operator name as it appears in JSON (e.g., ``"my_diffusion"``).
    handler : callable
        Function with signature ``(field: ScalarField, bc) -> ScalarField``.
    min_dim : int
        Minimum spatial grid dimension required for this operator.

    Raises
    ------
    ValueError
        If *name* would shadow a built-in operator.
    """
    if name in _OPERATOR_REGISTRY:
        msg = (
            f"Cannot register operator '{name}': "
            f"it shadows a built-in operator. Choose a different name."
        )
        raise ValueError(msg)
    _OPERATOR_REGISTRY[name] = (handler, min_dim)
    _CUSTOM_OPERATORS.add(name)


# ---------------------------------------------------------------------------
# Operator matrix builders: sparse matrix representations for constraint solving
#
# Each builder returns (matrix, vector) where:
#   matrix: scipy.sparse matrix (N x N) representing the linear operator
#   vector: scipy.sparse matrix (N x 1) for boundary condition offsets
#
# The convention matches py-pde's _get_laplace_matrix:
#   operator(field) = matrix @ field_flat + vector_flat
#
# To solve a constraint equation with arbitrary self-referencing operators:
#   0 = sum_i(c_i * Op_i(field)) + source
# we assemble A = sum_i(c_i * M_i) and solve A @ field = -source.
#
# This unified approach handles Poisson, Helmholtz, algebraic, anisotropic,
# and any future linear constraint type with a single code path.
# ---------------------------------------------------------------------------

# Type alias for matrix builder return type
_MatrixAndVector = tuple[sparse.spmatrix, sparse.spmatrix]


def _build_identity_matrix(
    grid: Any,  # noqa: ANN401
    bcs: Any,  # noqa: ANN401, ARG001
) -> _MatrixAndVector:
    """Build sparse identity matrix for the ``identity`` operator.

    The identity operator I satisfies I(field) = field, so its matrix
    representation is simply the identity matrix with zero BC offset.
    """
    n = math.prod(grid.shape)
    return sparse.eye(n, format="dok"), sparse.dok_matrix((n, 1))  # type: ignore[reportReturnType]


def _build_laplacian_matrix(
    grid: Any,  # noqa: ANN401, ARG001
    bcs: Any,  # noqa: ANN401
) -> _MatrixAndVector:
    """Build sparse Laplacian matrix by wrapping py-pde's internal builder.

    Reuses py-pde's ``_get_laplace_matrix`` which handles all grid dimensions
    (1D, 2D, 3D) and boundary condition types (periodic, Dirichlet, Neumann).

    Raises
    ------
    ImportError
        If py-pde's internal ``_get_laplace_matrix`` is not available.
    """
    try:
        from pde.grids.operators.cartesian import (  # noqa: PLC0415
            _get_laplace_matrix,  # noqa: PLC2701  # type: ignore[reportPrivateUsage]
        )
    except ImportError as e:
        msg = (
            "Cannot import py-pde's _get_laplace_matrix. "
            "The constraint matrix solver requires py-pde >= 0.30."
        )
        raise ImportError(msg) from e
    return _get_laplace_matrix(bcs)  # type: ignore[reportReturnType]


def _build_directional_laplacian_matrix(
    grid: Any,  # noqa: ANN401
    bcs: Any,  # noqa: ANN401
    *,
    axis: int,
) -> _MatrixAndVector:
    """Build sparse matrix for directional Laplacian ∂²/∂x_i².

    Computed as G_i @ G_i (composition of two gradient matrices along the
    same axis). This matches the function-based operator
    ``_op_directional_laplacian`` which applies ``field.gradient(bc)[axis]``
    twice, using py-pde's central-difference gradient stencil.

    .. note::

       The resulting wide stencil [1/(4dx²), 0, -2/(4dx²), 0, 1/(4dx²)]
       differs from the compact 3-point Laplacian stencil [1/dx², -2/dx²,
       1/dx²] used by py-pde's ``field.laplace(bc)``.  Consequently,
       ``laplacian_x_matrix + laplacian_y_matrix != laplacian_matrix``.
       This is intentional — each matrix matches its corresponding
       function-based operator exactly.
    """
    g_mat, g_vec = _build_gradient_matrix(grid, bcs, axis=axis)
    # dir_lap = G @ G: apply gradient twice along the same axis
    # G @ (G @ field + vec) + vec = G² @ field + G @ vec + vec
    dir_lap_matrix = g_mat @ g_mat  # type: ignore[reportOperatorIssue]
    dir_lap_vector = g_mat @ g_vec + g_vec  # type: ignore[reportOperatorIssue]
    return dir_lap_matrix, dir_lap_vector  # type: ignore[reportUnknownVariableType]


def _build_gradient_matrix(
    grid: Any,  # noqa: ANN401
    bcs: Any,  # noqa: ANN401
    *,
    axis: int,
) -> _MatrixAndVector:
    """Build sparse matrix for gradient ∂/∂x_i using central differences.

    Uses the standard central difference stencil [-1/(2dx), 0, 1/(2dx)]
    along the target axis. Boundary conditions modify boundary rows.
    """
    shape = grid.shape
    n = math.prod(shape)
    dx = grid.discretization[axis]
    scale = 1.0 / (2.0 * dx)

    matrix = sparse.dok_matrix((n, n))  # type: ignore[reportUnknownArgumentType]
    vector = sparse.dok_matrix((n, 1))  # type: ignore[reportUnknownArgumentType]

    bc_axis = bcs[axis]

    for flat_idx in range(n):
        multi_idx = list(np.unravel_index(flat_idx, shape))

        # Left neighbor: -1/(2dx)
        if multi_idx[axis] == 0:
            bc_idx = list(multi_idx)
            bc_idx[axis] = -1  # type: ignore[reportCallIssue,reportArgumentType]
            const, entries = bc_axis.get_sparse_matrix_data(tuple(bc_idx))
            vector[flat_idx, 0] += -const * scale
            for k, v in entries.items():
                neighbor_idx = list(multi_idx)
                neighbor_idx[axis] = k
                flat_neighbor = int(np.ravel_multi_index(neighbor_idx, shape))
                matrix[flat_idx, flat_neighbor] += -v * scale
        else:
            neighbor_idx = list(multi_idx)
            neighbor_idx[axis] -= 1
            flat_neighbor = int(np.ravel_multi_index(neighbor_idx, shape))
            matrix[flat_idx, flat_neighbor] += -scale

        # Right neighbor: +1/(2dx)
        if multi_idx[axis] == shape[axis] - 1:
            bc_idx = list(multi_idx)
            bc_idx[axis] = shape[axis]
            const, entries = bc_axis.get_sparse_matrix_data(tuple(bc_idx))
            vector[flat_idx, 0] += const * scale
            for k, v in entries.items():
                neighbor_idx = list(multi_idx)
                neighbor_idx[axis] = k
                flat_neighbor = int(np.ravel_multi_index(neighbor_idx, shape))
                matrix[flat_idx, flat_neighbor] += v * scale
        else:
            neighbor_idx = list(multi_idx)
            neighbor_idx[axis] += 1
            flat_neighbor = int(np.ravel_multi_index(neighbor_idx, shape))
            matrix[flat_idx, flat_neighbor] += scale

    return matrix, vector


def _build_cross_derivative_matrix(
    grid: Any,  # noqa: ANN401
    bcs: Any,  # noqa: ANN401
    *,
    axis1: int,
    axis2: int,
) -> _MatrixAndVector:
    """Build sparse matrix for cross derivative ∂²/(∂x_i ∂x_j).

    Computed as the product of two gradient matrices: G_i @ G_j.
    This naturally handles boundary conditions since each gradient
    matrix already incorporates BCs.
    """
    g1_mat, g1_vec = _build_gradient_matrix(grid, bcs, axis=axis1)
    g2_mat, g2_vec = _build_gradient_matrix(grid, bcs, axis=axis2)

    # Cross derivative = G1 @ G2 (apply axis2 gradient first, then axis1)
    # For the vector: G1 @ (G2 @ field + vec2) + vec1
    #   = (G1 @ G2) @ field + G1 @ vec2 + vec1
    cross_matrix = g1_mat @ g2_mat  # type: ignore[reportOperatorIssue]
    cross_vector = g1_mat @ g2_vec + g1_vec  # type: ignore[reportOperatorIssue]

    return cross_matrix, cross_vector  # type: ignore[reportUnknownVariableType]


def _build_biharmonic_matrix(
    grid: Any,  # noqa: ANN401
    bcs: Any,  # noqa: ANN401
) -> _MatrixAndVector:
    """Build sparse matrix for biharmonic operator ∇⁴ = ∇²(∇²).

    Computed as the square of the Laplacian matrix: L @ L.
    """
    lap_mat, lap_vec = _build_laplacian_matrix(grid, bcs)

    # L @ (L @ field + vec) + vec = L² @ field + L @ vec + vec
    biharm_matrix = lap_mat @ lap_mat  # type: ignore[reportOperatorIssue]
    biharm_vector = lap_mat @ lap_vec + lap_vec  # type: ignore[reportOperatorIssue]

    return biharm_matrix, biharm_vector  # type: ignore[reportUnknownVariableType]


#: Registry mapping operator names to sparse matrix builder functions.
#: Each builder takes (grid, bcs) and returns (sparse_matrix, sparse_vector).
#: To add constraint-solving support for a new operator, register it here.
_OPERATOR_MATRIX_REGISTRY: dict[str, Any] = {
    "identity": _build_identity_matrix,
    "laplacian": _build_laplacian_matrix,
    "laplacian_x": partial(_build_directional_laplacian_matrix, axis=0),
    "laplacian_y": partial(_build_directional_laplacian_matrix, axis=1),
    "laplacian_z": partial(_build_directional_laplacian_matrix, axis=2),
    "gradient_x": partial(_build_gradient_matrix, axis=0),
    "gradient_y": partial(_build_gradient_matrix, axis=1),
    "gradient_z": partial(_build_gradient_matrix, axis=2),
    "cross_derivative_xy": partial(_build_cross_derivative_matrix, axis1=0, axis2=1),
    "cross_derivative_xz": partial(_build_cross_derivative_matrix, axis1=0, axis2=2),
    "cross_derivative_yz": partial(_build_cross_derivative_matrix, axis1=1, axis2=2),
    "biharmonic": _build_biharmonic_matrix,
}


def _coalesce_directional_laplacians(
    terms: list[OperatorTerm],
    spatial_dimension: int,
) -> list[OperatorTerm]:
    """Replace directional Laplacians with a compact Laplacian when safe.

    When all spatial axes have ``laplacian_{axis}`` terms targeting the same
    field with the same constant coefficient, they are mathematically equivalent
    to a single ``laplacian`` operator.  The compact stencil couples all
    adjacent grid points, avoiding the checkerboard decoupling that arises
    from the wide stencil produced by composing gradient matrices (G @ G).

    Guards (all must pass):
    - All spatial axes present (e.g. ``laplacian_x`` AND ``laplacian_y`` in 2D)
    - All target the same field
    - All have the same numeric coefficient
    - All have the same ``coefficient_symbolic`` (different symbolic names may
      evaluate differently at runtime even if current numerics match)
    - None are time-dependent or position-dependent (the ``position_dependent``
      check also catches coordinate-dependent symbolic coefficients like
      ``"f(x)"`` that happen to have the same string across axes)

    Returns the original list unchanged if any guard fails.
    """
    expected_axes = set(AXIS_LETTERS[:spatial_dimension])
    expected_ops = {f"laplacian_{a}" for a in expected_axes}

    # Collect directional-laplacian terms
    dir_lap_terms: list[OperatorTerm] = []
    other_terms: list[OperatorTerm] = []
    for term in terms:
        if term.operator in expected_ops:
            dir_lap_terms.append(term)
        else:
            other_terms.append(term)

    # Guard: need exactly the right number of directional laplacians
    if len(dir_lap_terms) != spatial_dimension:
        return terms

    # Guard: all axes present
    found_ops = {t.operator for t in dir_lap_terms}
    if found_ops != expected_ops:
        return terms

    # Guard: all target the same field
    fields = {t.field for t in dir_lap_terms}
    if len(fields) != 1:
        return terms

    # Guard: all have the same coefficient AND coefficient_symbolic (different
    # symbolic names may evaluate to different values at runtime even if current
    # numerics match)
    coeffs = {t.coefficient for t in dir_lap_terms}
    symbolic_coeffs = {t.coefficient_symbolic for t in dir_lap_terms}
    if len(coeffs) != 1 or len(symbolic_coeffs) != 1:
        return terms

    # Guard: none are time-dependent or position-dependent
    if any(t.time_dependent or t.position_dependent for t in dir_lap_terms):
        return terms

    # All guards passed — coalesce into compact laplacian
    representative = dir_lap_terms[0]
    coalesced = OperatorTerm(
        coefficient=representative.coefficient,
        operator="laplacian",
        field=representative.field,
        coefficient_symbolic=representative.coefficient_symbolic,
        time_dependent=False,
        coordinate_dependent=(),
    )
    logger.debug(
        "Coalesced %d directional laplacians into compact laplacian "
        "(field=%s, coeff=%s)",
        len(dir_lap_terms),
        representative.field,
        representative.coefficient,
    )
    return [coalesced, *other_terms]


# ---------------------------------------------------------------------------
# FFT operator multiplier registry: discrete Fourier-space representations
#
# For fully periodic grids, solving constraint equations in Fourier space
# is O(N log N) instead of O(N²) for sparse matrix solve.
#
# Each multiplier function takes (k_grids, dx_array) and returns the
# discrete Fourier-space symbol of the operator. These are eigenvalues
# of the finite-difference stencils used by py-pde, NOT the continuous
# symbols (e.g., the Laplacian uses (2cos(k·dx)-2)/dx², not -k²).
# This ensures the FFT path produces identical results to the matrix path.
#
# k_grids: list of wavenumber arrays, one per axis (angular frequency)
# dx_array: grid discretization per axis
#
# The combined operator symbol is: sum_i(c_i * symbol_i(k, dx))
# and the solution is: field_hat = source_hat / combined_symbol.
# ---------------------------------------------------------------------------


def _fft_identity(
    k_grids: list[np.ndarray],
    dx_array: np.ndarray,  # noqa: ARG001
) -> np.ndarray:
    return np.ones_like(k_grids[0])


def _fft_laplacian(k_grids: list[np.ndarray], dx_array: np.ndarray) -> np.ndarray:
    """Discrete Laplacian: sum_i (2cos(k_i dx_i) - 2) / dx_i²."""
    result = np.zeros_like(k_grids[0], dtype=complex)
    for k, dx in zip(k_grids, dx_array, strict=True):
        result += (2.0 * np.cos(k * dx) - 2.0) / dx**2
    return result


def _fft_directional_laplacian(
    k_grids: list[np.ndarray], dx_array: np.ndarray, *, axis: int
) -> np.ndarray:
    """Discrete directional laplacian = (gradient)² along one axis.

    Uses the wide stencil matching _build_directional_laplacian_matrix:
    G² eigenvalue = (i·sin(k·dx)/dx)² = -sin²(k·dx)/dx².
    """
    k = k_grids[axis]
    dx = dx_array[axis]
    return -(np.sin(k * dx) ** 2) / dx**2


def _fft_gradient(
    k_grids: list[np.ndarray], dx_array: np.ndarray, *, axis: int
) -> np.ndarray:
    """Discrete gradient: i·sin(k·dx) / dx (central difference)."""
    k = k_grids[axis]
    dx = dx_array[axis]
    return 1j * np.sin(k * dx) / dx


def _fft_cross_derivative(
    k_grids: list[np.ndarray], dx_array: np.ndarray, *, axis1: int, axis2: int
) -> np.ndarray:
    """Discrete cross derivative: product of two gradient eigenvalues."""
    g1 = _fft_gradient(k_grids, dx_array, axis=axis1)
    g2 = _fft_gradient(k_grids, dx_array, axis=axis2)
    return g1 * g2


def _fft_biharmonic(k_grids: list[np.ndarray], dx_array: np.ndarray) -> np.ndarray:
    """Discrete biharmonic = (discrete laplacian)²."""
    lap = _fft_laplacian(k_grids, dx_array)
    return lap**2


_OPERATOR_FFT_MULTIPLIERS: dict[str, Any] = {
    "identity": _fft_identity,
    "laplacian": _fft_laplacian,
    "laplacian_x": partial(_fft_directional_laplacian, axis=0),
    "laplacian_y": partial(_fft_directional_laplacian, axis=1),
    "laplacian_z": partial(_fft_directional_laplacian, axis=2),
    "gradient_x": partial(_fft_gradient, axis=0),
    "gradient_y": partial(_fft_gradient, axis=1),
    "gradient_z": partial(_fft_gradient, axis=2),
    "cross_derivative_xy": partial(_fft_cross_derivative, axis1=0, axis2=1),
    "cross_derivative_xz": partial(_fft_cross_derivative, axis1=0, axis2=2),
    "cross_derivative_yz": partial(_fft_cross_derivative, axis1=1, axis2=2),
    "biharmonic": _fft_biharmonic,
}


@dataclass(frozen=True)
class ParsedFieldName:
    """Parsed field name components.

    Supports multiple field naming conventions:
    - standard: A_0, phi_1 (base_index)
    - tensor: stress_xy_0, u_x_1 (base_component_index)
    - compact: phi0, A1 (base+digits)
    - simple: phi, psi (letters only, index defaults to 0)
    """

    base: str
    index: int
    format: str

    @classmethod
    def parse(cls, name: str) -> ParsedFieldName:
        """Parse field name auto-detecting format.

        Parameters
        ----------
        name : str
            Field name to parse.

        Returns
        -------
        ParsedFieldName
            Parsed components with base, index, and format.
        """
        # Standard format: A_0, phi_1
        match = _STANDARD_FORMAT_RE.match(name)
        if match:
            return cls(
                base=match.group(1), index=int(match.group(2)), format="standard"
            )

        # Tensor format: stress_xy_0, u_x_1 (greedy match for base)
        match = _TENSOR_FORMAT_RE.match(name)
        if match:
            return cls(base=match.group(1), index=int(match.group(2)), format="tensor")

        # Compact format: phi0, A1
        match = _COMPACT_FORMAT_RE.match(name)
        if match:
            return cls(base=match.group(1), index=int(match.group(2)), format="compact")

        # Simple format: phi, psi (no index, defaults to 0)
        if _SIMPLE_FORMAT_RE.match(name):
            return cls(base=name, index=0, format="simple")

        # Fallback
        return cls(base=name, index=0, format="unknown")

    def to_momentum_name(self) -> str:
        """Convert to momentum field name."""
        return f"pi_{self.index}"


def parse_momentum_field_name(field_name: str) -> int | None:
    """Parse momentum field name and return index.

    Supports both pi_N and piN formats.

    Parameters
    ----------
    field_name : str
        Momentum field name like "pi_0", "pi0", "pi_1", "pi1".

    Returns
    -------
    int | None
        Index if valid momentum field name, None otherwise.
    """
    # Standard format: pi_N
    match = _MOMENTUM_STANDARD_RE.match(field_name)
    if match:
        return int(match.group(1))

    # Compact format: piN
    match = _MOMENTUM_COMPACT_RE.match(field_name)
    if match:
        return int(match.group(1))

    return None


class PDEFromSpec(PDEBase):
    """Generic PDE class built from JSON equation specification.

    This class dynamically constructs the evolution equations from a parsed
    specification. NO physics is hardcoded - all equation structure comes
    from the EquationSystem that was derived from a Lagrangian.

    The state layout is determined by each component's time derivative order:
    - Second-order (wave, time_order>=2): [field_i, momentum_i] pair
    - First-order (heat/diffusion, time_order=1): [field_i] only
    - Constraint (elliptic, time_order=0): [field_i] only (no evolution)

    For second-order components:
        d/dt field_i = momentum_i
        d/dt momentum_i = RHS from specification
    For first-order components:
        d/dt field_i = RHS from specification
    For constraint components:
        d/dt field_i = 0

    Parameters
    ----------
    spec : EquationSystem
        The equation specification loaded from JSON.

    Attributes
    ----------
    spec : EquationSystem
        The equation specification.
    n_components : int
        Number of field components.
    explicit_time_dependence : bool
        True to support time-dependent coefficients (e.g., de Sitter spacetime).

    Examples
    --------
    >>> from tidal.symbolic import load_equation_system
    >>> from tidal.symbolic.pde_builder import PDEFromSpec
    >>> spec = load_equation_system("examples/data/em_1d.json")
    >>> pde = PDEFromSpec(spec)
    >>> # pde can now be used with py-pde solvers
    """

    # Enable time-dependent coefficients for curved spacetime (e.g., de Sitter)
    # This allows evolution_rate to receive the current time t
    explicit_time_dependence = True

    def __init__(
        self,
        spec: EquationSystem,
        parameters: dict[str, float] | None = None,
        *,
        constraint_eps: float = 1e-14,
        coupled_svd_rcond: float = 0.01,
    ) -> None:
        """Initialize PDE from equation specification.

        Parameters
        ----------
        spec : EquationSystem
            The equation specification loaded from JSON.
        parameters : dict[str, float] | None
            Optional parameter values to override symbolic coefficients.
            Keys are symbolic names (e.g., "dSH", "dSm2", "kappa"), values are numeric.
            When a term has a coefficient_symbolic that matches a key in this
            dict, the parameter value is used instead of the numeric coefficient.

            For time-dependent coefficients in curved spacetime, all symbols
            appearing in the coefficient expression must be provided here.
            The expressions are evaluated by substituting these values.
        constraint_eps : float
            Tolerance for determining whether a Laplacian coefficient is
            effectively zero in the constraint solver.  Default ``1e-14``.
        coupled_svd_rcond : float
            Relative singular-value threshold for Tikhonov regularization
            in the coupled FFT constraint solver.  Singular values smaller
            than ``rcond * max(S)`` are attenuated instead of inverted
            directly, preventing noise amplification at near-singular
            wavenumbers (e.g., Helmholtz resonance at k² ≈ m²).
            Default ``0.01``.
        """
        super().__init__()
        self.spec = spec
        self.n_components = spec.n_components
        self._component_name_to_index = {
            name: i for i, name in enumerate(spec.component_names)
        }
        self._parameters = parameters or {}
        self._constraint_eps = constraint_eps
        self._coupled_svd_rcond = coupled_svd_rcond

        # Build slot maps from state_layout for mixed time-order support
        self._field_slot_map: dict[str, int] = {}
        self._momentum_slot_map: dict[str, int] = {}
        for slot_idx, (name, slot_type) in enumerate(spec.state_layout):
            if slot_type == "field":
                self._field_slot_map[name] = slot_idx
            else:
                self._momentum_slot_map[name] = slot_idx

        # B1: Cache for _mathematica_to_python() results (same symbolic expr → same output)
        self._expr_cache: dict[str, str] = {}

        # B2: Pre-build the static part of the coefficient namespace (math functions,
        # parameters). Only t and grid coordinates change per-call.
        self._base_namespace: dict[str, Any] = self._build_base_namespace()

        # B4: Pre-resolve constant coefficients (not position- or time-dependent).
        # Only pre-resolve when the symbolic coefficient is absent or resolvable
        # from the provided parameters — avoids emitting premature warnings for
        # terms whose symbolic cannot be resolved (the warning fires at runtime).
        self._preresolved: dict[tuple[int, int], float] = {}
        for eq_idx, eq in enumerate(spec.equations):
            for term_idx, term in enumerate(eq.rhs_terms):
                if (
                    not term.time_dependent
                    and not term.position_dependent
                    and self._is_resolvable(term)
                ):
                    self._preresolved[eq_idx, term_idx] = self._resolve_coefficient(
                        term
                    )

        # B5: Cache boundary conditions and grid coordinates (populated on first call)
        self._cached_bc: BCDescriptor | None = None
        self._cached_grid_id: int | None = None

        # Validate operator dimension requirements against spec at construction time
        self._validate_operator_dimensions()

    @staticmethod
    def _operator_min_dim(operator_name: str) -> int:
        """Return minimum spatial grid dimension required by an operator.

        Uses the same resolution logic as ``_get_operator`` (static registry
        then dynamic regex patterns) but without needing a grid or field.

        Parameters
        ----------
        operator_name : str
            Operator name (e.g. ``"gradient_y"``, ``"derivative_3_z"``).

        Returns
        -------
        int
            Minimum spatial dimension (1, 2, or 3).

        Raises
        ------
        ValueError
            If the operator name is not recognized.
        """
        entry = _OPERATOR_REGISTRY.get(operator_name)
        if entry is not None:
            return entry[1]
        m = _GENERIC_SINGLE_RE.match(operator_name)
        if m:
            return _AXIS_MIN_DIM[m.group(2)]
        m_multi = _GENERIC_MULTI_RE.match(operator_name)
        if m_multi:
            axes_and_orders = _parse_multi_axis_spec(m_multi.group(1))
            max_axis = max(axis for axis, _ in axes_and_orders)
            return _AXIS_MIN_DIM[_AXIS_LETTER[max_axis]]
        msg = (
            f"Unknown operator: '{operator_name}'. "
            f"Known operators: {sorted(_OPERATOR_REGISTRY.keys())}. "
            f"Dynamic patterns: derivative_N_x, derivative_Nx_My."
        )
        raise ValueError(msg)

    def _validate_operator_dimensions(self) -> None:
        """Validate all operators are compatible with the spec's spatial dimension.

        Called during ``__init__`` to fail fast when a JSON spec contains
        operators that require more spatial dimensions than the spec provides
        (e.g. ``gradient_z`` in a 2+1D spacetime).

        Raises
        ------
        ValueError
            If any operator requires more dimensions than ``spec.spatial_dimension``.
        """
        spatial_dim = self.spec.spatial_dimension
        for eq in self.spec.equations:
            for term in eq.rhs_terms:
                min_dim = self._operator_min_dim(term.operator)
                if min_dim > spatial_dim:
                    msg = (
                        f"Operator '{term.operator}' in equation for "
                        f"'{eq.field_name}' requires at least {min_dim}D "
                        f"spatial grid, but the spec has "
                        f"spatial_dimension={spatial_dim} "
                        f"(from {self.spec.dimension}D spacetime)."
                    )
                    raise ValueError(msg)

    def _is_resolvable(self, term: OperatorTerm) -> bool:
        """Check whether a term's coefficient can be resolved without warnings.

        Returns True if the term has no symbolic coefficient, or if its symbolic
        coefficient (possibly negated or compound like ``-2*m2``) can be
        evaluated from the parameters dict.
        """
        sym = term.coefficient_symbolic
        if sym is None:
            return True
        if sym.startswith("-") and sym[1:] in self._parameters:
            return True
        if sym in self._parameters:
            return True
        # Try compound expression evaluation (e.g., "-2*m2")
        try:
            py_expr = self._mathematica_to_python(sym)
            eval(py_expr, {"__builtins__": {}}, dict(self._base_namespace))  # noqa: S307
        except (NameError, SyntaxError, TypeError, ValueError, ZeroDivisionError):
            return False
        else:
            return True

    def _resolve_coefficient(self, term: OperatorTerm) -> float:
        """Resolve the effective coefficient for a term.

        If the term has a symbolic coefficient name and that name (or its
        negation) is in the parameters dict, use the parameter value.
        Otherwise use the numeric coefficient from the JSON.

        Parameters
        ----------
        term : OperatorTerm
            The term whose coefficient to resolve.

        Returns
        -------
        float
            The effective coefficient value.

        Raises
        ------
        ValueError
            If a symbolic coefficient cannot be resolved from parameters.
        """
        if term.coefficient_symbolic is not None:
            sym = term.coefficient_symbolic

            # Fast path: simple parameter name or negated name
            if sym.startswith("-") and sym[1:] in self._parameters:
                return -self._parameters[sym[1:]]
            if sym in self._parameters:
                return self._parameters[sym]

            # Compound expression: e.g., "-2*m2", "3*lambda"
            try:
                py_expr = self._mathematica_to_python(sym)
                result = eval(py_expr, {"__builtins__": {}}, dict(self._base_namespace))  # noqa: S307
            except (NameError, SyntaxError, TypeError, ValueError, ZeroDivisionError):
                pass
            else:
                return float(result)

            # Symbolic coefficient present but not resolvable from parameters.
            # Fail fast: wrong physics that looks right is worse than a crash.
            msg = (
                f"Symbolic coefficient '{sym[:80]}' could not be resolved from "
                f"parameters {sorted(self._parameters.keys())}. "
                f"Pass the required parameter via "
                f"parameters={{'{sym.lstrip('-')}': <value>}} to "
                f"PDEFromSpec or build_pde_from_json."
            )
            raise ValueError(msg)

        # Default: use numeric coefficient from JSON (no symbolic → purely numeric term)
        return term.coefficient

    def _build_base_namespace(self) -> dict[str, Any]:
        """Build the static part of the coefficient evaluation namespace.

        Contains numpy/scipy math functions and user-provided parameters.
        The time variable ``t`` and spatial grid coordinates are injected
        per-call in ``_resolve_coefficient_at_point``.
        """
        ns: dict[str, Any] = dict(self._parameters)
        ns["exp"] = np.exp
        # Basic trig
        ns["sin"] = np.sin
        ns["cos"] = np.cos
        ns["tan"] = np.tan
        # Reciprocal trig (no direct numpy equivalents)
        ns["cot"] = lambda x: np.cos(x) / np.sin(x)  # type: ignore[reportUnknownLambdaType]
        ns["sec"] = lambda x: 1.0 / np.cos(x)  # type: ignore[reportUnknownLambdaType]
        ns["csc"] = lambda x: 1.0 / np.sin(x)  # type: ignore[reportUnknownLambdaType]
        # Inverse trig
        ns["arcsin"] = np.arcsin
        ns["arccos"] = np.arccos
        ns["arctan"] = np.arctan
        ns["arctan2"] = np.arctan2
        # Hyperbolic
        ns["sinh"] = np.sinh
        ns["cosh"] = np.cosh
        ns["tanh"] = np.tanh
        # Inverse hyperbolic
        ns["arcsinh"] = np.arcsinh
        ns["arccosh"] = np.arccosh
        ns["arctanh"] = np.arctanh
        # Other
        ns["log"] = np.log
        ns["sqrt"] = np.sqrt
        ns["abs"] = np.abs
        ns["sign"] = np.sign
        ns["maximum"] = np.maximum
        ns["minimum"] = np.minimum
        # Step functions (UnitStep/HeavisideTheta)
        ns["heaviside"] = lambda x: np.heaviside(x, 0.5)  # type: ignore[reportUnknownLambdaType]
        # Piecewise (from Mathematica Simplify converting UnitStep products)
        ns["piecewise"] = np.where
        # Special functions (scipy.special)
        ns["erf"] = special.erf
        ns["jv"] = special.jv  # BesselJ
        ns["yv"] = special.yv  # BesselY
        # numpy module — needed for np.pi in position-dependent expressions
        ns["np"] = np
        # Mathematica booleans — eval uses __builtins__={} which strips Python builtins
        ns["True"] = True
        ns["False"] = False
        return ns

    @staticmethod
    def _convert_power_function(expr: str) -> str:
        """Convert Power[base, exponent] to (base)**(exponent).

        Handles nested expressions in both arguments. Must run before bracket
        conversion ([] → ()).

        Parameters
        ----------
        expr : str
            Expression potentially containing Power[...] syntax.

        Returns
        -------
        str
            Expression with Power[...] converted to (...)**(...).
        """
        # Pattern: Power[<arg1>, <arg2>] where args may contain nested brackets
        # (?:[^[\]]|\[[^\]]*\]) matches either non-bracket chars or [...] pairs
        pattern = r"Power\[((?:[^[\]]|\[[^\]]*\])*),\s*((?:[^[\]]|\[[^\]]*\])*)\]"

        def replacer(match: re.Match[str]) -> str:
            base = match.group(1).strip()
            exp = match.group(2).strip()
            return f"({base})**({exp})"

        # Multiple passes for nested Power calls
        prev = None
        result = expr
        while prev != result:
            prev = result
            result = re.sub(pattern, replacer, result)
        return result

    @staticmethod
    def _convert_arctan2(expr: str) -> str:
        """Convert ArcTan[x, y] to arctan2(y, x) with argument swap.

        Mathematica's ArcTan[x, y] computes atan2(y, x), so arguments must
        be swapped during conversion to match NumPy's arctan2(y, x) signature.

        Parameters
        ----------
        expr : str
            Expression potentially containing ArcTan[x, y] syntax.

        Returns
        -------
        str
            Expression with ArcTan[x, y] converted to arctan2(y, x).
        """
        # Pattern: ArcTan[<arg1>, <arg2>] - must handle before generic function conversion
        pattern = r"ArcTan\[((?:[^[\],]|\[[^\]]*\])*),\s*((?:[^[\]]|\[[^\]]*\])*)\]"

        def replacer(match: re.Match[str]) -> str:
            x = match.group(1).strip()
            y = match.group(2).strip()
            return f"arctan2({y}, {x})"  # Swap x and y!

        return re.sub(pattern, replacer, expr)

    # Comparison operator names used inside Inequality[...]
    _COMPARISON_OPS: ClassVar[dict[str, str]] = {
        "LessEqual": "<=",
        "Less": "<",
        "GreaterEqual": ">=",
        "Greater": ">",
        "Equal": "==",
    }

    @staticmethod
    def _split_bracket_aware(s: str) -> list[str]:
        """Split string on commas, respecting ``[...]`` bracket nesting.

        Parameters
        ----------
        s : str
            Content string (typically inside a Mathematica function call).

        Returns
        -------
        list[str]
            Comma-separated parts, preserving nested bracket structure.
        """
        parts: list[str] = []
        depth = 0
        current: list[str] = []
        for ch in s:
            if ch == "[":
                depth += 1
                current.append(ch)
            elif ch == "]":
                depth -= 1
                current.append(ch)
            elif ch == "," and depth == 0:
                parts.append("".join(current))
                current = []
            else:
                current.append(ch)
        if current:
            parts.append("".join(current))
        return parts

    @staticmethod
    def _convert_inequality(expr: str) -> str:
        """Convert ``Inequality[a, op, b, op, c, ...]`` to chained comparisons.

        Mathematica's ``Inequality`` represents compound comparisons like
        ``30 <= x <= 70``.  The function form is::

            Inequality[30, LessEqual, x[], LessEqual, 70]

        This is converted to::

            ((30) <= (x[])) & ((x[]) <= (70))

        Must run before bracket conversion (Step 5) and before Piecewise
        conversion (Inequality appears inside Piecewise conditions).

        Parameters
        ----------
        expr : str
            Expression potentially containing ``Inequality[...]`` syntax.

        Returns
        -------
        str
            Expression with ``Inequality`` converted to Python comparisons.
        """
        pattern = r"Inequality\[((?:[^[\]]|\[[^\]]*\])*)\]"

        def replacer(match: re.Match[str]) -> str:
            inner = match.group(1)
            args = PDEFromSpec._split_bracket_aware(inner)
            min_args = 3  # Inequality[a, op, b] minimum
            if len(args) < min_args or len(args) % 2 == 0:
                return match.group(0)  # Malformed, leave unchanged
            parts: list[str] = []
            for i in range(0, len(args) - 2, 2):
                left = args[i].strip()
                op_name = args[i + 1].strip()
                right = args[i + 2].strip()
                op_sym = PDEFromSpec._COMPARISON_OPS.get(op_name, op_name)
                parts.append(f"(({left}) {op_sym} ({right}))")
            return " & ".join(parts)

        prev = None
        result = expr
        while prev != result:
            prev = result
            result = re.sub(pattern, replacer, result)
        return result

    @staticmethod
    def _find_top_level_comma(s: str) -> int:
        """Find index of first comma not inside ``()[]{}``."""
        depth = 0
        for i, ch in enumerate(s):
            if ch in "([{":
                depth += 1
            elif ch in ")]}":
                depth -= 1
            elif ch == "," and depth == 0:
                return i
        return -1

    @staticmethod
    def _find_matching_brace(s: str) -> int:
        """Find index of the closing ``}`` that matches the opening ``{`` at s[0]."""
        depth = 0
        for i, ch in enumerate(s):
            if ch == "{":
                depth += 1
            elif ch == "}":
                depth -= 1
                if depth == 0:
                    return i
        return -1

    @staticmethod
    def _extract_brace_pairs(cases_str: str) -> list[tuple[str, str]]:
        """Extract ``{value, condition}`` pairs from the inner list string."""
        cases: list[tuple[str, str]] = []
        i = 0
        while i < len(cases_str):
            if cases_str[i] == "{":
                end = PDEFromSpec._find_matching_brace(cases_str[i:])
                if end < 0:
                    break
                case_content = cases_str[i + 1 : i + end]
                comma_idx = PDEFromSpec._find_top_level_comma(case_content)
                if comma_idx >= 0:
                    val = case_content[:comma_idx].strip()
                    cond = case_content[comma_idx + 1 :].strip()
                    cases.append((val, cond))
                i += end + 1
            else:
                i += 1
        return cases

    @staticmethod
    def _parse_piecewise_content(
        content: str,
    ) -> tuple[list[tuple[str, str]], str]:
        """Parse the inner content of a ``Piecewise[...]`` expression.

        Parameters
        ----------
        content : str
            Everything between ``Piecewise[`` and the closing ``]``.
            Expected form: ``{{val1, cond1}, {val2, cond2}, ...}, default``

        Returns
        -------
        tuple[list[tuple[str, str]], str]
            List of ``(value, condition)`` pairs and the default value string.
        """
        content = content.strip()
        if not content.startswith("{{"):
            return [], content

        end_idx = PDEFromSpec._find_matching_brace(content)
        if end_idx < 0:
            return [], content  # Malformed

        cases_str = content[1:end_idx]  # Strip outer { ... }
        default_str = content[end_idx + 1 :].strip()
        if default_str.startswith(","):
            default_str = default_str[1:].strip()
        if not default_str:
            default_str = "0"

        return PDEFromSpec._extract_brace_pairs(cases_str), default_str

    @staticmethod
    def _convert_piecewise(expr: str) -> str:
        """Convert ``Piecewise[{{val, cond}, ...}, default]`` to ``piecewise()`` calls.

        Mathematica's ``Simplify`` often converts products of ``UnitStep``
        functions into ``Piecewise`` form.  This converts back into calls to
        ``piecewise(cond, val_true, val_false)`` (a thin ``numpy.where``
        wrapper provided in the eval namespace).

        Must run AFTER ``_convert_inequality`` (conditions are then
        bracket-free comparison chains) and BEFORE bracket conversion.

        Parameters
        ----------
        expr : str
            Expression potentially containing ``Piecewise[...]`` syntax.

        Returns
        -------
        str
            Expression with ``Piecewise`` converted to ``piecewise()`` calls.
        """
        # After _convert_inequality, conditions contain comparisons with &
        # and coordinate symbols like x[].  Square-bracket nesting inside
        # Piecewise[...] is at most one level deep (e.g. x[]).
        pattern = r"Piecewise\[((?:[^[\]]|\[[^\]]*\])*)\]"

        def replacer(match: re.Match[str]) -> str:
            content = match.group(1)
            cases, default = PDEFromSpec._parse_piecewise_content(content)
            if not cases:
                return default
            result = default
            for val, cond in reversed(cases):
                result = f"piecewise({cond}, {val}, {result})"
            return result

        prev = None
        result = expr
        while prev != result:
            prev = result
            result = re.sub(pattern, replacer, result)
        return result

    def _mathematica_to_python(self, expr: str) -> str:
        """Convert Mathematica InputForm expression to evaluable Python.

        Handles common Mathematica syntax:
        - ``E^(...)`` to ``exp(...)`` (Euler's number)
        - ``Power[x,y]`` to ``(x)**(y)`` (function form of exponentiation)
        - ``Sin[x]`` to ``sin(x)``, ``Cos[x]`` to ``cos(x)``, ``Tan[x]`` to ``tan(x)``
        - ``Cot[x]`` to ``cot(x)``, ``Sec[x]`` to ``sec(x)``, ``Csc[x]`` to ``csc(x)``
        - ``ArcSin[x]`` to ``arcsin(x)``, ``ArcCos[x]`` to ``arccos(x)``, etc.
        - ``ArcTan[x, y]`` to ``arctan2(y, x)`` (note argument order swap!)
        - ``Sinh[x]`` to ``sinh(x)``, ``Cosh[x]`` to ``cosh(x)``, etc.
        - ``ArcSinh[x]`` to ``arcsinh(x)``, ``ArcCosh[x]`` to ``arccosh(x)``, etc.
        - ``Erf[x]`` to ``erf(x)`` (scipy.special)
        - ``BesselJ[n, x]`` to ``jv(n, x)``, ``BesselY[n, x]`` to ``yv(n, x)``
        - ``Rational[p, q]`` to ``(p)/(q)`` (exact fractions from xAct InputForm)
        - ``Pi`` to ``np.pi`` (mathematical constant)
        - ``Sign[x]`` to ``sign(x)``, ``Max[a,b]`` to ``maximum(a,b)``
        - ``UnitStep[x]`` to ``heaviside(x)``, ``HeavisideTheta[x]`` to ``heaviside(x)``
        - ``Inequality[a, LessEqual, b, LessEqual, c]`` to ``((a) <= (b)) & ((b) <= (c))``
        - ``Piecewise[{{val, cond}}, default]`` to ``piecewise(cond, val, default)``
        - ``t[]`` to ``t`` (xCoba coordinate symbols, using actual coordinate names)
        - Mathematica brackets ``[``, ``]`` to Python parens ``(``, ``)``
        - Mathematica ``^`` to Python ``**``

        Parameters
        ----------
        expr : str
            Mathematica InputForm expression string.

        Returns
        -------
        str
            Python-evaluable expression string.
        """
        result = expr

        # Step 1: E^(...) → exp(...) — Mathematica's Euler number
        result = re.sub(r"\bE\^", "exp", result)

        # Step 2: Power[x, y] → (x)**(y) — must handle before bracket conversion
        result = self._convert_power_function(result)

        # Step 3: ArcTan[x, y] → arctan2(y, x) — special 2-arg case with swap
        # Must handle before generic function conversion
        result = self._convert_arctan2(result)

        # Step 3.5: Rational[p, q] → (p)/(q) — Mathematica exact fractions
        # xAct outputs Rational[1,2] for 1/2 in InputForm. Must handle before
        # bracket conversion since Rational uses Mathematica brackets.
        result = re.sub(r"Rational\[([^,\]]+),\s*([^,\]]+)\]", r"(\1)/(\2)", result)

        # Step 3.7: Inequality[a, op, b, op, c] → chained comparisons
        # Must run before Piecewise (Inequality appears inside Piecewise conditions)
        result = self._convert_inequality(result)

        # Step 3.8: Piecewise[{{val, cond}, ...}, default] → piecewise() calls
        # Must run after Inequality (conditions are now comparison chains)
        result = self._convert_piecewise(result)

        # Step 4: Function name conversions (batch)
        function_map = [
            # Basic trig
            ("Sin", "sin"),
            ("Cos", "cos"),
            ("Tan", "tan"),
            # Reciprocal trig
            ("Cot", "cot"),
            ("Sec", "sec"),
            ("Csc", "csc"),
            # Inverse trig (1-arg)
            ("ArcSin", "arcsin"),
            ("ArcCos", "arccos"),
            ("ArcTan", "arctan"),  # 1-arg version only (2-arg handled above)
            # Hyperbolic
            ("Sinh", "sinh"),
            ("Cosh", "cosh"),
            ("Tanh", "tanh"),
            # Inverse hyperbolic
            ("ArcSinh", "arcsinh"),
            ("ArcCosh", "arccosh"),
            ("ArcTanh", "arctanh"),
            # Other
            ("Exp", "exp"),
            ("Log", "log"),
            ("Sqrt", "sqrt"),
            ("Abs", "abs"),
            ("Sign", "sign"),
            ("Max", "maximum"),
            ("Min", "minimum"),
            # Step functions
            ("UnitStep", "heaviside"),
            ("HeavisideTheta", "heaviside"),
            # Special functions (scipy.special)
            ("Erf", "erf"),
            ("BesselJ", "jv"),
            ("BesselY", "yv"),
        ]
        for mma_func, py_func in function_map:
            result = re.sub(rf"\b{mma_func}\b", py_func, result)

        # Step 4.5: Pi → np.pi — Mathematica constant (before bracket conversion)
        result = re.sub(r"\bPi\b", "np.pi", result)

        # Step 5: Mathematica brackets to Python parens (after function renaming)
        result = result.replace("[", "(").replace("]", ")")

        # Step 6: Mathematica ^ to Python ** (AFTER E^ → exp to avoid double-conversion)
        result = result.replace("^", "**")

        # Step 7: xCoba coordinate symbols: t() → t, x() → x, etc.
        # Uses actual coordinate names from equation system (not hardcoded x/y/z).
        for coord in self.spec.effective_coordinates:
            result = result.replace(f"{coord}()", coord)

        return result

    def _resolve_coefficient_at_point(
        self,
        term: OperatorTerm,
        t: float,
        grid: GridBase | None = None,
        coord_arrays: dict[str, NumericArray] | None = None,
    ) -> float | NumericArray:
        """Resolve a potentially coordinate-dependent coefficient.

        Handles constant, time-dependent, and position-dependent coefficients
        in a single unified code path. Uses numpy functions throughout so that
        the same evaluation works for both scalar (time-only) and array
        (position-dependent) results.

        Returns a scalar ``float`` for constant or time-only coefficients, or a
        ``numpy.ndarray`` (same shape as the grid) for position-dependent ones.

        Parameters
        ----------
        term : OperatorTerm
            The term whose coefficient to resolve.
        t : float
            Current simulation time.
        grid : GridBase | None
            Simulation grid. Required when the term is position-dependent.

        Returns
        -------
        float | NumericArray
            Scalar or grid-shaped array of coefficient values.

        Raises
        ------
        ValueError
            If required parameters/grid are missing or expression cannot be evaluated.
        """
        # Fast path: no coordinate dependence → use simple parameter lookup
        if not term.time_dependent and not term.position_dependent:
            return self._resolve_coefficient(term)

        # Position-dependent: need grid for spatial coordinates
        if term.position_dependent and grid is None:
            msg = (
                f"Position-dependent coefficient '{term.coefficient_symbolic}' "
                f"requires grid info but no grid was provided."
            )
            raise ValueError(msg)

        sym = term.coefficient_symbolic or ""

        # B1: Use cached Mathematica→Python conversion
        if sym not in self._expr_cache:
            self._expr_cache[sym] = self._mathematica_to_python(sym)
        py_expr = self._expr_cache[sym]

        # B2: Clone pre-built base namespace and inject dynamic variables
        namespace: dict[str, Any] = dict(self._base_namespace)
        namespace["t"] = t

        # C2: Inject spatial coordinates — use pre-extracted arrays if available
        if term.position_dependent:
            if coord_arrays is not None:
                namespace.update(coord_arrays)
            else:
                spatial_coords = self.spec.spatial_coordinates
                coords = grid.cell_coords  # type: ignore[union-attr]
                for i, name in enumerate(spatial_coords[: grid.num_axes]):  # type: ignore[union-attr]
                    namespace[name] = np.asarray(coords[..., i], dtype=np.float64)  # pyright: ignore[reportUnknownArgumentType]

        # Validate all symbols can be resolved
        identifiers = set(re.findall(r"\b[a-zA-Z_]\w*\b", py_expr))
        # Derive builtin names from namespace (all math functions, excluding parameters and t)
        builtin_names = set(namespace.keys()) - set(self._parameters.keys()) - {"t"}
        # Exclude coordinate variables if position-dependent
        if term.position_dependent:
            builtin_names -= set(self.spec.spatial_coordinates)
        coord_vars = set(self.spec.effective_coordinates)
        identifiers -= builtin_names | coord_vars
        missing = identifiers - set(self._parameters.keys())
        if missing:
            msg = (
                f"Parameters {sorted(missing)} are required for "
                f"coordinate-dependent coefficient '{sym}'. "
                f"Pass them via parameters={{...}} to PDEFromSpec or build_pde_from_json."
            )
            raise ValueError(msg)

        try:
            result = eval(py_expr, {"__builtins__": {}}, namespace)  # noqa: S307
        except Exception as e:
            msg = (
                f"Cannot evaluate coordinate-dependent coefficient '{sym}' "
                f"(Python form: '{py_expr}') at t={t}: {e}"
            )
            raise ValueError(msg) from e
        return self._validate_eval_result(result, sym, py_expr)

    @staticmethod
    def _validate_eval_result(
        result: object, sym: str, py_expr: str
    ) -> float | NumericArray:
        """Validate and coerce an eval() result to float or ndarray.

        Raises ValueError for complex, NaN, or Inf results with clear
        diagnostic messages pointing to the source expression.

        Raises
        ------
        TypeError
            If the result is complex.
        ValueError
            If the result is NaN or Inf.
        """
        if isinstance(result, complex):
            msg = (
                f"Coefficient '{sym}' evaluated to complex number {result} "
                f"(from '{py_expr}'). Only real-valued coefficients are supported."
            )
            raise TypeError(msg)

        if isinstance(result, np.ndarray):
            arr = np.asarray(result, dtype=np.float64)
            if np.isnan(arr).any():
                msg = (
                    f"Coefficient '{sym}' produced NaN values "
                    f"(from '{py_expr}'). Check for 0/0 or invalid operations."
                )
                raise ValueError(msg)
            if np.isinf(arr).any():
                msg = (
                    f"Coefficient '{sym}' produced Inf values "
                    f"(from '{py_expr}'). Check for division by zero."
                )
                raise ValueError(msg)
            return arr

        scalar = float(cast("SupportsFloat", result))
        if math.isnan(scalar):
            msg = (
                f"Coefficient '{sym}' evaluated to NaN "
                f"(from '{py_expr}'). Check for 0/0 or invalid operations."
            )
            raise ValueError(msg)
        if math.isinf(scalar):
            msg = (
                f"Coefficient '{sym}' evaluated to Inf "
                f"(from '{py_expr}'). Check for division by zero."
            )
            raise ValueError(msg)
        return scalar

    @staticmethod
    def _get_operator(
        operator_name: str, field: ScalarField, bc: BCDescriptor
    ) -> ScalarField:
        """Apply a named operator to a field.

        Uses the module-level ``_OPERATOR_REGISTRY`` for dispatch.
        Each operator specifies a handler function and minimum grid dimension.

        Parameters
        ----------
        operator_name : str
            Name of the operator ("laplacian", "identity", "gradient_x", etc.)
        field : ScalarField
            The field to operate on.
        bc : BCDescriptor
            Boundary condition specification.

        Returns
        -------
        ScalarField
            Result of applying the operator.

        Raises
        ------
        ValueError
            If the operator is not recognized or the grid dimension is too low.
        RuntimeError
            If the operator is in the registry but cannot be applied as a spatial operator.
        """
        entry = _OPERATOR_REGISTRY.get(operator_name)
        if entry is None:
            # Try dynamic resolution for generic Nth-order derivatives
            m = _GENERIC_SINGLE_RE.match(operator_name)
            if m:
                order = int(m.group(1))
                axis_letter = m.group(2)
                axis = _AXIS_INDEX[axis_letter]
                min_dim = _AXIS_MIN_DIM[axis_letter]
                entry = (_op_nth_derivative(axis, order), min_dim)
            else:
                # Try multi-axis pattern: derivative_2x_1y, derivative_1x_1y_1z
                m_multi = _GENERIC_MULTI_RE.match(operator_name)
                if m_multi:
                    axes_and_orders = _parse_multi_axis_spec(m_multi.group(1))
                    max_axis = max(axis for axis, _ in axes_and_orders)
                    min_dim = _AXIS_MIN_DIM[_AXIS_LETTER[max_axis]]
                    entry = (_op_multi_axis_derivative(axes_and_orders), min_dim)
                else:
                    msg = (
                        f"Unknown operator: '{operator_name}'. "
                        f"Known operators: {sorted(_OPERATOR_REGISTRY.keys())}. "
                        f"Dynamic patterns: derivative_N_x, derivative_Nx_My."
                    )
                    raise ValueError(msg)

        handler, min_dim = entry

        # Sentinel check: some operators (first_derivative_t) are in the registry
        # but handled specially elsewhere — they cannot be applied as spatial operators.
        if handler is None:
            msg = (
                f"Operator '{operator_name}' cannot be applied as a spatial operator. "
                f"It is handled specially in _compute_rhs_for_component."
            )
            raise RuntimeError(msg)

        if field.grid.dim < min_dim:
            msg = (
                f"Operator '{operator_name}' requires at least {min_dim}D grid, "
                f"but got {field.grid.dim}D grid."
            )
            raise ValueError(msg)

        return handler(field, bc)

    def _get_field_from_state(
        self,
        state: FieldCollection,
        field_name: str,
        virtual_momenta: dict[str, ScalarField] | None = None,
    ) -> ScalarField:
        """Get a field from state by name, supporting both field and momentum names.

        For mixed time-space derivatives like d_t d_x A, the Wolfram pipeline
        expresses these as gradients of momentum fields (e.g., gradient_x(pi_0)).
        This is valid because d_t A = pi, so d_x(d_t A) = d_x(pi).

        State layout is determined by ``spec.state_layout``:
        - Second-order components have [field, momentum] pairs
        - First-order components have [field] only (no momentum in state)

        For first-order components, ``pi_N`` references are resolved via the
        ``virtual_momenta`` dict, which contains the computed RHS of those
        first-order equations (since d_t(field) = RHS for first-order PDEs).

        Parameters
        ----------
        state : FieldCollection
            Current state with fields ordered by state_layout.
        field_name : str
            Name of the field to retrieve. Can be:
            - Regular field name like "A_0", "phi_0", "phi0"
            - Momentum field name like "pi_0", "pi_1", "pi0", "pi1"
        virtual_momenta : dict[str, ScalarField] | None
            Pre-computed RHS fields for first-order components, keyed by
            component name. Used when ``pi_N`` references a first-order
            component that has no momentum state variable.

        Returns
        -------
        ScalarField
            The requested field from the state.

        Raises
        ------
        TypeError
            If a state element is not a ``ScalarField``.
        ValueError
            If the field name is not recognized or momentum is unavailable.
        """
        # Check if this is a momentum field reference (pi_0, pi0, etc.)
        if field_name.startswith("pi"):
            momentum_idx = parse_momentum_field_name(field_name)
            if momentum_idx is not None:
                if not (0 <= momentum_idx < self.n_components):
                    msg = (
                        f"Momentum field index {momentum_idx} out of range. "
                        f"This system has {self.n_components} components "
                        f"(valid indices: 0 to {self.n_components - 1}). "
                        f"Field reference: '{field_name}'."
                    )
                    raise ValueError(msg)

                # Look up the component name for this index, then find its momentum slot
                comp_name = self.spec.component_names[momentum_idx]
                slot = self._momentum_slot_map.get(comp_name)
                if slot is not None:
                    # Second-order component: momentum is a state variable
                    momentum = state[slot]
                    if not isinstance(momentum, ScalarField):
                        msg = f"Expected ScalarField for momentum, got {type(momentum).__name__}"
                        raise TypeError(msg)
                    return momentum

                # First-order component: check virtual_momenta
                if virtual_momenta is not None and comp_name in virtual_momenta:
                    return virtual_momenta[comp_name]

                eq = self.spec.equations[momentum_idx]
                msg = (
                    f"Momentum field '{field_name}' referenced but component "
                    f"'{comp_name}' has time_derivative_order={eq.time_derivative_order} "
                    f"(no momentum in state, and no virtual momentum computed). "
                    f"This may indicate a circular dependency between first-order components."
                )
                raise ValueError(msg)
            # If it looks like momentum but couldn't be parsed, raise clear error
            msg = (
                f"Invalid momentum field format: '{field_name}'. "
                f"Expected 'pi_N' or 'piN' where N is a numeric index (e.g., 'pi_0', 'pi0')."
            )
            raise ValueError(msg)

        # Regular field lookup via slot map
        slot = self._field_slot_map.get(field_name)
        if slot is not None:
            field = state[slot]
            if not isinstance(field, ScalarField):
                msg = f"Expected ScalarField, got {type(field).__name__}"
                raise TypeError(msg)
            return field

        msg = f"Unknown field name: {field_name}"
        raise ValueError(msg)

    def _compute_rhs_for_component(  # noqa: C901, PLR0912, PLR0914
        self,
        component_idx: int,
        state: FieldCollection,
        bc: BCDescriptor,
        t: float = 0.0,
        virtual_momenta: dict[str, ScalarField] | None = None,
    ) -> ScalarField:
        """Compute the RHS for a single component's equation.

        This method evaluates all terms in the component's equation specification
        and sums them together.

        Parameters
        ----------
        component_idx : int
            Index of the component.
        state : FieldCollection
            Current state (all fields and momenta).
        bc : BCDescriptor
            Boundary condition specification.
        t : float
            Current simulation time (for time-dependent coefficients in curved spacetime).
        virtual_momenta : dict[str, ScalarField] | None
            Pre-computed RHS fields for first-order components. When a term
            references ``first_derivative_t`` or ``pi_N`` for a first-order
            component, the virtual momentum is used instead of a state variable.

        Returns
        -------
        ScalarField
            The computed RHS.

        Raises
        ------
        TypeError
            If a state element is not a ``ScalarField``.
        ValueError
            If a field or operator cannot be resolved.
        """
        eq = self.spec.equations[component_idx]
        grid = state.grid

        # C2: Pre-extract spatial coordinate arrays once for all position-dependent terms
        coord_arrays: dict[str, NumericArray] | None = None
        if any(term.position_dependent for term in eq.rhs_terms):
            spatial_coords = self.spec.spatial_coordinates
            raw_coords = grid.cell_coords  # pyright: ignore[reportUnknownVariableType]
            coord_arrays = {
                name: np.asarray(raw_coords[..., i], dtype=np.float64)  # pyright: ignore[reportUnknownArgumentType]
                for i, name in enumerate(spatial_coords[: grid.num_axes])
            }

        # C1: Per-timestep coefficient cache — same symbolic expression at
        # same (t, grid) produces the same value, so deduplicate eval() calls
        coeff_cache: dict[str, float | NumericArray] = {}

        # Start with zero field
        result = ScalarField(grid, data=0.0)

        # Sum all terms from the specification
        for term_idx, term in enumerate(eq.rhs_terms):
            target_field_name = term.field

            # Handle first_derivative_t operator specially
            # d_t(field) = momentum (second-order) or RHS (first-order via virtual_momenta)
            if term.operator == "first_derivative_t":
                target_idx = self._component_name_to_index.get(target_field_name)
                if target_idx is None:
                    msg = f"Unknown field for first_derivative_t: {target_field_name}"
                    raise ValueError(msg)

                comp_name = self.spec.component_names[target_idx]

                # Try state momentum slot first (second-order fields)
                slot = self._momentum_slot_map.get(comp_name)
                if slot is not None:
                    momentum = state[slot]
                    if not isinstance(momentum, ScalarField):
                        msg = f"Expected ScalarField for momentum, got {type(momentum).__name__}"
                        raise TypeError(msg)
                    operated = momentum.copy()
                elif virtual_momenta is not None and comp_name in virtual_momenta:
                    operated = virtual_momenta[comp_name].copy()
                else:
                    eq_target = self.spec.equations[target_idx]
                    msg = (
                        f"first_derivative_t references '{target_field_name}' "
                        f"but it has time_derivative_order={eq_target.time_derivative_order} "
                        f"and no virtual momentum was computed."
                    )
                    raise ValueError(msg)
            else:
                # Standard operator handling (passes virtual_momenta for pi_N resolution)
                target_field = self._get_field_from_state(
                    state, target_field_name, virtual_momenta
                )
                operated = self._get_operator(term.operator, target_field, bc)

            # Resolve coefficient: B4 preresolved → C1 timestep cache → full eval
            preresolved = self._preresolved.get((component_idx, term_idx))
            if preresolved is not None:
                coefficient: float | NumericArray = preresolved
            else:
                cache_key = term.coefficient_symbolic
                if cache_key is not None and cache_key in coeff_cache:
                    coefficient = coeff_cache[cache_key]
                else:
                    coefficient = self._resolve_coefficient_at_point(
                        term, t, grid, coord_arrays=coord_arrays
                    )
                    if cache_key is not None:
                        coeff_cache[cache_key] = coefficient

            # Add coefficient * operated to result
            if isinstance(coefficient, np.ndarray):
                contribution = ScalarField(grid, data=coefficient * operated.data)
            else:
                contribution = coefficient * operated
            result += contribution

        return result

    def _build_constraint_bc(
        self,
        config: ConstraintSolverConfig,
        grid: GridBase,
    ) -> Any:  # noqa: ANN401
        """Convert a ConstraintSolverConfig to a py-pde boundary condition.

        Parameters
        ----------
        config : ConstraintSolverConfig
            Constraint solver configuration with per-axis BCs.
        grid : GridBase
            The simulation grid (used to check periodicity).

        Returns
        -------
        str | dict
            py-pde boundary condition descriptor.
        """
        spatial_coords = self.spec.spatial_coordinates

        # If all BCs are periodic and grid is periodic, use shorthand
        all_periodic = all(
            config.boundary_conditions.get(coord, BoundaryCondition("periodic")).type
            == "periodic"
            for coord in spatial_coords
        )
        if all_periodic and hasattr(grid, "periodic") and all(grid.periodic):
            return "auto_periodic_neumann"

        # Build explicit per-axis BC dict using coordinate names
        # py-pde expects: {"x": bc_x, "y": bc_y} where bc is "periodic"
        # or {"value": V} or {"derivative": D}
        bc_dict: dict[str, Any] = {}
        for i, coord in enumerate(spatial_coords):
            bc_config = config.boundary_conditions.get(coord)
            if bc_config is None:
                # Default: periodic if grid is periodic on this axis, else Neumann
                if hasattr(grid, "periodic") and grid.periodic[i]:
                    bc_dict[coord] = "periodic"
                else:
                    bc_dict[coord] = {"derivative": 0.0}
            elif bc_config.type == "periodic":
                bc_dict[coord] = "periodic"
            elif bc_config.type == "dirichlet":
                bc_dict[coord] = {
                    "value": bc_config.value if bc_config.value is not None else 0.0
                }
            elif bc_config.type == "neumann":
                bc_dict[coord] = {
                    "derivative": bc_config.derivative
                    if bc_config.derivative is not None
                    else 0.0
                }

        return bc_dict

    def _partition_constraint_terms(
        self,
        component_idx: int,
    ) -> tuple[list[OperatorTerm], list[OperatorTerm]]:
        """Separate an equation's RHS into self-terms and source-terms.

        A self-term is any term where ``term.field == eq.field_name`` (the
        operator acts on the field being solved for). All other terms are
        source-terms (cross-field contributions or momentum references).

        Parameters
        ----------
        component_idx : int
            Index of the constraint equation in ``spec.equations``.

        Returns
        -------
        tuple[list[OperatorTerm], list[OperatorTerm]]
            (self_terms, source_terms) partitioning of the RHS.
        """
        eq = self.spec.equations[component_idx]
        self_terms: list[OperatorTerm] = []
        source_terms: list[OperatorTerm] = []
        for term in eq.rhs_terms:
            if term.field == eq.field_name:
                self_terms.append(term)
            else:
                source_terms.append(term)
        return self_terms, source_terms

    def _compute_constraint_source(  # noqa: PLR0913, PLR0917
        self,
        component_idx: int,  # noqa: ARG002
        state: FieldCollection,
        bc: BCDescriptor,
        t: float,
        source_terms: list[OperatorTerm],
        virtual_momenta: dict[str, ScalarField] | None = None,
    ) -> ScalarField:
        """Evaluate the cross-field source for a constraint equation.

        Computes S = sum_i(coeff_i * operator_i(field_i)) for all
        source-terms (cross-field and momentum references).

        Parameters
        ----------
        component_idx : int
            Index of the constraint equation in ``spec.equations``.
        state : FieldCollection
            Current state (fields + momenta).
        bc : BCDescriptor
            Boundary conditions for operator evaluation.
        t : float
            Current time (for time-dependent coefficients).
        source_terms : list[OperatorTerm]
            The cross-field terms to evaluate.
        virtual_momenta : dict[str, ScalarField] | None
            Virtual momenta for constraint/first-order fields (always zero
            for constraint fields). Needed to resolve ``pi_N`` references
            where ``h_N`` is a constraint field.

        Returns
        -------
        ScalarField
            The evaluated source field.
        """
        grid = state.grid
        rhs_source = ScalarField(grid, data=0.0)
        for term in source_terms:
            target_field = self._get_field_from_state(
                state, term.field, virtual_momenta
            )
            operated = self._get_operator(term.operator, target_field, bc)
            coefficient = self._resolve_coefficient_at_point(term, t, grid)

            if isinstance(coefficient, np.ndarray):
                contribution = ScalarField(grid, data=coefficient * operated.data)
            else:
                contribution = coefficient * operated
            rhs_source += contribution
        return rhs_source

    def _solve_constraint_equation(
        self,
        component_idx: int,
        state: FieldCollection,
        bc: BCDescriptor,
        t: float,
        virtual_momenta: dict[str, ScalarField] | None = None,
    ) -> FieldCollection:
        """Solve a constraint equation and update the state.

        Dispatches to the appropriate solver based on the constraint's
        ``method`` configuration:

        - **"poisson"**: Original py-pde Poisson solver. Requires exactly
          one ``laplacian(self_field)`` self-term. Backward compatible.

        - **"auto"**: Automatically selects the best solver:
          - FFT for fully periodic grids (O(N log N))
          - Sparse matrix for non-periodic grids (O(N) via LU)

        - **"fft"**: Force FFT solver (requires periodic grid).

        - **"matrix"**: Force sparse matrix solver.

        The unified (auto/fft/matrix) solver handles ANY linear constraint:
        Poisson, Helmholtz, algebraic, anisotropic, partial Helmholtz, or
        any linear combination of spatial operators acting on the field.

        Parameters
        ----------
        component_idx : int
            Index of the constraint equation in ``spec.equations``.
        state : FieldCollection
            Current state. Updated in-place with the solved field.
        bc : BCDescriptor
            Boundary conditions for evaluating source-term operators.
        t : float
            Current time (for time-dependent coefficients).
        virtual_momenta : dict[str, ScalarField] | None
            Virtual momenta for constraint/first-order fields. Needed to
            resolve ``pi_N`` references where field N is a constraint.

        Returns
        -------
        FieldCollection
            Updated state with the constraint field solved.

        Raises
        ------
        ValueError
            If the equation has no self-referencing terms, if the operator
            is singular, or if the solver method is unsupported.
        """
        eq = self.spec.equations[component_idx]
        method = eq.constraint_solver.method

        if method == "poisson":
            return self._solve_constraint_poisson(
                component_idx, state, bc, t, virtual_momenta
            )
        if method in {"auto", "fft", "matrix"}:
            return self._solve_constraint_unified(
                component_idx, state, bc, t, virtual_momenta
            )

        msg = (
            f"Unknown constraint solver method '{method}' for {eq.field_name}. "
            f"Expected 'auto', 'fft', 'matrix', or 'poisson'."
        )
        raise ValueError(msg)

    def _solve_constraint_poisson(
        self,
        component_idx: int,
        state: FieldCollection,
        bc: BCDescriptor,
        t: float,
        virtual_momenta: dict[str, ScalarField] | None = None,
    ) -> FieldCollection:
        """Solve a Poisson-type constraint using py-pde's built-in solver.

        This is the original solver, preserved for backward compatibility.
        Requires the equation to have exactly one laplacian(self_field) term.

        Equation form: 0 = c * laplacian(field) + source_terms
        Rearranged to: laplacian(field) = -source_terms / c

        Raises
        ------
        ValueError
            If the equation lacks a laplacian self-term, has zero coefficient,
            or the Poisson solver fails.
        TypeError
            If the computed RHS is not a ScalarField.
        """
        from pde import solve_poisson_equation  # noqa: PLC0415, I001  # type: ignore[reportUnknownVariableType]

        eq = self.spec.equations[component_idx]
        grid = state.grid
        field_slot = self._field_slot_map[eq.field_name]

        # Separate RHS into the laplacian-of-self term and source terms
        laplacian_coeff: float | None = None
        source_terms: list[OperatorTerm] = []

        for term in eq.rhs_terms:
            if term.operator == "laplacian" and term.field == eq.field_name:
                if laplacian_coeff is not None:
                    msg = (
                        f"Multiple laplacian({eq.field_name}) terms in constraint "
                        f"equation. Expected exactly one."
                    )
                    raise ValueError(msg)
                laplacian_coeff = self._resolve_coefficient(term)
            else:
                source_terms.append(term)

        # Warn if non-laplacian self-terms exist (Helmholtz-type equation)
        non_lap_self = [
            t
            for t in eq.rhs_terms
            if t.field == eq.field_name and t.operator != "laplacian"
        ]
        if non_lap_self:
            ops = [t.operator for t in non_lap_self]
            warnings.warn(
                f"Constraint {eq.field_name} has non-laplacian self-terms "
                f"({ops}). The 'poisson' method only handles pure Poisson "
                f"equations correctly. Use method='auto' for Helmholtz or "
                f"other mixed-operator constraints.",
                stacklevel=2,
            )

        # Validate equation structure
        if laplacian_coeff is None:
            msg = (
                f"Constraint equation for {eq.field_name} lacks a "
                f"laplacian({eq.field_name}) term. "
                f"The Poisson solver requires the form: "
                f"laplacian(field) + source = 0. "
                f"Use method='auto' for general constraint types."
            )
            raise ValueError(msg)

        if abs(laplacian_coeff) < self._constraint_eps:
            msg = (
                f"Laplacian coefficient for {eq.field_name} is effectively "
                f"zero ({laplacian_coeff}). Cannot solve elliptic equation."
            )
            raise ValueError(msg)

        # Compute source
        rhs_source = self._compute_constraint_source(
            component_idx, state, bc, t, source_terms, virtual_momenta
        )

        # Rearrange to Poisson form: nabla^2(phi) = -S / laplacian_coeff
        poisson_rhs = -rhs_source / laplacian_coeff
        if not isinstance(poisson_rhs, ScalarField):
            msg = f"Expected ScalarField for Poisson RHS, got {type(poisson_rhs).__name__}"
            raise TypeError(msg)

        solver_bc = self._build_constraint_bc(eq.constraint_solver, grid)

        try:
            solution = solve_poisson_equation(
                rhs=poisson_rhs,
                bc=solver_bc,
                label=eq.field_name,
            )
        except Exception as e:
            rhs_max = float(np.max(np.abs(poisson_rhs.data)))
            msg = (
                f"Poisson solver failed for constraint {eq.field_name}:\n"
                f"  RHS max |f|: {rhs_max:.3e}\n"
                f"  BC: {solver_bc}\n"
                f"  Error: {e}"
            )
            raise ValueError(msg) from e

        state[field_slot].data[:] = solution.data
        return state

    def _solve_constraint_unified(
        self,
        component_idx: int,
        state: FieldCollection,
        bc: BCDescriptor,
        t: float,
        virtual_momenta: dict[str, ScalarField] | None = None,
    ) -> FieldCollection:
        """Solve a general linear constraint via operator matrix assembly.

        For a constraint equation of the form::

            0 = sum_i(c_i * op_i(self_field)) + source_terms

        This assembles the combined operator A = sum_i(c_i * M_i) where M_i
        is the sparse matrix for op_i, and solves A @ field = -source.

        This handles ANY linear constraint type: Poisson, Helmholtz,
        algebraic, anisotropic, partial Helmholtz, and arbitrary linear
        combinations of spatial operators.

        Parameters
        ----------
        component_idx : int
            Index of the constraint equation in ``spec.equations``.
        state : FieldCollection
            Current state (updated in-place).
        bc : BCDescriptor
            Boundary conditions for operator evaluation.
        t : float
            Current time (for time-dependent coefficients).
        virtual_momenta : dict[str, ScalarField] | None
            Virtual momenta for constraint/first-order fields.

        Returns
        -------
        FieldCollection
            Updated state with the constraint field solved.

        Raises
        ------
        ValueError
            If no self-terms, unknown operator, or singular operator.
        """
        eq = self.spec.equations[component_idx]
        grid = state.grid
        field_slot = self._field_slot_map[eq.field_name]
        method = eq.constraint_solver.method

        # 1. Partition terms
        self_terms, source_terms = self._partition_constraint_terms(component_idx)

        if not self_terms:
            msg = (
                f"Constraint equation for {eq.field_name} has no "
                f"self-referencing terms (no operator acts on {eq.field_name}). "
                f"Cannot solve: equation does not depend on the unknown field."
            )
            raise ValueError(msg)

        # 2. Compute cross-field source
        source = self._compute_constraint_source(
            component_idx, state, bc, t, source_terms, virtual_momenta
        )

        # 3. Choose solver path
        all_periodic = hasattr(grid, "periodic") and all(grid.periodic)
        use_fft = (method == "fft") or (method == "auto" and all_periodic)

        if method == "fft" and not all_periodic:
            msg = (
                f"FFT constraint solver for {eq.field_name} requires a "
                f"fully periodic grid, but the grid is not periodic."
            )
            raise ValueError(msg)

        if use_fft:
            solution_data = self._solve_constraint_fft(eq, grid, self_terms, source, t)
        else:
            solution_data = self._solve_constraint_matrix(
                eq, grid, self_terms, source, bc, t
            )

        # 4. Validate and update state
        if not np.isfinite(solution_data).all():
            msg = (
                f"Constraint solver for {eq.field_name} produced non-finite "
                f"values (NaN or Inf). This typically indicates a singular or "
                f"near-singular operator. Check equation structure and "
                f"coefficient values."
            )
            raise ValueError(msg)
        state[field_slot].data[:] = solution_data
        return state

    def _solve_constraint_fft(  # noqa: PLR0914
        self,
        eq: ComponentEquation,
        grid: GridBase,
        self_terms: list[OperatorTerm],
        source: ScalarField,
        t: float,
    ) -> np.ndarray:
        """Solve constraint in Fourier space (periodic grids only).

        Assembles the combined operator symbol in Fourier space as
        sum_i(c_i * symbol_i(k)), then solves field_hat = -source_hat / symbol.

        Uses discrete eigenvalues of the finite-difference stencils to
        match the matrix-based operators exactly.

        Returns
        -------
        np.ndarray
            The solved field data, shaped to match the grid.

        Raises
        ------
        ValueError
            If the combined symbol is singular (zero at any wavenumber)
            or if an operator has no FFT multiplier.
        """
        dx_array = np.array(grid.discretization)

        # Build wavenumber grids
        k_arrays: list[np.ndarray] = [
            np.fft.fftfreq(n, d=dx) * 2 * np.pi  # type: ignore[reportUnknownMemberType]
            for n, dx in zip(grid.shape, dx_array, strict=True)
        ]
        k_grids = list(np.meshgrid(*k_arrays, indexing="ij"))

        # Coalesce directional laplacians into compact laplacian when possible.
        # This uses the compact FFT multiplier (2cos(k*dx)-2)/dx² which matches
        # the standard 3-point stencil, avoiding wide-stencil truncation artifacts.
        self_terms = _coalesce_directional_laplacians(
            self_terms, self.spec.spatial_dimension
        )

        # Sum Fourier multipliers for all self-terms
        combined_multiplier = np.zeros(grid.shape, dtype=complex)
        for term in self_terms:
            coeff = self._resolve_coefficient_at_point(term, t, grid)
            multiplier_fn = _OPERATOR_FFT_MULTIPLIERS.get(term.operator)
            if multiplier_fn is None:
                msg = (
                    f"No FFT multiplier for operator '{term.operator}'. "
                    f"Use method='matrix' for this constraint."
                )
                raise ValueError(msg)
            term_mult = multiplier_fn(k_grids, dx_array)

            if isinstance(coeff, np.ndarray):
                # Position-dependent coefficient: FFT not applicable
                msg = (
                    f"Position-dependent coefficient on self-term "
                    f"'{term.operator}({eq.field_name})' is not compatible "
                    f"with FFT solver. Use method='matrix'."
                )
                raise ValueError(msg)  # noqa: TRY004

            combined_multiplier += coeff * term_mult

        # Handle null-space at singular modes (e.g., k=0 for pure Laplacian).
        # For Poisson-type equations on periodic grids, the zero-wavenumber
        # mode is in the null space (Laplacian eigenvalue = 0). If the source
        # is compatible (zero mean), we can set the k=0 solution mode to zero,
        # giving the unique zero-mean solution. If incompatible, we raise.
        source_hat = np.fft.fftn(-source.data)
        # Use relative threshold scaled by the operator's maximum magnitude
        # (with floor of 1.0 to handle near-zero operators gracefully)
        mult_scale = max(float(np.max(np.abs(combined_multiplier))), 1.0)
        singular_mask = np.abs(combined_multiplier) < self._constraint_eps * mult_scale
        n_singular = int(singular_mask.sum())

        if n_singular > 0:
            # Check if source is compatible with the null space:
            # source_hat must be zero at singular wavenumbers.
            source_at_singular = np.abs(source_hat[singular_mask])
            max_source_at_null = float(np.max(source_at_singular))
            if max_source_at_null > self._constraint_eps * float(
                np.max(np.abs(source_hat))
            ):
                msg = (
                    f"Constraint operator for {eq.field_name} is singular in "
                    f"Fourier space ({n_singular} null-space mode(s)). "
                    f"The source has non-zero projection onto the null space "
                    f"(max|source_hat|={max_source_at_null:.2e}). "
                    f"No solution exists."
                )
                raise ValueError(msg)
            # Set solution to zero at singular modes (zero-mean solution)
            safe_multiplier = np.where(singular_mask, 1.0, combined_multiplier)
            solution_hat = source_hat / safe_multiplier
            solution_hat[singular_mask] = 0.0
        else:
            solution_hat = source_hat / combined_multiplier

        return np.fft.ifftn(solution_hat).real  # type: ignore[return-value]

    def _solve_constraint_matrix(  # noqa: PLR0913, PLR0917
        self,
        eq: ComponentEquation,
        grid: GridBase,
        self_terms: list[OperatorTerm],
        source: ScalarField,
        bc: BCDescriptor,  # noqa: ARG002
        t: float,
    ) -> np.ndarray:
        """Solve constraint via sparse matrix assembly and direct solve.

        Assembles A = sum_i(c_i * M_i) where M_i is the sparse matrix for
        operator_i, then solves A @ field = -source using scipy.sparse.linalg.spsolve.

        Works with any boundary conditions (periodic, Dirichlet, Neumann).

        Returns
        -------
        np.ndarray
            The solved field data, shaped to match the grid.

        Raises
        ------
        ValueError
            If an operator has no matrix builder, or the solve fails.
        """
        solver_bc = self._build_constraint_bc(eq.constraint_solver, grid)
        bcs = grid.get_boundary_conditions(solver_bc)

        self_terms = _coalesce_directional_laplacians(
            self_terms, self.spec.spatial_dimension
        )

        n = math.prod(grid.shape)
        combined_matrix = sparse.dok_matrix((n, n))
        combined_vector = sparse.dok_matrix((n, 1))

        for term in self_terms:
            coeff = self._resolve_coefficient_at_point(term, t, grid)
            builder = _OPERATOR_MATRIX_REGISTRY.get(term.operator)
            if builder is None:
                msg = (
                    f"No matrix builder for operator '{term.operator}'. "
                    f"Register it in _OPERATOR_MATRIX_REGISTRY to enable "
                    f"constraint solving for this operator type."
                )
                raise ValueError(msg)
            op_matrix, op_vector = builder(grid, bcs)

            if isinstance(coeff, np.ndarray):
                # Position-dependent coefficient: multiply each row by coeff
                diag = sparse.diags(coeff.ravel())
                combined_matrix += diag @ op_matrix
                combined_vector += diag @ op_vector
            else:
                combined_matrix += coeff * op_matrix
                combined_vector += coeff * op_vector

        # Convert to CSC for efficient direct solve
        a_csc = combined_matrix.tocsc()
        rhs = (
            np.ravel(-source.data) - np.asarray(combined_vector.toarray()).ravel()  # type: ignore[reportUnknownArgumentType]
        )

        try:
            solution: np.ndarray = spsolve(a_csc, rhs)  # type: ignore[reportUnknownVariableType]
        except (RuntimeError, ValueError) as e:
            msg = (
                f"Matrix solver failed for constraint {eq.field_name}:\n"
                f"  Matrix shape: {a_csc.shape}, nnz: {a_csc.nnz}\n"
                f"  RHS max: {float(np.max(np.abs(rhs))):.3e}\n"
                f"  Error: {e}"
            )
            raise ValueError(msg) from e

        return solution.reshape(grid.shape)  # type: ignore[return-value]

    def _has_coupled_constraints(
        self,
        enabled_indices: list[int],
    ) -> bool:
        """Check if any enabled constraints reference each other's fields.

        Two constraints are coupled if constraint A's source terms reference
        constraint B's field, or vice versa. When constraints are coupled,
        a single pass may not converge — Gauss-Seidel iteration is needed.

        Parameters
        ----------
        enabled_indices : list[int]
            Indices of enabled constraint equations.

        Returns
        -------
        bool
            True if mutual coupling is detected.
        """
        enabled_fields = {self.spec.equations[i].field_name for i in enabled_indices}
        for i in enabled_indices:
            eq = self.spec.equations[i]
            _, source_terms = self._partition_constraint_terms(i)
            for term in source_terms:
                # Check if this source term references another enabled constraint
                ref_field = term.field
                # Handle momentum references (pi_N → field N)
                if ref_field.startswith("pi"):
                    momentum_idx = parse_momentum_field_name(ref_field)
                    if (
                        momentum_idx is not None
                        and 0 <= momentum_idx < self.n_components
                    ):
                        ref_field = self.spec.component_names[momentum_idx]
                if ref_field in enabled_fields and ref_field != eq.field_name:
                    return True
        return False

    def _evolve_constraints(
        self,
        state: TState,
        bc: BCDescriptor,
        t: float,
    ) -> tuple[TState, dict[str, ScalarField]]:
        """Solve constraint equations (time_derivative_order=0).

        Constraints are solved elliptically when enabled, updating the state
        in-place so that dynamical equations see the resolved fields.

        When multiple enabled constraints reference each other's fields
        (coupled constraints), Gauss-Seidel iteration is used: constraints
        are solved sequentially in a loop until convergence or the maximum
        iteration count is reached.

        Returns
        -------
        tuple[TState, dict[str, ScalarField]]
            Updated state and virtual momenta dict with zero entries for
            each constraint field.

        Raises
        ------
        TypeError
            If state is not a FieldCollection.
        """
        # Type narrowing: constraints require FieldCollection
        if not isinstance(state, FieldCollection):
            msg = "Constraint solving requires FieldCollection state"
            raise TypeError(msg)

        grid = state.grid
        virtual_momenta: dict[str, ScalarField] = {}

        # Collect enabled constraint indices
        enabled_indices: list[int] = []
        for i, eq in enumerate(self.spec.equations):
            if eq.time_derivative_order == 0:
                if eq.constraint_solver.enabled:
                    enabled_indices.append(i)
                virtual_momenta[eq.field_name] = ScalarField(grid, data=0.0)

        if not enabled_indices:
            return state, virtual_momenta

        # Check for coupled constraints
        if self._has_coupled_constraints(enabled_indices):
            state = self._solve_coupled_constraints(
                enabled_indices, state, bc, t, virtual_momenta
            )
        else:
            # Uncoupled: single pass suffices
            for i in enabled_indices:
                state = self._solve_constraint_equation(
                    i, state, bc, t, virtual_momenta
                )

        return state, virtual_momenta

    def _solve_coupled_constraints(
        self,
        enabled_indices: list[int],
        state: FieldCollection,
        bc: BCDescriptor,
        t: float,
        virtual_momenta: dict[str, ScalarField] | None = None,
    ) -> FieldCollection:
        """Solve mutually-coupled constraints via block FFT or Gauss-Seidel.

        For fully periodic grids, uses Fourier-space block solve: at each
        wavenumber k, assembles a small NxN dense matrix (N = number of
        coupled constraints) and solves exactly. This handles arbitrary
        coupling strength including Laplacian cross-terms.

        For non-periodic grids, falls back to Gauss-Seidel iteration
        (sequential single-field solves until convergence).

        Parameters
        ----------
        enabled_indices : list[int]
            Indices of enabled constraint equations.
        state : FieldCollection
            Current state (updated in-place during iteration).
        bc : BCDescriptor
            Boundary conditions.
        t : float
            Current time.
        virtual_momenta : dict[str, ScalarField] | None
            Virtual momenta for constraint/first-order fields.

        Returns
        -------
        FieldCollection
            State with solved constraint fields.
        """
        grid = state.grid
        all_periodic = hasattr(grid, "periodic") and all(grid.periodic)

        if all_periodic:
            return self._solve_coupled_constraints_fft(
                enabled_indices, state, bc, t, virtual_momenta
            )
        return self._solve_coupled_constraints_gauss_seidel(
            enabled_indices, state, bc, t, virtual_momenta
        )

    def _solve_coupled_constraints_fft(  # noqa: C901, PLR0912, PLR0914, PLR0915
        self,
        enabled_indices: list[int],
        state: FieldCollection,
        bc: BCDescriptor,
        t: float,
        virtual_momenta: dict[str, ScalarField] | None = None,
    ) -> FieldCollection:
        """Solve coupled constraints exactly via Fourier-space block solve.

        At each wavenumber k, assembles a dense NxN coupling matrix M(k)
        where M[i,j] = sum of (coefficient * Fourier_multiplier(operator))
        for all terms in constraint_i that reference constraint_j. Then
        solves M(k) @ f_hat(k) = -s_hat(k) where s is the external source
        (non-constraint field references).

        This is exact (no iteration needed) and handles arbitrary coupling
        strength. Complexity: O(N_grid * N_fields^3 * log(N_grid)).

        Parameters
        ----------
        enabled_indices : list[int]
            Indices of coupled constraint equations.
        state : FieldCollection
            Current state (updated in-place).
        bc : BCDescriptor
            Boundary conditions for operator evaluation.
        t : float
            Current time.
        virtual_momenta : dict[str, ScalarField] | None
            Virtual momenta for constraint/first-order fields.

        Returns
        -------
        FieldCollection
            State with solved constraint fields.

        Raises
        ------
        ValueError
            If an operator has no FFT multiplier or the system is singular.
        """
        grid = state.grid
        n_constraints = len(enabled_indices)
        enabled_names = {self.spec.equations[i].field_name for i in enabled_indices}

        # Build wavenumber grids
        dx_array = np.array(grid.discretization)
        k_arrays: list[np.ndarray] = [
            np.fft.fftfreq(n, d=dx) * 2 * np.pi  # type: ignore[reportUnknownMemberType]
            for n, dx in zip(grid.shape, dx_array, strict=True)
        ]
        k_grids = np.meshgrid(*k_arrays, indexing="ij")

        # Map constraint field names to local indices (0..n_constraints-1)
        name_to_local: dict[str, int] = {}
        for local_idx, comp_idx in enumerate(enabled_indices):
            name_to_local[self.spec.equations[comp_idx].field_name] = local_idx

        # Build M(k) coupling matrix and s(k) source vector in Fourier space
        grid_shape = tuple(grid.shape)
        # M_hat[i, j, ...grid_shape] = coupling multiplier at each k
        m_hat = np.zeros((n_constraints, n_constraints, *grid_shape), dtype=complex)
        s_hat = np.zeros((n_constraints, *grid_shape), dtype=complex)

        for local_i, comp_idx in enumerate(enabled_indices):
            eq = self.spec.equations[comp_idx]

            # Separate terms into intra-cluster (go into M) and external (go into s)
            external_terms: list[OperatorTerm] = []

            for term in eq.rhs_terms:
                ref_field = term.field

                # Check if this is a momentum reference (pi_N → field h_N)
                momentum_idx = parse_momentum_field_name(ref_field)
                if momentum_idx is not None:
                    if 0 <= momentum_idx < self.n_components:
                        ref_comp_name = self.spec.component_names[momentum_idx]
                        if ref_comp_name in enabled_names:
                            # Momentum of a constraint field is always zero → skip
                            continue
                    # Non-constraint momentum → external source
                    external_terms.append(term)
                    continue

                # Check if this references another constraint in the cluster
                if ref_field in enabled_names:
                    local_j = name_to_local[ref_field]
                    # Get Fourier multiplier for this operator
                    multiplier_fn = _OPERATOR_FFT_MULTIPLIERS.get(term.operator)
                    if multiplier_fn is None:
                        msg = (
                            f"No FFT multiplier for operator '{term.operator}' "
                            f"in coupled constraint {eq.field_name}."
                        )
                        raise ValueError(msg)
                    coeff = self._resolve_coefficient_at_point(term, t, grid)
                    if isinstance(coeff, np.ndarray):
                        msg = (
                            f"Position-dependent coefficient on self-term "
                            f"'{term.operator}({ref_field})' in coupled "
                            f"constraint {eq.field_name} is not compatible "
                            f"with FFT block solver. Use method='matrix'."
                        )
                        raise ValueError(msg)
                    multiplier = multiplier_fn(k_grids, dx_array)
                    m_hat[local_i, local_j] += coeff * multiplier
                else:
                    external_terms.append(term)

            # Compute external source in physical space, then FFT
            if external_terms:
                source = self._compute_constraint_source(
                    comp_idx, state, bc, t, external_terms, virtual_momenta
                )
                s_hat[local_i] = np.fft.fftn(source.data)

        # Solve M(k) @ f_hat(k) = -s_hat(k) via SVD with Tikhonov regularization.
        # This handles near-singular wavenumbers (e.g., Helmholtz resonance at
        # k² ≈ m²) by attenuating rather than amplifying near-null modes.
        m_hat_transposed = np.moveaxis(m_hat, [0, 1], [-2, -1])  # (..., n, n)
        s_hat_transposed = np.moveaxis(s_hat, 0, -1)  # (..., n)
        rhs = -s_hat_transposed  # (..., n)

        # Batched SVD: M = U @ diag(S) @ Vh at each grid point
        try:
            u, s, vh = np.linalg.svd(m_hat_transposed, full_matrices=False)
        except np.linalg.LinAlgError as e:
            field_names = [self.spec.equations[i].field_name for i in enabled_indices]
            msg = (
                f"SVD failed for coupled constraint system in Fourier "
                f"space. Fields: {field_names}. Error: {e}"
            )
            raise ValueError(msg) from e

        # Tikhonov regularization: S_reg_inv = S / (S^2 + alpha^2)
        # alpha = rcond * max(S) ensures well-conditioned modes are unaffected
        # while near-singular modes are smoothly attenuated.
        s_max = float(np.max(s))
        alpha = self._coupled_svd_rcond * max(s_max, 1e-30)
        alpha_sq = alpha * alpha
        s_reg_inv = s / (s * s + alpha_sq)

        n_regularized = int((s < alpha).sum())
        total_svs = math.prod(s.shape)

        if n_regularized > 0:
            field_names = [self.spec.equations[i].field_name for i in enabled_indices]
            logger.debug(
                "Coupled FFT solver: %d/%d singular values (%.1f%%) below "
                "Tikhonov threshold %.2e (max SV %.2e). Fields: %s.",
                n_regularized,
                total_svs,
                100.0 * n_regularized / max(total_svs, 1),
                alpha,
                s_max,
                field_names,
            )

            # Only emit visible warning for truly pathological systems
            # where the majority of singular values need regularization.
            if n_regularized > total_svs // 2:
                warnings.warn(
                    f"Coupled FFT solver: {n_regularized}/{total_svs} "
                    f"singular values "
                    f"({100.0 * n_regularized / total_svs:.0f}%) required "
                    f"Tikhonov regularization, indicating a severely "
                    f"ill-conditioned constraint system. "
                    f"Fields: {field_names}.",
                    stacklevel=2,
                )

        # Compute f = Vh^H @ diag(S_reg_inv) @ U^H @ rhs
        u_h_rhs = np.einsum("...ji,...j->...i", u.conj(), rhs)
        scaled = s_reg_inv * u_h_rhs
        f_hat_transposed = np.einsum("...ji,...j->...i", vh.conj(), scaled)

        # Move back to (n_constraints, grid...) layout
        f_hat = np.moveaxis(f_hat_transposed, -1, 0)

        for local_i, comp_idx in enumerate(enabled_indices):
            eq = self.spec.equations[comp_idx]
            field_slot = self._field_slot_map[eq.field_name]
            solution = np.fft.ifftn(f_hat[local_i]).real
            if not np.isfinite(solution).all():
                msg = (
                    f"Coupled FFT solver for {eq.field_name} produced "
                    f"non-finite values. The coupled system may be "
                    f"near-singular at some wavenumber."
                )
                raise ValueError(msg)
            state[field_slot].data[:] = solution

        return state

    def _solve_coupled_constraints_gauss_seidel(
        self,
        enabled_indices: list[int],
        state: FieldCollection,
        bc: BCDescriptor,
        t: float,
        virtual_momenta: dict[str, ScalarField] | None = None,
    ) -> FieldCollection:
        """Solve coupled constraints via Gauss-Seidel iteration.

        Each iteration solves all enabled constraints sequentially, using
        the most recent field values (Gauss-Seidel, not Jacobi). Iteration
        stops when:

        - All fields change by less than ``tolerance`` (convergence), or
        - ``max_iterations`` is reached (issues a warning).

        Note: Gauss-Seidel may not converge for strongly-coupled systems
        where cross-terms include Laplacian operators. For periodic grids,
        the FFT block solver is preferred.

        Parameters
        ----------
        enabled_indices : list[int]
            Indices of enabled constraint equations.
        state : FieldCollection
            Current state (updated in-place during iteration).
        bc : BCDescriptor
            Boundary conditions.
        t : float
            Current time.
        virtual_momenta : dict[str, ScalarField] | None
            Virtual momenta for constraint/first-order fields.

        Returns
        -------
        FieldCollection
            State with solved constraint fields.
        """
        # Use convergence parameters from the first enabled constraint
        first_config = self.spec.equations[enabled_indices[0]].constraint_solver
        max_iter = first_config.max_iterations
        tol = first_config.tolerance
        max_change = 0.0

        for _iteration in range(max_iter):
            # Save current field values for convergence check
            prev_data: dict[str, np.ndarray] = {}
            for i in enabled_indices:
                eq = self.spec.equations[i]
                field_slot = self._field_slot_map[eq.field_name]
                prev_data[eq.field_name] = state[field_slot].data.copy()

            # Solve all constraints sequentially (Gauss-Seidel)
            for i in enabled_indices:
                state = self._solve_constraint_equation(
                    i, state, bc, t, virtual_momenta
                )

            # Check convergence: max|field_new - field_old| relative to field scale
            max_change = 0.0
            max_magnitude = 0.0
            for i in enabled_indices:
                eq = self.spec.equations[i]
                field_slot = self._field_slot_map[eq.field_name]
                change = float(
                    np.max(np.abs(state[field_slot].data - prev_data[eq.field_name]))
                )
                max_change = max(max_change, change)
                max_magnitude = max(
                    max_magnitude,
                    float(np.max(np.abs(state[field_slot].data))),
                )

            # Use relative tolerance: scale threshold by field magnitude
            # (with floor of 1.0 to avoid issues with near-zero fields)
            effective_tol = tol * max(1.0, max_magnitude)
            if max_change < effective_tol:
                return state

        # Did not converge — warn but don't error
        field_names = [self.spec.equations[i].field_name for i in enabled_indices]
        warnings.warn(
            f"Coupled constraint iteration did not converge after "
            f"{max_iter} iterations (max_change={max_change:.2e}, "
            f"tolerance={tol:.2e}). Fields: {field_names}. "
            f"Increase max_iterations or tolerance in ConstraintSolverConfig.",
            stacklevel=2,
        )
        return state

    def _evolve_first_order(
        self,
        state: TState,
        bc: BCDescriptor,
        t: float,
        virtual_momenta: dict[str, ScalarField],
    ) -> dict[str, ScalarField]:
        """Compute virtual momenta for first-order components.

        The RHS of each first-order equation becomes a "virtual momentum"
        that second-order equations can reference via ``pi_N`` or
        ``first_derivative_t``.  This runs after constraints are solved
        so that first-order equations see updated constraint fields.

        Returns
        -------
        dict[str, ScalarField]
            Updated virtual momenta dict including first-order entries.

        Raises
        ------
        TypeError
            If state is not a FieldCollection.
        """
        # Type narrowing: RHS computation requires FieldCollection
        if not isinstance(state, FieldCollection):
            msg = "First-order evolution requires FieldCollection state"
            raise TypeError(msg)

        for i, eq in enumerate(self.spec.equations):
            if eq.time_derivative_order == 1:
                virtual_momenta[eq.field_name] = self._compute_rhs_for_component(
                    i, state, bc, t
                )
        return virtual_momenta

    def _evolve_second_order(
        self,
        state: TState,
        bc: BCDescriptor,
        t: float,
        virtual_momenta: dict[str, ScalarField],
        rates: list[ScalarField | None],
    ) -> list[ScalarField | None]:
        """Compute rates for all components from the virtual momenta.

        - Second-order: d/dt field = momentum, d/dt momentum = RHS
        - First-order: d/dt field = virtual momentum (already computed)
        - Constraint: d/dt field = 0

        Returns
        -------
        list[ScalarField | None]
            Populated rates array.

        Raises
        ------
        TypeError
            If momentum field is not a ScalarField or state is not FieldCollection.
        """
        # Type narrowing: second-order evolution requires FieldCollection
        if not isinstance(state, FieldCollection):
            msg = "Second-order evolution requires FieldCollection state"
            raise TypeError(msg)

        grid = state.grid
        for i, eq in enumerate(self.spec.equations):
            field_slot = self._field_slot_map[eq.field_name]

            if eq.time_derivative_order >= 2:  # noqa: PLR2004
                momentum_slot = self._momentum_slot_map[eq.field_name]
                momentum = state[momentum_slot]
                if not isinstance(momentum, ScalarField):
                    msg = f"Expected ScalarField for momentum, got {type(momentum).__name__}"
                    raise TypeError(msg)
                rates[field_slot] = momentum.copy()
                rates[momentum_slot] = self._compute_rhs_for_component(
                    i, state, bc, t, virtual_momenta
                )
            elif eq.time_derivative_order == 1:
                rates[field_slot] = virtual_momenta[eq.field_name]
            else:
                rates[field_slot] = ScalarField(grid, data=0.0)

        return rates

    @override
    def evolution_rate(
        self,
        state: TState,
        t: float = 0.0,
    ) -> FieldCollection:
        """Compute the time derivatives for all fields.

        Supports mixed time-derivative orders:
        - Second-order (wave): d/dt field = momentum, d/dt momentum = RHS
        - First-order (heat/diffusion): d/dt field = RHS
        - Constraint (elliptic, order=0): d/dt field = 0

        For first-order components, the computed RHS is stored as a "virtual
        momentum" so that second-order equations referencing ``pi_N`` or
        ``first_derivative_t`` of a first-order component can access it.

        Parameters
        ----------
        state : FieldCollection
            Current state with ``spec.state_size`` fields.
        t : float
            Current time. Used for time-dependent coefficients in curved spacetime
            (e.g., Hubble friction in de Sitter expansion).

        Returns
        -------
        FieldCollection
            Time derivatives for all fields.

        Raises
        ------
        TypeError
            If ``state`` is not a ``FieldCollection``.
        ValueError
            If the state size or grid dimension does not match the spec.
        """
        if not isinstance(state, FieldCollection):
            msg = f"Expected FieldCollection, got {type(state).__name__}"
            raise TypeError(msg)
        expected_fields = self.spec.state_size
        if len(state) != expected_fields:
            msg = f"Expected {expected_fields} fields, got {len(state)}"
            raise ValueError(msg)

        # Validate grid dimension matches spec
        grid_dim = state.grid.dim
        expected_dim = self.spec.spatial_dimension
        if grid_dim != expected_dim:
            msg = (
                f"Grid dimension {grid_dim} does not match spec "
                f"spatial_dimension {expected_dim}. "
                f"The equation system expects a {expected_dim}D spatial grid "
                f"(from {self.spec.dimension}D spacetime)."
            )
            raise ValueError(msg)

        # B5: Cache boundary conditions (same grid always produces same BCs)
        grid = state.grid
        grid_id = id(grid)
        if self._cached_bc is None or self._cached_grid_id != grid_id:
            self._cached_bc = infer_bc_from_grid(grid)
            self._cached_grid_id = grid_id
        bc = self._cached_bc

        # Pass 1a: Solve constraints
        state, virtual_momenta = self._evolve_constraints(state, bc, t)

        # Pass 1b: Compute first-order virtual momenta
        virtual_momenta = self._evolve_first_order(state, bc, t, virtual_momenta)

        # Pass 2: Build the full rates array
        rates: list[ScalarField | None] = [None] * expected_fields
        rates = self._evolve_second_order(state, bc, t, virtual_momenta, rates)

        return FieldCollection(rates)  # type: ignore[arg-type]

    def check_stability(self, dt: float, grid: GridBase) -> list[str]:
        """Check CFL / stability conditions for explicit time-stepping.

        Estimates maximum wave speeds and diffusivities from the equation
        coefficients and checks whether ``dt`` satisfies the CFL condition.

        Parameters
        ----------
        dt : float
            Proposed time step.
        grid : GridBase
            The spatial grid (used for cell spacing).

        Returns
        -------
        list[str]
            Warning messages for any violated conditions. Empty if stable.
        """
        warnings: list[str] = []
        dx_min = min(grid.discretization)

        for eq in self.spec.equations:
            if eq.time_derivative_order < 2:  # noqa: PLR2004
                continue

            # Estimate max wave speed from laplacian coefficients
            max_laplacian_coeff = 0.0
            max_diffusive_coeff = 0.0
            for term in eq.rhs_terms:
                coeff_abs = abs(term.coefficient)
                if (
                    term.operator == "laplacian" and term.field == eq.field_name
                ) or term.operator.startswith("laplacian_"):
                    max_laplacian_coeff = max(max_laplacian_coeff, coeff_abs)
                elif term.operator == "biharmonic":
                    max_diffusive_coeff = max(max_diffusive_coeff, coeff_abs)

            # CFL for wave equation: dt < dx / c where c = sqrt(laplacian_coeff)
            if max_laplacian_coeff > 0:
                c_max = math.sqrt(max_laplacian_coeff)
                cfl_limit = dx_min / c_max
                if dt > cfl_limit:
                    warnings.append(
                        f"CFL violated for {eq.field_name}: "
                        f"dt={dt:.3e} > dx/c={cfl_limit:.3e} "
                        f"(c={c_max:.3e}, dx={dx_min:.3e})"
                    )

            # Stability for biharmonic: dt < dx^4 / (2 * coeff)
            if max_diffusive_coeff > 0:
                biharm_limit = dx_min**4 / (2 * max_diffusive_coeff)
                if dt > biharm_limit:
                    warnings.append(
                        f"Biharmonic stability violated for {eq.field_name}: "
                        f"dt={dt:.3e} > dx^4/(2D)={biharm_limit:.3e}"
                    )

        return warnings


def build_pde_from_json(
    json_path: Path | str,
    parameters: dict[str, float] | None = None,
    *,
    constraint_eps: float = 1e-14,
    coupled_svd_rcond: float = 0.01,
) -> PDEFromSpec:
    """Build a PDE from a JSON equation specification file.

    This is the main entry point for the Lagrangian-to-PDE pipeline on the
    Python side. Given a JSON file exported from Mathematica/xAct, this
    function creates a py-pde compatible PDE class.

    Parameters
    ----------
    json_path : Path | str
        Path to the JSON file containing the equation specification.
    parameters : dict[str, float] | None
        Optional parameter values to override symbolic coefficients.
        Keys are symbolic names (e.g., "m2", "kappa"), values are numeric.
        Example: {"m2": 0.5, "kappa": 1.0}
    constraint_eps : float
        Tolerance for the constraint solver's Laplacian coefficient check.
        Default ``1e-14``.
    coupled_svd_rcond : float
        Relative singular-value threshold for Tikhonov regularization in
        the coupled FFT constraint solver.  Default ``0.01``.

    Returns
    -------
    PDEFromSpec
        A PDE instance ready for use with py-pde solvers.

    Examples
    --------
    >>> pde = build_pde_from_json("examples/data/em_1d.json")
    >>> # Create initial state and run simulation
    >>> from pde import CartesianGrid, ScalarField, FieldCollection
    >>> grid = CartesianGrid([(0, 100)], 256, periodic=True)
    >>> # ... create initial conditions and solve

    >>> # With custom parameter values:
    >>> pde = build_pde_from_json("examples/data/proca_1d.json", parameters={"m2": 2.0})
    """
    spec = load_equation_system(json_path)
    # Auto-load metadata parameter defaults when none explicitly provided
    if parameters is None:
        meta_params: dict[str, object] | None = spec.metadata.get("parameters")
        if isinstance(meta_params, dict):
            parameters = {
                k: float(v)
                for k, v in meta_params.items()
                if isinstance(v, (int, float))
            }
    return PDEFromSpec(
        spec,
        parameters=parameters,
        constraint_eps=constraint_eps,
        coupled_svd_rcond=coupled_svd_rcond,
    )


def create_initial_state(
    grid: GridBase,
    spec: EquationSystem,
    field_data: dict[str, NDArray[np.float64]] | None = None,
    momentum_data: dict[str, NDArray[np.float64]] | None = None,
) -> FieldCollection:
    """Create initial state for a PDEFromSpec simulation.

    The state layout is determined by ``spec.state_layout``:
    - Second-order components get [field, momentum] pairs
    - First-order/constraint components get [field] only

    Parameters
    ----------
    grid : GridBase
        The simulation grid.
    spec : EquationSystem
        The equation specification.
    field_data : dict[str, NDArray] | None
        Initial data for each field component. Keys are component names.
        Components not specified default to zero.
    momentum_data : dict[str, NDArray] | None
        Initial data for each momentum component. Keys are component names.
        Components not specified default to zero. Ignored for first-order
        and constraint components (which have no momentum slot).

    Returns
    -------
    FieldCollection
        Initial state with ``spec.state_size`` fields.
    """
    field_data = field_data or {}
    momentum_data = momentum_data or {}

    fields: list[ScalarField] = []

    for name, slot_type in spec.state_layout:
        if slot_type == "field":
            if name in field_data:
                fields.append(ScalarField(grid, data=field_data[name]))
            else:
                fields.append(ScalarField(grid, data=0.0))
        # Momentum slot (only for second-order components)
        elif name in momentum_data:
            fields.append(ScalarField(grid, data=momentum_data[name]))
        else:
            fields.append(ScalarField(grid, data=0.0))

    return FieldCollection(fields)
