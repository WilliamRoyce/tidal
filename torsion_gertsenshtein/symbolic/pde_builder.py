"""Build py-pde PDEBase subclasses from equation specifications.

This module provides the core functionality for converting symbolically-derived
field equations (loaded from JSON) into executable py-pde PDE classes.

The key principle is that NO physics is hardcoded here - all equation structure
comes from the specification that was derived from the Lagrangian.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from pde import FieldCollection, PDEBase, ScalarField
from typing_extensions import override

from torsion_gertsenshtein.kgsim.utils import infer_bc_from_grid
from torsion_gertsenshtein.symbolic.json_loader import (
    EquationSystem,
    load_equation_system,
)

if TYPE_CHECKING:
    from pathlib import Path

    import numpy as np
    from numpy.typing import NDArray
    from pde.grids.base import GridBase
    from pde.pdes.base import TState

    from torsion_gertsenshtein.utils import BCDescriptor

    NumericArray = NDArray[np.float64]


class PDEFromSpec(PDEBase):
    """Generic PDE class built from JSON equation specification.

    This class dynamically constructs the evolution equations from a parsed
    specification. NO physics is hardcoded - all equation structure comes
    from the EquationSystem that was derived from a Lagrangian.

    The state is organized as a FieldCollection with 2 * n_components fields:
        [field_0, momentum_0, field_1, momentum_1, ...]

    For each component i, the evolution equations are:
        d/dt field_i = momentum_i
        d/dt momentum_i = sum of RHS terms from specification

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
        Always False (autonomous system).

    Examples
    --------
    >>> from torsion_gertsenshtein.symbolic import load_equation_system
    >>> from torsion_gertsenshtein.symbolic.pde_builder import PDEFromSpec
    >>> spec = load_equation_system("examples/data/em_1d.json")
    >>> pde = PDEFromSpec(spec)
    >>> # pde can now be used with py-pde solvers
    """

    explicit_time_dependence = False

    def __init__(self, spec: EquationSystem) -> None:
        """Initialize PDE from equation specification.

        Parameters
        ----------
        spec : EquationSystem
            The equation specification loaded from JSON.
        """
        super().__init__()
        self.spec = spec
        self.n_components = spec.n_components
        self._component_name_to_index = {
            name: i for i, name in enumerate(spec.component_names)
        }

    @staticmethod
    def _get_operator(  # noqa: C901, PLR0911
        operator_name: str, field: ScalarField, bc: BCDescriptor
    ) -> ScalarField:
        """Apply a named operator to a field.

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
            If the operator is not recognized.
        """
        if operator_name == "laplacian":
            # laplace() returns ScalarField or DataFieldBase depending on context
            return field.laplace(bc=bc)
        if operator_name == "identity":
            return field.copy()
        if operator_name == "gradient_x":
            # gradient() returns FieldCollection, index to get ScalarField
            grad = field.gradient(bc=bc)
            component = grad[0]
            assert isinstance(component, ScalarField)
            return component
        if operator_name == "gradient_y":
            grad = field.gradient(bc=bc)
            component = grad[1]
            assert isinstance(component, ScalarField)
            return component
        if operator_name == "gradient_z":
            grad = field.gradient(bc=bc)
            component = grad[2]
            assert isinstance(component, ScalarField)
            return component
        if operator_name == "cross_derivative_xy":
            # d/dx d/dy f = d/dx (d/dy f)
            # This is NOT a laplacian - it's a cross spatial derivative
            if field.grid.dim < 2:  # noqa: PLR2004
                msg = "cross_derivative_xy requires at least 2D grid"
                raise ValueError(msg)
            # First compute d/dy
            grad_y = field.gradient(bc=bc)[1]
            assert isinstance(grad_y, ScalarField)
            # Then compute d/dx of that
            grad_xy = grad_y.gradient(bc=bc)[0]
            assert isinstance(grad_xy, ScalarField)
            return grad_xy
        if operator_name == "cross_derivative_xz":
            # d/dx d/dz f = d/dx (d/dz f)
            # Cross spatial derivative for 3+1D
            if field.grid.dim < 3:  # noqa: PLR2004
                msg = "cross_derivative_xz requires at least 3D grid"
                raise ValueError(msg)
            grad_z = field.gradient(bc=bc)[2]
            assert isinstance(grad_z, ScalarField)
            grad_xz = grad_z.gradient(bc=bc)[0]
            assert isinstance(grad_xz, ScalarField)
            return grad_xz
        if operator_name == "cross_derivative_yz":
            # d/dy d/dz f = d/dy (d/dz f)
            # Cross spatial derivative for 3+1D
            if field.grid.dim < 3:  # noqa: PLR2004
                msg = "cross_derivative_yz requires at least 3D grid"
                raise ValueError(msg)
            grad_z = field.gradient(bc=bc)[2]
            assert isinstance(grad_z, ScalarField)
            grad_yz = grad_z.gradient(bc=bc)[1]
            assert isinstance(grad_yz, ScalarField)
            return grad_yz
        if operator_name == "laplacian_x":
            # Pure ∂²/∂x² - second derivative in x only
            # For anisotropic equations like Navier-Cauchy elasticity
            grad_x = field.gradient(bc=bc)[0]
            assert isinstance(grad_x, ScalarField)
            d2_x = grad_x.gradient(bc=bc)[0]
            assert isinstance(d2_x, ScalarField)
            return d2_x
        if operator_name == "laplacian_y":
            # Pure ∂²/∂y² - second derivative in y only
            if field.grid.dim < 2:  # noqa: PLR2004
                msg = "laplacian_y requires at least 2D grid"
                raise ValueError(msg)
            grad_y = field.gradient(bc=bc)[1]
            assert isinstance(grad_y, ScalarField)
            d2_y = grad_y.gradient(bc=bc)[1]
            assert isinstance(d2_y, ScalarField)
            return d2_y
        if operator_name == "laplacian_z":
            # Pure ∂²/∂z² - second derivative in z only
            if field.grid.dim < 3:  # noqa: PLR2004
                msg = "laplacian_z requires at least 3D grid"
                raise ValueError(msg)
            grad_z = field.gradient(bc=bc)[2]
            assert isinstance(grad_z, ScalarField)
            d2_z = grad_z.gradient(bc=bc)[2]
            assert isinstance(d2_z, ScalarField)
            return d2_z

        msg = f"Unknown operator: {operator_name}"
        raise ValueError(msg)

    def _get_field_from_state(
        self, state: FieldCollection, field_name: str
    ) -> ScalarField:
        """Get a field from state by name, supporting both field and momentum names.

        For mixed time-space derivatives like d_t d_x A, the Wolfram pipeline
        expresses these as gradients of momentum fields (e.g., gradient_x(pi_0)).
        This is valid because d_t A = pi, so d_x(d_t A) = d_x(pi).

        The state is organized as: [field_0, pi_0, field_1, pi_1, ...]
        - Regular fields (A_0, A_1) are at even indices: state[2*i]
        - Momentum fields (pi_0, pi_1) are at odd indices: state[2*i + 1]

        Parameters
        ----------
        state : FieldCollection
            Current state with interleaved field/momentum pairs.
        field_name : str
            Name of the field to retrieve. Can be:
            - Regular field name like "A_0", "phi_0"
            - Momentum field name like "pi_0", "pi_1"

        Returns
        -------
        ScalarField
            The requested field from the state.

        Raises
        ------
        ValueError
            If the field name is not recognized.
        """
        # Check if this is a momentum field reference (pi_0, pi_1, etc.)
        if field_name.startswith("pi_"):
            parts = field_name.split("_")
            # Validate format: exactly "pi_N" where N is numeric
            if len(parts) != 2:  # noqa: PLR2004
                msg = (
                    f"Invalid momentum field format: '{field_name}'. "
                    f"Expected 'pi_N' where N is a numeric index (e.g., 'pi_0', 'pi_1')."
                )
                raise ValueError(msg)

            idx_str = parts[1]
            if not idx_str.isdigit():
                msg = (
                    f"Invalid momentum field index in '{field_name}'. "
                    f"Expected numeric index, got '{idx_str}'."
                )
                raise ValueError(msg)

            idx = int(idx_str)
            if not (0 <= idx < self.n_components):
                msg = (
                    f"Momentum field index {idx} out of range. "
                    f"This system has {self.n_components} components "
                    f"(valid indices: 0 to {self.n_components - 1}). "
                    f"Field reference: '{field_name}'."
                )
                raise ValueError(msg)

            # Momentum is at odd indices: state[2*i + 1]
            momentum = state[2 * idx + 1]
            assert isinstance(momentum, ScalarField)
            return momentum

        # Regular field lookup
        target_idx = self._component_name_to_index.get(field_name)
        if target_idx is not None:
            # Field is at even indices: state[2*i]
            field = state[2 * target_idx]
            assert isinstance(field, ScalarField)
            return field

        msg = f"Unknown field name: {field_name}"
        raise ValueError(msg)

    def _compute_rhs_for_component(
        self,
        component_idx: int,
        state: FieldCollection,
        bc: BCDescriptor,
    ) -> ScalarField:
        """Compute the RHS for a single component's momentum equation.

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

        Returns
        -------
        ScalarField
            The computed RHS for d/dt momentum_i.
        """
        eq = self.spec.equations[component_idx]
        grid = state.grid

        # Start with zero field
        result = ScalarField(grid, data=0.0)

        # Sum all terms from the specification
        for term in eq.rhs_terms:
            # Find which field this term operates on
            # Supports both regular fields (A_0, phi_0) and momentum fields (pi_0, pi_1)
            # Momentum field support enables mixed time-space derivative handling
            target_field_name = term.field
            target_field = self._get_field_from_state(state, target_field_name)

            # Apply the operator
            operated = self._get_operator(term.operator, target_field, bc)

            # Add coefficient * operated to result
            contribution = term.coefficient * operated
            result += contribution

        return result

    @override
    def evolution_rate(
        self,
        state: TState,
        t: float = 0.0,
    ) -> FieldCollection:
        """Compute the time derivatives for all fields.

        For each component i:
            d/dt field_i = momentum_i
            d/dt momentum_i = RHS from specification

        Parameters
        ----------
        state : FieldCollection
            Current state with 2 * n_components fields.
        t : float
            Current time (unused for autonomous systems).

        Returns
        -------
        FieldCollection
            Time derivatives for all fields.

        Raises
        ------
        ValueError
            If the state does not have exactly 2 * n_components fields.
        """
        assert isinstance(state, FieldCollection)
        expected_fields = 2 * self.n_components
        if len(state) != expected_fields:
            msg = f"Expected {expected_fields} fields, got {len(state)}"
            raise ValueError(msg)

        bc = infer_bc_from_grid(state.grid)
        rates: list[ScalarField] = []

        for i in range(self.n_components):
            field_i = state[2 * i]
            momentum_i = state[2 * i + 1]

            assert isinstance(field_i, ScalarField)
            assert isinstance(momentum_i, ScalarField)

            # d/dt field_i = momentum_i
            d_field_dt = momentum_i.copy()

            # d/dt momentum_i = RHS from specification
            d_momentum_dt = self._compute_rhs_for_component(i, state, bc)

            rates.extend((d_field_dt, d_momentum_dt))

        return FieldCollection(rates)

    def _cache_key(self) -> dict[str, Any]:
        """Return a cache key for this PDE.

        The key includes the specification metadata to ensure different
        equation systems don't share cached operators.
        """
        return {
            "n_components": self.n_components,
            "component_names": self.spec.component_names,
            "metadata_hash": hash(frozenset(self.spec.metadata.items())),
        }


def build_pde_from_json(json_path: Path | str) -> PDEFromSpec:
    """Build a PDE from a JSON equation specification file.

    This is the main entry point for the Lagrangian-to-PDE pipeline on the
    Python side. Given a JSON file exported from Mathematica/xAct, this
    function creates a py-pde compatible PDE class.

    Parameters
    ----------
    json_path : Path | str
        Path to the JSON file containing the equation specification.

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
    """
    spec = load_equation_system(json_path)
    return PDEFromSpec(spec)


def create_initial_state(
    grid: GridBase,
    spec: EquationSystem,
    field_data: dict[str, NDArray[np.float64]] | None = None,
    momentum_data: dict[str, NDArray[np.float64]] | None = None,
) -> FieldCollection:
    """Create initial state for a PDEFromSpec simulation.

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
        Components not specified default to zero.

    Returns
    -------
    FieldCollection
        Initial state with 2 * n_components fields.
    """
    field_data = field_data or {}
    momentum_data = momentum_data or {}

    fields: list[ScalarField] = []

    for name in spec.component_names:
        # Field
        if name in field_data:
            field = ScalarField(grid, data=field_data[name])
        else:
            field = ScalarField(grid, data=0.0)

        # Momentum
        if name in momentum_data:
            momentum = ScalarField(grid, data=momentum_data[name])
        else:
            momentum = ScalarField(grid, data=0.0)

        fields.extend((field, momentum))

    return FieldCollection(fields)
