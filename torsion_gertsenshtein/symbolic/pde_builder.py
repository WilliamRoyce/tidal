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

    def _get_operator(
        self, operator_name: str, field: ScalarField, bc: Any
    ) -> ScalarField:
        """Apply a named operator to a field.

        Parameters
        ----------
        operator_name : str
            Name of the operator ("laplacian", "identity", "gradient_x", etc.)
        field : ScalarField
            The field to operate on.
        bc : Any
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
            result = field.laplace(bc=bc)
            # Ensure we return a ScalarField
            if not isinstance(result, ScalarField):
                result = ScalarField(field.grid, data=result.data)
            return result
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

        msg = f"Unknown operator: {operator_name}"
        raise ValueError(msg)

    def _compute_rhs_for_component(
        self,
        component_idx: int,
        state: FieldCollection,
        bc: Any,
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
        bc : Any
            Boundary condition specification.

        Returns
        -------
        ScalarField
            The computed RHS for d/dt momentum_i.

        Raises
        ------
        ValueError
            If a term references an unknown field not in the component names.
        """
        eq = self.spec.equations[component_idx]
        grid = state.grid

        # Start with zero field
        result = ScalarField(grid, data=0.0)

        # Sum all terms from the specification
        for term in eq.rhs_terms:
            # Find which field this term operates on
            target_field_name = term.field
            target_idx = self._component_name_to_index.get(target_field_name)

            if target_idx is None:
                msg = f"Unknown field in equation: {target_field_name}"
                raise ValueError(msg)

            # Get the field (not momentum) for this component
            target_field = state[2 * target_idx]
            assert isinstance(target_field, ScalarField)

            # Apply the operator
            operated = self._get_operator(term.operator, target_field, bc)

            # Add coefficient * operated to result
            contribution = term.coefficient * operated
            result += contribution

        return result

    @override
    def evolution_rate(
        self,
        state: FieldCollection,
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
