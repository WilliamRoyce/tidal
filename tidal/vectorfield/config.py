"""Configuration classes for multi-component field simulations.

This module provides configuration dataclasses that work with equation
specifications loaded from the symbolic pipeline.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from tidal.symbolic.json_loader import EquationSystem


@dataclass(frozen=True)
class ComponentFieldParams:
    """Parameters for a multi-component field system.

    This class encapsulates the field configuration derived from an
    equation system specification.

    Attributes
    ----------
    n_components : int
        Number of field components.
    component_names : tuple[str, ...]
        Names of the field components (e.g., ("A_0", "A_1")).
    spatial_dimension : int
        Number of spatial dimensions.
    is_massless : bool
        Whether all components have zero mass.

    Examples
    --------
    >>> from tidal.symbolic import load_equation_system
    >>> spec = load_equation_system("examples/data/em_1d.json")
    >>> params = ComponentFieldParams.from_equation_system(spec)
    >>> params.n_components
    2
    """

    n_components: int
    component_names: tuple[str, ...]
    spatial_dimension: int
    is_massless: bool = True

    def __post_init__(self) -> None:
        """Validate parameters.

        Raises
        ------
        ValueError
            If n_components < 1, component_names length doesn't match
            n_components, or spatial_dimension < 1.
        """
        if self.n_components < 1:
            msg = "n_components must be at least 1"
            raise ValueError(msg)

        if len(self.component_names) != self.n_components:
            msg = (
                f"component_names length {len(self.component_names)} "
                f"!= n_components {self.n_components}"
            )
            raise ValueError(msg)

        if self.spatial_dimension < 1:
            msg = "spatial_dimension must be at least 1"
            raise ValueError(msg)

    @classmethod
    def from_equation_system(cls, spec: EquationSystem) -> ComponentFieldParams:
        """Create parameters from an equation system specification.

        Parameters
        ----------
        spec : EquationSystem
            The equation specification loaded from JSON.

        Returns
        -------
        ComponentFieldParams
            Configuration extracted from the specification.
        """
        # Check if massless by examining mass matrix diagonal.
        # Also check symbolic matrix: a symbolic mass entry (e.g., for
        # coordinate-dependent mass) means the field is NOT massless even
        # if the numeric placeholder is ambiguous.
        is_massless = all(
            spec.mass_matrix[i][i] == 0.0
            and (
                not spec.mass_matrix_symbolic
                or spec.mass_matrix_symbolic[i][i] is None
            )
            for i in range(spec.n_components)
        )

        return cls(
            n_components=spec.n_components,
            component_names=spec.component_names,
            spatial_dimension=spec.spatial_dimension,
            is_massless=is_massless,
        )

    def get_component_index(self, name: str) -> int:
        """Get the index of a component by name.

        Parameters
        ----------
        name : str
            Component name (e.g., "A_0" or "phi").

        Returns
        -------
        int
            Index of the component.

        Raises
        ------
        ValueError
            If the component name is not found.
        """
        try:
            return self.component_names.index(name)
        except ValueError:
            msg = f"Unknown component: {name}. Valid names: {self.component_names}"
            raise ValueError(msg) from None

    def get_state_field_indices(self) -> dict[str, tuple[int, int]]:
        """Get state array indices for each component.

        The state is organized as [field_0, momentum_0, field_1, momentum_1, ...].

        Returns
        -------
        dict[str, tuple[int, int]]
            Mapping from component name to (field_index, momentum_index).
        """
        return {name: (2 * i, 2 * i + 1) for i, name in enumerate(self.component_names)}
