"""Initial conditions for multi-component field simulations.

This module provides classes and functions for creating initial conditions
for field systems derived from Lagrangians.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, cast

import numpy as np
from pde import FieldCollection, ScalarField

if TYPE_CHECKING:
    from numpy.typing import NDArray
    from pde.grids.base import GridBase

    from tidal.symbolic.json_loader import EquationSystem


@dataclass
class ComponentGaussianPulse:
    """Gaussian pulse initial condition for specified field components.

    Creates a Gaussian pulse in one or more field components:
        field(x) = amplitude * exp(-|x - center|^2 / (2 * width^2))

    Attributes
    ----------
    center : tuple[float, ...]
        Center position of the Gaussian in each spatial dimension.
    width : float
        Width (standard deviation) of the Gaussian.
    amplitude : float
        Peak amplitude of the Gaussian.
    active_components : dict[str, float] | None
        Mapping from component names to their relative amplitudes.
        If None, all components get the pulse with amplitude 1.0.

    Examples
    --------
    >>> pulse = ComponentGaussianPulse(
    ...     center=(50.0,),
    ...     width=5.0,
    ...     amplitude=1.0,
    ...     active_components={"A_1": 1.0}
    ... )
    >>> state = pulse.create(grid, spec)
    """

    center: tuple[float, ...]
    width: float
    amplitude: float = 1.0
    active_components: dict[str, float] | None = None

    def __post_init__(self) -> None:
        """Validate parameters.

        Raises
        ------
        ValueError
            If width is not positive.
        """
        if self.width <= 0:
            msg = "width must be positive"
            raise ValueError(msg)

    def compute_gaussian(self, grid: GridBase) -> NDArray[np.float64]:
        """Compute Gaussian profile on the grid.

        Parameters
        ----------
        grid : GridBase
            The simulation grid.

        Returns
        -------
        NDArray[np.float64]
            Gaussian profile values at each grid point.
        """
        # Get grid coordinates
        coords = cast("np.ndarray", grid.cell_coords)

        # Compute squared distance from center
        dist_sq = np.zeros(grid.shape, dtype=np.float64)
        for dim in range(grid.num_axes):
            if dim < len(self.center):
                dist_sq += (coords[..., dim] - self.center[dim]) ** 2
            # If center doesn't specify this dimension, treat as centered at 0

        # Compute Gaussian
        return self.amplitude * np.exp(-dist_sq / (2 * self.width**2))

    def create(self, grid: GridBase, spec: EquationSystem) -> FieldCollection:
        """Create initial state with Gaussian pulses in specified components.

        Parameters
        ----------
        grid : GridBase
            The simulation grid.
        spec : EquationSystem
            The equation specification (determines component count and names).

        Returns
        -------
        FieldCollection
            Initial state with 2 * n_components fields (field + momentum pairs).
        """
        gaussian = self.compute_gaussian(grid)

        # Determine which components get the pulse
        if self.active_components is None:
            # Default: all components get the pulse
            active = dict.fromkeys(spec.component_names, 1.0)
        else:
            active = self.active_components

        fields: list[ScalarField] = []

        for name in spec.component_names:
            # Field component
            if name in active:
                field_data = active[name] * gaussian
            else:
                field_data = np.zeros(grid.shape, dtype=np.float64)

            # Use extend to add both field and momentum components at once
            fields.extend(
                (
                    ScalarField(grid, data=field_data),
                    ScalarField(grid, data=0.0),  # Momentum component (always zero)
                )
            )

        return FieldCollection(fields)


@dataclass
class ComponentPlaneWave:
    """Plane wave initial condition for specified field components.

    Creates a plane wave:
        field(x) = amplitude * cos(k · x + phase)

    With corresponding momentum for a propagating wave:
        momentum(x) = -amplitude * omega * sin(k · x + phase)

    where omega = |k| for massless fields.

    Attributes
    ----------
    wavevector : tuple[float, ...]
        Wavevector k in each spatial dimension.
    amplitude : float
        Amplitude of the wave.
    phase : float
        Initial phase in radians.
    active_components : dict[str, float] | None
        Mapping from component names to their relative amplitudes.

    Examples
    --------
    >>> wave = ComponentPlaneWave(
    ...     wavevector=(0.1,),  # k_x = 0.1
    ...     amplitude=1.0,
    ...     phase=0.0,
    ...     active_components={"A_1": 1.0}
    ... )
    >>> state = wave.create(grid, spec)
    """

    wavevector: tuple[float, ...]
    amplitude: float = 1.0
    phase: float = 0.0
    active_components: dict[str, float] | None = None

    def compute_plane_wave(
        self, grid: GridBase
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Compute plane wave field and momentum on the grid.

        Returns
        -------
        tuple[NDArray, NDArray]
            (field_data, momentum_data) arrays.
        """
        coords = cast("np.ndarray", grid.cell_coords)

        # Compute k · x
        k_dot_x = np.zeros(grid.shape, dtype=np.float64)
        for dim in range(min(grid.num_axes, len(self.wavevector))):
            k_dot_x += self.wavevector[dim] * coords[..., dim]

        # Compute |k| (omega for massless)
        k_mag = np.sqrt(sum(k**2 for k in self.wavevector))

        # Field: amplitude · cos(k·x + phase)
        field_data = self.amplitude * np.cos(k_dot_x + self.phase)

        # Momentum: -amplitude · omega · sin(k·x + phase) for right-moving wave
        momentum_data = -self.amplitude * k_mag * np.sin(k_dot_x + self.phase)

        return field_data, momentum_data

    def create(self, grid: GridBase, spec: EquationSystem) -> FieldCollection:
        """Create initial state with plane waves in specified components.

        Parameters
        ----------
        grid : GridBase
            The simulation grid.
        spec : EquationSystem
            The equation specification.

        Returns
        -------
        FieldCollection
            Initial state with propagating plane waves.
        """
        field_template, momentum_template = self.compute_plane_wave(grid)

        # Determine which components get the wave
        if self.active_components is None:
            active = dict.fromkeys(spec.component_names, 1.0)
        else:
            active = self.active_components

        fields: list[ScalarField] = []

        for name in spec.component_names:
            if name in active:
                rel_amp = active[name]
                field_data = rel_amp * field_template
                momentum_data = rel_amp * momentum_template
            else:
                field_data = np.zeros(grid.shape, dtype=np.float64)
                momentum_data = np.zeros(grid.shape, dtype=np.float64)

            fields.extend(
                (
                    ScalarField(grid, data=field_data),
                    ScalarField(grid, data=momentum_data),
                )
            )

        return FieldCollection(fields)


def create_gaussian_pulse_state(  # noqa: PLR0913
    grid: GridBase,
    spec: EquationSystem,
    center: tuple[float, ...],
    width: float,
    amplitude: float = 1.0,
    *,
    component: str | None = None,
) -> FieldCollection:
    """Return a Gaussian pulse initial state.

    Parameters
    ----------
    grid : GridBase
        The simulation grid.
    spec : EquationSystem
        The equation specification.
    center : tuple[float, ...]
        Center position of the Gaussian.
    width : float
        Width of the Gaussian.
    amplitude : float
        Peak amplitude.
    component : str | None
        If specified, only this component gets the pulse.
        If None, all components get the pulse.

    Returns
    -------
    FieldCollection
        Initial state with Gaussian pulse.

    Examples
    --------
    >>> state = create_gaussian_pulse_state(
    ...     grid, spec,
    ...     center=(50.0,), width=5.0,
    ...     component="A_1"
    ... )
    """
    if component is not None:
        active_components = {component: amplitude}
        amplitude = 1.0  # Amplitude handled in active_components
    else:
        active_components = None

    pulse = ComponentGaussianPulse(
        center=center,
        width=width,
        amplitude=amplitude,
        active_components=active_components,
    )

    return pulse.create(grid, spec)
