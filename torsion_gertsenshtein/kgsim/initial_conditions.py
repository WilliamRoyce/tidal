"""Initial condition classes and factory functions for Klein-Gordon fields.

This module provides both class-based (InitialCondition, GaussianPulse, RingPulse2D)
and function-based (gaussian_pulse, ring_pulse_2d, plane_wave) APIs for creating
initial field configurations.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, cast

import numpy as np
from pde import CartesianGrid, FieldCollection, ScalarField
from typing_extensions import override

from torsion_gertsenshtein.kgsim.utils import extract_grid_coordinates, natural_center

if TYPE_CHECKING:
    from collections.abc import Sequence


class InitialCondition(ABC):
    """Abstract base class for Klein-Gordon initial conditions.

    This class provides common functionality for creating initial conditions
    on CartesianGrid, eliminating coordinate handling duplication across
    different IC types.

    Subclasses should implement the `_compute_phi()` and optionally
    `_compute_pi()` methods to define the specific initial field profiles.

    Serialization
    -------------
    Instances can be pickled/deepcopied without issues.
    """

    def _get_coordinates(self, grid: CartesianGrid) -> np.ndarray:  # noqa: PLR6301
        """Extract and normalize grid coordinates.

        Parameters
        ----------
        grid : CartesianGrid
            The computational grid

        Returns
        -------
        np.ndarray
            Coordinates array in shape (N_cells, dim) for uniform handling
        """
        return extract_grid_coordinates(grid, flatten=True)

    def _compute_distances_from_center(
        self, grid: CartesianGrid, center: Sequence[float] | None = None
    ) -> np.ndarray:
        """Compute Euclidean distances from a center point.

        Parameters
        ----------
        grid : CartesianGrid
            The computational grid
        center : Sequence[float] | None
            Center coordinates. If None, uses natural_center(grid.axes_bounds)

        Returns
        -------
        np.ndarray
            Flattened array of distances, shape (N_cells,)

        Raises
        ------
        ValueError
            If center dimension does not match grid dimension
        """
        if center is None:
            center = natural_center(grid.axes_bounds)

        if len(center) != grid.dim:
            msg = f"center must have length {grid.dim}, got {len(center)}"
            raise ValueError(msg)

        coords = self._get_coordinates(grid)
        center_arr = np.asarray(center)

        # Compute Euclidean distance
        return np.sqrt(np.sum((coords - center_arr) ** 2, axis=1))

    @abstractmethod
    def _compute_phi(self, grid: CartesianGrid) -> np.ndarray:
        """Compute the phi field data.

        Parameters
        ----------
        grid : CartesianGrid
            The computational grid

        Returns
        -------
        np.ndarray
            Flattened phi field values

        Raises
        ------
        ValueError
            If grid dimensions are invalid
        """

    def _compute_pi(self, grid: CartesianGrid, phi_data: np.ndarray) -> np.ndarray:  # noqa: ARG002, PLR6301
        """Compute the pi (time derivative) field data.

        Default implementation returns zeros. Override for non-zero initial velocity.

        Parameters
        ----------
        grid : CartesianGrid
            The computational grid
        phi_data : np.ndarray
            Flattened phi field values

        Returns
        -------
        np.ndarray
            Flattened pi field values
        """
        return np.zeros_like(phi_data)

    def build(self, grid: CartesianGrid) -> FieldCollection:
        """Build the initial condition FieldCollection.

        Parameters
        ----------
        grid : CartesianGrid
            The computational grid

        Returns
        -------
        FieldCollection
            Initial state containing [phi, pi] ScalarFields
        """
        phi_data = self._compute_phi(grid)
        pi_data = self._compute_pi(grid, phi_data)

        # Reshape to grid shape
        phi = ScalarField(grid, data=phi_data.reshape(grid.shape), label="phi")
        pi = ScalarField(grid, data=pi_data.reshape(grid.shape), label="pi")

        return FieldCollection([phi, pi], labels=["phi", "pi"])


class GaussianPulse(InitialCondition):
    """Gaussian pulse initial condition using the base class pattern.

    Creates a Gaussian pulse:
        phi(x) = amplitude * exp(-|x - center|^2 / (2 * width^2))
        pi(x) = initial_velocity * phi(x)

    Examples
    --------
    >>> from torsion_gertsenshtein.kgsim import make_grid, GridConfig
    >>> grid = make_grid(GridConfig(dim=1, shape=(64,), bounds=((0, 50),)))
    >>> ic = GaussianPulse(amplitude=1.0, width=2.0, center=[25.0])
    >>> state = ic.build(grid)
    """

    def __init__(
        self,
        *,
        amplitude: float = 1.0,
        width: float | None = None,
        center: list[float] | tuple[float, ...] | None = None,
        initial_velocity: float = 0.0,
        sigma: float | None = None,
    ) -> None:
        """Initialize Gaussian pulse parameters.

        Parameters
        ----------
        amplitude : float
            Peak amplitude (default: 1.0)
        width : float | None
            Gaussian width (standard deviation), must be positive (default: 2.0)
        center : list[float] | tuple[float, ...] | None
            Center coordinates. If None, uses grid center (default: None)
        initial_velocity : float
            Initial time derivative scaling factor (default: 0.0)
        sigma : float | None
            Alias for width (for consistency with RingPulse2D). Cannot be
            specified together with width.

        Raises
        ------
        ValueError
            If both width and sigma are specified, or if width <= 0, or if values are not finite

        Notes
        -----
        Both `width` and `sigma` refer to the Gaussian standard deviation.
        The `sigma` parameter is provided for consistency with RingPulse2D naming.
        If neither is specified, defaults to width=2.0.
        """
        # Reject if both parameters are specified
        if width is not None and sigma is not None:
            msg = "Cannot specify both 'width' and 'sigma' parameters. Use one or the other."
            raise ValueError(msg)

        # Use sigma if provided, otherwise use width, otherwise default
        if sigma is not None:
            final_width = sigma
        elif width is not None:
            final_width = width
        else:
            final_width = 2.0  # Default value

        # Validate finiteness BEFORE range checks
        if not np.isfinite(amplitude):
            msg = f"amplitude must be finite, got {amplitude}"
            raise ValueError(msg)

        if not np.isfinite(final_width):
            msg = f"width must be finite, got {final_width}"
            raise ValueError(msg)

        if not np.isfinite(initial_velocity):
            msg = f"initial_velocity must be finite, got {initial_velocity}"
            raise ValueError(msg)

        # Now check ranges
        if final_width <= 0:
            msg = f"width must be positive, got {final_width}"
            raise ValueError(msg)

        # Validate and copy center to prevent mutation
        center_copy = list(center) if center is not None else None

        self.amplitude = amplitude
        self.width = final_width
        self.center = center_copy
        self.initial_velocity = initial_velocity

    @override
    def _compute_phi(self, grid: CartesianGrid) -> np.ndarray:
        """Compute Gaussian phi field."""
        r = self._compute_distances_from_center(grid, self.center)
        return self.amplitude * np.exp(-(r**2) / (2.0 * self.width**2))

    @override
    def _compute_pi(self, grid: CartesianGrid, phi_data: np.ndarray) -> np.ndarray:
        """Compute pi field with initial velocity."""
        return self.initial_velocity * phi_data


class RingPulse2D(InitialCondition):
    """2D ring-shaped Gaussian pulse initial condition.

    Creates a radially symmetric ring:
        phi(r) = amplitude * exp(-(r - initial_radius)^2 / (2 * sigma^2))

    where r is the radial distance from the grid center (midpoint of bounds).

    Examples
    --------
    >>> grid = make_grid(GridConfig(dim=2, shape=(64, 64),
    ...                             bounds=((-20, 20), (-20, 20))))
    >>> ic = RingPulse2D(amplitude=1.0, initial_radius=5.0, sigma=1.0)
    >>> state = ic.build(grid)
    """

    def __init__(
        self,
        *,
        amplitude: float = 1.0,
        initial_radius: float = 5.0,
        sigma: float | None = None,
        width: float | None = None,
    ) -> None:
        """Initialize ring pulse parameters.

        Parameters
        ----------
        amplitude : float
            Peak amplitude (default: 1.0)
        initial_radius : float
            Radial position of the ring (default: 5.0)
        sigma : float | None
            Gaussian width of the ring (default: 1.0)
        width : float | None
            Alias for sigma (for consistency with GaussianPulse). Cannot be
            specified together with sigma.

        Raises
        ------
        ValueError
            If both sigma and width are specified, or if sigma <= 0,
            or if initial_radius < 0

        Notes
        -----
        Both `sigma` and `width` refer to the Gaussian standard deviation.
        The `width` parameter is provided for consistency with GaussianPulse naming.
        If neither is specified, defaults to sigma=1.0.
        """
        # Reject if both parameters are specified
        if width is not None and sigma is not None:
            msg = "Cannot specify both 'sigma' and 'width' parameters. Use one or the other."
            raise ValueError(msg)

        # Use width if provided, otherwise use sigma, otherwise default
        if width is not None:
            final_sigma = width
        elif sigma is not None:
            final_sigma = sigma
        else:
            final_sigma = 1.0  # Default value

        # Validate finiteness BEFORE range checks
        if not np.isfinite(amplitude):
            msg = f"amplitude must be finite, got {amplitude}"
            raise ValueError(msg)

        if not np.isfinite(final_sigma):
            msg = f"sigma must be finite, got {final_sigma}"
            raise ValueError(msg)

        if not np.isfinite(initial_radius):
            msg = f"initial_radius must be finite, got {initial_radius}"
            raise ValueError(msg)

        # Now check ranges
        if final_sigma <= 0:
            msg = f"sigma must be positive, got {final_sigma}"
            raise ValueError(msg)
        if initial_radius < 0:
            msg = f"initial_radius must be non-negative, got {initial_radius}"
            raise ValueError(msg)

        self.amplitude = amplitude
        self.initial_radius = initial_radius
        self.width = final_sigma  # Store as .width for consistency with GaussianPulse

    @override
    def _compute_phi(self, grid: CartesianGrid) -> np.ndarray:
        """Compute ring-shaped phi field centered at grid midpoint.

        Raises
        ------
        ValueError
            If grid is not 2D
        """
        if grid.dim != 2:  # noqa: PLR2004
            msg = f"RingPulse2D requires 2D grid, got {grid.dim}D"
            raise ValueError(msg)

        # Use grid center instead of origin for consistency with original behavior
        center = natural_center(grid.axes_bounds)
        r = self._compute_distances_from_center(grid, center)
        return self.amplitude * np.exp(
            -((r - self.initial_radius) ** 2) / (2.0 * self.width**2)
        )


# Keep original function-based API for backward compatibility


def gaussian_pulse(
    grid: CartesianGrid,
    *,
    amplitude: float = 1.0,
    center: list[float] | None = None,
    width: float = 2.0,
    initial_velocity: float = 0.0,
) -> FieldCollection:
    """
    Create a Gaussian pulse on a CartesianGrid and its conjugate momentum.

    The scalar field phi is defined as
        phi(x) = amplitude * exp(-|x - center|^2 / (2 * width^2))
    and the conjugate momentum pi is taken to be
        pi(x) = initial_velocity * phi(x).

    Parameters
    ----------
    grid : CartesianGrid
        The computational grid. The function uses the grid's cell coordinates and
        shape to evaluate the Gaussian and to construct ScalarField objects that
        conform to the grid.
    amplitude : float, optional
        Peak amplitude A of the Gaussian (default: 1.0).
    center : list[float] or None, optional
        Coordinates of the Gaussian center x0. If None (default), the center is
        chosen via natural_center(grid.axes_bounds). The length of the list must
        match the dimensionality of the grid.
    width : float, optional
        Gaussian width (standard deviation sigma). Must be positive (default: 2.0).
    initial_velocity : float, optional
        Scalar factor used to set the initial momentum pi = initial_velocity * phi
        (default: 0.0).

    Returns
    -------
    FieldCollection
        A FieldCollection containing two ScalarField objects, labeled "phi" and
        "pi", whose data arrays have the same shape as the provided grid.

    Notes
    -----
    - Distances are computed using Euclidean norm on the grid cell coordinates.
    - This function is a convenience wrapper around the GaussianPulse class.
      For more flexibility, consider using GaussianPulse directly.

    See Also
    --------
    GaussianPulse : Class-based API for more extensibility and customization.
    """
    ic = GaussianPulse(
        amplitude=amplitude,
        width=width,
        center=center,
        initial_velocity=initial_velocity,
    )
    return ic.build(grid)


def ring_pulse_2d(
    grid: CartesianGrid,
    *,
    amplitude: float = 1.0,
    initial_radius: float = 5.0,
    sigma: float = 1.0,
) -> FieldCollection:
    """
    Create a 2D ring-shaped Gaussian pulse on the provided CartesianGrid.

    The function constructs a scalar field phi defined by a circular Gaussian
    envelope centered at the origin and an accompanying momentum field pi
    initialized to zero. The radial profile is

        phi(r) = amplitude * exp(- (r - initial_radius)**2 / (2 * sigma**2) )

    where r is the distance from the origin. The fields are returned as a
    FieldCollection containing two ScalarField instances labeled "phi" and "pi".

    Parameters
    ----------
    grid : CartesianGrid
        A 2-dimensional grid object. If grid.dim != 2 a ValueError is raised.
    amplitude : float, optional
        Peak amplitude of the Gaussian ring (default: 1.0).
    initial_radius : float, optional
        Radius of the ring (distance from the origin where the Gaussian is
        centered) (default: 5.0).
    sigma : float, optional
        Standard deviation (width) of the Gaussian envelope (default: 1.0).

    Returns
    -------
    FieldCollection
        A FieldCollection([phi, pi], labels=["phi", "pi"]) where:
          - phi is a ScalarField on `grid` containing the ring Gaussian values
            with shape == grid.shape, and
          - pi is a ScalarField of zeros with the same shape.

    Notes
    -----
    - The returned phi array has the same layout/shape as grid.shape.
    - Ring is centered at the grid midpoint, not the origin.
    - This function is a convenience wrapper around the RingPulse2D class.
      For more flexibility, consider using RingPulse2D directly.

    See Also
    --------
    RingPulse2D : Class-based API for more extensibility and customization.
    """
    ic = RingPulse2D(
        amplitude=amplitude,
        initial_radius=initial_radius,
        sigma=sigma,
    )
    return ic.build(grid)


def plane_wave(
    grid: CartesianGrid,
    *,
    amplitude: float = 1.0,
    k_vec: list[float] | None = None,
    mass: float = 1.0,
    phase: float = 0.0,
) -> FieldCollection:
    """
    Create a Klein-Gordon plane-wave initial condition on a CartesianGrid.

    The returned FieldCollection contains two ScalarField objects:
    - "phi":    φ(x) = A * cos(k · x + phase)
    - "pi":     π(x) = -A * ω * sin(k · x + phase)

    The angular frequency ω satisfies the dispersion relation
        ω^2 = |k|^2 + m^2

    Parameters
    ----------
    grid : CartesianGrid
        Grid object providing geometry and cell coordinates. The function uses
        grid.cell_coords (flattened coordinates) for evaluating the plane wave
        and reshapes results to grid.shape when constructing the ScalarField
        instances.
    amplitude : float, optional
        Amplitude A of the plane wave. Default is 1.0.
    k_vec : list[float] or None, optional
        Wavevector k. If None, the function chooses the fundamental mode along
        the first axis with magnitude 2π / L_x (where L_x is the axis length),
        and zeros for the remaining components. If provided, its length should
        match grid.dim.
    mass : float, optional
        Mass parameter m appearing in the dispersion relation. Default is 1.0.
    phase : float, optional
        Global phase offset (in radians) added to k · x. Default is 0.0.

    Returns
    -------
    FieldCollection
        A FieldCollection containing two ScalarField objects, labeled
        ["phi", "pi"], each with data shaped to grid.shape corresponding to the
        cell-centered evaluation of φ and π over the grid.

    Notes
    -----
    - The function computes k · x using the flattened cell coordinates returned
      by grid.cell_coords; ensure those coordinates are in the desired units.
    - If k_vec has a different length than grid.dim, behavior is undefined and
      a mismatch will typically raise an exception when performing the dot
      product or when reshaping the arrays.
    - Omega is computed from k_vec and mass via the dispersion relation.
    """
    if k_vec is None:
        # fundamental mode along x
        l_x = grid.axes_bounds[0][1] - grid.axes_bounds[0][0]
        k_vec = [2 * np.pi / l_x] + [0] * (grid.dim - 1)
    k = np.array(k_vec, dtype=float)

    omega = np.sqrt(np.dot(k, k) + mass**2)

    x_coordinate = cast("np.ndarray", grid.cell_coords)
    phase_array = x_coordinate @ k + phase
    phi_array = amplitude * np.cos(phase_array)
    pi_array = -amplitude * omega * np.sin(phase_array)

    phi = ScalarField(grid, data=phi_array.reshape(grid.shape))
    pi = ScalarField(grid, data=pi_array.reshape(grid.shape))
    return FieldCollection([phi, pi], labels=["phi", "pi"])


def _expand_param(
    name: str, values: Sequence[float] | float | None, n: int
) -> np.ndarray:
    """Normalize scalar-or-sequence parameter to a 1D numpy array of length n.

    Parameters
    ----------
    name : str
        Name of the parameter (used for error messages).
    values : Sequence[float] | float | None
        A scalar or sequence to normalize to a 1D numpy array.
        If a scalar is provided, it will be broadcast to all fields.
        If a sequence is provided, it must have exactly `n` elements.
    n : int
        Target length of the returned array.

    Returns
    -------
    np.ndarray
        A 1D numpy array of length `n` containing the parameter values.

    Raises
    ------
    ValueError
        If `values` is None, or if a sequence is provided with length != `n`.
    """
    if values is None:
        msg = f"{name} must be provided"
        raise ValueError(msg)

    # Check if input is a scalar (not a sequence)
    # np.ndim returns 0 for scalars, >0 for arrays/sequences
    is_scalar = np.ndim(values) == 0

    arr = np.atleast_1d(np.asarray(values, dtype=float))

    # Only broadcast if the original input was a scalar
    if is_scalar and n > 1:
        arr = np.full((n,), float(arr.item()), dtype=float)

    if arr.size != n:
        msg = f"{name} must have length {n} (got {arr.size}), or be a scalar for broadcasting"
        raise ValueError(msg)

    return arr


def multi_gaussian(
    grid: CartesianGrid,
    amplitudes: Sequence[float],
    widths: Sequence[float],
    centers: Sequence[float] | None = None,
    velocities: Sequence[float] | None = None,
) -> FieldCollection:
    """1D-only multi-gaussian initializer (safe handling of shapes and broadcasting).

    Parameters
    ----------
    grid : CartesianGrid
        A 1D CartesianGrid.
    amplitudes : Sequence[float]
        Sequence of amplitudes, one per Gaussian.
    widths : Sequence[float]
        Sequence of widths (standard deviations), one per Gaussian.
    centers : Sequence[float] or None, optional
        Sequence of center positions (or None to use domain midpoint for each).
    velocities : Sequence[float] or None, optional
        Sequence of initial velocities (or None to use zeros).

    Returns
    -------
    FieldCollection
        FieldCollection containing (phi, pi) pairs for each Gaussian.

    Raises
    ------
    ValueError
        If grid.dim != 1, if any width is non-positive, or if any of the
        parameter sequences does not match the expected number of fields.
    """
    if grid.dim != 1:
        msg = "multi_gaussian currently supports 1D grids only"
        raise ValueError(msg)

    field_count = len(amplitudes)

    # validate non-empty amplitudes (prevent silent no-op with zero fields)
    if field_count == 0:
        msg = "amplitudes must be a non-empty sequence"
        raise ValueError(msg)

    amp_arr = _expand_param("amplitudes", amplitudes, field_count)
    width_arr = _expand_param("widths", widths, field_count)

    # coordinates: ensure concrete numpy array for typing and numeric ops
    coords = cast("np.ndarray", grid.cell_coords)
    x = coords[:, 0] if coords.ndim > 1 else coords

    # default center: midpoint of domain
    mid = float((x.min() + x.max()) / 2.0)
    centers_arr = (
        _expand_param("centers", centers, field_count)
        if centers is not None
        else np.full((field_count,), mid, dtype=float)
    )

    vel_arr = (
        _expand_param("velocities", velocities, field_count)
        if velocities is not None
        else np.zeros((field_count,), dtype=float)
    )

    fields: list[ScalarField] = []
    labels: list[str | None] = []
    # Build fields using vectorized expression per field
    for i in range(field_count):
        sigma = width_arr[i]
        if sigma <= 0:
            msg = "widths must be positive"
            raise ValueError(msg)
        dx = x - centers_arr[i]
        phi_arr = amp_arr[i] * np.exp(-(dx * dx) / (2.0 * sigma * sigma))
        pi_arr = vel_arr[i] * phi_arr

        fields.extend(
            (
                ScalarField(grid, data=phi_arr.reshape(grid.shape)),
                ScalarField(grid, data=pi_arr.reshape(grid.shape)),
            )
        )
        labels.extend([f"phi{i}", f"pi{i}"])

    return FieldCollection(fields, labels=labels)


def multi_gaussian_2d(  # noqa: PLR0914
    grid: CartesianGrid,
    amplitudes: Sequence[float],
    widths: Sequence[float],
    centers: Sequence[tuple[float, float]] | None = None,
    velocities: Sequence[float] | None = None,
) -> FieldCollection:
    """2D multi-gaussian initializer for N coupled Klein-Gordon fields.

    This function creates a FieldCollection with N field pairs (phi_i, pi_i),
    each initialized as a Gaussian centered at a specified 2D position. The
    Gaussian for field i has the form:

        phi_i(x, y) = A_i * exp(-r_i^2 / (2 * sigma_i^2))

    where r_i = sqrt((x - x_i)^2 + (y - y_i)^2) is the Euclidean distance from
    the center (x_i, y_i), A_i is the amplitude, and sigma_i is the width.

    The conjugate momentum is:
        pi_i(x, y) = v_i * phi_i(x, y)

    where v_i is the initial velocity for field i.

    Parameters
    ----------
    grid : CartesianGrid
        A 2D CartesianGrid. The function will raise ValueError if grid.dim != 2.
    amplitudes : Sequence[float]
        Sequence of amplitudes [A_0, A_1, ..., A_{N-1}], one per field.
        Must be non-empty.
    widths : Sequence[float]
        Sequence of widths (standard deviations) [sigma_0, sigma_1, ..., sigma_{N-1}],
        one per field. All widths must be positive.
    centers : Sequence[tuple[float, float]] or None, optional
        Sequence of center positions [(x_0, y_0), (x_1, y_1), ...], one per
        field. If None, all fields are centered at the midpoint of the domain.
    velocities : Sequence[float] or None, optional
        Sequence of initial velocities [v_0, v_1, ...], one per field.
        If None, all velocities are set to zero.

    Returns
    -------
    FieldCollection
        A FieldCollection containing 2N ScalarField objects with labels
        ["phi0", "pi0", "phi1", "pi1", ..., "phi{N-1}", "pi{N-1}"].
        Each phi_i and pi_i has shape == grid.shape.

    Raises
    ------
    ValueError
        - If grid.dim != 2.
        - If amplitudes is empty.
        - If any width is non-positive.
        - If centers is provided but has length != len(amplitudes).
        - If velocities is provided but has length != len(amplitudes).

    Notes
    -----
    - This function is designed for initializing coupled Klein-Gordon systems
      with spatially separated or overlapping Gaussian pulses.
    - For single-field 2D Gaussians, use `gaussian_pulse` instead.
    - For radial ring-shaped pulses in 2D, use `ring_pulse_2d`.

    Examples
    --------
    Create two Gaussians separated along the x-axis::

        grid = make_grid(GridConfig(dim=2, shape=(128, 128),
                                    bounds=((-50, 50), (-50, 50)),
                                    periodic=True))
        state = multi_gaussian_2d(
            grid,
            amplitudes=[1.0, 0.0],
            widths=[5.0, 5.0],
            centers=[(-10.0, 0.0), (10.0, 0.0)],
        )
    """
    if grid.dim != 2:  # noqa: PLR2004
        msg = "multi_gaussian_2d requires a 2D grid (grid.dim == 2)"
        raise ValueError(msg)

    field_count = len(amplitudes)

    # Validate non-empty amplitudes
    if field_count == 0:
        msg = "amplitudes must be a non-empty sequence"
        raise ValueError(msg)

    amp_arr = _expand_param("amplitudes", amplitudes, field_count)
    width_arr = _expand_param("widths", widths, field_count)

    # Get 2D coordinates - cell_coords may be either (nx, ny, 2) or (N_cells, 2) depending on py-pde version
    coords = cast("np.ndarray", grid.cell_coords)
    if coords.ndim == 3 and coords.shape[-1] == grid.dim:  # noqa: PLR2004
        # Coords are in shape (nx, ny, 2) - extract x and y arrays
        x = coords[:, :, 0].ravel()
        y = coords[:, :, 1].ravel()
    elif coords.ndim == 2 and coords.shape[1] == grid.dim:  # noqa: PLR2004
        # Coords are already flattened (N_cells, 2)
        x = coords[:, 0]
        y = coords[:, 1]
    else:
        msg = f"Unexpected cell_coords shape: {coords.shape}"
        raise ValueError(msg)

    # Default center: midpoint of domain in both dimensions
    mid_x = float((grid.axes_bounds[0][0] + grid.axes_bounds[0][1]) / 2.0)
    mid_y = float((grid.axes_bounds[1][0] + grid.axes_bounds[1][1]) / 2.0)

    # Handle centers parameter
    if centers is None:
        centers_arr = np.array(
            [(mid_x, mid_y) for _ in range(field_count)], dtype=float
        )
    else:
        if len(centers) != field_count:
            msg = f"centers must have length {field_count} (one per field)"
            raise ValueError(msg)
        centers_arr = np.array(centers, dtype=float)
        if centers_arr.shape != (field_count, 2):
            msg = f"centers must be a sequence of (x, y) tuples, got shape {centers_arr.shape}"
            raise ValueError(msg)

    # Handle velocities parameter
    vel_arr = (
        _expand_param("velocities", velocities, field_count)
        if velocities is not None
        else np.zeros((field_count,), dtype=float)
    )

    fields: list[ScalarField] = []
    labels: list[str | None] = []

    # Build fields using vectorized expressions per field
    for i in range(field_count):
        sigma = width_arr[i]
        if sigma <= 0:
            msg = "all widths must be positive"
            raise ValueError(msg)

        # Compute Euclidean distance from center
        cx, cy = centers_arr[i]
        dx = x - cx
        dy = y - cy
        r_squared = dx * dx + dy * dy

        # Gaussian profile
        phi_arr = amp_arr[i] * np.exp(-r_squared / (2.0 * sigma * sigma))
        pi_arr = vel_arr[i] * phi_arr

        # Reshape to grid shape and create ScalarFields
        fields.extend(
            (
                ScalarField(grid, data=phi_arr.reshape(grid.shape)),
                ScalarField(grid, data=pi_arr.reshape(grid.shape)),
            )
        )
        labels.extend([f"phi{i}", f"pi{i}"])

    return FieldCollection(fields, labels=labels)
