from __future__ import annotations

from typing import cast

import numpy as np
from pde import CartesianGrid, FieldCollection, ScalarField

from torsion_gertsenshtein.kgsim.utils import natural_center


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

    Raises
    ------
    ValueError
        If `center` is provided but its dimensionality does not match the grid.

    Notes
    -----
    - Distances are computed using Euclidean norm on the grid cell coordinates.
    """
    coordinates = cast("np.ndarray", grid.cell_coords)  # (N, dim)
    if center is None:
        center = natural_center(grid.axes_bounds)
    # validate center dimensionality and raise explicitly to match the docstring
    if len(center) != grid.dim:
        msg = f"center must have length {grid.dim}, got {len(center)}"
        raise ValueError(msg)
    r2 = np.sum((coordinates - np.array(center)) ** 2, axis=1)
    phi_arr = amplitude * np.exp(-r2 / (2.0 * width**2))
    pi_arr = initial_velocity * phi_arr

    phi = ScalarField(grid, data=phi_arr.reshape(grid.shape))
    pi = ScalarField(grid, data=pi_arr.reshape(grid.shape))
    return FieldCollection([phi, pi], labels=["phi", "pi"])


def ring_pulse_2d(
    grid: CartesianGrid,
    *,
    amplitude: float = 1.0,
    initial_radius: float = 5.0,
    sigma: float = 1.0,
) -> FieldCollection:
    """
    Create a 2D ring-shaped Gaussian initial condition on a CartesianGrid.

    Parameters
    ----------
    grid : CartesianGrid
        A 2D CartesianGrid on which to evaluate the ring pulse. The function
        requires grid.dim == 2. The grid is expected to provide:
        - grid.cell_coords: an (N_cells, 2) array of cell center coordinates,
        - grid.axes_bounds: sequence of two (min, max) bounds for each axis,
        - grid.shape: the shape used to reshape 1D data arrays into the grid layout.
    amplitude : float, optional
        Peak amplitude of the Gaussian ring (default: 1.0).
    initial_radius : float, optional
        Radial distance from the grid center where the ring is centered
        (default: 5.0).
    sigma : float, optional
        Standard deviation of the radial Gaussian profile (controls ring width)
        (default: 1.0).

    Returns
    -------
    FieldCollection
        A collection containing two ScalarField objects with labels ["phi", "pi"]:
        - phi: the ring-shaped Gaussian field with data shaped to grid.shape,
          computed as amplitude * exp(-((radius - initial_radius)**2) / (2*sigma**2)),
          where radius is the distance from the grid center for each cell.
        - pi: a ScalarField of zeros with the same shape as phi.

    Raises
    ------
    ValueError
        If the provided grid is not two-dimensional (grid.dim != 2).

    Notes
    -----
    - The grid center is computed as the midpoint of each axis bounds:
      center_i = (axes_bounds[i][0] + axes_bounds[i][1]) / 2.
    - radius is computed per cell via Euclidean distance from the center using
      grid.cell_coords.
    - The returned phi data is created from a flattened radial profile and then
      reshaped to grid.shape to match the grid layout.
    """
    if grid.dim != 2:  # noqa: PLR2004
        msg = "ring_pulse_2d requires a 2D grid."
        raise ValueError(msg)
    coordinates = cast("np.ndarray", grid.cell_coords)
    # center is mid of bounds
    center_x = sum(grid.axes_bounds[0]) / 2
    center_y = sum(grid.axes_bounds[1]) / 2
    radius = np.sqrt(
        (coordinates[:, 0] - center_x) ** 2 + (coordinates[:, 1] - center_y) ** 2
    )
    phi_array = amplitude * np.exp(-((radius - initial_radius) ** 2) / (2 * sigma**2))
    phi = ScalarField(grid, data=phi_array.reshape(grid.shape))
    pi = ScalarField(grid, data=np.zeros_like(phi.data))
    return FieldCollection([phi, pi], labels=["phi", "pi"])


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
