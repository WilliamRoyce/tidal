"""Observer callbacks for tracking simulation quantities.

This module provides callback factories for use with py-pde's tracker system,
including energy conservation monitoring.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, cast

import numpy as np
from pde import ScalarField

from tidal.kgsim.utils import infer_bc_from_grid

if TYPE_CHECKING:
    from collections.abc import Callable

    from pde import FieldCollection


def total_energy_observer(
    mass: float,
) -> Callable[[FieldCollection, float], dict[str, Any]]:
    """Create an observer callback that computes total field energy.

    The returned callback computes the spatial integral of the canonical energy
    density for a Klein-Gordon-like field:

        E = ∫ dV [ 0.5 π² + 0.5 |∇φ|² + 0.5 m² φ² ]

    where φ is the scalar field, π its conjugate momentum, and m is the mass.

    Parameters
    ----------
    mass : float
        The scalar field mass m. The callback uses m² in the energy density.

    Returns
    -------
    Callable[[FieldCollection, float], dict[str, Any]]
        A function f(state, t) -> {"time": t, "total_energy": E} that:

        - Expects `state` to be a FieldCollection where state[0] is φ and
          state[1] is π.
        - Validates that both φ and π are ScalarField instances.
        - Computes the gradient of φ, forms |∇φ|², evaluates the pointwise
          energy density, multiplies by grid cell volumes, and returns the
          summed total energy.
        - Raises TypeError if state does not contain ScalarField instances.

    Notes
    -----
    The callback relies on `phi.grid.cell_volumes` for the integration measure.
    Numerical precision and conservation depend on the discretization and
    boundary conditions used when evolving the fields.
    """
    m2 = mass**2

    def _callback(state: FieldCollection, t: float) -> dict[str, Any]:
        phi = state[0]
        pi = state[1]
        if not isinstance(phi, ScalarField) or not isinstance(pi, ScalarField):
            msg = "evolution_rate requires ScalarField instances for (phi, pi)"
            raise TypeError(msg)

        bc = infer_bc_from_grid(phi.grid)
        if bc is None:
            grad_phi = phi.gradient(bc="periodic")
        else:
            # Cast to Any to satisfy the type checker; runtime accepts the mapping form.
            grad_phi = phi.gradient(bc=cast("Any", bc))
        grad_phi = cast("FieldCollection", grad_phi)

        grad_sq = sum(g.data**2 for g in grad_phi)  # |∇φ|^2
        density = 0.5 * (pi.data**2 + grad_sq + m2 * phi.data**2)

        # cell_volume accounts for grid measure
        volume_element = cast("np.ndarray", phi.grid.cell_volumes)

        total_energy = float(np.sum(density * volume_element.reshape(phi.data.shape)))
        return {"time": t, "total_energy": total_energy}

    return _callback
