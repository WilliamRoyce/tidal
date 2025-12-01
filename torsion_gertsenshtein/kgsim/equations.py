from __future__ import annotations

from typing import TYPE_CHECKING, Any, cast

import numpy as np
from pde import PDE, DataFieldBase, FieldCollection, PDEBase, ScalarField
from typing_extensions import override

from torsion_gertsenshtein.kgsim.utils import infer_bc_from_grid

if TYPE_CHECKING:
    from collections.abc import Sequence

    from torsion_gertsenshtein.kgsim.config import MultiFieldParams

    from .config import KGParameters
import logging

_logger = logging.getLogger(__name__)


class KleinGordonPDE(PDE):
    """
    Klein-Gordon equation represented in first-order (in time) form.

    This class implements the following system of equations:
        d/dt phi = pi
        d/dt pi  = laplace(phi) - m^2 * phi

    It is intended to be used with the PDE solver framework and operates on
    FieldCollection objects containing two ScalarField instances: (phi, pi).

    Parameters
    ----------
    params : KGParameters
        Parameter bundle for the Klein-Gordon equation. Only the `mass` attribute
        is used; its square is stored as the internal `m2` attribute.

    Attributes
    ----------
    m2 : float
        Square of the Klein-Gordon mass (mass**2), computed from `params.mass`.
    explicit_time_dependence : bool
        Set to False to indicate the PDE has no explicit dependence on time.

    """

    explicit_time_dependence = False

    def __init__(self, params: KGParameters) -> None:
        m2 = float(params.mass) ** 2
        self.m2: float = m2  # retained for energy observers / metadata

        # Expression engine: two fields named "phi" and "pi"
        # d_t phi = pi
        # d_t pi  = laplace(phi) - m2 * phi
        rhs = {
            "phi": "pi",
            "pi": "laplace(phi) - m2 * phi",
        }
        consts = {"m2": m2}
        super().__init__(rhs, consts=cast("Any", consts))

    def _cache_key(self) -> dict[str, Any]:
        return {"m2": self.m2}


class InhomogeneousKGPDE(PDEBase):
    """
    Inhomogeneous Klein-Gordon with spatial coefficients.

        d_t phi = pi
        d_t pi  = laplace(phi) - m2(x) * phi + V(x) * phi

    Notes
    -----
    - `m2_field` and `V_field` are ScalarField on the same grid as `state`.
    - backend="numba" is NOT supported for this custom RHS; use backend="numpy".
      (Your runner can auto-fallback to numpy.)
    """

    explicit_time_dependence = False

    def __init__(
        self, m2_field: ScalarField, potential_field: ScalarField | None = None
    ) -> None:
        super().__init__()
        if potential_field is None:
            potential_field = ScalarField(m2_field.grid, data=0.0)
        if m2_field.grid != potential_field.grid:
            msg = "m2_field and V_field must live on the same grid"
            raise ValueError(msg)
        # Normalize numeric dtype to float64 to avoid unexpected dtype/broadcast behavior
        m2_field = ScalarField(
            m2_field.grid, data=np.asarray(m2_field.data, dtype=float)
        )
        potential_field = ScalarField(
            potential_field.grid, data=np.asarray(potential_field.data, dtype=float)
        )
        self.m2_field = m2_field
        self.potential_field = potential_field

    @override
    def evolution_rate(
        self, state: FieldCollection | DataFieldBase, t: float = 0.0
    ) -> FieldCollection:
        # runtime sanity checks to prevent silent broadcasting which leads to blow-ups
        if not isinstance(state, FieldCollection):
            msg = "state must be FieldCollection(phi, pi)"
            raise TypeError(msg)
        phi = state[0]
        pi = state[1]
        if not isinstance(phi, ScalarField) or not isinstance(pi, ScalarField):
            msg = "state must contain ScalarField(phi) and ScalarField(pi)"
            raise TypeError(msg)

        # Ensure coefficient fields live on the same grid as phi
        if phi.grid != self.m2_field.grid or phi.grid != self.potential_field.grid:
            msg = "state fields and coefficient fields must live on the same grid"
            raise ValueError(msg)

        # For periodic grids: pass bc="periodic" to let py-pde handle it
        # For non-periodic: use inferred BC or periodic as fallback
        periodic = getattr(phi.grid, "periodic", None)
        # Narrow the type for the type checker before calling any()
        if periodic is True:
            is_periodic = True
        elif isinstance(periodic, (list, tuple)):
            periodic_seq = cast("Sequence[bool]", periodic)
            is_periodic = any(periodic_seq)
        else:
            is_periodic = False

        if is_periodic:
            # Periodic grid: use default boundary handling (NO bc argument)
            _logger.debug("Using periodic boundary handling (no explicit bc)")
            lap_phi = phi.laplace(bc="periodic")
        else:
            # Non-periodic: infer appropriate BC
            bc = infer_bc_from_grid(phi.grid)
            lap_phi = (
                phi.laplace(bc="periodic")
                if bc is None
                else phi.laplace(bc=cast("Any", bc))
            )

        dphi_dt = pi.copy()

        dpi_dt = lap_phi - self.m2_field * phi + self.potential_field * phi
        # Ensure runtime types are ScalarField so we can safely cast for the type checker.
        if not isinstance(dpi_dt, ScalarField):
            msg = "dpi_dt computed non-ScalarField results"
            raise TypeError(msg)

        return FieldCollection([dphi_dt, dpi_dt], labels=["phi", "pi"])

    def _cache_key(self) -> dict[str, Any]:
        # Hash by checksums to avoid serializing full arrays
        return {
            "m2_checksum": float(self.m2_field.data.sum()),
            "V_checksum": float(self.potential_field.data.sum()),
        }


def make_coupled_kg_pde(params: MultiFieldParams) -> PDE:
    """
    Build an expression-based PDE for N coupled Klein-Gordon fields.

    This class implements the following system of equations:
        d_t phi_i = pi_i
        d_t pi_i  = laplace(phi_i) - sum_j M2_ij * phi_j
        M2 = diag(m_i^2) + coupling (assumed symmetric).

    Raises
    ------
    ValueError
        If `coupling` is not an NxN array matching the number of masses.
    """
    masses = np.asarray(params.masses, dtype=float)
    field_count = masses.size
    coup = np.asarray(params.coupling, dtype=float)
    if coup.shape != (field_count, field_count):
        msg = "coupling must be NxN matching len(masses)"
        raise ValueError(msg)
    mass_matrix = np.diag(masses**2) + coup

    # Field names
    phi = [f"phi{i}" for i in range(field_count)]
    pi = [f"pi{i}" for i in range(field_count)]

    # Build RHS dict
    rhs: dict[str, str] = {}
    for i in range(field_count):
        rhs[phi[i]] = pi[i]
        # laplace(phi_i) - sum_j M2_ij * phi_j
        mix = " + ".join([f"(-M2_{i}_{j})*{phi[j]}" for j in range(field_count)]) or "0"
        rhs[pi[i]] = f"laplace({phi[i]}) + {mix}"

    # Constants for all M2_ij
    consts: dict[str, float] = {
        f"M2_{i}_{j}": float(mass_matrix[i, j])
        for i in range(field_count)
        for j in range(field_count)
    }

    return PDE(rhs, consts=cast("Any", consts))
