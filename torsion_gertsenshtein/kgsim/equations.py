from __future__ import annotations

from typing import TYPE_CHECKING, Any

from pde import FieldCollection, PDEBase, ScalarField
from typing_extensions import override

from torsion_gertsenshtein.kgsim.utils import mul_scalar_field, sub_scalar_fields

if TYPE_CHECKING:
    from pde.fields.datafield_base import DataFieldBase

    from .config import KGParameters


class KleinGordonPDE(PDEBase):
    """
    Klein-Gordon equation represented in first-order (in time) form.

    This class implements the system:
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

    Methods
    -------
    evolution_rate(state, t=0)
        Compute the time derivatives of the field pair (phi, pi).
        - Expects `state` to be a FieldCollection with exactly two entries:
          the ScalarField `phi` and the ScalarField `pi`.
        - Validates runtime types and raises TypeError with a clear message if
          the input is not a FieldCollection or the entries are not ScalarField
          instances.
        - Computes laplace(phi) using natural boundary conditions via
          `phi.laplace(bc="natural")`.
        - Returns a FieldCollection([dphi_dt, dpi_dt], labels=["phi", "pi"]) where
          dphi_dt is a copy of `pi` and dpi_dt = laplace(phi) - m2 * phi.
        - The method signature includes an optional time argument `t` for
          compatibility with time integrators; the PDE itself does not depend on t.

    _cache_key()
        Returns a dictionary describing the equation configuration used for
        caching / hashing (currently {"m2": self.m2}).

    Exceptions
    ----------
    TypeError
        Raised from evolution_rate when `state` is not a FieldCollection or when
        the contained fields are not ScalarField instances.

    Notes
    -----
    - Arithmetic between fields uses helpers (e.g. mul_scalar_field, sub_scalar_fields)
      to preserve the concrete ScalarField types expected by the framework.
    - Boundary conditions for the spatial Laplacian are set to "natural" by
      default in evolution_rate; change this call if different BCs are required.
    """

    explicit_time_dependence = False

    def __init__(self, params: KGParameters) -> None:
        super().__init__()
        self.m2 = float(params.mass) ** 2

    @override
    def evolution_rate(
        self, state: FieldCollection | DataFieldBase, t: float = 0
    ) -> FieldCollection:
        # Narrow the union to the concrete type we require.
        if not isinstance(state, FieldCollection):
            msg = "evolution_rate requires a FieldCollection with fields (phi, pi)"
            raise TypeError(msg)

        # Extract fields and validate runtime types so the type checker can accept the assignment.
        phi = state[0]
        pi = state[1]
        if not isinstance(phi, ScalarField) or not isinstance(pi, ScalarField):
            msg = "evolution_rate requires ScalarField instances for (phi, pi)"
            raise TypeError(msg)

        # Compute laplacian(φ)
        lap_phi = phi.laplace(bc="natural")

        dphi_dt = pi.copy()
        # ensure the subtraction is typed as a ScalarField
        scaled_phi = mul_scalar_field(self.m2, phi)
        dpi_dt = sub_scalar_fields(lap_phi, scaled_phi)

        return FieldCollection([dphi_dt, dpi_dt], labels=["phi", "pi"])

    def _cache_key(self) -> dict[str, Any]:
        return {"m2": self.m2}
