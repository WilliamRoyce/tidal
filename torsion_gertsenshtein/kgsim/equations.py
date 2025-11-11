from __future__ import annotations

from typing import TYPE_CHECKING, Any

from pde import PDE

if TYPE_CHECKING:
    from .config import KGParameters


class KleinGordonPDE(PDE):
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

    """

    explicit_time_dependence = False

    def __init__(self, params: KGParameters) -> None:
        m2 = float(params.mass) ** 2
        # store for observers, metadata, etc.
        self.m2: float = m2

        # Expression system with constants: two fields, 'phi' and 'pi'
        super().__init__(
            {
                "phi": "pi",
                "pi": "laplace(phi) - m2 * phi",
            },
            consts={"m2": m2},
        )

    # Optional: keep a cache key for reproducibility / memoization if you use it elsewhere
    def _cache_key(self) -> dict[str, Any]:
        return {"m2": self.m2}
