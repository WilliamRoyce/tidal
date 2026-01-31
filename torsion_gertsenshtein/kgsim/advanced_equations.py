"""Advanced Klein-Gordon equations with anisotropy and higher-order terms.

This module provides PDE classes for Klein-Gordon equations beyond the standard
isotropic case, including:
- Anisotropic wave propagation (directional wave speeds)
- Higher-order derivative terms (fourth-order dispersion, biharmonic)
- Directional-only dynamics (evolution in selected spatial dimensions)
- Base classes for custom operator combinations

All classes follow the first-order (in time) formulation:
    d/dt phi = pi
    d/dt pi  = [spatial operators] - m^2 * phi

Technical Notes
---------------
True anisotropic operators are implemented using py-pde's gradient() method:
- phi.gradient(bc) returns FieldCollection [∂φ/∂x, ∂φ/∂y, ∂φ/∂z, ...]
- Each component is a ScalarField that can be differentiated again
- grad_phi[i].gradient(bc)[i] gives ∂²φ/∂x_i² (second derivative in direction i)
- This approach works for arbitrary order: chain N gradients for 2N-th derivative

References
----------
- py-pde documentation: https://py-pde.readthedocs.io/
- ScalarField.gradient(): Returns VectorField/FieldCollection of partial derivatives
- ScalarField.laplace(): Computes isotropic Laplacian Σ_i ∂²φ/∂x_i²
- Custom operators: Subclass PDEBase, override evolution_rate(), use numpy backend
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from numba import (  # pyright: ignore[reportMissingTypeStubs,reportUnknownVariableType]
    jit,  # pyright: ignore[reportMissingTypeStubs,reportUnknownVariableType]
)
from pde import FieldCollection, PDEBase, ScalarField
from typing_extensions import override

from torsion_gertsenshtein.kgsim.utils import infer_bc_from_grid

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence

    from numpy.typing import NDArray
    from pde import FieldBase

    NumericArray = NDArray[np.float64]

import logging

_logger = logging.getLogger(__name__)


class AnisotropicKGPDE(PDEBase):
    """Klein-Gordon equation with anisotropic (direction-dependent) wave speeds.

    This implements the anisotropic wave equation:
        d/dt phi = pi
        d/dt pi  = c_x^2 * d^2 phi/dx^2 + c_y^2 * d^2 phi/dy^2 + ... - m^2 * phi

    where c_i are directional wave speeds. When all c_i are equal, this reduces
    to the standard isotropic Klein-Gordon equation.

    Implementation Details
    ----------------------
    Uses py-pde's gradient() chaining to compute true directional second derivatives:

    1. grad_phi = phi.gradient(bc) returns [∂φ/∂x, ∂φ/∂y, ...]
    2. For each direction i: grad_phi[i].gradient(bc)[i] gives ∂²φ/∂x_i²
    3. Sum with directional weights: Σ_i c_i² * ∂²φ/∂x_i²

    This is NOT an approximation - it computes true anisotropic Laplacian components.

    Parameters
    ----------
    mass : float
        Klein-Gordon mass parameter.
    speeds : Sequence[float]
        Wave speeds in each spatial direction [c_x, c_y, c_z]. Length must match
        the grid dimension. For isotropic case, use speeds=[1.0, 1.0, ...].

    Attributes
    ----------
    m2 : float
        Square of the mass parameter.
    speeds : np.ndarray
        Array of wave speeds in each direction.
    speeds_squared : np.ndarray
        Pre-computed squares of wave speeds for efficiency.

    Examples
    --------
    >>> from pde import CartesianGrid, FieldCollection, ScalarField
    >>> grid = CartesianGrid([[0, 10], [0, 10]], 64)
    >>> pde = AnisotropicKGPDE(mass=0.5, speeds=[2.0, 0.5])  # Fast in x, slow in y
    >>> phi = ScalarField.random_uniform(grid, -0.1, 0.1)
    >>> pi = ScalarField.random_uniform(grid, -0.1, 0.1)
    >>> state = FieldCollection([phi, pi])
    >>> result = pde.solve(state, t_range=10, dt=0.01, backend="numpy")

    Notes
    -----
    - Requires numpy backend (no numba support for custom evolution_rate)
    - Boundary conditions are inferred from the grid
    - For highly anisotropic cases, ensure adequate grid resolution in all directions
    - Grid spacing should satisfy CFL condition for fastest wave speed
    >>> pde = AnisotropicKGPDE(mass=0.5, speeds=[2.0, 0.5])  # Fast in x, slow in y
    >>> phi = ScalarField.random_uniform(grid, -0.1, 0.1)
    >>> pi = ScalarField.random_uniform(grid, -0.1, 0.1)
    >>> state = FieldCollection([phi, pi])
    >>> result = pde.solve(state, t_range=10, dt=0.01)

    Notes
    -----
    - Uses custom evolution_rate() with numpy backend only (no numba support)
    - Boundary conditions are inferred from the grid
    - For highly anisotropic cases, ensure adequate grid resolution in all directions
    """

    explicit_time_dependence = False

    def __init__(self, mass: float, speeds: Sequence[float]) -> None:
        """Initialize anisotropic Klein-Gordon PDE.

        Parameters
        ----------
        mass : float
            Klein-Gordon mass parameter (m in the equation).
        speeds : Sequence[float]
            Wave speeds in each spatial direction. Length must match grid dimension.

        Raises
        ------
        ValueError
            If speeds is empty or contains non-positive values.
        """
        super().__init__()
        self.m2 = mass**2
        self.speeds = np.array(speeds, dtype=float)
        self.speeds_squared = self.speeds**2

        if len(self.speeds) == 0:
            msg = "speeds must be a non-empty sequence"
            raise ValueError(msg)
        if np.any(self.speeds <= 0):
            msg = "All wave speeds must be positive"
            raise ValueError(msg)

        _logger.info(
            "Initialized AnisotropicKGPDE with mass=%.3f, speeds=%s",
            mass,
            self.speeds.tolist(),
        )

    @override
    def evolution_rate(
        self,
        state: FieldBase,
        t: float = 0.0,
    ) -> FieldCollection:
        """Compute the right-hand side of the PDE system.

        This implementation uses py-pde's gradient() method to compute true
        anisotropic second derivatives. The gradient returns a VectorField
        (or FieldCollection) where grad_phi[i] is ∂φ/∂x_i. We then compute
        grad_phi[i].gradient()[i] to get ∂²φ/∂x_i² for each direction i.

        The anisotropic wave equation is:
            ∂²φ/∂t² = Σ_i c_i² ∂²φ/∂x_i² - m²φ

        References
        ----------
        - py-pde ScalarField.gradient() returns FieldCollection of partial derivatives
        - Each component can be differentiated again to get second derivatives
        - This pattern is demonstrated in AnisotropicHigherOrderKGPDE (lines 536-550)

        Parameters
        ----------
        state : FieldBase
            Current state as FieldCollection([phi, pi]).
        t : float, optional
            Current time (unused for time-independent PDE).

        Returns
        -------
        FieldCollection
            Time derivatives [d_phi/dt, d_pi/dt] = [pi, Σ c_i²∂²φ/∂x_i² - m²φ].

        Raises
        ------
        ValueError
            If speeds length does not match grid dimension.
        """
        assert isinstance(state, FieldCollection)
        assert len(state) == 2, "State must contain [phi, pi]"  # noqa: PLR2004

        phi = state[0]
        pi = state[1]
        assert isinstance(phi, ScalarField)
        assert isinstance(pi, ScalarField)
        # Validate dimension
        grid_dim = phi.grid.dim
        if len(self.speeds) != grid_dim:
            msg = f"speeds length ({len(self.speeds)}) must match grid dimension ({grid_dim})"
            raise ValueError(msg)

        # Infer boundary conditions
        bc = infer_bc_from_grid(phi.grid)

        # Compute anisotropic second derivatives: sum_i c_i^2 * d^2phi/dx_i^2
        # Using gradient chains (same pattern as AnisotropicHigherOrderKGPDE)
        grad_phi = phi.gradient(bc=bc)  # Returns FieldCollection: [∂φ/∂x, ∂φ/∂y, ...]

        spatial_term: ScalarField = ScalarField(phi.grid, data=np.zeros_like(phi.data))
        for i, c_squared in enumerate(self.speeds_squared):
            # ∂²φ/∂x_i² = ∂/∂x_i (∂φ/∂x_i) = (gradient of i-th component)[i]
            grad_component = grad_phi[i]
            assert isinstance(grad_component, ScalarField)
            d2_phi_dxi2: ScalarField = grad_component.gradient(bc=bc)[i]
            spatial_term += c_squared * d2_phi_dxi2

        # Klein-Gordon evolution: d_pi/dt = anisotropic_laplacian - m^2 * phi
        dpi_dt: ScalarField = spatial_term - self.m2 * phi  # type: ignore[assignment]
        # First-order system: d_phi/dt = pi
        # Note: Must create new ScalarField, not return reference to input pi
        dphi_dt = ScalarField(phi.grid, data=pi.data.copy())
        return FieldCollection([dphi_dt, dpi_dt])

    def make_pde_rhs_numba(
        self, state: FieldCollection
    ) -> Callable[[NumericArray, float], NumericArray]:
        """Create a compiled function evaluating the right hand side of the PDE.

        This implements anisotropic wave propagation using gradient chaining:
        For each direction i: compute c_i² * ∂²φ/∂x_i² using two gradient applications.

        Args:
            state (:class:`~pde.fields.FieldCollection`):
                An example for the state defining the grid and data types

        Returns
        -------
            A function with signature `(state_data, t)`, which can be called with an
            instance of :class:`~numpy.ndarray` of the state data and the time to
            obtained an instance of :class:`~numpy.ndarray` giving the evolution rate.

        Raises
        ------
        ValueError
            If speeds length does not match grid dimension.
        """
        # Copy attributes locally (frozen during compilation)
        m2 = self.m2
        speeds_squared = self.speeds_squared.copy()
        dim = len(self.speeds)

        # Validate dimension matches grid
        if dim != state.grid.dim:
            msg = f"speeds length ({dim}) must match grid dimension ({state.grid.dim})"
            raise ValueError(msg)

        # Get boundary conditions
        bc = infer_bc_from_grid(state.grid)

        # Create compiled gradient operator
        gradient = state.grid.make_operator("gradient", bc=bc, backend="numba")  # type: ignore[arg-type]

        @jit
        def pde_rhs(state_data: NumericArray, t: float = 0) -> NumericArray:
            """Compiled helper function evaluating right hand side."""
            # Extract phi and pi from state array
            phi_data = state_data[0]
            pi_data = state_data[1]

            # Prepare output array
            rate = np.empty_like(state_data)

            # dphi/dt = pi
            rate[0] = pi_data

            # dpi/dt = Σ_i c_i² * ∂²φ/∂x_i² - m²φ
            # Start with mass term
            rate[1] = -m2 * phi_data

            # Anisotropic Laplacian via gradient chaining
            # grad_phi has shape (dim, *grid_shape)
            grad_phi = gradient(phi_data, args={"t": t})

            for i in range(dim):
                # Extract i-th component: ∂φ/∂x_i
                grad_i = grad_phi[i]

                # Compute gradient of this component
                # grad_grad_i has shape (dim, *grid_shape)
                grad_grad_i = gradient(grad_i, args={"t": t})

                # Extract i-th component: ∂²φ/∂x_i²
                d2_phi_dxi2 = grad_grad_i[i]

                # Add anisotropic contribution
                rate[1] += speeds_squared[i] * d2_phi_dxi2

            return rate

        return pde_rhs


class HigherOrderKGPDE(PDEBase):
    """Klein-Gordon equation with higher-order derivative terms.

    This implements Klein-Gordon with fourth-order (and optionally sixth-order)
    spatial derivatives:
        d/dt phi = pi
        d/dt pi  = alpha_2 * laplace(phi) - alpha_4 * laplace^2(phi)
                   + alpha_6 * laplace^3(phi) - m^2 * phi

    Higher-order terms appear in:
    - Quantum corrections to classical field theories
    - Beam equations and plate theories
    - Regularization of singular solutions
    - Modified dispersion relations

    Parameters
    ----------
    mass : float
        Klein-Gordon mass parameter.
    alpha_2 : float, optional
        Coefficient for standard Laplacian term (default: 1.0).
    alpha_4 : float, optional
        Coefficient for bi-Laplacian (fourth-order) term (default: 0.0).
    alpha_6 : float, optional
        Coefficient for tri-Laplacian (sixth-order) term (default: 0.0).

    Attributes
    ----------
    m2 : float
        Square of the mass parameter.
    alpha_2, alpha_4, alpha_6 : float
        Coefficients for Laplacian, bi-Laplacian, and tri-Laplacian terms.

    Examples
    --------
    >>> from pde import UnitGrid, FieldCollection, ScalarField
    >>> grid = UnitGrid([128], periodic=True)  # Periodic BC recommended
    >>> pde = HigherOrderKGPDE(mass=0.5, alpha_2=1.0, alpha_4=0.01)
    >>> phi = ScalarField.random_uniform(grid, -0.1, 0.1)
    >>> pi = ScalarField.random_uniform(grid, -0.1, 0.1)
    >>> state = FieldCollection([phi, pi])
    >>> result = pde.solve(state, t_range=10, dt=0.001)  # Smaller dt needed!

    Notes
    -----
    - Higher-order terms require smaller time steps for stability
    - Periodic boundary conditions strongly recommended
    - Fourth-order: dt ~ O(dx^4), Sixth-order: dt ~ O(dx^6)
    - Nested Laplacians may amplify boundary artifacts with non-periodic BCs
    """

    explicit_time_dependence = False

    def __init__(
        self,
        mass: float,
        alpha_2: float = 1.0,
        alpha_4: float = 0.0,
        alpha_6: float = 0.0,
    ) -> None:
        """Initialize higher-order Klein-Gordon PDE.

        Parameters
        ----------
        mass : float
            Klein-Gordon mass parameter.
        alpha_2 : float, optional
            Coefficient for laplace(phi) term (default: 1.0 for standard KG).
        alpha_4 : float, optional
            Coefficient for laplace^2(phi) term (default: 0.0).
        alpha_6 : float, optional
            Coefficient for laplace^3(phi) term (default: 0.0).
        """
        super().__init__()
        self.m2 = mass**2
        self.alpha_2 = alpha_2
        self.alpha_4 = alpha_4
        self.alpha_6 = alpha_6

        _logger.info(
            "Initialized HigherOrderKGPDE with mass=%.3f, "
            "alpha_2=%.3e, alpha_4=%.3e, alpha_6=%.3e",
            mass,
            alpha_2,
            alpha_4,
            alpha_6,
        )

        if alpha_4 != 0 or alpha_6 != 0:
            _logger.warning(
                "Higher-order terms require smaller time steps for stability. "
                "Consider dt ~ O(dx^4) for fourth-order, dt ~ O(dx^6) for sixth-order."
            )

    @override
    def evolution_rate(
        self,
        state: FieldBase,
        t: float = 0.0,
    ) -> FieldCollection:
        """Compute the right-hand side of the PDE system.

        Parameters
        ----------
        state : FieldBase
            Current state as FieldCollection([phi, pi]).
        t : float, optional
            Current time (unused for time-independent PDE).

        Returns
        -------
        FieldCollection
            Time derivatives [d_phi/dt, d_pi/dt].
        """
        assert isinstance(state, FieldCollection)
        assert len(state) == 2, "State must contain [phi, pi]"  # noqa: PLR2004

        phi = state[0]
        pi = state[1]
        assert isinstance(phi, ScalarField)
        assert isinstance(pi, ScalarField)
        # Infer boundary conditions
        bc = infer_bc_from_grid(phi.grid)

        # Build spatial terms progressively
        spatial_term: ScalarField = ScalarField(phi.grid, data=np.zeros_like(phi.data))

        # Second-order: alpha_2 * laplace(phi)
        lap_phi: ScalarField | None
        if self.alpha_2 != 0:
            lap_phi = phi.laplace(bc=bc)
            spatial_term += self.alpha_2 * lap_phi
        else:
            lap_phi = None

        # Fourth-order: -alpha_4 * laplace^2(phi)
        lap2_phi: ScalarField | None
        if self.alpha_4 != 0:
            if lap_phi is None:
                lap_phi = phi.laplace(bc=bc)
            lap2_phi = lap_phi.laplace(bc=bc)
            spatial_term -= self.alpha_4 * lap2_phi
        else:
            lap2_phi = None

        # Sixth-order: +alpha_6 * laplace^3(phi)
        if self.alpha_6 != 0:
            if lap2_phi is None:
                if lap_phi is None:
                    lap_phi = phi.laplace(bc=bc)
                lap2_phi = lap_phi.laplace(bc=bc)
            lap3_phi: ScalarField = lap2_phi.laplace(bc=bc)
            spatial_term += self.alpha_6 * lap3_phi

        # Klein-Gordon evolution: d_pi/dt = spatial_terms - m^2 * phi
        dpi_dt = spatial_term - self.m2 * phi

        # First-order system: d_phi/dt = pi
        dphi_dt = ScalarField(phi.grid, data=pi.data.copy())
        return FieldCollection([dphi_dt, dpi_dt])  # type: ignore[arg-type]

    def make_pde_rhs_numba(
        self, state: FieldCollection
    ) -> Callable[[NumericArray, float], NumericArray]:
        """Create a compiled function evaluating the right hand side of the PDE.

        Args:
            state (:class:`~pde.fields.FieldCollection`):
                An example for the state defining the grid and data types

        Returns
        -------
            A function with signature `(state_data, t)`, which can be called with an
            instance of :class:`~numpy.ndarray` of the state data and the time to
            obtained an instance of :class:`~numpy.ndarray` giving the evolution rate.
        """
        # Copy attributes locally (frozen during compilation)
        m2 = self.m2
        alpha_2 = self.alpha_2
        alpha_4 = self.alpha_4
        alpha_6 = self.alpha_6

        # Get boundary conditions
        bc = infer_bc_from_grid(state.grid)

        # Create compiled operators
        laplace = state.grid.make_operator("laplace", bc=bc, backend="numba")  # type: ignore[arg-type]

        @jit
        def pde_rhs(state_data: NumericArray, t: float = 0) -> NumericArray:
            """Compiled helper function evaluating right hand side."""
            # Extract phi and pi from state array
            phi_data = state_data[0]
            pi_data = state_data[1]

            # Prepare output array
            rate = np.empty_like(state_data)

            # dphi/dt = pi
            rate[0] = pi_data

            # dpi/dt = alpha_2 * ∇²φ - alpha_4 * ∇⁴φ + alpha_6 * ∇⁶φ - m²φ
            rate[1] = -m2 * phi_data

            # Nested Laplacians for higher-order terms
            if alpha_2 != 0 or alpha_4 != 0 or alpha_6 != 0:
                lap1 = laplace(phi_data, args={"t": t})
                if alpha_2 != 0:
                    rate[1] += alpha_2 * lap1

                if alpha_4 != 0 or alpha_6 != 0:
                    lap2 = laplace(lap1, args={"t": t})
                    if alpha_4 != 0:
                        rate[1] -= alpha_4 * lap2

                    if alpha_6 != 0:
                        lap3 = laplace(lap2, args={"t": t})
                        rate[1] += alpha_6 * lap3

            return rate

        return pde_rhs


class DirectionalKGPDE(PDEBase):
    """Klein-Gordon equation with evolution only in selected spatial directions.

    This implements a "directional waveguide" where the field evolves only along
    specified axes:
        d/dt phi = pi
        d/dt pi  = sum_{i in active} d^2 phi/dx_i^2 - m^2 * phi

    Useful for:
    - Studying dimensional reduction
    - Waveguide physics (1D propagation in 2D/3D space)
    - Testing anisotropic limit cases
    - Validating multi-dimensional codes against 1D results

    Implementation Details
    ----------------------
    Uses gradient chaining to compute only active directional second derivatives:
    - For each active direction i: compute ∂²φ/∂x_i²
    - Inactive directions contribute zero to spatial evolution
    - This is exact, not an approximation

    Parameters
    ----------
    mass : float
        Klein-Gordon mass parameter.
    active_directions : Sequence[bool]
        Boolean flags indicating which spatial directions are active.
        Length must match grid dimension. E.g., (True, False) means
        evolution only in x-direction for a 2D grid.

    Attributes
    ----------
    m2 : float
        Square of the mass parameter.
    active_directions : tuple[bool, ...]
        Tuple of boolean flags for active directions.

    Examples
    --------
    >>> from pde import CartesianGrid, FieldCollection, ScalarField
    >>> grid = CartesianGrid([[0, 10], [0, 10]], 64)
    >>> # Evolution only in x-direction (waveguide along x-axis)
    >>> pde = DirectionalKGPDE(mass=0.5, active_directions=(True, False))
    >>> phi = ScalarField.random_uniform(grid, -0.1, 0.1)
    >>> pi = ScalarField.random_uniform(grid, -0.1, 0.1)
    >>> state = FieldCollection([phi, pi])
    >>> result = pde.solve(state, t_range=10, dt=0.01, backend="numpy")

    Notes
    -----
    - At least one direction must be active
    - For all directions active, equivalent to standard isotropic KG
    - For one direction active, effective 1D propagation in higher-D space
    >>> # Evolution only in x-direction (waveguide along x-axis)
    >>> pde = DirectionalKGPDE(mass=0.5, active_directions=(True, False))
    >>> phi = ScalarField.random_uniform(grid, -0.1, 0.1)
    >>> pi = ScalarField.random_uniform(grid, -0.1, 0.1)
    >>> state = FieldCollection([phi, pi])
    >>> result = pde.solve(state, t_range=10, dt=0.01)

    Notes
    -----
    - At least one direction must be active
    - For all directions active, equivalent to standard isotropic KG
    - For one direction active, effective 1D propagation in higher-D space
    """

    explicit_time_dependence = False

    def __init__(self, mass: float, active_directions: Sequence[bool]) -> None:
        """Initialize directional Klein-Gordon PDE.

        Parameters
        ----------
        mass : float
            Klein-Gordon mass parameter.
        active_directions : Sequence[bool]
            Boolean flags for active spatial directions. Length must match grid dimension.

        Raises
        ------
        ValueError
            If active_directions is empty or if no directions are active.
        """
        super().__init__()
        self.m2 = mass**2
        self.active_directions = tuple(active_directions)

        if len(self.active_directions) == 0:
            msg = "active_directions must be non-empty"
            raise ValueError(msg)
        if not any(self.active_directions):
            msg = "At least one direction must be active"
            raise ValueError(msg)

        active_count = sum(self.active_directions)
        _logger.info(
            "Initialized DirectionalKGPDE with mass=%.3f, "
            "active_directions=%s (%d/%d active)",
            mass,
            self.active_directions,
            active_count,
            len(self.active_directions),
        )

    @override
    def evolution_rate(
        self,
        state: FieldBase,
        t: float = 0.0,
    ) -> FieldCollection:
        """Compute the right-hand side of the PDE system.

        Uses gradient chaining to compute second derivatives only in active directions.

        Parameters
        ----------
        state : FieldBase
            Current state as FieldCollection([phi, pi]).
        t : float, optional
            Current time (unused for time-independent PDE).

        Returns
        -------
        FieldCollection
            Time derivatives [d_phi/dt, d_pi/dt].

        Raises
        ------
        ValueError
            If active_directions length does not match grid dimension.
        """
        assert isinstance(state, FieldCollection)
        assert len(state) == 2, "State must contain [phi, pi]"  # noqa: PLR2004

        phi = state[0]
        pi = state[1]

        # Validate dimension
        grid_dim = phi.grid.dim
        if len(self.active_directions) != grid_dim:
            msg = (
                f"active_directions length ({len(self.active_directions)}) "
                f"must match grid dimension ({grid_dim})"
            )
            raise ValueError(msg)

        # Infer boundary conditions
        bc = infer_bc_from_grid(phi.grid)

        # Compute gradient (all directions)
        grad_phi = phi.gradient(bc=bc)  # type: ignore[attr-defined]

        # Sum second derivatives only for active directions
        spatial_term: ScalarField = ScalarField(phi.grid, data=np.zeros_like(phi.data))
        for i, is_active in enumerate(self.active_directions):
            if is_active:
                # d^2phi/dx_i^2 = d/dx_i (d phi/dx_i)
                grad_component: ScalarField = grad_phi[i]  # pyright: ignore[reportUnknownVariableType]
                assert isinstance(grad_component, ScalarField)
                second_deriv: ScalarField = grad_component.gradient(bc=bc)[i]
                spatial_term += second_deriv

        # Klein-Gordon evolution: d_pi/dt = directional_laplacian - m^2 * phi
        dpi_dt = spatial_term - self.m2 * phi

        # First-order system: d_phi/dt = pi
        dphi_dt = ScalarField(phi.grid, data=pi.data.copy())
        return FieldCollection([dphi_dt, dpi_dt])  # type: ignore[arg-type]

    def make_pde_rhs_numba(
        self, state: FieldCollection
    ) -> Callable[[NumericArray, float], NumericArray]:
        """Create a compiled function evaluating the right hand side of the PDE.

        This implements directional evolution with selective gradient chaining:
        Only compute ∂²φ/∂x_i² for active directions.

        Args:
            state (:class:`~pde.fields.FieldCollection`):
                An example for the state defining the grid and data types

        Returns
        -------
            A function with signature `(state_data, t)`, which can be called with an
            instance of :class:`~numpy.ndarray` of the state data and the time to
            obtained an instance of :class:`~numpy.ndarray` giving the evolution rate.

        Raises
        ------
        ValueError
            If active_directions length does not match grid dimension.
        """
        # Copy attributes locally (frozen during compilation)
        m2 = self.m2
        active_directions = tuple(self.active_directions)
        dim = len(active_directions)

        # Validate dimension matches grid
        if dim != state.grid.dim:
            msg = f"active_directions length ({dim}) must match grid dimension ({state.grid.dim})"
            raise ValueError(msg)

        # Get boundary conditions
        bc = infer_bc_from_grid(state.grid)

        # Create compiled gradient operator
        gradient = state.grid.make_operator("gradient", bc=bc, backend="numba")  # type: ignore[arg-type]

        @jit
        def pde_rhs(state_data: NumericArray, t: float = 0) -> NumericArray:
            """Compiled helper function evaluating right hand side."""
            # Extract phi and pi from state array
            phi_data = state_data[0]
            pi_data = state_data[1]

            # Prepare output array
            rate = np.empty_like(state_data)

            # dphi/dt = pi
            rate[0] = pi_data

            # dpi/dt = sum over active directions of d2phi/dxi2 - m2*phi
            rate[1] = -m2 * phi_data

            # Directional Laplacian via gradient chaining
            grad_phi = gradient(phi_data, args={"t": t})

            for i in range(dim):
                if active_directions[i]:
                    # Extract i-th component: ∂φ/∂x_i
                    grad_i = grad_phi[i]

                    # Compute gradient of this component
                    grad_grad_i = gradient(grad_i, args={"t": t})

                    # Extract i-th component: ∂²φ/∂x_i²
                    d2_phi_dxi2 = grad_grad_i[i]

                    # Add contribution from active direction
                    rate[1] += d2_phi_dxi2

            return rate

        return pde_rhs


class AnisotropicHigherOrderKGPDE(PDEBase):
    """Combined anisotropic and higher-order Klein-Gordon equation.

    This combines directional wave speeds with higher-order dispersion:
        d/dt phi = pi
        d/dt pi  = sum_i [alpha_2 * c_i^2 * d^2/dx_i^2
                          - alpha_4 * c_i^4 * d^4/dx_i^4 + ...] phi - m^2 * phi

    Parameters
    ----------
    mass : float
        Klein-Gordon mass parameter.
    speeds : Sequence[float]
        Wave speeds in each spatial direction.
    alpha_2 : float, optional
        Coefficient for second-order terms (default: 1.0).
    alpha_4 : float, optional
        Coefficient for fourth-order terms (default: 0.0).

    Examples
    --------
    >>> pde = AnisotropicHigherOrderKGPDE(
    ...     mass=0.5,
    ...     speeds=[2.0, 0.5],
    ...     alpha_2=1.0,
    ...     alpha_4=0.01
    ... )

    Notes
    -----
    - Combines both anisotropy and dispersion effects
    - Requires even smaller time steps than either effect alone
    - Useful for realistic material models with direction-dependent dispersion
    """

    explicit_time_dependence = False

    def __init__(
        self,
        mass: float,
        speeds: Sequence[float],
        alpha_2: float = 1.0,
        alpha_4: float = 0.0,
    ) -> None:
        """Initialize combined anisotropic higher-order KG PDE.

        Parameters
        ----------
        mass : float
            Klein-Gordon mass parameter.
        speeds : Sequence[float]
            Wave speeds in each spatial direction.
        alpha_2 : float, optional
            Coefficient for second-order terms (default: 1.0).
        alpha_4 : float, optional
            Coefficient for fourth-order terms (default: 0.0).

        Raises
        ------
        ValueError
            If speeds is empty or contains non-positive values.
        """
        super().__init__()
        self.m2 = mass**2
        self.speeds = np.array(speeds, dtype=float)
        self.alpha_2 = alpha_2
        self.alpha_4 = alpha_4

        if len(self.speeds) == 0:
            msg = "speeds must be non-empty"
            raise ValueError(msg)
        if np.any(self.speeds <= 0):
            msg = "All wave speeds must be positive"
            raise ValueError(msg)

        _logger.info(
            "Initialized AnisotropicHigherOrderKGPDE with mass=%.3f, "
            "speeds=%s, alpha_2=%.3e, alpha_4=%.3e",
            mass,
            self.speeds.tolist(),
            alpha_2,
            alpha_4,
        )

    @override
    def evolution_rate(
        self,
        state: FieldBase,
        t: float = 0.0,
    ) -> FieldCollection:
        """Compute the right-hand side of the PDE system.

        Parameters
        ----------
        state : FieldBase
            Current state as FieldCollection([phi, pi]).
        t : float, optional
            Current time (unused for time-independent PDE).

        Returns
        -------
        FieldCollection
            Time derivatives [d_phi/dt, d_pi/dt].

        Raises
        ------
        ValueError
            If speeds length does not match grid dimension.
        """
        assert isinstance(state, FieldCollection)
        assert len(state) == 2, "State must contain [phi, pi]"  # noqa: PLR2004

        phi = state[0]
        pi = state[1]
        assert isinstance(phi, ScalarField)
        assert isinstance(pi, ScalarField)
        grid_dim = phi.grid.dim
        if len(self.speeds) != grid_dim:
            msg = f"speeds length ({len(self.speeds)}) must match grid dimension ({grid_dim})"
            raise ValueError(msg)

        bc = infer_bc_from_grid(phi.grid)

        # Compute directional derivatives
        grad_phi = phi.gradient(bc=bc)
        spatial_term: ScalarField = ScalarField(phi.grid, data=np.zeros_like(phi.data))

        for i, c_i in enumerate(self.speeds):
            # Second-order: alpha_2 * c_i^2 * d^2phi/dx_i^2
            if self.alpha_2 != 0:
                grad_component = grad_phi[i]
                assert isinstance(grad_component, ScalarField)
                d2_phi: ScalarField = grad_component.gradient(bc=bc)[i]
                spatial_term += self.alpha_2 * c_i**2 * d2_phi

                # Fourth-order: -alpha_4 * c_i^4 * d^4phi/dx_i^4
                if self.alpha_4 != 0:
                    # d^4/dx^4 = d^2/dx^2(d^2/dx^2)
                    d2_phi_grad = d2_phi.gradient(bc=bc)
                    d2_phi_grad_component = d2_phi_grad[i]
                    assert isinstance(d2_phi_grad_component, ScalarField)
                    d4_phi: ScalarField = d2_phi_grad_component.gradient(bc=bc)[i]
                    spatial_term -= self.alpha_4 * c_i**4 * d4_phi

        dpi_dt = spatial_term - self.m2 * phi
        dphi_dt = ScalarField(phi.grid, data=pi.data.copy())
        return FieldCollection([dphi_dt, dpi_dt])  # type: ignore[arg-type]

    def make_pde_rhs_numba(
        self, state: FieldCollection
    ) -> Callable[[NumericArray, float], NumericArray]:
        """Create a compiled function evaluating the right hand side of the PDE.

        This implements combined anisotropic and higher-order dispersion using
        nested gradient chaining: up to 4x gradient applications per direction.

        Args:
            state (:class:`~pde.fields.FieldCollection`):
                An example for the state defining the grid and data types

        Returns
        -------
            A function with signature `(state_data, t)`, which can be called with an
            instance of :class:`~numpy.ndarray` of the state data and the time to
            obtained an instance of :class:`~numpy.ndarray` giving the evolution rate.

        Raises
        ------
        ValueError
            If speeds length does not match grid dimension.
        """
        # Copy attributes locally (frozen during compilation)
        m2 = self.m2
        speeds = self.speeds.copy()
        alpha_2 = self.alpha_2
        alpha_4 = self.alpha_4
        dim = len(self.speeds)

        # Validate dimension matches grid
        if dim != state.grid.dim:
            msg = f"speeds length ({dim}) must match grid dimension ({state.grid.dim})"
            raise ValueError(msg)

        # Get boundary conditions
        bc = infer_bc_from_grid(state.grid)

        # Create compiled gradient operator
        gradient = state.grid.make_operator("gradient", bc=bc, backend="numba")  # type: ignore[arg-type]

        @jit
        def pde_rhs(state_data: NumericArray, t: float = 0) -> NumericArray:
            """Compiled helper function evaluating right hand side."""
            # Extract phi and pi from state array
            phi_data = state_data[0]
            pi_data = state_data[1]

            # Prepare output array
            rate = np.empty_like(state_data)

            # dphi/dt = pi
            rate[0] = pi_data

            # dpi/dt = anisotropic second and fourth order terms minus mass term
            rate[1] = -m2 * phi_data

            # Anisotropic higher-order terms via gradient chaining
            grad_phi = gradient(phi_data, args={"t": t})

            for i in range(dim):
                c_i = speeds[i]

                # Second-order: α₂ * c_i² * ∂²φ/∂x_i²
                if alpha_2 != 0:
                    # First gradient: ∂φ/∂x_i
                    grad_i = grad_phi[i]

                    # Second gradient
                    grad_grad_i = gradient(grad_i, args={"t": t})
                    d2_phi_dxi2 = grad_grad_i[i]  # ∂²φ/∂x_i²

                    rate[1] += alpha_2 * c_i**2 * d2_phi_dxi2

                    # Fourth-order: -α₄ * c_i⁴ * ∂⁴φ/∂x_i⁴
                    if alpha_4 != 0:
                        # Third gradient: ∂/∂x_j (∂²φ/∂x_i²) for all j
                        grad_d2phi = gradient(d2_phi_dxi2, args={"t": t})

                        # Fourth gradient: ∂/∂x_i (∂³φ/∂x_i²∂x_i)
                        grad_grad_d2phi = gradient(grad_d2phi[i], args={"t": t})
                        d4_phi_dxi4 = grad_grad_d2phi[i]  # ∂⁴φ/∂x_i⁴

                        rate[1] -= alpha_4 * c_i**4 * d4_phi_dxi4

            return rate

        return pde_rhs
