"""Build py-pde PDEBase subclasses from equation specifications.

This module provides the core functionality for converting symbolically-derived
field equations (loaded from JSON) into executable py-pde PDE classes.

The key principle is that NO physics is hardcoded here - all equation structure
comes from the specification that was derived from the Lagrangian.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from pde import FieldCollection, PDEBase, ScalarField
from typing_extensions import override

from torsion_gertsenshtein.kgsim.utils import infer_bc_from_grid
from torsion_gertsenshtein.symbolic.json_loader import (
    EquationSystem,
    OperatorTerm,
    load_equation_system,
)

if TYPE_CHECKING:
    from pathlib import Path

    import numpy as np
    from numpy.typing import NDArray
    from pde.grids.base import GridBase
    from pde.pdes.base import TState

    from torsion_gertsenshtein.utils import BCDescriptor

    NumericArray = NDArray[np.float64]


@dataclass(frozen=True)
class ParsedFieldName:
    """Parsed field name components.

    Supports multiple field naming conventions:
    - standard: A_0, phi_1 (base_index)
    - tensor: stress_xy_0, u_x_1 (base_component_index)
    - compact: phi0, A1 (base+digits)
    - simple: phi, psi (letters only, index defaults to 0)
    """

    base: str
    index: int
    format: str

    @classmethod
    def parse(cls, name: str) -> ParsedFieldName:
        """Parse field name auto-detecting format.

        Parameters
        ----------
        name : str
            Field name to parse.

        Returns
        -------
        ParsedFieldName
            Parsed components with base, index, and format.
        """
        # Standard format: A_0, phi_1
        match = re.match(r"^([a-zA-Z]+)_([0-9]+)$", name)
        if match:
            return cls(base=match.group(1), index=int(match.group(2)), format="standard")

        # Tensor format: stress_xy_0, u_x_1 (greedy match for base)
        match = re.match(r"^(.+)_([0-9]+)$", name)
        if match:
            return cls(base=match.group(1), index=int(match.group(2)), format="tensor")

        # Compact format: phi0, A1
        match = re.match(r"^([a-zA-Z]+)([0-9]+)$", name)
        if match:
            return cls(base=match.group(1), index=int(match.group(2)), format="compact")

        # Simple format: phi, psi (no index, defaults to 0)
        if re.match(r"^[a-zA-Z]+$", name):
            return cls(base=name, index=0, format="simple")

        # Fallback
        return cls(base=name, index=0, format="unknown")

    def to_momentum_name(self) -> str:
        """Convert to momentum field name."""
        return f"pi_{self.index}"


def parse_momentum_field_name(field_name: str) -> int | None:
    """Parse momentum field name and return index.

    Supports both pi_N and piN formats.

    Parameters
    ----------
    field_name : str
        Momentum field name like "pi_0", "pi0", "pi_1", "pi1".

    Returns
    -------
    int | None
        Index if valid momentum field name, None otherwise.
    """
    # Standard format: pi_N
    match = re.match(r"^pi_([0-9]+)$", field_name)
    if match:
        return int(match.group(1))

    # Compact format: piN
    match = re.match(r"^pi([0-9]+)$", field_name)
    if match:
        return int(match.group(1))

    return None


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
        True to support time-dependent coefficients (e.g., de Sitter spacetime).

    Examples
    --------
    >>> from torsion_gertsenshtein.symbolic import load_equation_system
    >>> from torsion_gertsenshtein.symbolic.pde_builder import PDEFromSpec
    >>> spec = load_equation_system("examples/data/em_1d.json")
    >>> pde = PDEFromSpec(spec)
    >>> # pde can now be used with py-pde solvers
    """

    # Enable time-dependent coefficients for curved spacetime (e.g., de Sitter)
    # This allows evolution_rate to receive the current time t
    explicit_time_dependence = True

    def __init__(
        self,
        spec: EquationSystem,
        parameters: dict[str, float] | None = None,
    ) -> None:
        """Initialize PDE from equation specification.

        Parameters
        ----------
        spec : EquationSystem
            The equation specification loaded from JSON.
        parameters : dict[str, float] | None
            Optional parameter values to override symbolic coefficients.
            Keys are symbolic names (e.g., "m2", "kappa"), values are numeric.
            When a term has a coefficient_symbolic that matches a key in this
            dict, the parameter value is used instead of the numeric coefficient.

            For time-dependent coefficients in curved spacetime:
            - "H": Hubble parameter for de Sitter expansion (e.g., 0.1)
            - "m2": Mass squared for Klein-Gordon field

            Time-dependent terms like -2H∂_t φ are evaluated using these parameters.
        """
        super().__init__()
        self.spec = spec
        self.n_components = spec.n_components
        self._component_name_to_index = {
            name: i for i, name in enumerate(spec.component_names)
        }
        self._parameters = parameters or {}

    def _resolve_coefficient(self, term: OperatorTerm) -> float:
        """Resolve the effective coefficient for a term.

        If the term has a symbolic coefficient name and that name (or its
        negation) is in the parameters dict, use the parameter value.
        Otherwise use the numeric coefficient from the JSON.

        Parameters
        ----------
        term : OperatorTerm
            The term whose coefficient to resolve.

        Returns
        -------
        float
            The effective coefficient value.
        """
        if term.coefficient_symbolic is not None:
            sym = term.coefficient_symbolic

            # Check for negated symbol like "-m2"
            if sym.startswith("-") and sym[1:] in self._parameters:
                return -self._parameters[sym[1:]]
            if sym in self._parameters:
                return self._parameters[sym]

        # Default: use numeric coefficient from JSON
        return term.coefficient

    def _resolve_coefficient_at_time(  # noqa: C901
        self, term: OperatorTerm, t: float
    ) -> float:
        """Resolve coefficient that may depend on time.

        For curved spacetime with time-dependent metrics (e.g., de Sitter expansion),
        coefficients may depend on time. This method handles common patterns:

        - "-2*H" (Hubble friction): Returns -2*H where H is from parameters
        - "m2*exp(2*H*t)" (time-dependent mass): Returns -m2*exp(2*H*t)
        - Other symbolic expressions: Attempts to evaluate or falls back to base coefficient

        Parameters
        ----------
        term : OperatorTerm
            The term whose coefficient to resolve.
        t : float
            Current simulation time.

        Returns
        -------
        float
            The effective coefficient value at time t.

        Raises
        ------
        ValueError
            If required parameters (H, m2) are missing or coefficient cannot be parsed.
        """
        # For non-time-dependent terms, use standard resolution
        if not term.time_dependent:
            return self._resolve_coefficient(term)

        # Get symbolic expression
        sym = term.coefficient_symbolic or ""

        # Pattern: n*H (Hubble friction, where n = spatial dimensions)
        # In 1+1D (n=1): -H∂_t φ, in 2+1D (n=2): -2H∂_t φ, in 3+1D (n=3): -3H∂_t φ
        if "H" in sym or "h" in sym.lower():
            # Require explicit Hubble parameter - no defaults
            if "H" not in self._parameters:
                msg = (
                    f"Hubble parameter 'H' is required for time-dependent coefficient '{sym}'. "
                    f"Pass parameters={{'H': value}} to PDEFromSpec or build_pde_from_json."
                )
                raise ValueError(
                    msg
                )
            hubble_param = self._parameters["H"]

            # Match patterns like "-2*H", "2H", "-H", "H", "-2 H", etc.
            match = re.match(r"^(-?\d*\.?\d*)\s*\*?\s*[Hh]$", sym.strip())
            if not match:
                msg = (
                    f"Cannot parse Hubble friction coefficient '{sym}'. "
                    f"Expected format like '-H', '-2*H', '3H', etc."
                )
                raise ValueError(
                    msg
                )

            multiplier_str = match.group(1)
            # Validate the multiplier string before conversion
            if multiplier_str in {"", None}:
                multiplier = 1.0
            elif multiplier_str == "-":
                multiplier = -1.0
            else:
                try:
                    multiplier = float(multiplier_str)
                except ValueError as e:
                    msg = f"Invalid numeric multiplier '{multiplier_str}' in coefficient '{sym}'"
                    raise ValueError(
                        msg
                    ) from e
            return multiplier * hubble_param

        # Pattern: exp(2*H*t) or similar exponential
        # This appears in time-dependent mass terms
        if "exp" in sym.lower():
            # Require explicit parameters - no defaults
            if "m2" not in self._parameters:
                msg = (
                    f"Mass parameter 'm2' is required for time-dependent coefficient '{sym}'. "
                    f"Pass parameters={{'m2': value}} to PDEFromSpec or build_pde_from_json."
                )
                raise ValueError(
                    msg
                )
            if "H" not in self._parameters:
                msg = (
                    f"Hubble parameter 'H' is required for time-dependent coefficient '{sym}'. "
                    f"Pass parameters={{'H': value}} to PDEFromSpec or build_pde_from_json."
                )
                raise ValueError(
                    msg
                )
            m2 = self._parameters["m2"]
            hubble_param = self._parameters["H"]
            # De Sitter mass term: m² Ω² = m² e^{2Ht}
            return -m2 * math.exp(2 * hubble_param * t)

        # No fallbacks - if we can't parse it, fail explicitly
        msg = (
            f"Cannot resolve time-dependent coefficient '{sym}' for term with "
            f"operator '{term.operator}'. Supported patterns: '-n*H' for Hubble friction, "
            f"'exp(...)' for exponential mass terms."
        )
        raise ValueError(
            msg
        )

    @staticmethod
    def _get_operator(  # noqa: C901, PLR0911, PLR0912, PLR0915
        operator_name: str, field: ScalarField, bc: BCDescriptor
    ) -> ScalarField:
        """Apply a named operator to a field.

        Parameters
        ----------
        operator_name : str
            Name of the operator ("laplacian", "identity", "gradient_x", etc.)
        field : ScalarField
            The field to operate on.
        bc : BCDescriptor
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
            return field.laplace(bc=bc)
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
        if operator_name == "cross_derivative_xy":
            # d/dx d/dy f = d/dx (d/dy f)
            # This is NOT a laplacian - it's a cross spatial derivative
            if field.grid.dim < 2:  # noqa: PLR2004
                msg = "cross_derivative_xy requires at least 2D grid"
                raise ValueError(msg)
            # First compute d/dy
            grad_y = field.gradient(bc=bc)[1]
            assert isinstance(grad_y, ScalarField)
            # Then compute d/dx of that
            grad_xy = grad_y.gradient(bc=bc)[0]
            assert isinstance(grad_xy, ScalarField)
            return grad_xy
        if operator_name == "cross_derivative_xz":
            # d/dx d/dz f = d/dx (d/dz f)
            # Cross spatial derivative for 3+1D
            if field.grid.dim < 3:  # noqa: PLR2004
                msg = "cross_derivative_xz requires at least 3D grid"
                raise ValueError(msg)
            grad_z = field.gradient(bc=bc)[2]
            assert isinstance(grad_z, ScalarField)
            grad_xz = grad_z.gradient(bc=bc)[0]
            assert isinstance(grad_xz, ScalarField)
            return grad_xz
        if operator_name == "cross_derivative_yz":
            # d/dy d/dz f = d/dy (d/dz f)
            # Cross spatial derivative for 3+1D
            if field.grid.dim < 3:  # noqa: PLR2004
                msg = "cross_derivative_yz requires at least 3D grid"
                raise ValueError(msg)
            grad_z = field.gradient(bc=bc)[2]
            assert isinstance(grad_z, ScalarField)
            grad_yz = grad_z.gradient(bc=bc)[1]
            assert isinstance(grad_yz, ScalarField)
            return grad_yz
        if operator_name == "laplacian_x":
            # Pure ∂²/∂x² - second derivative in x only
            # For anisotropic equations like Navier-Cauchy elasticity
            grad_x = field.gradient(bc=bc)[0]
            assert isinstance(grad_x, ScalarField)
            d2_x = grad_x.gradient(bc=bc)[0]
            assert isinstance(d2_x, ScalarField)
            return d2_x
        if operator_name == "laplacian_y":
            # Pure ∂²/∂y² - second derivative in y only
            if field.grid.dim < 2:  # noqa: PLR2004
                msg = "laplacian_y requires at least 2D grid"
                raise ValueError(msg)
            grad_y = field.gradient(bc=bc)[1]
            assert isinstance(grad_y, ScalarField)
            d2_y = grad_y.gradient(bc=bc)[1]
            assert isinstance(d2_y, ScalarField)
            return d2_y
        if operator_name == "laplacian_z":
            # Pure ∂²/∂z² - second derivative in z only
            if field.grid.dim < 3:  # noqa: PLR2004
                msg = "laplacian_z requires at least 3D grid"
                raise ValueError(msg)
            grad_z = field.gradient(bc=bc)[2]
            assert isinstance(grad_z, ScalarField)
            d2_z = grad_z.gradient(bc=bc)[2]
            assert isinstance(d2_z, ScalarField)
            return d2_z
        if operator_name == "biharmonic":
            # ∇⁴f = ∇²(∇²f) = laplacian(laplacian(f))
            # Used for plate theory and thin-shell mechanics
            lap = field.laplace(bc=bc)
            assert isinstance(lap, ScalarField)
            bilap = lap.laplace(bc=bc)
            assert isinstance(bilap, ScalarField)
            return bilap

        msg = f"Unknown operator: {operator_name}"
        raise ValueError(msg)

    def _get_field_from_state(
        self, state: FieldCollection, field_name: str
    ) -> ScalarField:
        """Get a field from state by name, supporting both field and momentum names.

        For mixed time-space derivatives like d_t d_x A, the Wolfram pipeline
        expresses these as gradients of momentum fields (e.g., gradient_x(pi_0)).
        This is valid because d_t A = pi, so d_x(d_t A) = d_x(pi).

        The state is organized as: [field_0, pi_0, field_1, pi_1, ...]
        - Regular fields (A_0, phi_0, phi0, etc.) are at even indices: state[2*i]
        - Momentum fields (pi_0, pi0, pi_1, pi1, etc.) are at odd indices: state[2*i + 1]

        Supports flexible field naming conventions:
        - Standard: A_0, phi_1
        - Compact: A0, phi1, pi0, pi1
        - Tensor: stress_xy_0
        - Simple: phi (index defaults to 0)

        Parameters
        ----------
        state : FieldCollection
            Current state with interleaved field/momentum pairs.
        field_name : str
            Name of the field to retrieve. Can be:
            - Regular field name like "A_0", "phi_0", "phi0"
            - Momentum field name like "pi_0", "pi_1", "pi0", "pi1"

        Returns
        -------
        ScalarField
            The requested field from the state.

        Raises
        ------
        ValueError
            If the field name is not recognized.
        """
        # Check if this is a momentum field reference (pi_0, pi0, etc.)
        # First check if it looks like a momentum field (starts with "pi")
        if field_name.startswith("pi"):
            momentum_idx = parse_momentum_field_name(field_name)
            if momentum_idx is not None:
                if not (0 <= momentum_idx < self.n_components):
                    msg = (
                        f"Momentum field index {momentum_idx} out of range. "
                        f"This system has {self.n_components} components "
                        f"(valid indices: 0 to {self.n_components - 1}). "
                        f"Field reference: '{field_name}'."
                    )
                    raise ValueError(msg)

                # Momentum is at odd indices: state[2*i + 1]
                momentum = state[2 * momentum_idx + 1]
                assert isinstance(momentum, ScalarField)
                return momentum
            # If it looks like momentum but couldn't be parsed, raise clear error
            msg = (
                f"Invalid momentum field format: '{field_name}'. "
                f"Expected 'pi_N' or 'piN' where N is a numeric index (e.g., 'pi_0', 'pi0')."
            )
            raise ValueError(msg)

        # Regular field lookup - first try direct name match
        target_idx = self._component_name_to_index.get(field_name)
        if target_idx is not None:
            # Field is at even indices: state[2*i]
            field = state[2 * target_idx]
            assert isinstance(field, ScalarField)
            return field

        msg = f"Unknown field name: {field_name}"
        raise ValueError(msg)

    def _compute_rhs_for_component(
        self,
        component_idx: int,
        state: FieldCollection,
        bc: BCDescriptor,
        t: float = 0.0,
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
        bc : BCDescriptor
            Boundary condition specification.
        t : float
            Current simulation time (for time-dependent coefficients in curved spacetime).

        Returns
        -------
        ScalarField
            The computed RHS for d/dt momentum_i.

        Raises
        ------
        ValueError
            If a field or operator cannot be resolved.
        """
        eq = self.spec.equations[component_idx]
        grid = state.grid

        # Start with zero field
        result = ScalarField(grid, data=0.0)

        # Sum all terms from the specification
        for term in eq.rhs_terms:
            # Find which field this term operates on
            # Supports both regular fields (A_0, phi_0) and momentum fields (pi_0, pi_1)
            # Momentum field support enables mixed time-space derivative handling
            target_field_name = term.field

            # Handle first_derivative_t operator specially
            # For d_t(field_i) = momentum_i, so first_derivative_t(field) returns momentum
            if term.operator == "first_derivative_t":
                # Get the field index for this term's target field
                target_idx = self._component_name_to_index.get(target_field_name)
                if target_idx is not None:
                    # Momentum is at odd indices: state[2*i + 1]
                    momentum = state[2 * target_idx + 1]
                    assert isinstance(momentum, ScalarField)
                    operated = momentum.copy()
                else:
                    msg = f"Unknown field for first_derivative_t: {target_field_name}"
                    raise ValueError(msg)
            else:
                # Standard operator handling
                target_field = self._get_field_from_state(state, target_field_name)
                operated = self._get_operator(term.operator, target_field, bc)

            # Resolve coefficient (use time-dependent resolution if needed)
            coefficient = self._resolve_coefficient_at_time(term, t)

            # Add coefficient * operated to result
            contribution = coefficient * operated
            result += contribution

        return result

    @override
    def evolution_rate(
        self,
        state: TState,
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
            Current time. Used for time-dependent coefficients in curved spacetime
            (e.g., Hubble friction in de Sitter expansion).

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
            # Pass current time for time-dependent coefficients (curved spacetime)
            d_momentum_dt = self._compute_rhs_for_component(i, state, bc, t)

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


def build_pde_from_json(
    json_path: Path | str,
    parameters: dict[str, float] | None = None,
) -> PDEFromSpec:
    """Build a PDE from a JSON equation specification file.

    This is the main entry point for the Lagrangian-to-PDE pipeline on the
    Python side. Given a JSON file exported from Mathematica/xAct, this
    function creates a py-pde compatible PDE class.

    Parameters
    ----------
    json_path : Path | str
        Path to the JSON file containing the equation specification.
    parameters : dict[str, float] | None
        Optional parameter values to override symbolic coefficients.
        Keys are symbolic names (e.g., "m2", "kappa"), values are numeric.
        Example: {"m2": 0.5, "kappa": 1.0}

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

    >>> # With custom parameter values:
    >>> pde = build_pde_from_json("examples/data/proca_1d.json", parameters={"m2": 2.0})
    """
    spec = load_equation_system(json_path)
    return PDEFromSpec(spec, parameters=parameters)


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
