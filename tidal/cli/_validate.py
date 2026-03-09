"""``tidal validate`` — Validate a JSON equation specification."""

from __future__ import annotations

import contextlib
import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from argparse import Namespace
    from pathlib import Path

    from tidal.symbolic.json_loader import EquationSystem

# Tolerance matching check_pointwise_mass_stability in validation.py
_STABILITY_TOLERANCE: float = 1e-10


def _check_file_exists(json_path: Path) -> list[str]:
    """Check that the JSON file exists and is readable."""
    if not json_path.exists():
        return [f"File not found: {json_path}"]
    if not json_path.is_file():
        return [f"Not a file: {json_path}"]
    return []


def _check_json_parse(json_path: Path) -> tuple[object | None, list[str]]:
    """Attempt to parse the JSON and load as EquationSystem."""
    import json

    errors: list[str] = []

    try:
        with json_path.open(encoding="utf-8") as f:
            json.load(f)
    except json.JSONDecodeError as exc:
        errors.append(f"Invalid JSON: {exc}")
        return None, errors

    try:
        from tidal.symbolic.json_loader import load_equation_system

        spec = load_equation_system(json_path)
    except Exception as exc:  # noqa: BLE001
        errors.append(f"Failed to load equation system: {exc}")
        return None, errors

    return spec, errors


def _check_operators(spec: object) -> list[str]:
    """Check that all operators in the spec are recognized."""
    from tidal.symbolic.json_loader import (
        EquationSystem,
        is_known_operator,
    )

    if not isinstance(spec, EquationSystem):
        return []

    errors: list[str] = []
    for eq in spec.equations:
        errors.extend(
            f"Unknown operator '{term.operator}' in equation for {eq.field_name}"
            for term in eq.rhs_terms
            if not is_known_operator(term.operator)
        )
    return errors


def _check_field_references(spec: object) -> list[str]:
    """Check that all field references in terms point to existing fields."""
    from tidal.symbolic.json_loader import EquationSystem

    if not isinstance(spec, EquationSystem):
        return []

    errors: list[str] = []
    # Build set of all valid field names (components + velocities)
    valid_names = set(spec.component_names)
    valid_names.update(f"v_{name}" for name in spec.component_names)

    for eq in spec.equations:
        errors.extend(
            f"Unknown field reference '{term.field}' in equation for {eq.field_name}"
            for term in eq.rhs_terms
            if term.field not in valid_names
        )
    return errors


def _check_parameters(spec: object) -> list[str]:
    """Check for symbolic parameters that have no default values."""
    from tidal.cli._inspect import discover_parameters
    from tidal.symbolic.json_loader import EquationSystem

    if not isinstance(spec, EquationSystem):
        return []

    warnings: list[str] = []
    param_map = discover_parameters(spec)
    defaults = spec.metadata.get("parameters", {})

    for param, fields in sorted(param_map.items()):
        if isinstance(defaults, dict) and param not in defaults:
            warnings.append(
                f"Parameter '{param}' (used in: {', '.join(fields)}) has no default value"
            )
    return warnings


def _parse_validate_params(raw: list[str], spec: EquationSystem) -> dict[str, float]:
    """Parse --param KEY=VAL arguments, merging metadata defaults.

    Raises
    ------
    ValueError
        If parameter format is invalid or value is non-numeric.
    """
    params: dict[str, float] = {}

    # Start with metadata defaults
    meta_params = spec.metadata.get("parameters", {})
    if isinstance(meta_params, dict):
        for key, val in meta_params.items():  # type: ignore[union-attr]
            with contextlib.suppress(ValueError, TypeError):
                params[str(key)] = float(val)  # type: ignore[arg-type]

    # Override with CLI params
    for item in raw:
        if "=" not in item:
            msg = f"Invalid --param format: '{item}'. Expected KEY=VALUE (e.g. --param m2=1.0)"
            raise ValueError(msg)
        key, val_str = item.split("=", 1)
        key = key.strip()
        try:
            params[key] = float(val_str.strip())
        except ValueError:
            msg = f"Invalid parameter value: '{val_str}' for key '{key}'. Must be a number."
            raise ValueError(msg) from None

    return params


def _check_tachyons(
    spec: EquationSystem,
    params: dict[str, float],
) -> tuple[list[str], list[str]]:
    """Check mass matrix stability (tachyon detection).

    Builds the potential matrix M[i,j] from identity-operator RHS terms for
    dynamical fields and verifies all eigenvalues are non-negative. A negative
    eigenvalue indicates an exponentially growing (tachyonic) mode.

    Uses a minimal 1-point grid; for constant-coefficient systems (the common
    case) this is exact. For position-dependent coefficients it evaluates at
    the origin, which may miss instabilities elsewhere.

    Parameters
    ----------
    spec : EquationSystem
    params : dict[str, float]
        Parameter values (merged from metadata defaults + CLI overrides).

    Returns
    -------
    errors : list[str]
        Tachyon instability errors.
    notes : list[str]
        Informational notes (e.g. asymmetric matrix detected).
    """
    from tidal.solver.coefficients import CoefficientEvaluator
    from tidal.solver.grid import GridInfo
    from tidal.solver.validation import check_pointwise_mass_stability

    spatial_dim = spec.spatial_dimension
    bounds = tuple((0.0, 1.0) for _ in range(spatial_dim))
    shape = tuple(2 for _ in range(spatial_dim))  # minimum 2 cells required by GridInfo
    periodic = tuple(True for _ in range(spatial_dim))
    grid = GridInfo(bounds=bounds, shape=shape, periodic=periodic)

    coeff_eval = CoefficientEvaluator(spec, grid, params)
    result = check_pointwise_mass_stability(coeff_eval, spec, grid)
    return result.errors, result.notes


def _check_ghosts(spec: EquationSystem) -> list[str]:
    """Check for ghost modes (wrong-sign kinetic terms) from Hamiltonian terms.

    A ghost is a field with a negative coefficient in front of ``(d_t phi)^2``.
    Detects diagonal kinetic terms (the common case) by scanning
    ``canonical.hamiltonian_terms`` for ``"time_derivative"`` x ``"time_derivative"``
    bilinears with negative coefficients.

    Off-diagonal kinetic mixing (e.g. ``d_t phi * d_t chi`` terms) would require the
    full kinetic matrix K[i,j] to be exported from Wolfram — this is deferred.

    Returns
    -------
    list[str]
        Ghost mode error messages (one per ghost field found).
    """
    if spec.canonical is None:
        return []

    ghosts: list[str] = []
    for term in spec.canonical.hamiltonian_terms:
        if (
            term.factor_a.operator == "time_derivative"
            and term.factor_b.operator == "time_derivative"
            and term.factor_a.field == term.factor_b.field
            and term.coefficient < -_STABILITY_TOLERANCE
        ):
            sym = f" ({term.coefficient_symbolic!r})" if term.coefficient_symbolic else ""
            ghosts.append(
                f"Ghost mode: field '{term.factor_a.field}' has negative kinetic "
                f"coefficient {term.coefficient:.6g}{sym}. "
                f"Wrong-sign kinetic term — the theory has a ghost degree of freedom."
            )
    return ghosts


def _run_stability_checks(
    spec: EquationSystem,
    raw_params: list[str],
) -> tuple[list[str], list[str]]:
    """Run tachyon and ghost checks, returning (errors, warnings).

    Tachyon detection (negative mass-matrix eigenvalue) is reported as an
    **error** — the system will have exponentially growing modes at runtime.

    Ghost detection (wrong-sign kinetic term) is reported as a **warning**
    because the check is heuristic: in GR and tensor theories with Lorentzian
    signature, the Hamiltonian kinetic coefficients are naturally negative
    (e.g. ``H = -½(∂_t h)²`` for metric perturbations, from ``L = -½(∂_t h)²``
    and ``H = p·v - L`` with ``p = -v``).  False positives are common for fields
    derived from linearized GR Lagrangians.  True ghost instabilities require
    a full analysis of the physical kinetic matrix K[i,j] in the TT sector.

    Parameters
    ----------
    spec : EquationSystem
    raw_params : list[str]
        Raw --param KEY=VAL strings from CLI.

    Returns
    -------
    errors : list[str]
        Tachyon instability errors (fatal — exponentially growing modes).
    warnings : list[str]
        Ghost candidate warnings (informational — heuristic, GR false-positives
        common).
    """
    try:
        params = _parse_validate_params(raw_params, spec)
    except ValueError as exc:
        return [str(exc)], []

    tachyon_errors, notes = _check_tachyons(spec, params)
    ghost_warnings = _check_ghosts(spec)
    return tachyon_errors, notes + ghost_warnings


def validate_command(args: Namespace) -> int:  # noqa: C901
    """Execute the validate command.

    Parameters
    ----------
    args : Namespace
        Parsed CLI arguments.

    Returns
    -------
    int
        Exit code (0 = valid, 1 = errors found).
    """
    from pathlib import Path

    json_path = Path(args.json_path)

    errors: list[str] = []
    warnings: list[str] = []

    # Check 1: File exists
    errors.extend(_check_file_exists(json_path))
    if errors:
        for err in errors:
            print(f"ERROR: {err}", file=sys.stderr)
        return 1

    # Check 2: Parse JSON + load EquationSystem
    spec, parse_errors = _check_json_parse(json_path)
    errors.extend(parse_errors)
    if errors:
        for err in errors:
            print(f"ERROR: {err}", file=sys.stderr)
        return 1

    if spec is None:
        return 1

    # Check 3: Operators
    errors.extend(_check_operators(spec))

    # Check 4: Field references
    errors.extend(_check_field_references(spec))

    # Check 5: Parameters (warnings, not errors)
    warnings.extend(_check_parameters(spec))

    # Check 6: Stability — tachyon and ghost mode detection
    if getattr(args, "stability", False):
        from tidal.symbolic.json_loader import EquationSystem

        if isinstance(spec, EquationSystem):
            raw_params: list[str] = getattr(args, "param", []) or []
            stab_errors, stab_notes = _run_stability_checks(spec, raw_params)
            errors.extend(stab_errors)
            warnings.extend(stab_notes)

    # Report results
    if errors:
        for err in errors:
            print(f"ERROR: {err}", file=sys.stderr)
        return 1

    for warn in warnings:
        print(f"WARNING: {warn}", file=sys.stderr)

    suffix = " [stable]" if getattr(args, "stability", False) else ""
    print(f"OK: {json_path.name} is valid{suffix}")
    return 0
