"""``tg validate`` — Validate a JSON equation specification."""

from __future__ import annotations

import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from argparse import Namespace
    from pathlib import Path


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
        from torsion_gertsenshtein.symbolic.json_loader import load_equation_system

        spec = load_equation_system(json_path)
    except Exception as exc:  # noqa: BLE001
        errors.append(f"Failed to load equation system: {exc}")
        return None, errors

    return spec, errors


def _check_operators(spec: object) -> list[str]:
    """Check that all operators in the spec are recognized."""
    from torsion_gertsenshtein.symbolic.json_loader import (
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
    from torsion_gertsenshtein.symbolic.json_loader import EquationSystem

    if not isinstance(spec, EquationSystem):
        return []

    errors: list[str] = []
    # Build set of all valid field names (components + momenta)
    valid_names = set(spec.component_names)
    valid_names.update(f"pi_{name}" for name in spec.component_names)

    for eq in spec.equations:
        errors.extend(
            f"Unknown field reference '{term.field}' in equation for {eq.field_name}"
            for term in eq.rhs_terms
            if term.field not in valid_names
        )
    return errors


def _check_parameters(spec: object) -> list[str]:
    """Check for symbolic parameters that have no default values."""
    from torsion_gertsenshtein.cli._inspect import discover_parameters
    from torsion_gertsenshtein.symbolic.json_loader import EquationSystem

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


def validate_command(args: Namespace) -> int:
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

    # Report results
    if errors:
        for err in errors:
            print(f"ERROR: {err}", file=sys.stderr)
        return 1

    for warn in warnings:
        print(f"WARNING: {warn}", file=sys.stderr)

    print(f"OK: {json_path.name} is valid")
    return 0
