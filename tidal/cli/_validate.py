"""``tidal validate`` — Validate a JSON equation specification."""

from __future__ import annotations

import contextlib
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
                f"Parameter '{param}' (used in: {', '.join(fields)}) has no default value",
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


def _leading_sign(value: object) -> int:
    """Sign of a coefficient that may be numeric or a symbolic expression.

    Accepts ``None`` (treated as ``+1``, i.e. an implicit unit coefficient),
    a float, or a string such as ``"-1 + 2*B0^2*rho"``.  Only the leading
    sign is inspected -- enough to compare orientations, and robust against
    expressions we cannot evaluate without parameter values.
    """
    if value is None:
        return 1
    if isinstance(value, (int, float)):
        return 1 if value > 0 else (-1 if value < 0 else 0)
    text = str(value).strip().lstrip("(").strip()
    return -1 if text.startswith("-") else 1


# NOTE(GH #401): duplicates logic that belongs on the model.  Refactor onto
# the shared semantic accessors on ComponentEquation once they land.
def _effective_self_sign(equation: object, operator_prefix: str) -> int | None:
    """Orientation of a field's own ``operator_prefix`` terms in its equation.

    Returns ``sign(sum of matching self-terms) * sign(kinetic coefficient)``,
    or ``None`` when the field has no such term.

    Three details matter, each of which produced a wrong answer during the
    GH #397 investigation when omitted:

    * **Normalise by the kinetic coefficient.** Euler-Heisenberg photons carry
      ``lap = -1`` with ``kin = -1 + 2*B0^2*rho``, i.e. orientation ``+1``.
      Comparing raw coefficients marks them inconsistent when they are correct.
    * **Sum every matching term.** EH components have a base term *and* a
      ``B0^2*rho`` correction; taking a single term gives the wrong sign.
    * **Symbolic coefficients count.** Torsion self-laplacians appear as
      ``-xi``; skipping non-numeric coefficients misses them entirely.
    """
    total = 0
    found = False
    for term in equation.rhs_terms:  # type: ignore[attr-defined]
        if term.field != equation.field_name:  # type: ignore[attr-defined]
            continue
        if not term.operator.startswith(operator_prefix):
            continue
        found = True
        raw = term.coefficient_symbolic
        total += _leading_sign(raw if raw is not None else term.coefficient)
    if not found:
        return None
    kin = _leading_sign(equation.kinetic_coefficient_symbolic)  # type: ignore[attr-defined]
    return (1 if total > 0 else (-1 if total < 0 else 0)) * kin


def _check_volume_element_consistency(spec: object) -> list[str]:
    """Require Christoffel gradient terms when the volume element is non-constant (GH #394).

    ``ComponentEulerLagrange`` varies a Lagrangian *density*.  When the metric
    has a non-constant ``sqrt|g|`` the equations must carry the first-derivative
    terms it generates.  If ``canonical.volume_element`` is non-constant but no
    equation has a spatial gradient term, the measure was dropped during
    derivation and the energy integral no longer matches the operator that
    evolves the fields -- exactly the #393/#394 failure.
    """
    from tidal.symbolic.json_loader import EquationSystem

    if not isinstance(spec, EquationSystem):
        return []
    canonical = spec.canonical
    if canonical is None or canonical.volume_element is None:
        return []

    ve = str(canonical.volume_element)
    coords = set(spec.coordinates)
    if not any(f"{c}[]" in ve for c in coords):
        return []  # constant measure: factors out, nothing to check

    has_gradient = any(
        term.operator.startswith("gradient")
        for eq in spec.equations
        for term in eq.rhs_terms
    )
    if has_gradient:
        return []
    return [
        f"Non-constant volume_element '{ve}' but no spatial gradient term in any "
        "equation. The sqrt|g| measure was likely dropped during Euler-Lagrange "
        "variation, so the energy integral will not match the evolution operator "
        "(see GH #394).",
    ]


def _check_temporal_component_sign(spec: object) -> list[str]:
    """Check that a vector field's temporal and spatial components agree in sign (GH #397).

    ``Box A_mu = 0`` acts identically on every component, so a rank-1 field's
    temporal component must carry the same effective self-laplacian orientation
    as its spatial siblings.  Pre-GH-#381 exports left the temporal component
    un-normalised (``lhsCoeff = -1`` never divided through), producing
    ``d2_t a_0 = -laplacian a_0`` -- a temporal-only tachyon.

    Deliberately narrow, because each restriction prevents a real false positive:

    * **Bare ``d2_t`` LHS only.** ``h_5`` carries ``laplacian = -kappa^(-2)``
      matched by its stored kinetic coefficient and is correctly normalised;
      constraint equations have a conventional overall sign.
    * **Rank-1 families only.** Torsion is rank 3 and its irreducible components
      legitimately carry different normalisations -- that non-uniformity is
      physical and must not be reported.
    """
    from tidal.symbolic.json_loader import EquationSystem

    if not isinstance(spec, EquationSystem):
        return []

    tensor_meta = spec.metadata.get("tensor_metadata")
    has_meta = isinstance(tensor_meta, dict) and bool(tensor_meta)

    by_head: dict[str, dict[int, object]] = {}
    for eq in spec.equations:
        if eq.time_derivative_order != 2 or eq.kinetic_coefficient_symbolic is not None:
            continue

        if has_meta:
            meta = tensor_meta.get(eq.field_name)  # type: ignore[union-attr]
            if not isinstance(meta, dict) or meta.get("tensor_rank") != 1:
                continue
            indices = meta.get("tensor_indices") or []
            if len(indices) != 1:
                continue
            head, component = str(meta.get("tensor_head")), int(indices[0])
        else:
            # Pre-78374c1 exports carry no tensor metadata, and those are the
            # OLDEST specs -- precisely the ones most likely to predate the
            # GH #381 fix.  Skipping them would blind the guard to its main
            # target, so fall back to the "<head>_<index>" naming convention.
            #
            # Restricted to families whose indices are exactly 0..dim-1, which
            # identifies a rank-1 field: a rank-2 or rank-3 field in the same
            # notation spans far more components (h_0..h_9, t_0..t_23), so it
            # cannot be mistaken for a vector.
            name = eq.field_name
            if "_" not in name:
                continue
            head, _, suffix = name.rpartition("_")
            if not suffix.isdigit():
                continue
            component = int(suffix)

        by_head.setdefault(head, {})[component] = eq

    if not has_meta:
        # Discriminate rank by index range against the UNREDUCED dimension: a
        # rank-1 field spans 0..D-1, while rank-2 spans 0..D^2-1 and rank-3
        # far more.  Using the max index rather than the full set keeps this
        # correct when some components were filtered out above.
        #
        # `spacetime.dimension` is the *reduced* dimension for plane-wave
        # specs (2), while the field indices still run over the original 4 --
        # so prefer metadata.reduction.original_dimension where present.
        reduction = spec.metadata.get("reduction")
        dim = spec.dimension
        if isinstance(reduction, dict):
            dim = int(reduction.get("original_dimension", dim) or dim)
        by_head = {h: c for h, c in by_head.items() if c and max(c) < dim}

    errors: list[str] = []
    for head, comps in sorted(by_head.items()):
        temporal = comps.get(0)
        if temporal is None or len(comps) < 2:
            continue
        t_sign = _effective_self_sign(temporal, "laplacian")
        if not t_sign:
            continue
        spatial = {
            i: _effective_self_sign(eq, "laplacian")
            for i, eq in comps.items()
            if i != 0
        }
        opposed = sorted(i for i, s in spatial.items() if s and s * t_sign < 0)
        if opposed and len(opposed) == len([s for s in spatial.values() if s]):
            errors.append(
                f"Field '{head}' temporal component "
                f"'{temporal.field_name}' has self-laplacian sign opposite to its "  # type: ignore[attr-defined]
                f"spatial components {[comps[i].field_name for i in opposed]}. "  # type: ignore[attr-defined]
                "Box A_mu = 0 acts identically on all components, so this is a "
                "temporal-only tachyon -- the spec predates the GH #381 LHS "
                "normalisation fix and needs re-deriving (see GH #397).",
            )
    return errors


def _check_perturbative_consistency(spec: object) -> tuple[list[str], list[str]]:
    """Check perturbative-reduction scope for theories with [perturbation].

    Returns ``(errors, warnings)``.  Active when the spec has non-empty
    ``metadata.perturbation.small_parameters``:

    1. Term-vs-coefficient consistency: each ``order_in_eps > 0`` term has
       a ``coefficient_symbolic`` that mentions a declared small parameter.
    2. No HamiltonianTerm operator matches the LPS-invariant violation
       pattern (``mixed_T_*`` with T >= 2) on theories where LPS has run
       successfully.
    3. Constraint-promotion detection: warn (not error) when the spec has
       both ``[perturbation]`` AND any equation with
       ``time_derivative_order > 2`` whose order-0 kinetic coefficient
       would vanish (signature of constraint promotion).  Names the
       affected fields and links to issue #321 / the TeX write-up.
    """
    from tidal.symbolic.json_loader import EquationSystem

    if not isinstance(spec, EquationSystem):
        return [], []

    errors: list[str] = []
    warnings: list[str] = []

    pert = (spec.metadata or {}).get("perturbation", {})
    small_params = list(pert.get("small_parameters", []) or [])
    if not small_params:
        return [], []

    # Check 1: term-vs-coefficient consistency on Hamiltonian terms
    if spec.canonical is not None:
        for idx, term in enumerate(spec.canonical.hamiltonian_terms):
            if term.order_in_eps and term.order_in_eps > 0:
                coeff_sym = term.coefficient_symbolic or ""
                if not any(p in coeff_sym for p in small_params):
                    warnings.append(
                        f"hamiltonian_terms[{idx}] has order_in_eps="
                        f"{term.order_in_eps} but coefficient_symbolic "
                        f"({coeff_sym!r}) does not mention any declared "
                        f"small_parameter {small_params}.  Re-derive the "
                        f"theory if the JSON predates v0.36."
                    )

    # Check 2: LPS invariant — no mixed_T_* with T >= 2
    if spec.canonical is not None:
        import re as _re

        bad_op = _re.compile(r"^mixed_(\d+)_")
        for idx, term in enumerate(spec.canonical.hamiltonian_terms):
            for label, factor in (("a", term.factor_a), ("b", term.factor_b)):
                m = bad_op.match(factor.operator or "")
                if m and int(m.group(1)) >= 2:
                    warnings.append(
                        f"hamiltonian_terms[{idx}].factor_{label}.operator="
                        f"{factor.operator!r} encodes time-derivative order "
                        f"{m.group(1)}.  This is the irreducible LPS residue "
                        f"of the constraint-promotion case (the small parameter "
                        f"promotes an algebraic constraint to a higher-derivative "
                        f"dynamical field).  See "
                        f"docs/tex/perturbative_reduction_constraint_barrier.tex "
                        f"and issue #321 for the architectural barrier."
                    )

    # Check 3: constraint-promotion detection on equations
    constraint_promoted_fields: list[str] = []
    for eq in spec.equations:
        if eq.time_derivative_order <= 2:
            continue
        # Field is dynamical at full theory but its kinetic coefficient
        # might vanish at small_param=0 (constraint promotion).  Heuristic:
        # if kinetic_coefficient_symbolic contains a declared small_parameter,
        # treat the field as constraint-promoted.
        kc = eq.kinetic_coefficient_symbolic or ""
        if any(p in kc for p in small_params):
            constraint_promoted_fields.append(eq.field_name)
    if constraint_promoted_fields:
        warnings.append(
            f"Constraint-promotion case detected for fields "
            f"{constraint_promoted_fields}.  Lagrangian Perturbative "
            f"Substitution (LPS, v6 Phase 2) cannot reduce these theories' "
            f"Hamiltonians at the standard substitution level.  The JSON's "
            f"existing hamiltonian_terms are the pre-LPS form (containing "
            f"mixed_T_* with T>=2 residues).  See issue #321 and "
            f"docs/tex/perturbative_reduction_constraint_barrier.tex for "
            f"the architectural barrier and three open research directions."
        )

    return errors, warnings


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


def _run_stability_checks(
    spec: EquationSystem,
    raw_params: list[str],
) -> tuple[list[str], list[str]]:
    """Run tachyon check, returning (errors, notes).

    Tachyon detection (negative mass-matrix eigenvalue) is reported as an
    **error** — the system will have exponentially growing modes at runtime.

    Ghost detection is not implemented here: determining whether a theory has
    ghost modes from ``hamiltonian_terms`` alone is unreliable because, in
    linearized GR and other gauge theories, the naive Hamiltonian kinetic
    coefficients are negative even in ghost-free theories (gauge-structure
    artifact).  Use a dedicated tool (e.g. xAct/Mathematica) to verify
    ghost-freeness for a specific theory.

    Parameters
    ----------
    spec : EquationSystem
    raw_params : list[str]
        Raw --param KEY=VAL strings from CLI.

    Returns
    -------
    errors : list[str]
        Tachyon instability errors (fatal — exponentially growing modes).
    notes : list[str]
        Informational notes from the tachyon check (e.g. asymmetric matrix).
    """
    try:
        params = _parse_validate_params(raw_params, spec)
    except ValueError as exc:
        return [str(exc)], []

    tachyon_errors, notes = _check_tachyons(spec, params)
    return tachyon_errors, notes


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

    from tidal.cli._console import error as _error
    from tidal.cli._console import success as _success
    from tidal.cli._console import warn as _warn

    json_path = Path(args.json_path)

    errors: list[str] = []
    warnings: list[str] = []

    # Check 1: File exists
    errors.extend(_check_file_exists(json_path))
    if errors:
        from tidal.cli._console import error_with_hint

        for err in errors:
            error_with_hint(
                err,
                hints=["Use `tidal list` to find specs, or `tidal derive` to generate"],
            )
        return 1

    # Check 2: Parse JSON + load EquationSystem
    spec, parse_errors = _check_json_parse(json_path)
    errors.extend(parse_errors)
    if errors:
        from tidal.cli._console import error_with_hint

        for err in errors:
            if err.startswith("Invalid JSON"):
                error_with_hint(
                    err,
                    hints=[
                        "Check JSON syntax. Validate with: `python -m json.tool <file>`",
                    ],
                )
            elif err.startswith("Failed to load"):
                error_with_hint(
                    err,
                    hints=["Ensure JSON was generated by `tidal derive`"],
                )
            else:
                _error(err)
        return 1

    if spec is None:
        return 1

    # Check 3: Operators
    errors.extend(_check_operators(spec))

    # Check 4: Field references
    errors.extend(_check_field_references(spec))

    # Check 5: Parameters (warnings, not errors)
    warnings.extend(_check_parameters(spec))

    # Check 5c: Derivation-integrity guards.  Both catch defects that are
    # silent at derivation time and only surface as wrong physics later --
    # a dropped sqrt|g| measure (#394) and an un-normalised temporal
    # component (#397).
    errors.extend(_check_volume_element_consistency(spec))
    errors.extend(_check_temporal_component_sign(spec))

    # Check 5b: Perturbative-reduction scope (warnings).  Active when the
    # spec has [perturbation]; flags constraint-promotion theories where
    # LPS does not apply, with pointers to issue #321 and the TeX
    # write-up.
    pert_errors, pert_warnings = _check_perturbative_consistency(spec)
    errors.extend(pert_errors)
    warnings.extend(pert_warnings)

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
        from tidal.cli._console import error_with_hint

        for err in errors:
            if "Unknown operator" in err or "Unknown field reference" in err:
                error_with_hint(
                    err,
                    hints=[
                        "Run `tidal inspect <json>` to check field and operator names",
                    ],
                )
            else:
                _error(err)
        return 1

    for w in warnings:
        _warn(w)

    suffix = " [stable]" if getattr(args, "stability", False) else ""
    _success(f"{json_path.name} is valid{suffix}")
    return 0
