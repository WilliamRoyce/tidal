"""Query surface for ``tidal inspect`` — read parts of a spec and check them.

Backs ``--equation``, ``--coefficient``, ``--families`` and ``--diff``, all of
which render :mod:`tidal.symbolic.spec_query` results.  No physics reasoning
lives here; this module only presents what the accessors decided, always
including *which tactic* decided it so a reader can audit a verdict without
re-deriving it.

Output conventions, chosen so the same commands serve a person at a terminal
and a program parsing the result:

* ``--json`` emits an ``{"items": [...]}`` envelope, never a bare array, so
  metadata can be added later without breaking consumers.
* Data goes to stdout and diagnostics to stderr, so a ``| jq`` pipeline never
  breaks on a warning.
* ``--diff`` exits ``1`` when it finds real changes.  That is a *declared
  outcome* in the ``diff(1)`` sense, not an error, which makes the command
  usable directly in a shell conditional.  Genuine failures use ``2``.
* ``--fields`` projects the JSON output.  These specs are large — a torsion
  spec has 32 components — so a full dump is thousands of tokens when the
  question needs one verdict.
"""

from __future__ import annotations

import json
import re
import sys
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from argparse import Namespace
    from collections.abc import Iterable, Mapping, Sequence

    from tidal.symbolic.json_loader import EquationSystem
    from tidal.symbolic.sign_algebra import SignResult

# Exit codes. 0 success, 1 "differences found" (an outcome, not a failure),
# 2 usage/lookup error. Kept distinct so a caller can branch on them.
EXIT_OK = 0
EXIT_DIFFERENCES = 1
EXIT_ERROR = 2

# "h_5:laplacian_x(a_1)"  ->  equation h_5, operator laplacian_x, field a_1
_COEFFICIENT_RE = re.compile(
    r"^\s*(?P<equation>\w+)\s*:\s*(?P<operator>\w+)\s*\(\s*(?P<field>\w+)\s*\)\s*$",
)


def parse_coefficient_key(spec_text: str) -> tuple[str, str, str]:
    """Parse a ``EQUATION:operator(field)`` selector.

    Parameters
    ----------
    spec_text : str
        Selector such as ``"h_5:laplacian_x(a_1)"``.

    Returns
    -------
    tuple[str, str, str]
        ``(equation_field, operator, field)``.

    Raises
    ------
    ValueError
        If the selector is malformed.
    """
    match = _COEFFICIENT_RE.match(spec_text)
    if match is None:
        msg = (
            f"Invalid --coefficient selector: {spec_text!r}. "
            "Expected EQUATION:operator(field), e.g. 'h_5:laplacian_x(a_1)'"
        )
        raise ValueError(msg)
    return match["equation"], match["operator"], match["field"]


def _split_names(raw: str | None) -> tuple[str, ...]:
    """Split a comma-separated option value into names."""
    if not raw:
        return ()
    return tuple(part.strip() for part in raw.split(",") if part.strip())


def _project(payload: dict[str, Any], fields: tuple[str, ...]) -> dict[str, Any]:
    """Restrict a record to *fields*, leaving it untouched when none are given."""
    if not fields:
        return payload
    return {k: v for k, v in payload.items() if k in fields}


def _sign_payload(result: SignResult) -> dict[str, Any]:
    """Render a :class:`~tidal.symbolic.sign_algebra.SignResult` as a record."""
    return {
        "sign": result.sign.value,
        "symbol": result.sign.symbol,
        "definite": result.is_definite,
        "tactic": result.tactic,
        "assumptions": list(result.assumptions),
        "free_names": list(result.free_names),
        "value": str(result.value) if result.value is not None else None,
        "numeric": result.numeric,
    }


# --- --equation ---


def show_equation(  # noqa: PLR0913
    spec: EquationSystem,
    field: str,
    *,
    as_json: bool,
    fields: tuple[str, ...],
    parameters: Mapping[str, float] | None,
    assume_positive: Iterable[str] | None,
    assume_nonzero: Iterable[str] | None,
) -> int:
    """Print one or more equations with their effective coefficients.

    Parameters
    ----------
    field : str
        A component name, a comma-separated list, or ``all``.

    Returns
    -------
    int
        ``EXIT_OK``, or ``EXIT_ERROR`` if any component is unknown.
    """
    requested = _resolve_fields(spec, field)
    if requested is None:
        return _fail(
            "unknown_component",
            f"unknown component in {field!r}",
            f"Use --families to list components, or --equation all. "
            f"Components: {', '.join(spec.component_names)}",
            as_json=as_json,
        )

    records: list[dict[str, Any]] = []
    text: list[str] = []
    for name in requested:
        text.extend(
            _equation_lines(
                spec,
                name,
                records,
                as_json=as_json,
                parameters=parameters,
                assume_positive=assume_positive,
                assume_nonzero=assume_nonzero,
            )
        )

    if as_json:
        print(json.dumps({"items": [_project(r, fields) for r in records]}, indent=2))
        return EXIT_OK
    for line in text:
        print(line)
    return EXIT_OK


def _resolve_fields(spec: EquationSystem, requested: str) -> list[str] | None:
    """Expand an ``--equation`` value into component names.

    Accepts a single name, a comma-separated list, or ``all``. Returns ``None``
    if any name is unknown, so the caller reports one clear error rather than
    printing part of the answer.
    """
    if requested.strip() == "all":
        return list(spec.component_names)
    names = [n.strip() for n in requested.split(",") if n.strip()]
    if any(n not in spec.equation_map for n in names):
        return None
    return names


def _equation_lines(  # noqa: PLR0913
    spec: EquationSystem,
    field: str,
    records: list[dict[str, Any]],
    *,
    as_json: bool,  # noqa: ARG001
    parameters: Mapping[str, float] | None,
    assume_positive: Iterable[str] | None,
    assume_nonzero: Iterable[str] | None,
) -> list[str]:
    """Render one component, appending its records for the JSON form."""
    from tidal.symbolic.spec_query import effective_coefficient

    index = spec.equation_map[field]
    equation = spec.equations[index]
    keys = sorted({(t.operator, t.field) for t in equation.rhs_terms})

    if equation.time_derivative_order:
        lhs = f"d{equation.time_derivative_order}_t({field})"
    elif field in spec.second_order_sector.promoted:
        # Algebraic LHS carrying inter-constraint time derivatives: the
        # row belongs to the second-order sector (GH #457).
        lhs = f"{field} (algebraic LHS — promoted to second-order sector, GH #457)"
    else:
        lhs = f"{field} (constraint)"
    text = [
        f"{lhs}   kinetic = {equation.kinetic_coefficient_symbolic or '1'}",
        "",
        "  effective coefficients (all matching terms summed, kinetic divided out):",
    ]
    for operator, target in keys:
        eff = effective_coefficient(equation, target, operator)
        result = eff.sign(
            assume_positive=assume_positive,
            assume_nonzero=assume_nonzero,
            parameters=parameters,
        )
        note = f"  [{eff.term_count} terms summed]" if eff.term_count > 1 else ""
        text.extend(
            (
                f"    {operator}({target}): {eff.expression}{note}",
                f"        sign {result.describe()}",
            )
        )
        records.append(
            {
                "equation": field,
                "operator": operator,
                "field": target,
                "expression": eff.expression,
                "numerator": eff.numerator,
                "kinetic": eff.kinetic,
                "term_count": eff.term_count,
                **_sign_payload(result),
            },
        )

    return text


# --- --coefficient ---


def show_coefficient(  # noqa: PLR0913, C901
    spec: EquationSystem,
    selector: str,
    *,
    as_json: bool,
    fields: tuple[str, ...],
    parameters: Mapping[str, float] | None,
    assume_positive: Iterable[str] | None,
    assume_nonzero: Iterable[str] | None,
) -> int:
    """Print one coefficient and every place it is recorded.

    Returns
    -------
    int
        ``EXIT_OK``, or ``EXIT_ERROR`` on a bad selector or unknown component.
    """
    from tidal.symbolic.spec_query import coefficient_provenance

    try:
        equation_field, operator, field = parse_coefficient_key(selector)
    except ValueError as exc:
        return _fail(
            "invalid_selector",
            str(exc),
            "Example: --coefficient 'h_5:laplacian_x(h_5)'. "
            "Use --equation FIELD to list the operator(field) keys available.",
            as_json=as_json,
        )

    try:
        prov = coefficient_provenance(spec, equation_field, field, operator)
    except KeyError:
        return _fail(
            "unknown_component",
            f"unknown component: {equation_field!r}",
            f"Use --families to see components grouped by tensor family. "
            f"Components: {', '.join(spec.component_names)}",
            as_json=as_json,
        )

    eff = prov.effective
    result = eff.sign(
        assume_positive=assume_positive,
        assume_nonzero=assume_nonzero,
        parameters=parameters,
    )

    record: dict[str, Any] = {
        "equation": equation_field,
        "operator": operator,
        "field": field,
        "expression": eff.expression,
        "numerator": eff.numerator,
        "kinetic": eff.kinetic,
        "term_count": eff.term_count,
        "terms": [
            {
                "coefficient": t.coefficient,
                "coefficient_symbolic": t.coefficient_symbolic,
                "order_in_eps": t.order_in_eps,
            }
            for t in eff.terms
        ],
        "order_spread": {str(k): list(v) for k, v in prov.order_spread.items()},
        "matrix_entry": prov.matrix_entry,
        "hamiltonian_terms": [
            {
                "coefficient_symbolic": t.coefficient_symbolic,
                "factor_a": f"{t.factor_a.operator}({t.factor_a.field})",
                "factor_b": f"{t.factor_b.operator}({t.factor_b.field})",
                "term_class": t.term_class,
            }
            for t in prov.hamiltonian_terms
        ],
        "checks": [
            {"name": c.name, "status": c.status, "detail": c.detail}
            for c in prov.checks
        ],
        **_sign_payload(result),
    }

    if as_json:
        print(json.dumps({"items": [_project(record, fields)]}, indent=2))
        return EXIT_OK

    print(f"{operator}({field})  in equation {equation_field}")
    print()
    print("  THIS COEFFICIENT")
    print(f"    effective     : {eff.expression}")
    print(f"    sign          : {result.describe()}")
    print(f"    from {eff.term_count} RHS term(s):")
    for term in eff.terms:
        print(
            f"      {term.coefficient_symbolic or term.coefficient}"
            f"   (order_in_eps={term.order_in_eps})",
        )
    print(f"    LHS kinetic   : {eff.kinetic or '1 (bare LHS)'}")

    if prov.matrix_entry is not None:
        print()
        print("  REDUNDANT RE-ENCODING (same quantity, different convention)")
        print(f"    mass/coupling matrix : {prov.matrix_entry!r}")
        print("      note: stored verbatim, un-normalized by the kinetic coefficient;")
        print("      the numeric matrix negates it (see GH #404)")

    if prov.hamiltonian_terms:
        print()
        print("  RELATED BUT DISTINCT (different formalism, its own factor convention)")
        for term in prov.hamiltonian_terms:
            print(
                f"    hamiltonian {term.coefficient_symbolic}  "
                f"{term.factor_a.operator}({term.factor_a.field})"
                f" * {term.factor_b.operator}({term.factor_b.field})",
            )
        print("      note: NOT the same number as the EOM coefficient above")

    if prov.checks:
        statuses = {c.status for c in prov.checks}
        print()
        print(f"  CONSISTENCY: {', '.join(sorted(statuses))}")
        for check in prov.checks:
            if check.status != "ok":
                print(f"    {check.status}: {check.detail}")
    return EXIT_OK


# --- --families ---


def show_families(
    spec: EquationSystem,
    *,
    as_json: bool,
    fields: tuple[str, ...],
) -> int:
    """Print the spec's tensor families and their index structure.

    Returns
    -------
    int
        Always ``EXIT_OK``.
    """
    from tidal.symbolic.spec_query import field_families

    records: list[dict[str, Any]] = []
    for family in field_families(spec):
        slots = family.group_by_temporal_slots()
        records.append(
            {
                "head": family.head,
                "rank": family.rank,
                "members": list(family.members),
                "exact": family.exact,
                "temporal_slots": {str(k): list(v) for k, v in slots.items()},
            },
        )

    if as_json:
        print(json.dumps({"items": [_project(r, fields) for r in records]}, indent=2))
        return EXIT_OK

    for record in records:
        marker = (
            "" if record["exact"] else "   (grouping GUESSED from name, no metadata)"
        )
        print(
            f"family {record['head']!r}  rank={record['rank']}  "
            f"n={len(record['members'])}{marker}",
        )
        for slots, members in record["temporal_slots"].items():
            print(f"    {slots} temporal index/indices: {', '.join(members)}")
        if not record["temporal_slots"]:
            print(f"    members: {', '.join(record['members'])}")
    return EXIT_OK


# --- --diff ---


def show_diff(
    left: EquationSystem,
    right: EquationSystem,
    *,
    as_json: bool,
    fields: tuple[str, ...],
) -> int:
    """Compare two specs, separating real changes from rescalings.

    Returns
    -------
    int
        ``EXIT_DIFFERENCES`` when real changes exist (a declared outcome, in
        the sense of ``diff(1)``), otherwise ``EXIT_OK``.
    """
    from tidal.symbolic.spec_query import diff_systems

    diff = diff_systems(left, right)
    records = [
        {
            "field": c.field,
            "verdict": c.verdict,
            "changed_keys": list(c.changed_keys),
            "undecided_keys": list(c.undecided_keys),
            "detail": c.detail,
        }
        for c in diff.comparisons
        if c.verdict != "identical"
    ]

    if as_json:
        print(
            json.dumps(
                {
                    "items": [_project(r, fields) for r in records],
                    "only_left": list(diff.only_left),
                    "only_right": list(diff.only_right),
                    "has_real_changes": diff.has_real_changes,
                },
                indent=2,
            ),
        )
        return EXIT_DIFFERENCES if diff.has_real_changes else EXIT_OK

    if diff.only_left or diff.only_right:
        print(f"components only in first : {', '.join(diff.only_left) or 'none'}")
        print(f"components only in second: {', '.join(diff.only_right) or 'none'}")
        print()

    if diff.real:
        print("REAL CHANGES (physics differs):")
        for comparison in diff.real:
            print(f"  {comparison.field}: {', '.join(comparison.changed_keys)}")
    if diff.representational:
        print()
        print("REPRESENTATIONAL ONLY (rescaled; physics unchanged):")
        for comparison in diff.representational:
            print(f"  {comparison.field}: {comparison.detail}")
    if diff.undecided:
        print()
        print("UNDECIDED (could not be settled either way):")
        for comparison in diff.undecided:
            print(f"  {comparison.field}: {', '.join(comparison.undecided_keys)}")
    if not records:
        print("no differences")

    return EXIT_DIFFERENCES if diff.has_real_changes else EXIT_OK


def parse_params(raw: Sequence[str] | None, spec: EquationSystem) -> dict[str, float]:
    """Parse ``--param KEY=VALUE`` arguments over the spec's metadata defaults.

    Parameters
    ----------
    raw : Sequence[str] | None
        Raw ``KEY=VALUE`` strings.
    spec : EquationSystem
        Spec supplying default parameter values.

    Returns
    -------
    dict[str, float]
        Merged parameter values.

    Raises
    ------
    ValueError
        If an entry is malformed or non-numeric.
    """
    params: dict[str, float] = {}
    meta = spec.metadata.get("parameters")
    if isinstance(meta, dict):
        for key, value in meta.items():
            try:
                params[str(key)] = float(value)
            except (TypeError, ValueError):
                continue

    for item in raw or ():
        if "=" not in item:
            msg = f"Invalid --param format: {item!r}. Expected KEY=VALUE"
            raise ValueError(msg)
        key, value_text = item.split("=", 1)
        try:
            params[key.strip()] = float(value_text.strip())
        except ValueError:
            msg = f"Invalid --param value {value_text!r} for {key.strip()!r}"
            raise ValueError(msg) from None
    return params


def run_query(args: Namespace, spec: EquationSystem) -> int | None:  # noqa: PLR0911
    """Dispatch whichever query flag was given, or ``None`` if none was.

    Returns
    -------
    int | None
        Exit code, or ``None`` when no query flag is active so the caller
        falls through to the default display.
    """
    as_json = bool(getattr(args, "json_output", False))
    fields = _split_names(getattr(args, "fields", None))
    assume_positive = _split_names(getattr(args, "assume_positive", None))
    assume_nonzero = _split_names(getattr(args, "assume_nonzero", None))

    raw_params = getattr(args, "param", None)
    try:
        parameters = parse_params(raw_params, spec) if raw_params else None
    except ValueError as exc:
        return _fail(
            "invalid_parameter",
            str(exc),
            "Example: --param B0=0.01",
            as_json=as_json,
        )

    if getattr(args, "families", False):
        return show_families(spec, as_json=as_json, fields=fields)

    if getattr(args, "equation", None):
        return show_equation(
            spec,
            args.equation,
            as_json=as_json,
            fields=fields,
            parameters=parameters,
            assume_positive=assume_positive,
            assume_nonzero=assume_nonzero,
        )

    if getattr(args, "coefficient", None):
        return show_coefficient(
            spec,
            args.coefficient,
            as_json=as_json,
            fields=fields,
            parameters=parameters,
            assume_positive=assume_positive,
            assume_nonzero=assume_nonzero,
        )

    if getattr(args, "diff", None):
        from pathlib import Path

        from tidal.symbolic.json_loader import load_equation_system

        other_path = Path(args.diff)
        if not other_path.exists():
            return _fail(
                "file_not_found",
                f"file not found: {other_path}",
                "Pass the JSON spec to compare against",
                as_json=as_json,
            )
        other = load_equation_system(other_path, strict_v6=False)
        return show_diff(spec, other, as_json=as_json, fields=fields)

    return None


def _fail(
    kind: str,
    message: str,
    hint: str,
    *,
    as_json: bool,
) -> int:
    """Report an error and return :data:`EXIT_ERROR`.

    Always prints the human form via ``error_with_hint``. Under ``--json`` it
    additionally writes a machine-readable envelope as the last line of stderr,
    so a caller can branch on a stable ``kind`` instead of matching on prose.
    Data stays on stdout, so a ``| jq`` pipeline is unaffected either way.

    Parameters
    ----------
    kind : str
        Stable, machine-readable error classification.
    message : str
        Human-readable description.
    hint : str
        Actionable suggestion.
    as_json : bool
        Whether structured output was requested.

    Returns
    -------
    int
        ``EXIT_ERROR``, so callers can ``return _fail(...)``.
    """
    from tidal.cli._console import error_with_hint

    error_with_hint(message, hints=[hint])
    if as_json:
        payload = {"error": {"kind": kind, "message": message, "hint": hint}}
        print(json.dumps(payload), file=sys.stderr)
    return EXIT_ERROR
