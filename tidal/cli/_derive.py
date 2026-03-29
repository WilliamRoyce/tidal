"""``tidal derive`` — Derive equations from Lagrangian via Wolfram/xAct."""

from __future__ import annotations

import hashlib
import json as _json_mod
import re
import shutil
import subprocess
import sys
import tempfile
import tomllib
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from argparse import Namespace

# Re-export validation API (moved to _derive_validate.py for readability)
from tidal.cli._derive_validate import (  # noqa: F401
    _GAUGE_PRESETS,  # type: ignore[reportPrivateUsage]
    _validate_config,  # type: ignore[reportPrivateUsage]
    _validate_reduction,  # type: ignore[reportPrivateUsage, reportUnusedImport]
)
from tidal.cli._wls_helpers import (
    validate_wls_brackets,
    wl_bg_rule_entry,
    wl_component_value,
    wl_diag_matrix,
    wl_flatten_list,
    wl_index,
    wl_list,
    wl_zero_component,
)

# Coordinates per spacetime dimension
_COORDS: dict[int, list[str]] = {
    2: ["t", "x"],
    3: ["t", "x", "y"],
    4: ["t", "x", "y", "z"],
    5: ["t", "x", "y", "z", "w"],
    6: ["t", "x", "y", "z", "w", "v"],
    7: ["t", "x", "y", "z", "w", "v", "u"],
}

_INDEX_LETTERS = list("abcdefghijklmnop")  # 16 letters, enough for dim 7 + 4 = 11

# Minkowski signatures per spacetime dimension
_MINKOWSKI_SIGNATURES: dict[int, list[int]] = {
    2: [-1, 1],
    3: [-1, 1, 1],
    4: [-1, 1, 1, 1],
}


_MIN_PREFIX_LEN = 2

# --- Code generation helpers ---


def _make_prefix(config: dict[str, Any]) -> str:
    """Generate a 2-3 letter symbol prefix from the theory name."""
    name = config.get("theory", {}).get("name", "")
    if not name:
        return "tidal"
    # Take initials of first 2-3 words (alpha only), lowercase
    words = name.split()
    initials = [w[0].lower() for w in words if w[0].isalpha()]
    prefix = "".join(initials[:3])
    return prefix if len(prefix) >= _MIN_PREFIX_LEN else "tidal"


def _generate_metric_code(config: dict[str, Any], prefix: str) -> str:
    """Generate Wolfram code for metric definition.

    Raises
    ------
    ValueError
        If metric type is unknown or dimension mismatch in diagonal/matrix.
    """
    dim = config["spacetime"]["dimension"]
    metric_type = config["spacetime"].get("metric", "minkowski")

    if metric_type == "minkowski":
        sig = _MINKOWSKI_SIGNATURES.get(dim)
        if sig is None:
            sig = [-1] + [1] * (dim - 1)
        diag_str = ", ".join(str(s) for s in sig)
        return f"{prefix}MetricMatrix = {wl_diag_matrix(diag_str)};"

    if metric_type == "diagonal":
        entries = config["spacetime"].get("diagonal", [])
        if len(entries) != dim:
            msg = f"[spacetime].diagonal must have {dim} entries, got {len(entries)}"
            raise ValueError(msg)
        # Entries can be numbers or strings (for coordinate-dependent)
        parts = [str(e) for e in entries]
        return f"{prefix}MetricMatrix = {wl_diag_matrix(', '.join(parts))};"

    if metric_type == "matrix":
        matrix = config["spacetime"].get("matrix", [])
        if len(matrix) != dim:
            msg = f"[spacetime].matrix must be {dim}x{dim}"
            raise ValueError(msg)
        rows: list[str] = []
        for row in matrix:
            if len(row) != dim:
                msg = f"[spacetime].matrix rows must have {dim} entries"
                raise ValueError(msg)
            rows.append(wl_list(*(str(e) for e in row)))
        return f"{prefix}MetricMatrix = {wl_list(*rows)};"

    msg = f"Unknown metric type: '{metric_type}'. Use: minkowski, diagonal, matrix"
    raise ValueError(msg)


def _generate_field_def(
    field: dict[str, Any],
    prefix: str,
    manifold: str,
    *,
    head_override: str | None = None,
) -> str:
    """Generate DefTensor code for a field."""
    name = field["name"]
    ftype = field["type"]
    prefixed = head_override or f"{prefix}{name.capitalize()}"

    if ftype == "scalar":
        return f"If[!xTensorQ[{prefixed}],\n  DefTensor[{prefixed}[], {manifold}]\n];"

    if ftype == "vector":
        return f"If[!xTensorQ[{prefixed}],\n  DefTensor[{prefixed}[-a], {manifold}]\n];"

    # Tensor with rank and symmetry
    rank = field["rank"]
    symmetry = field.get("symmetry", "none")
    indices = ", ".join(f"-{_INDEX_LETTERS[i]}" for i in range(rank))

    if symmetry == "antisymmetric":
        sym_spec = f", Antisymmetric[{wl_list(indices)}]"
    elif symmetry == "symmetric":
        sym_spec = f", Symmetric[{wl_list(indices)}]"
    elif symmetry.startswith("antisymmetric_"):
        # Partial antisymmetry: "antisymmetric_23" → antisymmetric in slots 2,3
        # xAct syntax: Antisymmetric[{-b, -c}] for rank-3 T[-a,-b,-c]
        slot_digits = symmetry.split("_", 1)[1]
        slot_indices = ", ".join(f"-{_INDEX_LETTERS[int(s) - 1]}" for s in slot_digits)
        sym_spec = f", Antisymmetric[{wl_list(slot_indices)}]"
    elif symmetry.startswith("symmetric_"):
        # Partial symmetry: "symmetric_12" → symmetric in slots 1,2
        slot_digits = symmetry.split("_", 1)[1]
        slot_indices = ", ".join(f"-{_INDEX_LETTERS[int(s) - 1]}" for s in slot_digits)
        sym_spec = f", Symmetric[{wl_list(slot_indices)}]"
    else:
        sym_spec = ""

    return (
        f"If[!xTensorQ[{prefixed}],\n"
        f"  DefTensor[{prefixed}[{indices}], {manifold}{sym_spec},\n"
        f'    PrintAs -> "{name}"]\n'
        f"];"
    )


def _field_expression(
    field: dict[str, Any],
    prefix: str,
    *,
    head_override: str | None = None,
) -> str:
    """Return the xAct expression for a field reference (e.g., 'phi[]' or 'C[-a,-b,-c]')."""
    name = field["name"]
    ftype = field["type"]
    prefixed = head_override or f"{prefix}{name.capitalize()}"

    if ftype == "scalar":
        return f"{prefixed}[]"
    if ftype == "vector":
        return f"{prefixed}[-a]"
    rank = field["rank"]
    indices = ", ".join(f"-{_INDEX_LETTERS[i]}" for i in range(rank))
    return f"{prefixed}[{indices}]"


def _substitute_field_names(
    expression: str,
    fields: list[dict[str, Any]],
    prefix: str,
    *,
    derived_fields: list[dict[str, Any]] | None = None,
    background_fields: list[dict[str, Any]] | None = None,
) -> str:
    """Replace user field names with prefixed xAct names in the Lagrangian."""
    result = expression

    # Substitute built-in names (eta, CD, CDT, bg) FIRST — before field names —
    # to prevent short field names (e.g., "a") from corrupting built-in
    # identifiers (e.g., "eta[" contains "a[", which would become "etgeA[").
    # CDT must be substituted BEFORE CD to avoid CD matching inside CDT.
    result = result.replace("eta[", f"{prefix}Eta[")
    result = result.replace("bg[", f"{prefix}Bg[")
    # Torsion CovD (CDT) and Levi-Civita CovD (CD).
    # Replace "CDT" in ALL positions (standalone CDT[-a], compound RicciScalarCDT[],
    # TorsionCDT[], etc.).  Single global replace avoids double-prefixing.
    result = result.replace("CDT", f"{prefix}CDT")
    # Levi-Civita CD — must NOT match the "CD" inside "{prefix}CDT".
    # Replace only bare "CD" not followed by "T" (negative lookahead).
    result = re.sub(r"CD(?!T)", f"{prefix}CD", result)
    # Substitute chart placeholder for component-derivative notation
    # e.g., CD[{0, -chart}][ux[]] → {prefix}CD[{0, -{prefix}Cart}][...]
    result = result.replace("-chart}", f"-{prefix}Cart}}")

    # Merge fundamental, derived, and background fields for substitution
    all_fields = list(fields)
    if derived_fields:
        all_fields.extend(derived_fields)
    if background_fields:
        all_fields.extend(background_fields)

    # Sort by name length descending to avoid partial replacements
    sorted_fields = sorted(all_fields, key=lambda f: len(f["name"]), reverse=True)

    for field in sorted_fields:
        name = field["name"]
        prefixed_name = f"{prefix}{name.capitalize()}"
        # Replace field name references using word-boundary-aware regex.
        # Field names appear after non-alphanumeric chars or at start of string.
        # Must match "name[" or "name " but NOT "eta[" matching "a[" inside it.
        # Lookbehind: (?<![a-zA-Z]) ensures we don't match inside longer identifiers.
        result = re.sub(
            rf"(?<![a-zA-Z]){re.escape(name)}\[",
            f"{prefixed_name}[",
            result,
        )
        result = re.sub(
            rf"(?<![a-zA-Z]){re.escape(name)} ",
            f"{prefixed_name} ",
            result,
        )

    return result


# --- WLS: Context & header ---


@dataclass(frozen=True)
class _WlsContext:
    """Shared context for WLS script generation helpers."""

    prefix: str
    dim: int
    fields: list[dict[str, Any]]
    constants: list[str]
    coords: list[str]
    manifold: str
    metric: str
    cd: str
    chart: str
    torsion: dict[str, Any] | None  # [torsion] TOML section (None = no torsion)
    cdt: str  # torsion-full CovD name (e.g. "geCDT")
    theory_name: str
    output_path: str
    lagrangian_expr: str
    is_multi: bool
    pipeline_path: str
    parameters: dict[str, float]
    derived_fields: list[dict[str, Any]]
    background_fields: list[dict[str, Any]]
    linearization: dict[str, Any] | None
    constraint_solver: dict[str, Any] | None
    gauge: list[dict[str, Any]]
    reduction: dict[str, Any] | None
    metric_diagonal: list[str]
    metric_type: str  # "minkowski", "diagonal", or "matrix"


def _wls_mem_print(label: str) -> str:
    """Generate a Print statement that includes current Wolfram memory usage."""
    safe_label = label.replace('"', '\\"')
    return f'Print["[", Round[MemoryInUse[]/1024.^2], " MB] {safe_label}"];'


def _wls_timing_start(timer_var: str) -> str:
    """Generate WLS code to start a named timer using AbsoluteTime[]."""
    return f"{timer_var} = AbsoluteTime[];"


def _wls_timing_end(timer_var: str, label: str) -> str:
    """Generate WLS code to print elapsed time since timer start."""
    safe_label = label.replace('"', '\\"')
    return (
        f'Print["[TIMING] {safe_label}: ",'
        f" Round[AbsoluteTime[] - {timer_var}, 0.1],"
        f' " seconds"];'
    )


def _wls_header(ctx: _WlsContext) -> list[str]:
    """Generate script header lines."""
    return [
        "#!/usr/bin/env wolframscript",
        f"(* Auto-generated by tidal derive: {ctx.theory_name} *)",
        "",
        "(* Prevent Wolfram from caching all In/Out expressions — reduces memory *)",
        "$HistoryLength = 0;",
        "",
        "(* Increase recursion limit for complex tensor theories — xAct internal   *)",
        "(* evaluation can exceed the default 1024 for expressions with many        *)",
        "(* contracted index pairs (4D R̃² torsion, Einstein-Maxwell, etc.).       *)",
        "$RecursionLimit = 16384;",
        "",
        f'Print["=== {ctx.theory_name} ==="];',
        'Print[""];',
        "",
    ]


def _wls_packages(
    pipeline_path: str,
    *,
    load_xpert: bool = False,
    load_gauge: bool = False,
) -> list[str]:
    """Generate xAct package loading and pipeline import lines.

    Parameters
    ----------
    pipeline_path : str
        Absolute path to the ``tidal/wolfram/`` directory.
    load_xpert : bool
        If *True*, also load ``xAct`xPert`` and ``Linearize.wl``.
    load_gauge : bool
        If *True*, also load ``GaugeFix.wl``.
    """
    escaped = pipeline_path.replace("\\", "\\\\").replace('"', '\\"')
    lines = [
        "(* Load xAct packages *)",
        "<< xAct`xTensor`;",
        "<< xAct`xCoba`;",
    ]
    if load_xpert:
        lines.append("<< xAct`xPert`;")
    lines.extend(
        (
            "",
            "(* Load pipeline modules *)",
            f'pipelinePath = "{escaped}";',
            'Get[FileNameJoin[{pipelinePath, "CommonUtilities.wl"}]];',
            'Get[FileNameJoin[{pipelinePath, "EulerLagrange.wl"}]];',
            'Get[FileNameJoin[{pipelinePath, "ComponentDecompose.wl"}]];',
            'Get[FileNameJoin[{pipelinePath, "ExportJSON.wl"}]];',
        )
    )
    if load_xpert:
        lines.append('Get[FileNameJoin[{pipelinePath, "Linearize.wl"}]];')
    if load_gauge:
        lines.append('Get[FileNameJoin[{pipelinePath, "GaugeFix.wl"}]];')
    lines.extend(
        (
            "",
            "$DefInfoQ = False;",
            "(* Suppress xAct internal messages that add cumulative overhead.         *)",
            "(* Ref: supervisor's Initialization.m lines 10-15.                      *)",
            "Off@PrintAsCharacter::argx;",
            "Off@CollectTensors::denominator;",
            "Unprotect@AutomaticRules;",
            "Options[AutomaticRules] = {Verbose -> False};",
            "Protect@AutomaticRules;",
            "",
        )
    )
    return lines


# --- Coordinate-values pre-evaluation ---


def _apply_coord_values(exprs: list[str], coord_values: dict[str, str]) -> list[str]:
    """Substitute killed coordinates into Wolfram expression strings at the Python level.

    Applies ``coord_values`` to all Python-generated Wolfram string expressions before
    they are emitted into the WLS. This ensures that:

    * the TT-traceless substitution weights (``metric_diagonal``) have no killed-coordinate
      factors (e.g. ``x[]^2*Sin[Pi/2]^2`` rather than ``x[]^2*Sin[y[]]^2``),
    * background-field component strings evaluate cleanly if they depend on killed coords.

    For example, with ``coord_values = {"y": "Pi/2"}`` the replacement is::

        "x[]^2*Sin[y[]]^2"  →  "x[]^2*Sin[Pi/2]^2"

    Wolfram then auto-evaluates ``Sin[Pi/2] → 1``, so the expensive ``Simplify``
    in the TT-weight expression immediately collapses to a trivial result.

    NOTE: This function does NOT apply to the Wolfram-level metric matrix used for
    ``MetricInBasis``/``SetMetricDownValues``.  Substituting killed coordinates into the
    background metric before xPert's L^(2) expansion creates a non-flat background and
    produces spurious Riemann coupling terms in tensor field equations.  The coordinate
    evaluation of ``fieldEquations`` is handled by ``_wls_plane_wave_coordinate_evaluation``.

    Parameters
    ----------
    exprs:
        List of Wolfram expression strings (e.g. ``metric_diagonal``,
        ``background_field.components``).
    coord_values:
        Mapping ``{coord_name: wolfram_value}`` from the ``[reduction]`` TOML section.
        E.g. ``{"y": "Pi/2"}`` replaces every occurrence of ``y[]`` with ``Pi/2``.
    """
    if not coord_values:
        return exprs
    result = list(exprs)
    for coord, val in coord_values.items():
        result = [e.replace(f"{coord}[]", val) for e in result]
    return result


# --- WLS: Spacetime, fields & Lagrangian ---


def _wls_spacetime(config: dict[str, Any], ctx: _WlsContext) -> list[str]:
    """Emit spacetime manifold, metric, chart, and coordinate definitions."""
    coord_funcs = ", ".join(f"{c}[]" for c in ctx.coords)
    indices = ", ".join(str(i) for i in range(ctx.dim))
    idx_str = ", ".join(_INDEX_LETTERS[: min(ctx.dim + 4, 8)])

    lines = [
        f"(* Step 1: Define {ctx.dim}D spacetime *)",
        f"If[!xTensorQ[{ctx.manifold}],",
        f"  DefManifold[{ctx.manifold}, {ctx.dim}, {wl_list(idx_str)}]",
        "];",
        "",
        f"If[!MetricQ[{ctx.metric}],",
        f"  DefMetric[-1, {ctx.metric}[-a, -b], {ctx.cd},",
        '    SymbolOfCovD -> {";", "\\[Del]"},',
        '    PrintAs -> "\\[Eta]"]',
        "];",
        "",
        f"If[!ChartQ[{ctx.chart}],",
        f"  DefChart[{ctx.chart}, {ctx.manifold}, {wl_list(indices)}, {wl_list(coord_funcs)}]",
        "];",
        "",
        _generate_metric_code(config, ctx.prefix),
    ]

    # NOTE: coordinate_values (e.g. {y = "Pi/2"}) must NOT be applied to the metric
    # matrix here. Substituting y→Pi/2 into the metric before xPert's L^(2) expansion
    # creates a non-flat background: diag(-1,1,r²,r²) has R^θ_{φφ}θ = -1 ≠ 0, while
    # the true spherical metric diag(-1,1,r²,r²sin²θ) is flat (Riemann=0 everywhere).
    # xPert's background Riemann coupling terms in L^(2) would then contribute spurious
    # curvature corrections to tensor component equations (h_θθ, h_θφ) that are NOT
    # zeroed by the {RicciCD→0} substitution (Riemann ≠ Ricci).
    # The coordinate evaluation is correctly deferred to _wls_plane_wave_coordinate_evaluation,
    # which applies "fieldEquations /. {y[] -> Pi/2}" after full symbolic derivation.
    # The Python-level _apply_coord_values is still applied to metric_diagonal strings
    # (for TT-traceless weights and background field components) since those are correct.

    lines += [
        f"MetricInBasis[{ctx.metric}, -{ctx.chart}, {ctx.prefix}MetricMatrix];",
        f"SetMetricDownValues[{ctx.metric}, {ctx.chart}, {ctx.prefix}MetricMatrix];",
        "",
        "(* Pre-compute ALL geometric quantities (Christoffels, Riemann, Ricci) as   *)",
        "(* xCoba ComponentValues.  ToValues[] then auto-substitutes these during    *)",
        "(* the staggered ToBasis pipeline, eliminating the need for manual          *)",
        "(* EvaluateChristoffelComponents/EvaluateCurvatureComponents post-steps.    *)",
        "(* Ref: supervisor's xCoba example (commit 4a89164) — MetricCompute[G,chart,All] *)",
        f"Block[{{Print = Null}}, MetricCompute[{ctx.metric}, {ctx.chart}, All]];",
        f'Print["MetricCompute: pre-computed all geometric quantities for {ctx.chart}"];',
        "",
    ]

    # Define torsion-full covariant derivative for PGT theories.
    # CDT is metric-compatible but has torsion: Gamma = Christoffel + contortion.
    # xAct auto-creates TorsionCDT[a, -b, -c], RicciScalarCDT[], etc.
    # The Lagrangian can then reference these operators directly.
    # Before perturbation, ChangeCurvature/ChangeTorsion decomposes to CD form
    # for xPert compatibility.
    # Ref: xAct DefCovD documentation; Blagojevic & Hehl (2013).
    if ctx.torsion:
        lines += [
            "(* === Torsion-full covariant derivative (Poincaré gauge theory) === *)",
            f"If[!CovDQ[{ctx.cdt}],",
            f"  DefCovD[{ctx.cdt}[-a], "
            f'{{"#", "DT"}}, '
            f"FromMetric -> {ctx.metric}, "
            f"Torsion -> True]",
            "];",
            f'Print["Torsion CovD defined: ", {ctx.cdt}];',
            "",
        ]

    return lines


def _needs_bg_tensor(config: dict[str, Any]) -> bool:
    """Check if any expression references the background metric ``bg[...]``.

    The background metric is a non-dynamical reference tensor used in theories
    like massive gravity (Fierz-Pauli) and bimetric gravity.  It is defined
    via ``DefTensor`` so that xPert treats it as unperturbed.
    """
    exprs: list[str] = []
    if "lagrangian" in config:
        exprs.append(config["lagrangian"].get("expression", ""))
    if "linearization" in config:
        exprs.append(config["linearization"].get("expression", ""))
    exprs.extend(df.get("definition", "") for df in config.get("derived_fields", []))
    return any("bg[" in expr for expr in exprs)


# --- WLS: Background field setup ---


def _wls_background_component_values(
    field: dict[str, Any],
    prefix: str,
    chart: str,
    dim: int,
) -> list[str]:
    """Emit xCoba ComponentValue declarations for background field components."""
    comps: list[int | float | str] = field.get("components", [])
    prefixed = f"{prefix}{field['name'].capitalize()}"
    ftype = field["type"]
    lines: list[str] = []

    if ftype == "scalar":
        lines.append(f"ComponentValue[{prefixed}[], {comps[0]}];")
    elif ftype == "vector":
        for idx, val in enumerate(comps):
            lines.append(f"ComponentValue[{prefixed}[{wl_index(idx, chart)}], {val}];")
    else:
        # Tensor rank 2+: iterate over all index tuples
        rank = field.get("rank", 2)
        for flat_idx, val in enumerate(comps):
            multi_idx: list[int] = []
            remaining = flat_idx
            for _ in range(rank):
                multi_idx.append(remaining % dim)
                remaining //= dim
            multi_idx.reverse()
            idx_str = ", ".join(wl_index(k, chart) for k in multi_idx)
            lines.append(f"ComponentValue[{prefixed}[{idx_str}], {val}];")

    return lines


def _wls_validate_backgrounds_after_decompose(
    ctx: _WlsContext,
    comp_var: str,
) -> list[str]:
    """Generate validation that vector/tensor backgrounds resolved after ToBasis.

    After DecomposeToComponents, all background tensor symbols should be
    fully resolved to coordinate expressions via ComponentValue + ToBasis.
    This emits a check that catches silent substitution failures.

    Only emitted when there are non-scalar background fields (scalars use
    explicit ReplaceAll, so they're always resolved).
    """
    non_scalar_bgs = [f for f in ctx.background_fields if f["type"] != "scalar"]
    if not non_scalar_bgs:
        return []

    bg_heads = ", ".join(
        f"{ctx.prefix}{f['name'].capitalize()}" for f in non_scalar_bgs
    )
    return [
        "(* Validate background field substitution succeeded *)",
        f"ValidateNoUnresolvedBackgrounds[{comp_var}, {wl_list(bg_heads)}];",
    ]


def _wls_scalar_background_substitution(
    ctx: _WlsContext,
    eom_var: str,
) -> list[str]:
    """Generate explicit ``ReplaceAll`` for scalar background fields.

    Scalar backgrounds (rank 0) have no indices, so ``ToBasis`` inside
    ``DecomposeToComponents`` does **not** trigger ``ComponentValue``
    substitution.  We must substitute explicitly before decomposition.

    Vector/tensor backgrounds need explicit substitution AFTER
    decomposition — see ``_wls_vector_background_substitution``.
    """
    lines: list[str] = []
    for field in ctx.background_fields:
        if field["type"] != "scalar":
            continue
        prefixed = f"{ctx.prefix}{field['name'].capitalize()}"
        comps = field.get("components", [])
        if not comps:
            continue
        value = comps[0]
        lines.extend(
            (
                f"(* Substitute scalar background {field['name']} -> {value} *)",
                f"{eom_var} = {eom_var} /. {wl_list(f'{prefixed}[] -> {value}')};",
            )
        )
    return lines


def _compute_contra_components(
    comps: list[int | float | str],
    metric_diagonal: list[str],
) -> list[str]:
    """Compute contravariant component values A^μ from covariant A_μ.

    For a diagonal metric g_{μμ}: A^μ = A_μ / g_{μμ}.
    For Minkowski (empty metric_diagonal): A^μ = A_μ (g^{ii} = 1 spatially;
    temporal component A^0 = -A_0 but Abar_0 is always 0 in our gauge).
    Returns a list of Wolfram expression strings.
    """
    contra: list[str] = []
    for idx, val in enumerate(comps):
        val_str = str(val)
        if val_str.strip() in {"0", "0.0"}:
            contra.append("0")
        elif metric_diagonal:
            g_diag = metric_diagonal[idx]
            contra.append(f"Simplify[({val_str}) / ({g_diag})]")
        else:
            # Minkowski: spatial A^i = A_i / 1 = A_i
            contra.append(val_str)
    return contra


def _wls_vector_background_substitution(  # noqa: PLR0914
    ctx: _WlsContext,
    comp_var: str,
) -> list[str]:
    """Generate explicit ``ReplaceAll`` for vector/tensor background fields.

    Unlike scalar backgrounds (substituted BEFORE decomposition because
    ``ToBasis`` doesn't trigger ``ComponentValue`` for rank-0 tensors),
    vector/tensor backgrounds are substituted AFTER decomposition so that
    xAct handles the index algebra (contractions, metric raising/lowering)
    before we inject numeric component values.

    For curved diagonal metrics, contravariant components ``A^μ`` differ
    from covariant components ``A_μ``:  ``A^μ = g^{μμ} A_μ = A_μ / g_{μμ}``.
    When ``ctx.metric_diagonal`` is non-empty, the correct contravariant
    value is computed via ``Simplify[(A_μ) / (g_{μμ})]`` so Wolfram can
    cancel coordinate factors algebraically.

    For Minkowski (``ctx.metric_diagonal`` empty), both orientations get the
    same value (``g^{μμ} = ±1`` so sign is already encoded in A_μ).
    """
    non_scalar_bgs = [f for f in ctx.background_fields if f["type"] != "scalar"]
    if not non_scalar_bgs:
        return []

    lines: list[str] = []
    for field in non_scalar_bgs:
        prefixed = f"{ctx.prefix}{field['name'].capitalize()}"
        comps: list[int | float | str] = field.get("components", [])
        if not comps:
            continue

        rules: list[str] = []
        if field["type"] == "vector":
            for idx, val in enumerate(comps):
                # Compute correct contravariant component value.
                # For curved diagonal metrics: A^μ = A_μ / g_{μμ}.
                # For Minkowski (metric_diagonal empty): A^μ = A_μ (flat).
                val_str = str(val)
                if ctx.metric_diagonal and val_str.strip() not in {"0", "0.0"}:
                    g_diag = ctx.metric_diagonal[idx]
                    contra_val = f"Simplify[({val_str}) / ({g_diag})]"
                else:
                    contra_val = val_str
                rules.extend(
                    (
                        f"{prefixed}[{wl_index(idx, ctx.chart)}] -> {val_str}",
                        f"{prefixed}[{wl_index(idx, ctx.chart, contra=True)}] -> {contra_val}",
                        # Component function form (after ReplaceTensorFieldComponents
                        # converts vbdB[{i, -chart}] -> vbdBi[t, x, y])
                        f"{prefixed}{idx}[__] -> {val_str}",
                    )
                )
        else:
            # Tensor rank 2+: iterate over all index tuples.
            # For curved diagonal metrics, the fully-contravariant value T^{uv}
            # needs two metric factors: T^{uv} = (T_{uv}) / (g_{uu} * g_{vv}).
            # This generalizes the vector fix above.
            rank = field.get("rank", 2)
            for flat_idx, val in enumerate(comps):
                multi_idx: list[int] = []
                remaining = flat_idx
                for _ in range(rank):
                    multi_idx.append(remaining % ctx.dim)
                    remaining //= ctx.dim
                multi_idx.reverse()
                idx_down = ", ".join(wl_index(k, ctx.chart) for k in multi_idx)
                idx_up = ", ".join(
                    wl_index(k, ctx.chart, contra=True) for k in multi_idx
                )
                # Component function name: head + concatenated index digits
                comp_name = "".join(str(k) for k in multi_idx)
                val_str = str(val)
                # Fully contravariant: divide by each metric diagonal entry
                if ctx.metric_diagonal and val_str.strip() not in {"0", "0.0"}:
                    g_factors = " * ".join(
                        f"({ctx.metric_diagonal[k]})" for k in multi_idx
                    )
                    contra_val = f"Simplify[({val_str}) / ({g_factors})]"
                else:
                    contra_val = val_str
                rules.extend(
                    (
                        f"{prefixed}[{idx_down}] -> {val_str}",
                        f"{prefixed}[{idx_up}] -> {contra_val}",
                        # Component function form (after ReplaceTensorFieldComponents)
                        f"{prefixed}{comp_name}[__] -> {val_str}",
                    )
                )

        rules_str = ", ".join(rules)
        lines.extend(
            [
                f"(* Substitute vector/tensor background {field['name']} *)",
                f"{comp_var} = {comp_var} /. {wl_list(rules_str)};",
            ]
        )

    return lines


def _wls_fields(ctx: _WlsContext, *, include_bg: bool = False) -> list[str]:
    """Generate field definitions and constant symbol declarations."""
    lines: list[str] = ["(* Step 2: Define fields *)"]
    for field in ctx.fields:
        lines.extend((_generate_field_def(field, ctx.prefix, ctx.manifold), ""))

    if ctx.constants:
        lines.append("(* Define constant symbols *)")
        lines.extend(
            f"If[!ConstantSymbolQ[{c}], DefConstantSymbol[{c}]];" for c in ctx.constants
        )
        lines.append("")

    # Background fields — non-dynamical tensors (not varied by VarD)
    if ctx.background_fields:
        lines.append(
            "(* Background fields — non-dynamical (not varied in Euler-Lagrange) *)"
        )
        for field in ctx.background_fields:
            lines.extend((_generate_field_def(field, ctx.prefix, ctx.manifold), ""))
            # Set component values via xCoba's ComponentValue mechanism.
            lines.extend(
                _wls_background_component_values(field, ctx.prefix, ctx.chart, ctx.dim)
            )
            # Also set direct Mathematica DownValues for auto-collapse during
            # TraceBasisDummy (xCoba ComponentValues don't auto-evaluate).
            bg_prefixed = f"{ctx.prefix}{field['name'].capitalize()}"
            bg_comps = field.get("components", [])
            if field["type"] == "vector" and bg_comps:
                comps_str = ", ".join(str(c) for c in bg_comps)
                contra_comps = _compute_contra_components(bg_comps, ctx.metric_diagonal)
                contra_str = ", ".join(contra_comps)
                lines.append(
                    f"SetBackgroundFieldDownValues[{bg_prefixed}, {ctx.chart},"
                    f" {wl_list(comps_str)}, {wl_list(contra_str)}];"
                )
            lines.append("")
        lines.append("")

    if include_bg:
        bg_name = f"{ctx.prefix}Bg"
        lines.extend(
            [
                "(* Background/reference metric — not perturbed by xPert *)",
                f"If[!xTensorQ[{bg_name}],",
                f'  DefTensor[{bg_name}[-a, -b], {ctx.manifold}, Symmetric[{{-a, -b}}], PrintAs -> "bg"]',
                "];",
                "(* Explicit zero perturbation: bg is non-dynamical *)",
                "Unprotect[Perturbation];",
                f"Perturbation[{bg_name}[__], ___] := 0;",
                "Protect[Perturbation];",
                "",
            ]
        )

    return lines


def _wls_derived_fields(ctx: _WlsContext) -> list[str]:
    """Generate DefTensor and MakeRule for each derived field."""
    if not ctx.derived_fields:
        return []

    lines: list[str] = ["(* Derived field definitions *)"]
    for field in ctx.derived_fields:
        # DefTensor (reuse existing helper)
        lines.extend((_generate_field_def(field, ctx.prefix, ctx.manifold), ""))

        # Build the MakeRule LHS: prefixed tensor with lowered indices
        fexpr = _field_expression(field, ctx.prefix)
        # Substitute field names in the definition
        defn = _substitute_field_names(
            field["definition"].strip(),
            ctx.fields,
            ctx.prefix,
            derived_fields=ctx.derived_fields,
            background_fields=ctx.background_fields,
        )
        rule_var = f"{ctx.prefix}{field['name'].capitalize()}Rules"
        lines.extend(
            (
                f"{rule_var} = MakeRule[{wl_list(f'{fexpr}, {defn}')}, MetricOn -> All, ContractMetrics -> True];",
                "",
            )
        )

    return lines


def _wls_lagrangian(ctx: _WlsContext) -> list[str]:
    """Generate Lagrangian definition lines."""
    prefixed = _substitute_field_names(
        ctx.lagrangian_expr,
        ctx.fields,
        ctx.prefix,
        derived_fields=ctx.derived_fields,
        background_fields=ctx.background_fields,
    )
    lines = [
        "(* Step 3: Lagrangian *)",
        f"{ctx.prefix}Lagrangian = (",
        f"  {prefixed}",
        ");",
        "",
    ]

    if ctx.derived_fields:
        lines.append("(* Expand derived field definitions *)")
        for field in ctx.derived_fields:
            rule_var = f"{ctx.prefix}{field['name'].capitalize()}Rules"
            lines.append(
                f"{ctx.prefix}Lagrangian = {ctx.prefix}Lagrangian /. {rule_var};"
            )
        lines.extend(
            (
                f"{ctx.prefix}Lagrangian = ToCanonical[{ctx.prefix}Lagrangian];",
                f"{ctx.prefix}Lagrangian = ContractMetric[{ctx.prefix}Lagrangian, {ctx.metric}];",
                "",
                f'Print["Lagrangian (expanded): ", {ctx.prefix}Lagrangian];',
                "",
            )
        )
    else:
        lines.extend(
            (
                f'Print["Lagrangian: ", {ctx.prefix}Lagrangian];',
                "",
            )
        )

    return lines


def _xpert_index_pattern(
    pert_sym: str,
    full_head: str,
    field: dict[str, Any],
) -> tuple[str, str]:
    """Return (pert_idx, full_idx) for DefTensorPerturbation."""
    ftype = field["type"]
    if ftype == "scalar":
        return f"{pert_sym}[LI[order]]", f"{full_head}[]"
    if ftype == "vector":
        return f"{pert_sym}[LI[order], -a]", f"{full_head}[-a]"
    rank = field.get("rank", 2)
    idx_str = ", ".join(f"-{_INDEX_LETTERS[i]}" for i in range(rank))
    return f"{pert_sym}[LI[order], {idx_str}]", f"{full_head}[{idx_str}]"


def _pert_field_dict(pert_name: str, source_field: dict[str, Any]) -> dict[str, Any]:
    """Build a field dict for a perturbation field (inherits type/rank/symmetry)."""
    d: dict[str, Any] = {"name": pert_name, "type": source_field["type"]}
    if source_field["type"] not in {"scalar", "vector"}:
        d["rank"] = source_field.get("rank", 2)
        d["symmetry"] = source_field.get("symmetry", "none")
    return d


# --- WLS: Linearization (xPert) ---


def _wls_matter_perturbation_setup(  # noqa: PLR0914
    ctx: _WlsContext,
    matter_perts: list[dict[str, Any]],
    eps_sym: str,
) -> tuple[list[dict[str, str]], list[str]]:
    """Generate DefTensorPerturbation + DefTensor for each matter perturbation.

    Returns (matter_pert_info, wls_lines) where matter_pert_info contains
    the symbol mapping for each perturbation field.
    """
    p = ctx.prefix
    info: list[dict[str, str]] = []
    lines: list[str] = []

    for mp in matter_perts:
        mf_name = mp["field"]
        mp_name = mp["perturbation_name"]
        mf = next(f for f in ctx.fields if f["name"] == mf_name)
        mf_type = mf["type"]

        mp_sym = f"{p}{mp_name}Pert"
        mf_head = f"{p}{mf_name.capitalize()}"
        # Avoid collision: if capitalize(perturbation_name) == capitalize(field_name),
        # use the perturbation name as-is (Mathematica is case-sensitive).
        # E.g., field "A" → geA, perturbation "a" → gea (not geA).
        mp_head_candidate = f"{p}{mp_name.capitalize()}"
        mp_head = f"{p}{mp_name}" if mp_head_candidate == mf_head else mp_head_candidate

        # Build DefTensorPerturbation arguments + perturbation field DefTensor
        pert_idx, full_idx = _xpert_index_pattern(mp_sym, mf_head, mf)
        mp_def = _generate_field_def(
            _pert_field_dict(mp_name, mf),
            ctx.prefix,
            ctx.manifold,
            head_override=mp_head,
        )

        lines.extend(
            [
                f"(* Matter perturbation: {mf_name} = {mp['background']} + {eps_sym} * {mp_name} *)",
                f"(* Define perturbation field tensor {mp_name} *)",
                mp_def,
                "",
                f"DefTensorPerturbation[{pert_idx}, {full_idx}, {ctx.manifold}];",
                f'Print["Matter perturbation: {mf_name} -> background + {eps_sym} * {mp_name}"];',
                "",
            ]
        )

        bg_name = mp.get("background", "")
        bg_head = f"{p}{bg_name.capitalize()}" if bg_name else ""
        info.append(
            {
                "field_name": mf_name,
                "pert_name": mp_name,
                "pert_sym": mp_sym,
                "field_head": mf_head,
                "pert_head": mp_head,
                "bg_head": bg_head,
                "field_type": mf_type,
                "field_rank": str(mf.get("rank", 2)),
            }
        )

    return info, lines


def _wls_component_metadata(
    field_name: str, fexpr: str, comp_var: str, dim: int
) -> list[str]:
    """Generate Wolfram code to build tensor component metadata for a field.

    Re-calls ``EnumerateComponentTuples`` to reconstruct the flat-index → tuple
    mapping that was used inside ``DecomposeToComponents`` but not returned.
    For scalars (rank 0), ``SlotsOfTensor`` is not defined, so we check
    ``xTensorQ`` on the head first and fall back to rank=0 with empty indices.
    """
    head_expr = f"Head[{fexpr}]"
    return [
        f"(* Build tensor component metadata for {field_name} *)",
        "Module[{fieldRank, fHead},",
        f"  fHead = {head_expr};",
        "  fieldRank = If[xTensorQ[fHead], Length[SlotsOfTensor[fHead]], 0];",
        "  If[fieldRank > 0,",
        "    Module[{allTuples, flatMap},",
        f"      allTuples = EnumerateComponentTuples[fHead, {dim}];",
        "      flatMap = Association[Table[allTuples[[j]] -> j - 1, {j, Length[allTuples]}]];",
        "      componentMetadata = Join[componentMetadata, Association[Table[",
        "        Module[{flatIdx, tuple},",
        f"          flatIdx = {comp_var}[[k, 1]];",
        "          tuple = SelectFirst[Keys[flatMap], flatMap[#] == flatIdx &, {}];",
        f'          "{field_name}_" <> ToString[flatIdx] -> <|',
        f'            "head" -> "{field_name}",',
        '            "rank" -> fieldRank,',
        '            "indices" -> tuple',
        "          |>",
        "        ],",
        f"        {{k, Length[{comp_var}]}}",
        "      ]]];",
        "    ],",
        "    (* Scalar field: rank 0, empty indices *)",
        "    componentMetadata = Join[componentMetadata, Association[Table[",
        f'      "{field_name}_" <> ToString[{comp_var}[[k, 1]]] -> <|',
        f'        "head" -> "{field_name}",',
        '        "rank" -> 0,',
        '        "indices" -> {}',
        "      |>,",
        f"      {{k, Length[{comp_var}]}}",
        "    ]]];",
        "  ];",
        "];",
    ]


def _wls_shorthand_cd_tensors(  # noqa: PLR0914, PLR0915
    ctx: _WlsContext,
    dyn_fields: list[dict[str, Any]],
) -> list[str]:
    """Generate shorthand tensor definitions + MakeRule substitutions for CD[field].

    For each dynamical field ``f``, defines a shorthand tensor ``CD1f`` that
    represents ``CD[-a]@f[indices]``, creates a ``MakeRule`` substitution rule,
    and applies the substitution to all EOM variables.

    This follows the supervisor's ``ApplyShorthand`` pattern (commit 4a89164):
    eliminate ALL raw ``CD[-a]@field[...]`` operators from the abstract EOM
    BEFORE component projection, so ExpandScalarWrappers and StaggeredToBasis
    can resolve everything cleanly via ComponentValues.

    The DefTensor and MakeRule are created entirely in Wolfram using xAct's
    introspection (``SlotsOfTensor``, ``SymmetryGroupOfTensor``). No Python
    string parsing of index expressions.
    """
    lines: list[str] = [
        "",
        "(* === Shorthand CD[field] tensors (supervisor's pattern) ===              *)",
        "(* Replace CD[-a]@field[...] → CDfield[-a,...] in abstract EOM before     *)",
        "(* component projection. Enables ToValues to resolve all tensor ops.      *)",
        "(* Suppress index validation: CD[-a]@field[-a,-b] has contracted -a       *)",
        "(* that triggers Validate::repeated → Throw[Null].                        *)",
        "Off[Validate::repeated]; Off[Validate::inhom];",
        "Unprotect[xAct`xTensor`Private`ValidateIndices];",
        "savedValidateIndicesShorthand = xAct`xTensor`Private`ValidateIndices;",
        "xAct`xTensor`Private`ValidateIndices = (True &);",
    ]

    # Collect all shorthand rule variable names
    rule_vars: list[str] = []

    for df in dyn_fields:
        head = df["head"]
        df["fexpr"]  # e.g. "tidalH[-a,-b]" or "tidalT[a,-b,-c]"
        cd1_head = f"CD1{head}"
        rule_var = f"toCD1{head}"
        rule_vars.append(rule_var)

        # Extract the index pattern from fexpr for the MakeRule LHS/RHS
        # The fexpr is like "tidalH[-a, -b]" — we need the indices part
        # Use Wolfram introspection rather than Python parsing
        lines.extend(
            [
                f"(* Shorthand for CD[{df['name']}] *)",
                f"Module[{{slots = SlotsOfTensor[{head}], rank}},",
                "  rank = Length[slots];",
                "  If[rank == 0,",
                "    (* Scalar field: CD[idx]@f[] → CD1f[idx] *)",
                f"    If[!xTensorQ[{cd1_head}], DefTensor[{cd1_head}[-a], {ctx.manifold}]];",
                f"    {rule_var} = {ctx.cd}[idx_][{head}[]] :> {cd1_head}[idx],",
                "    (* Tensor field: CD[idx]@f[args] → CD1f[idx, args] *)",
                f"    If[!xTensorQ[{cd1_head}],",
                "      Module[{dummyIdxs, cdIdxList},",
                f"        dummyIdxs = Table[DummyIn[Tangent{ctx.manifold}], {{n, rank}}];",
                "        cdIdxList = Join[{-a}, MapThread[If[#1 === 1, #2, -#2] &, {slots, dummyIdxs}]];",
                f"        DefTensor[{cd1_head} @@ cdIdxList, {ctx.manifold}]",
                "      ]",
                "    ];",
                f"    {rule_var} = {ctx.cd}[idx_][{head}[args__]] :> {cd1_head}[idx, args]",
                "  ];",
                "];",
                "",
            ]
        )

    # CD shorthand rules are stored in $CDShorthandRules for use by
    # DecomposeScalarExpression during canonical Lagrangian decomposition.
    # EOM application was removed — Component E-L builds equations from
    # the decomposed Lagrangian, not from abstract EOM variables.

    # --- Second-order CD shorthands: CD2field = CD@CD1field ---
    # For R̃² torsion theories: the Lagrangian has ∇∇field terms (second covariant
    # derivatives) wrapped in Scalar[]. After CD1 substitution, these become
    # CD[-a]@CD1field[-b,...] patterns. Define CD2field to absorb these.
    # Ref: supervisor's DDClockField = CD@DClockField (commit 4a89164).
    cd2_rule_vars: list[str] = []
    for df in dyn_fields:
        head = df["head"]
        cd1_head = f"CD1{head}"
        cd2_head = f"CD2{head}"
        cd2_rule_var = f"toCD2{head}"
        cd2_rule_vars.append(cd2_rule_var)

        lines.extend(
            [
                f"(* Second-order shorthand for CD[CD[{df['name']}]] *)",
                f"If[xTensorQ[{cd1_head}],",
                f"  Module[{{cd1Slots = SlotsOfTensor[{cd1_head}], cd1Rank, dummyIdxs, cd2IdxList, cd1Expr, cdIdx}},",
                "    cd1Rank = Length[cd1Slots];",
                f"    If[!xTensorQ[{cd2_head}],",
                f"      dummyIdxs = Table[DummyIn[Tangent{ctx.manifold}], {{n, cd1Rank}}];",
                "      cd2IdxList = Join[{-a}, MapThread[If[#1 === 1, #2, -#2] &, {cd1Slots, dummyIdxs}]];",
                f"      DefTensor[{cd2_head} @@ cd2IdxList, {ctx.manifold}]",
                "    ];",
                "    (* Direct structural rule: CD[idx]@CD1field[args] → CD2field[idx, args] *)",
                f"    {cd2_rule_var} = {ctx.cd}[idx_][{cd1_head}[args__]] :> {cd2_head}[idx, args]",
                "  ],",
                f"  {cd2_rule_var} = {{}};",
                "];",
                "",
            ]
        )

    # --- Third-order CD shorthands: CD3field = CD@CD2field ---
    # For R̃² terms with three nested covariant derivatives.
    cd3_rule_vars: list[str] = []
    for df in dyn_fields:
        head = df["head"]
        cd2_head = f"CD2{head}"
        cd3_head = f"CD3{head}"
        cd3_rule_var = f"toCD3{head}"
        cd3_rule_vars.append(cd3_rule_var)

        lines.extend(
            [
                f"(* Third-order shorthand for CD[CD[CD[{df['name']}]]] *)",
                f"If[xTensorQ[{cd2_head}],",
                f"  Module[{{cd2Slots = SlotsOfTensor[{cd2_head}], cd2Rank, dummyIdxs, cd3IdxList}},",
                "    cd2Rank = Length[cd2Slots];",
                f"    If[!xTensorQ[{cd3_head}],",
                f"      dummyIdxs = Table[DummyIn[Tangent{ctx.manifold}], {{n, cd2Rank}}];",
                "      cd3IdxList = Join[{-a}, MapThread[If[#1 === 1, #2, -#2] &, {cd2Slots, dummyIdxs}]];",
                f"      DefTensor[{cd3_head} @@ cd3IdxList, {ctx.manifold}]",
                "    ];",
                f"    {cd3_rule_var} = {ctx.cd}[idx_][{cd2_head}[args__]] :> {cd3_head}[idx, args]",
                "  ],",
                f"  {cd3_rule_var} = {{}};",
                "];",
                "",
            ]
        )

    # --- Fourth-order CD shorthands: CD4field = CD@CD3field ---
    # For R̃² terms with four nested covariant derivatives (b5 ∇⁴h).
    # Only natural-placement ComponentValues computed (flat Minkowski: SeparateMetric
    # handles raising via explicit eta factors, so ToValues only needs natural).
    cd4_rule_vars: list[str] = []
    for df in dyn_fields:
        head = df["head"]
        cd3_head = f"CD3{head}"
        cd4_head = f"CD4{head}"
        cd4_rule_var = f"toCD4{head}"
        cd4_rule_vars.append(cd4_rule_var)

        lines.extend(
            [
                f"(* Fourth-order shorthand for CD⁴[{df['name']}] *)",
                f"If[xTensorQ[{cd3_head}],",
                f"  Module[{{cd3Slots = SlotsOfTensor[{cd3_head}], cd3Rank, dummyIdxs, cd4IdxList}},",
                "    cd3Rank = Length[cd3Slots];",
                f"    If[!xTensorQ[{cd4_head}],",
                f"      dummyIdxs = Table[DummyIn[Tangent{ctx.manifold}], {{n, cd3Rank}}];",
                "      cd4IdxList = Join[{-a}, MapThread[If[#1 === 1, #2, -#2] &, {cd3Slots, dummyIdxs}]];",
                f"      DefTensor[{cd4_head} @@ cd4IdxList, {ctx.manifold}]",
                "    ];",
                f"    {cd4_rule_var} = {ctx.cd}[idx_][{cd3_head}[args__]] :> {cd4_head}[idx, args]",
                "  ],",
                f"  {cd4_rule_var} = {{}};",
                "];",
                "",
            ]
        )

    # Store ALL rules in global $CDShorthandRules for reuse in ComponentDecompose.wl
    # (after ExpandScalarWrappers introduces new CD operators from Scalar contents)
    all_rule_vars = rule_vars + cd2_rule_vars + cd3_rule_vars + cd4_rule_vars
    if all_rule_vars:
        rules_list = ", ".join(all_rule_vars)
        lines.append(f"$CDShorthandRules = {wl_flatten_list(rules_list)};")
    else:
        lines.append("$CDShorthandRules = {};")

    # Generate REVERSE rules: CD1field[idx, args] → CD[idx][field[args]]
    # These undo the shorthand substitution before DecomposeToComponents,
    # because StaggeredToBasis can handle raw CD operators (via ToBasis +
    # CD→PD→Derivative chain) but NOT CD shorthand tensors (which cause
    # recursion in ToBasis and bypass isCDlikeQ in the force-resolve loop).
    reverse_rules: list[str] = []
    for df in dyn_fields:
        head = df["head"]
        cd1_head = f"CD1{head}"
        cd2_head = f"CD2{head}"
        cd3_head = f"CD3{head}"
        cd4_head = f"CD4{head}"
        # Reverse rules: CDNfield → CD@CD(N-1)field (highest first for ReplaceRepeated)
        reverse_rules.extend(
            [
                f"If[xTensorQ[{cd4_head}], {cd4_head}[idx_, args__] :> {ctx.cd}[idx][{cd3_head}[args]], {{}}]",
                f"If[xTensorQ[{cd3_head}], {cd3_head}[idx_, args__] :> {ctx.cd}[idx][{cd2_head}[args]], {{}}]",
                f"If[xTensorQ[{cd2_head}], {cd2_head}[idx_, args__] :> {ctx.cd}[idx][{cd1_head}[args]], {{}}]",
            ]
        )
        # CD1field[idx, args] → CD[idx][field[args]]
        if df.get("fexpr", "").endswith("[]"):
            # Scalar field: CD1f[idx] → CD[idx][f[]]
            reverse_rules.append(
                f"If[xTensorQ[{cd1_head}], {cd1_head}[idx_] :> {ctx.cd}[idx][{head}[]], {{}}]"
            )
        else:
            reverse_rules.append(
                f"If[xTensorQ[{cd1_head}], {cd1_head}[idx_, args__] :> {ctx.cd}[idx][{head}[args]], {{}}]"
            )
    if reverse_rules:
        rev_list = ", ".join(reverse_rules)
        lines.append(f"$CDShorthandReverseRules = {wl_flatten_list(rev_list)};")
    else:
        lines.append("$CDShorthandReverseRules = {};")

    lines.extend(
        [
            "(* Restore index validation after shorthand substitution *)",
            "xAct`xTensor`Private`ValidateIndices = savedValidateIndicesShorthand;",
            "Protect[xAct`xTensor`Private`ValidateIndices];",
            "On[Validate::repeated]; On[Validate::inhom];",
            'Print["Shorthand CD substitution complete"];',
            "",
        ]
    )
    return lines


def _wls_precompute_cd_component_values(
    ctx: _WlsContext,
    dyn_fields: list[dict[str, Any]],
) -> list[str]:
    """Pre-compute ComponentValues for CD shorthand tensors (supervisor's pattern).

    Follows the approach in ``SphericalEuclidean.m`` (lines 265-395): for each
    CD shorthand tensor (CD1field through CD4field), compute component arrays
    via ``StaggeredToBasis`` and assign them as ``ComponentValue`` entries so
    that ``ToValues`` resolves them in O(1) during Phase A decomposition.

    Each level builds iteratively on the previous level's ComponentValues:
    CD1field from base field → CD2field from CD1field → etc.  This means the
    cost is O(dim^2) per level, not O(dim^{2K}) for K nested CDs.

    Uses xAct introspection (``SlotsOfTensor``, ``ComponentArray``) throughout,
    making it theory-agnostic: works for any tensor rank, any field count, any
    coupling structure.

    References
    ----------
    - Supervisor's ``SphericalEuclidean.m``, lines 265-395 (commit 4a89164)
    - Supervisor's ``VerifyComponents`` helper, line 253-263
    """
    chart = ctx.chart
    manifold = ctx.manifold
    cd = ctx.cd
    lines: list[str] = [
        "",
        "(* === Pre-compute CD shorthand ComponentValues (supervisor's pattern) ===  *)",
        "(* For each CDN-field tensor, compute component arrays via StaggeredToBasis *)",
        "(* and assign as ComponentValue entries.  ToValues then resolves them O(1)  *)",
        "(* during Phase A canonical decomposition, reducing TraceBasisDummy K from  *)",
        "(* ~4+ contracted pairs to ~0.  Each level builds on the previous.         *)",
        "(* Ref: supervisor's SphericalEuclidean.m lines 265-395.                   *)",
        "(* Theory-agnostic: uses SlotsOfTensor introspection for any field rank.   *)",
        "Off[Validate::repeated]; Off[Validate::inhom];",
        _wls_timing_start("tCDPrecomp"),
        "",
        "(* Helper: Splinter pipeline for rank-N tensors (supervisor's pattern).    *)",
        "(* Unlike StaggeredToBasis (designed for scalar expressions), this         *)",
        "(* inserts ComponentArray BETWEEN ToBasis and ToValues to produce a        *)",
        "(* component ARRAY for rank-N tensors with free indices.                  *)",
        "(* Ref: supervisor's EuclideanSplinter (SphericalEuclidean.m:239-250).    *)",
        "tidalSplinter[expr_, chart2_] := Module[{e = expr, metric2, covd2, christPD2},",
        "  Module[{manifold2, tb2, ms2},",
        "    manifold2 = ManifoldOfChart[chart2];",
        '    tb2 = Symbol["Tangent" <> ToString[manifold2]];',
        "    ms2 = MetricsOfVBundle[tb2];",
        "    metric2 = First[ms2]; covd2 = CovDOfMetric[metric2];",
        "  ];",
        '  christPD2 = Symbol["Christoffel" <> ToString[covd2] <> "PD" <> ToString[chart2]];',
        "  (* Staggered ToBasis: peel covariant derivative layers *)",
        "  e = ToBasis[chart2][e];",
        "  e = e /. christPD2 -> Zero;",
        "  e = ToBasis[chart2][e];",
        "  e = e /. christPD2 -> Zero;",
        "  e = SeparateMetric[metric2][e];",
        "  e = ToBasis[chart2][e];",
        "  e = e /. christPD2 -> Zero;",
        "  (* ComponentArray: expand free basis indices into array *)",
        "  e = ComponentArray[e];",
        "  (* ToValues: substitute pre-computed ComponentValues *)",
        "  e = ToValues[e];",
        "  e = TraceBasisDummy[e];",
        "  e = ToValues[e];",
        "  e",
        "];",
        "",
    ]

    # Only pre-compute ComponentValues for MULTI-FIELD theories.
    # Single-field theories (graviton_torsion, scalar_field) don't have
    # cross-field coupling → StaggeredToBasis handles them efficiently.
    # Multi-field theories (torsion+EM) have cross-coupling that creates
    # O(dim^{2K}) blowup in TraceBasisDummy — pre-computation prevents this.
    #
    # Limit to CD1 only: each CD level adds a rank to the tensor, making
    # the ComponentArray exponentially larger (4^rank). CD1 resolves first
    # covariant derivatives (the most common), giving the best cost/benefit.
    # Higher levels are handled by StaggeredToBasis fallback during Phase A.
    min_fields_for_precompute = 2
    if len(dyn_fields) < min_fields_for_precompute:
        lines.extend(
            [
                "",
                "(* Single-field theory — skip CD pre-computation (already fast) *)",
                'Print["CD pre-computation: skipped (single dynamical field)"];',
            ]
        )
        lines.extend(
            [
                "",
                "On[Validate::repeated]; On[Validate::inhom];",
                _wls_timing_end(
                    "tCDPrecomp", "CD shorthand ComponentValue pre-computation"
                ),
                "",
            ]
        )
        return lines

    max_cd_precompute = 2
    for df in dyn_fields:
        head = df["head"]
        field_name = df["name"]
        # Field rank determines CD shorthand tensor rank: CD_N has rank = base_rank + N
        # Limit pre-computation to avoid dim^rank blowup in ComponentArray:
        # rank ≤ 4 → ≤ 256 entries (4D), fast; rank 5 → 1024 entries, slow
        field_dict = df.get("field", {})
        ftype = field_dict.get("type", df.get("type", "scalar"))
        if ftype == "scalar":
            base_rank = 0
        elif ftype == "vector":
            base_rank = 1
        else:
            base_rank = field_dict.get("rank", 2)
        # For this field, limit CD level so CDN tensor rank ≤ 4 (≤ 256 entries in 4D).
        # CD2 for rank-3 torsion (rank 5 = 1024 entries) tested: costs 55s precompute
        # AND makes per-term SplinterToArray slower (41-67s vs 12-34s without).
        # The dense CD2T ComponentValues cause TraceBasisDummy to enumerate MORE
        # combinations, not fewer. Only CD1 is beneficial for torsion.
        field_max_cd = min(max_cd_precompute, max(1, 4 - base_rank))

        for cd_level in range(1, field_max_cd + 1):
            prev_head = f"CD{cd_level - 1}{head}" if cd_level > 1 else head
            cd_head = f"CD{cd_level}{head}"

            lines.extend(
                [
                    "",
                    f"(* Pre-compute CD{cd_level}{field_name} ComponentValues *)",
                    f"If[xTensorQ[{cd_head}],",
                    "  Module[{prevSlots, prevDownIdxs, cdExpr, compDown, compUp},",
                    f"    prevSlots = SlotsOfTensor[{prev_head}];",
                    "",
                    "    (* Compute natural-down via tidalSplinter (produces component ARRAY) *)",
                    f"    prevDownIdxs = Table[-DummyIn[Tangent{manifold}], {{Length[prevSlots]}}];",
                    f"    cdExpr = {cd}[-a][{prev_head} @@ prevDownIdxs];",
                    f"    compDown = tidalSplinter[cdExpr, {chart}];",
                    "    (* Simplify components immediately after projection (supervisor's pattern). *)",
                    "    (* Removes algebraic junk before components feed into every downstream term. *)",
                    f"    compDown = Map[Simplify, compDown, {{Length[SlotsOfTensor[{cd_head}]]}}];",
                ]
            )
            # Apply plane-wave reduction to CD pre-computed components:
            # Zero entries where the FIRST index (the CD derivative direction)
            # is a killed (transverse) coordinate. This makes transverse-derivative
            # terms resolve to 0 immediately via ToValues during decomposition.
            if ctx.reduction is not None:
                killed_slots = [
                    ctx.coords.index(c)
                    for c in ctx.coords[1:]
                    if c != ctx.reduction["propagation_axis"]
                ]
                lines.extend(
                    f"    compDown[[{ks + 1}]] = 0 * compDown[[{ks + 1}]];"
                    for ks in killed_slots
                )
                lines.append(
                    '    Print["    plane-wave: zeroed transverse CD components"];'
                )
            lines.extend(
                [
                    "    (* Assign natural-down ComponentValues *)",
                    f"    Module[{{downIdxs = Table[-{{IndexRange[a, z][[i]], {chart}}}, "
                    f"{{i, Length[SlotsOfTensor[{cd_head}]]}}]}},",
                    "      Block[{Print = Null},",
                    f"        ComponentValue[ComponentArray[{cd_head} @@ downIdxs], compDown]",
                    "      ]",
                    "    ];",
                    "",
                    "    (* Compute all-up by raising (uses metric ComponentValues) *)",
                    f"    Module[{{upIdxs = Table[{{IndexRange[a, z][[i]], {chart}}}, "
                    f"{{i, Length[SlotsOfTensor[{cd_head}]]}}]}},",
                    f"      If[DiagonalMatrixQ[{ctx.prefix}MetricMatrix],",
                    "        (* Fast path: diagonal metric — raise by element-wise multiplication *)",
                    "        (* with inverse metric diagonal.  Avoids full tidalSplinter call.    *)",
                    "        (* Ref: supervisor only computes needed placements (line 379).       *)",
                    f"        Module[{{invDiag = Diagonal[Inverse[{ctx.prefix}MetricMatrix]], rank = Length[SlotsOfTensor[{cd_head}]]}},",
                    "          compUp = compDown * Outer[Times, Sequence @@ Table[invDiag, rank]];",
                    "          compUp = Map[Simplify, compUp, {rank}];",
                    '          Print["    fast-raise: diagonal metric, skipped tidalSplinter"];',
                    "        ],",
                    "        (* General path: full tidalSplinter for non-diagonal metric *)",
                    f"        compUp = tidalSplinter[{cd_head} @@ upIdxs, {chart}];",
                    f"        compUp = Map[Simplify, compUp, {{Length[SlotsOfTensor[{cd_head}]]}}];",
                    "      ];",
                ]
            )
            if ctx.reduction is not None:
                killed_slots = [
                    ctx.coords.index(c)
                    for c in ctx.coords[1:]
                    if c != ctx.reduction["propagation_axis"]
                ]
                lines.extend(
                    f"      compUp[[{ks + 1}]] = 0 * compUp[[{ks + 1}]];"
                    for ks in killed_slots
                )
            lines.extend(
                [
                    "      Block[{Print = Null},",
                    f"        ComponentValue[ComponentArray[{cd_head} @@ upIdxs], compUp]",
                    "      ]",
                    "    ];",
                    "",
                    f'    Print["  {cd_head}: ", Length[SlotsOfTensor[{cd_head}]], '
                    f'" indices, ComponentValues assigned"];',
                    "  ],",
                    f"  Null (* {cd_head} not defined for this theory *)",
                    "];",
                ]
            )

    lines.extend(
        [
            "",
            "On[Validate::repeated]; On[Validate::inhom];",
            _wls_timing_end(
                "tCDPrecomp", "CD shorthand ComponentValue pre-computation"
            ),
            "",
        ]
    )
    return lines


def _wls_multi_field_eom(
    ctx: _WlsContext,
    dyn_fields: list[dict[str, Any]],
) -> list[str] | None:
    """Initialize empty fieldEquations and set up CD shorthand tensors.

    The actual equations are computed later by ComponentEulerLagrange
    (after Phase A decomposes the Lagrangian to component form).
    This function sets up the CD shorthand tensors (CD1field, CD2field, etc.)
    whose $CDShorthandRules are used by DecomposeScalarExpression.
    """
    lines: list[str] = []

    if ctx.lagrangian_expr:
        # CD shorthand setup: define CDN-field tensors and substitution rules.
        # Ref: supervisor's ApplyShorthand pattern (commit 4a89164).
        lines.extend(_wls_shorthand_cd_tensors(ctx, dyn_fields))
        # Pre-compute ComponentValues for all CD shorthand tensors so that
        # ToValues resolves them in O(1) during Phase A decomposition.
        # This is the key optimization from the supervisor's pattern.
        lines.extend(_wls_precompute_cd_component_values(ctx, dyn_fields))
        lines.extend(
            [
                'Print["[TIMING] CD shorthand setup + pre-computation: ", '
                'Round[AbsoluteTime[] - tSubPhase, 0.1], " seconds"];',
                "tSubPhase = AbsoluteTime[];",
            ]
        )
        lines.extend(
            [
                "(* Initialize empty — Component E-L builds equations after Phase A *)",
                "fieldEquations = {};",
                "componentMetadata = <||>;",
                'Print["Equations from component E-L (after Phase A decomposition)"];',
            ]
        )
        return lines
    return None


def _wls_matter_pert_truncation(mpi: dict[str, str]) -> list[str]:
    """Generate LI[2] drop + LI[1] replacement + background substitution.

    After xPert perturbation, L^(2) contains:
    - ``pertSym[LI[2], ...]``: 2nd-order perturbation (dropped)
    - ``pertSym[LI[1], ...]``: 1st-order perturbation (→ perturbation field)
    - Original field symbol (e.g., ``geA``): zeroth-order = background

    The original field symbol must be replaced with the background field
    so that ComponentValues from ``[[background_fields]]`` apply during
    component decomposition.
    """
    mp_sym, mp_head = mpi["pert_sym"], mpi["pert_head"]
    field_head, bg_head = mpi["field_head"], mpi.get("bg_head", "")
    pname = mpi["pert_name"]
    lines: list[str] = []
    if mpi["field_type"] == "scalar":
        lines.extend(
            [
                f"(* Drop 2nd-order matter perturbation {pname}^(2) *)",
                f"l2Raw = l2Raw /. {mp_sym}[LI[2]] :> 0;",
                f"(* Replace xPert notation -> perturbation field {pname} *)",
                f"l2Raw = l2Raw /. {mp_sym}[LI[1]] :> {mp_head}[];",
            ]
        )
    else:
        lines.extend(
            [
                f"(* Drop 2nd-order matter perturbation {pname}^(2) *)",
                f"l2Raw = l2Raw /. {mp_sym}[LI[2], idx__] :> 0;",
                f"(* Replace xPert notation -> perturbation field {pname} *)",
                f"l2Raw = l2Raw /. {mp_sym}[LI[1], idx__] :> {mp_head}[idx];",
            ]
        )
    # Replace original field -> background so ComponentValues evaluate
    if bg_head and bg_head != field_head:
        lines.extend(
            [
                f"(* Replace zeroth-order field {field_head} -> background {bg_head} *)",
                f"l2Raw = l2Raw /. {field_head} -> {bg_head};",
            ]
        )
    return lines


def _wls_torsion_curvature_decomposition(ctx: _WlsContext) -> list[str]:
    """Decompose Riemann-Cartan curvature tensors into LC curvature + torsion.

    When a ``[torsion]`` section is present and the Lagrangian contains CDT
    curvature tensors (``RicciScalarCDT[]``, ``RicciCDT[-a,-b]``, etc.),
    this function generates Wolfram code that:

    1. Applies ``ChangeCurvature[L, CDT, CD]`` — xAct's built-in identity that
       decomposes R̃ curvature into R^{LC} curvature + contortion (K) terms.

    2. Substitutes contortion → torsion using the algebraic identity:
       ``K^a_{bc} = ½(T^a_{bc} + T_b^{a}_{c} - T_{bc}^{a})``
       where K = ChristoffelCDCDT (the connection difference CDT - CD).

    The result contains only ``RicciCD``, ``CD[...][TorsionCDT[...]]``, and
    ``TorsionCDT[...]`` — all of which have xPert perturbation rules and can
    be linearized via ``Perturbation[L, 2]`` + ``ExpandPerturbation``.

    **This decomposition is an exact identity** (not an approximation).  It
    is unique by the fundamental theorem of Riemannian geometry: every metric
    determines a unique torsion-free, metric-compatible connection (the
    Levi-Civita connection CD).

    If no CDT curvature tensors are present in the Lagrangian, this function
    is a no-op — the expression passes through unchanged.

    Returns
    -------
    list[str]
        Wolfram code lines.  Empty if no torsion or no CDT curvature detected.

    References
    ----------
    - Shapiro (2002), hep-th/0103093, eq. 2.9 — curvature decomposition
    - Hehl et al. (1976), Rev. Mod. Phys. 48:393 — PGT foundations
    - Sezgin & van Nieuwenhuizen (1980), Phys. Rev. D 21:3269 — linearized PGT
    """
    if not ctx.torsion:
        return []

    p = ctx.prefix
    cdt = ctx.cdt
    cd = f"{p}CD"
    metric = ctx.metric

    # The Christoffel name is auto-generated by xAct as "Christoffel" + CD + CDT
    christoffel_name = f"Christoffel{cd}{cdt}"

    # The TorsionCDT name is auto-generated by xAct
    torsion_cdt = f"Torsion{cdt}"

    return [
        "",
        "(* === Torsion curvature decomposition: R̃ → R + T ===                   *)",
        "(* Poincaré gauge theory Lagrangians may use the Riemann-Cartan          *)",
        "(* curvature R̃ (from the torsionful connection CDT).  xPert cannot      *)",
        "(* perturb R̃ directly because the underlying contortion K = CDT - CD    *)",
        "(* is a derived xAct object without perturbation rules.                  *)",
        "(*                                                                        *)",
        "(* Solution: two-step exact decomposition.                                *)",
        "(*   Step 1: ChangeCurvature[L, CDT, CD] — xAct built-in identity       *)",
        "(*     R̃^a_{bcd} = R^a_{bcd} + ∇_c K^a_{db} - ∇_d K^a_{cb}            *)",
        "(*                  + K^a_{ce} K^e_{db} - K^a_{de} K^e_{cb}              *)",
        "(*   Step 2: Contortion → torsion substitution                           *)",
        "(*     K^a_{bc} = ½(T^a_{bc} + T_b^a_c - T_{bc}^a)                      *)",
        "(*                                                                        *)",
        "(* After decomposition, the Lagrangian contains only RicciCD (LC Ricci), *)",
        "(* CD[...][TorsionCDT[...]] (LC derivatives of torsion), TorsionCDT      *)",
        "(* (torsion tensor), and the metric — all perturbable by xPert.          *)",
        "(*                                                                        *)",
        "(* This identity is exact and unique (fundamental theorem of Riemannian   *)",
        "(* geometry).  It does not change the physics — only rewrites L in a form *)",
        "(* that the perturbation pipeline can handle.                             *)",
        "(*                                                                        *)",
        "(* Ref: Shapiro (2002, hep-th/0103093, eq. 2.9)                          *)",
        "(*      Hehl et al. (1976, Rev. Mod. Phys. 48:393)                       *)",
        f"If[!FreeQ[lOriginal, RicciScalar{cdt} | Ricci{cdt} | Riemann{cdt} | Einstein{cdt}],",
        '  Print["Decomposing Riemann-Cartan curvature R̃ → Levi-Civita R + torsion T..."];',
        f"  lOriginal = ChangeCurvature[lOriginal, {cdt}, {cd}];",
        "",
        "  (* Contortion → torsion: K^a_bc = ½(T^a_bc} + T_b^a_c - T_bc^a)      *)",
        "  (* The contortion is a derived xAct tensor that DefTensorPerturbation  *)",
        "  (* cannot perturb.  We must express it in terms of TorsionCDT (which  *)",
        "  (* IS perturbable) before xPert expansion.                             *)",
        "  (* Use xAct's WithDummies to generate proper dummy indices.            *)",
        f"  lOriginal = lOriginal /. {christoffel_name}[a1_, a2_, a3_] :>",
        f"    With[{{d1 = DummyIn[VBundleOfMetric[{metric}]],",
        f"           d2 = DummyIn[VBundleOfMetric[{metric}]]}},",
        f"      1/2 ({torsion_cdt}[a1, a2, a3] +",
        f"           {metric}[a2, d1] {metric}[a1, d2] {torsion_cdt}[-d1, -d2, a3] -",
        f"           {metric}[a3, d1] {metric}[a1, d2] {torsion_cdt}[a2, -d1, -d2])",
        "    ];",
        "  lOriginal = Expand[lOriginal];",
        "",
        f'  Print["R̃ decomposition complete: contains RicciCD=", !FreeQ[lOriginal, Ricci{cd}],',
        f'    ", TorsionCDT=", !FreeQ[lOriginal, {torsion_cdt}],',
        f'    ", ChristoffelCDCDT=", !FreeQ[lOriginal, {christoffel_name}]];',
        f"  If[!FreeQ[lOriginal, {christoffel_name}],",
        f'    Print["WARNING: residual contortion ({christoffel_name}) detected after decomposition"]',
        "  ];",
        '  Print["Decomposed Lagrangian: ", Short[lOriginal, 3]];',
        "];",
        'Print["[TIMING] R̃ decomposition: ", Round[AbsoluteTime[] - tSubPhase, 0.1], " seconds"];',
        "tSubPhase = AbsoluteTime[];",
        "",
    ]


def _wls_linearize_from_lagrangian(  # noqa: C901, PLR0912, PLR0914, PLR0915
    ctx: _WlsContext,
    *,
    include_bg: bool = False,
) -> list[str]:
    """Lagrangian-first linearization: L → L^(2) → EOM + canonical from L^(2).

    Single-path approach using xPert's 2nd-order perturbation:

    1. Multiply Lagrangian by ``√(-g)`` to form the action density, then
       perturb to second order: ``δ²(√|g| L)``.  xPert natively handles
       ``Perturbation[Sqrt[-Detg[]], n]`` via the metric determinant symbol
       created by ``DefMetric``.  After expansion, divide out the background
       ``√|g₀|`` (covariantly constant for Levi-Civita connection).

    2. ``L^(2) = δ²(√|g| L) / (2 √|g₀|)`` → quadratic Lagrangian.

    3. Expand ``Scalar[x]^n`` → ``∏ Scalar[RenameDummies[x]]`` so that VarD
       can vary through each copy independently (fixes the index-collision
       problem from Fierz-Pauli trace-squared terms).

    4. ``VarD[H[-a,-b], CD][L^(2)]`` → correct linearized EOM for each
       dynamical field.

    5. Same L^(2) serves the canonical pipeline (momenta π, Hamiltonian H).

    Supports multi-field perturbation via ``[[linearization.matter_perturbations]]``:
    when present, each matter field is split via ``DefTensorPerturbation`` into
    background + perturbation, and ``VarD`` is called for each dynamical field.
    Cross-coupling terms (e.g., graviton-photon via background F̄) emerge
    automatically from the algebra.

    Works for any background metric.  For flat Minkowski ``√|g₀| = 1`` and
    the volume element factor is trivial.  For curved backgrounds, xPert
    correctly captures all ``δ(√|g|)`` contributions.

    Raises
    ------
    ValueError
        If ``ctx.linearization`` is ``None``.
    RuntimeError
        If linearization produces an unexpected number of equations.
    """
    if ctx.linearization is None:
        msg = "_wls_linearize_from_lagrangian called without linearization config"
        raise ValueError(msg)
    lin = ctx.linearization
    pert_field_name = lin.get("perturbation_field")  # May be None for matter-only
    matter_perts: list[dict[str, Any]] = lin.get("matter_perturbations", [])

    # Resolve the metric perturbation field (if any)
    has_metric_pert = pert_field_name is not None
    pert_field: dict[str, Any] = {}
    fexpr: str = ""
    if has_metric_pert:
        pert_field = next(f for f in ctx.fields if f["name"] == pert_field_name)
        fexpr = _field_expression(pert_field, ctx.prefix)

    pert_sym = f"{ctx.prefix}hpert"
    eps_sym = f"{ctx.prefix}Epsilon"
    field_head = (
        f"{ctx.prefix}{pert_field_name.capitalize()}" if has_metric_pert else None
    )
    bg_name = f"{ctx.prefix}Bg"
    ricci_sym = f"Ricci{ctx.cd}"
    ricci_scalar_sym = f"RicciScalar{ctx.cd}"
    riemann_sym = f"Riemann{ctx.cd}"
    einstein_sym = f"Einstein{ctx.cd}"

    p = ctx.prefix
    det_sym = f"Det{ctx.metric}"

    # ------------------------------------------------------------------
    # Step 1: L^(2) = δ²(√|g| L) / (2 √|g₀|) via xPert
    # ------------------------------------------------------------------
    lines: list[str] = [
        "",
        "(* ============================================================ *)",
        "(* Lagrangian-first linearization (single-path via L^(2))       *)",
        "(* L^(2) = d^2(sqrt|g| L) / (2 sqrt|g0|)  via xPert           *)",
        "(* Then: VarD[field, CD][L^(2)] -> linearized EOM              *)",
        "(* Same L^(2) feeds canonical pipeline (pi, H)                  *)",
        "(* ============================================================ *)",
        "",
        "(* Save original nonlinear Lagrangian *)",
        _wls_timing_start("tLinearize"),
        "tSubPhase = AbsoluteTime[];",
        f"lOriginal = {p}Lagrangian;",
        "",
    ]

    # --- Torsion: automatic R̃ decomposition ---
    # For PGT theories, the user may write the Lagrangian using the
    # Riemann-Cartan curvature tensors (RicciScalarCDT[], RicciCDT[-a,-b],
    # etc.) which are the natural objects in Poincaré gauge theory.
    #
    # However, xPert cannot perturb R̃ directly — it treats RicciScalarCDT[]
    # as opaque because the underlying connection CDT involves contortion
    # (Christoffel[CD,CDT]), a derived xAct object without perturbation rules.
    #
    # Solution: automatic two-step decomposition:
    #   1. ChangeCurvature[L, CDT, CD]: decomposes R̃ curvature into
    #      R^{LC} + contortion (K) terms using xAct's built-in identity
    #   2. Contortion → torsion substitution: replaces ChristoffelCDCDT
    #      with explicit TorsionCDT expressions using the algebraic identity
    #      K^a_{bc} = ½(T^a_{bc} + T_b^a_c - T_{bc}^a)
    #
    # The result contains only RicciCD, CD[...][TorsionCDT[...]], and
    # TorsionCDT[...] — all of which are perturbable by xPert.
    #
    # This decomposition is an exact identity (not an approximation) and is
    # unique by the fundamental theorem of Riemannian geometry.
    #
    # Ref: Shapiro (2002, hep-th/0103093, eq. 2.9);
    #      Hehl et al. (1976, Rev. Mod. Phys. 48:393)
    lines.extend(_wls_torsion_curvature_decomposition(ctx))

    # --- Metric perturbation setup ---
    # xPert requires SetupMetricPerturbation to register the global perturbation
    # parameter, even for matter-only perturbation. Without it, Perturbation[]
    # and ExpandPerturbation[] cannot expand Det[g] or matter fields.
    lines.extend(
        [
            "(* Set up metric perturbation parameter via xPert *)",
            f"{pert_sym}Tensor = SetupMetricPerturbation[{ctx.metric}, {pert_sym}, {eps_sym}];",
        ]
    )
    if has_metric_pert:
        lines.append(
            f'Print["Perturbation: {ctx.metric} -> {ctx.metric} + {eps_sym} * {pert_field_name}"];',
        )
    else:
        lines.append(
            'Print["Metric perturbation parameter registered (matter-only; h terms will be dropped)"];',
        )
    lines.append("")

    # --- Matter field perturbation setup (multi-field xPert) ---
    matter_pert_info, mp_lines = _wls_matter_perturbation_setup(
        ctx,
        matter_perts,
        eps_sym,
    )
    lines.extend(mp_lines)

    # --- Torsion perturbation (automatic for PGT) ---
    # When torsion = true, auto-register DefTensorPerturbation for TorsionCDT.
    # TorsionCDT[a, -b, -c] is auto-created by DefCovD[CDT, Torsion->True]
    # with antisymmetry in (-b, -c).  Background torsion = 0 (flat Minkowski).
    # The perturbation tensor inherits the same index structure and symmetry.
    if ctx.torsion:
        torsion_xact_head = f"Torsion{ctx.cdt}"  # xAct auto-created
        torsion_pert_name = ctx.torsion["perturbation_name"]
        torsion_field_head = f"{p}{torsion_pert_name.capitalize()}"
        torsion_pert_label = f"{p}{torsion_pert_name}Pert"
        torsion_print_as = torsion_pert_name
        lines.extend(
            [
                "(* === Torsion perturbation (from [torsion] section) === *)",
                f"(* Perturbation field: {torsion_pert_name}[a,-b,-c] for VarD/decomposition *)",
                f"(* xPert label: {torsion_pert_name}Pert[LI[n],a,-b,-c]                     *)",
                "(* Connected to TorsionCDT via DefTensorPerturbation                       *)",
                f"If[!xTensorQ[{torsion_field_head}],",
                f"  DefTensor[{torsion_field_head}[a, -b, -c], {ctx.manifold}, "
                f'Antisymmetric[{{-b, -c}}], PrintAs -> "{torsion_print_as}"]',
                "];",
                f"DefTensorPerturbation[{torsion_pert_label}[LI[order], a, -b, -c], "
                f"{torsion_xact_head}[a, -b, -c], {ctx.manifold}];",
                f'Print["Torsion perturbation: {torsion_xact_head} = 0 + {eps_sym} * {torsion_pert_name}"];',
                "",
            ]
        )

    lines.extend(
        [
            "(* Include sqrt(-g) volume element: S = int sqrt(-g) L d^n x      *)",
            "(* xPert natively perturbs Sqrt[-Detg[]], handling all orders.     *)",
            f"lDensity = Sqrt[-{det_sym}[]] * lOriginal;",
            "",
            'Print["[TIMING] Pre-xPert setup: ", Round[AbsoluteTime[] - tSubPhase, 0.1], " seconds"];',
            "tSubPhase = AbsoluteTime[];",
            "(* 2nd-order perturbation of the full Lagrangian density *)",
            "l2Raw = Perturbation[lDensity, 2];",
            "l2Raw = ExpandPerturbation[l2Raw];",
            'Print["[TIMING] xPert Perturbation + ExpandPerturbation: ", Round[AbsoluteTime[] - tSubPhase, 0.1], " seconds"];',
            "tSubPhase = AbsoluteTime[];",
            _wls_mem_print("L^(2) density expanded"),
            "",
            "(* Validate that xPert fully expanded *)",
            "If[!FreeQ[l2Raw, Perturbation],",
            '  Throw["Linearization: ExpandPerturbation did not fully expand L^(2)."]',
            "];",
            "",
            "(* Divide out background volume element sqrt(-g0)                   *)",
            "(* (covariantly constant for Levi-Civita, passes through VarD).     *)",
            f"l2Raw = l2Raw / Sqrt[-{det_sym}[]];",
            "",
            "(* Replace background metric determinant with evaluated value       *)",
            f"l2Raw = l2Raw /. {det_sym}[] -> Det[{p}MetricMatrix];",
            "(* Expand to distribute volume element — deferred for large expressions.*)",
            "(* For theories with rank-2 R̃ coupling (non-minimal), l2Raw has     *)",
            "(* 5000+ deeply nested tensor products; Expand takes minutes.        *)",
            "(* Use LeafCount as a size proxy: large expressions skip Expand and  *)",
            "(* rely on downstream pattern-matching to propagate zeros.           *)",
            "Module[{lc = LeafCount[l2Raw]},",
            '  Print["L^(2) LeafCount: ", lc];',
            "  If[lc > 50000,",
            '    Print["L^(2) volume Expand: SKIPPED (LeafCount=", lc, " > 50000)"],',
            "    l2Raw = Expand[l2Raw];",
            '    Print["L^(2) volume Expand: ", If[Head[l2Raw]===Plus, Length[l2Raw], 1], " terms"]',
            "  ]",
            "];",
            "",
        ]
    )

    # --- Drop LI[2] and replace LI[1] for ALL perturbation fields ---
    #
    # xPert expands each field as:  Φ = Φ̄ + ε·δΦ^(1) + ε²·δΦ^(2) + ...
    #   LI[1] = δΦ^(1)  — the linear perturbation field (dynamical DOF)
    #   LI[2] = δΦ^(2)  — second-order correction (NOT a separate DOF)
    #
    # The second-order action L^(2) = ½ δ²(S) contains:
    #   (δΦ^(1))²  — quadratic products of 1st-order perturbations → KEPT
    #                 (these give the linearized EOM upon variation)
    #   Φ̄ · δΦ^(2) — background x 2nd-order correction → ZERO (Φ̄ = 0 for
    #                 flat background torsion/curvature)
    #   δΦ^(2)      — appears only in L^(1), not in L^(2)
    #
    # Therefore LI[2] → 0 is correct: the physical content (wave equations
    # from quadratic δΦ^(1) x δΦ^(1) products) is entirely in LI[1] terms.
    # This applies identically to metric (h), torsion (t), and matter (a)
    # perturbation fields.
    #
    # Ref: Brizuela, Martin-Garcia, Mena Marugan (2009, Gen.Rel.Grav. 41:2415)
    #      — xPert perturbation theory formalism
    if has_metric_pert:
        lines.extend(
            [
                "(* Drop 2nd-order metric perturbation h^(2) -- keep h^(1)*h^(1) only *)",
                f"l2Raw = l2Raw /. {pert_sym}[LI[2], idx__] :> 0;",
                "(* Expand eagerly for small expressions; skip for large (LeafCount>50000) *)",
                "If[LeafCount[l2Raw] <= 50000, l2Raw = Expand[l2Raw]];",
                'Print["After LI[2] drop: ", If[Head[l2Raw]===Plus, Length[l2Raw], 1], " terms"];',
                "",
                "(* Replace xPert metric notation with declared field tensor *)",
                f"l2Raw = l2Raw /. {pert_sym}[LI[1], idx__] :> {field_head}[idx];",
            ]
        )
    else:
        # Matter-only: SetupMetricPerturbation was called to register the
        # perturbation parameter, but we don't want any metric perturbation
        # terms. Drop ALL orders of h.
        lines.extend(
            [
                "(* Matter-only: drop all metric perturbation orders *)",
                f"l2Raw = l2Raw /. {pert_sym}[LI[_], idx___] :> 0;",
                "l2Raw = Expand[l2Raw];",
            ]
        )

    for mpi in matter_pert_info:
        lines.extend(_wls_matter_pert_truncation(mpi))

    # Torsion perturbation truncation: same pattern as matter perturbations.
    # torsionPert[LI[2], ...] → 0  (drop 2nd-order)
    # torsionPert[LI[1], ...] → TorsionField[...]  (replace with physical field)
    # TorsionCDT[...] → 0  (background torsion = 0 for flat Minkowski)
    if ctx.torsion:
        torsion_pert_name = ctx.torsion["perturbation_name"]
        torsion_pert_label = f"{p}{torsion_pert_name}Pert"
        torsion_field_head = f"{p}{torsion_pert_name.capitalize()}"
        torsion_xact_head = f"Torsion{ctx.cdt}"
        lines.extend(
            [
                f"(* Drop 2nd-order torsion perturbation ({torsion_pert_name}) *)",
                f"l2Raw = l2Raw /. {torsion_pert_label}[LI[2], idx__] :> 0;",
                f"(* Replace xPert notation → physical torsion field {torsion_pert_name} *)",
                f"l2Raw = l2Raw /. {torsion_pert_label}[LI[1], idx__] :> {torsion_field_head}[idx];",
                "(* Background torsion = 0 (flat Minkowski, no torsion source) *)",
                f"l2Raw = l2Raw /. {torsion_xact_head}[__] :> 0;",
                "(* Expand eagerly for small expressions; skip for large (LeafCount>50000) *)",
                "If[LeafCount[l2Raw] <= 50000, l2Raw = Expand[l2Raw]];",
                'Print["After torsion truncation: ", If[Head[l2Raw]===Plus, Length[l2Raw], 1], " terms"];',
            ]
        )

    # Replace bg → metric if bg tensor is used
    if include_bg:
        lines.extend(
            [
                "",
                "(* Replace reference metric bg -> background metric *)",
                f"l2Raw = l2Raw /. {bg_name} -> {ctx.metric};",
            ]
        )

    # Scalar BG substitutions
    lines.extend(_wls_scalar_background_substitution(ctx, "l2Raw"))

    lines.extend(
        [
            "",
            "(* Set ALL background curvature to zero (flat Minkowski in any coordinates) *)",
            "(* Riemann[__]:>0 catches curvature-coupling terms R^{abcd}h_ac h_bd from xPert. *)",
            "(* Without this, spherical metrics leave symbolic RiemannCD unevaluated, *)",
            "(* corrupting equations for tensor components (h_theta_theta etc.). *)",
            f"l2Raw = l2Raw /. {wl_list(f'{riemann_sym}[__] :> 0, {ricci_sym}[__] :> 0, {ricci_scalar_sym}[] :> 0, {einstein_sym}[__] :> 0')};",
            "(* Expand eagerly for small expressions; skip for large (>500 terms) *)",
            "If[!Head[l2Raw] === Plus || Length[l2Raw] <= 500, l2Raw = Expand[l2Raw]];",
            'Print["After background curvature zeroing: ", If[Head[l2Raw]===Plus, Length[l2Raw], 1], " terms"];',
        ]
    )

    # NOTE: Background torsion zeroing (TorsionCDT → 0) is done in the
    # torsion perturbation truncation block above (alongside LI[2]→0 and
    # LI[1]→TorsionField substitutions), following the same pattern as
    # the matter perturbation truncation.

    # Zero background contortion (connection difference Γ̃ - Γ).
    #
    # Physics: The contortion K^a_{bc} is algebraically determined by torsion:
    #   K^a_{bc} = ½(T^a_{bc} + T_b^{a}_c - T_{bc}^{a})
    # For zero background torsion (T̄ = 0), the background contortion K̄ = 0.
    #
    # In the xPert perturbation expansion:
    #   - Perturbation[K, 0] = K̄ = 0  ← ZEROED HERE (bare ChristoffelCDCDT)
    #   - Perturbation[K, 1] = δK(t)   ← PRESERVED (encodes torsion wave)
    #   - Perturbation[K, 2]           ← already dropped by LI[2]→0 truncation
    #
    # This is exactly analogous to zeroing RicciScalarCD[] (flat background)
    # while keeping Perturbation[RicciScalarCD[]] (graviton wave).
    #
    # Technical necessity: Without this, R̃²-decomposed Lagrangians retain
    # Scalar[K·K·η] products through VarD that ExpandScalarWrappers cannot
    # resolve (Christoffel symbols aren't standard xAct tensors — ToBasis
    # and RenameDummies don't work on them).
    if ctx.torsion is not None:
        cd = f"{p}CD"
        cdt = ctx.cdt
        christoffel_cdt_cd = f"Christoffel{cd}{cdt}"
        lines.extend(
            [
                "",
                "(* Zero background contortion: K = CDT Christoffel - CD Christoffel *)",
                "(* For zero background torsion, K^a_{bc} = 0 identically.           *)",
                "(* Perturbation[K] terms already replaced with torsion field T above.*)",
                f"l2Raw = l2Raw /. {christoffel_cdt_cd}[__] :> 0;",
                "l2Raw = Expand[l2Raw];",
                'Print["After zeroing background contortion: ", If[Head[l2Raw]===Plus, Length[l2Raw], 1], " terms"];',
            ]
        )

    lines.extend(
        [
            "",
            "(* === Adaptive canonicalization with cost-based bypass ===               *)",
            "(* ToCanonical + ContractMetric scale super-linearly with expression     *)",
            "(* size. Measure first-batch cost; if >= 5s, SKIP canonicalization       *)",
            "(* entirely — downstream DecomposeScalarExpression handles non-canonical *)",
            "(* expressions correctly via ToBasis + TraceBasisDummy at component      *)",
            "(* level where per-component expressions are much smaller.               *)",
            "(* Ref: issue #201 — non-minimal R̃[μν]F coupling: 5000+ terms, 467s.   *)",
            "Module[{l2Terms, nTerms, canonBatchSize, tFirstBatch,",
            "         firstBatchTime, canonResult, batchStart, batchEnd},",
            "  l2Terms = If[Head[l2Raw] === Plus, List @@ l2Raw, {l2Raw}];",
            "  nTerms = Length[l2Terms];",
            '  Print["Canonicalizing ", nTerms, " L^(2) terms..."];',
            "",
            "  (* Measure first batch cost to decide: canonicalize or bypass *)",
            "  canonBatchSize = Min[20, nTerms];",
            "  tFirstBatch = AbsoluteTime[];",
            f"  canonResult = ContractMetric[ToCanonical[Total[l2Terms[[1 ;; canonBatchSize]]]], {ctx.metric}];",
            "  firstBatchTime = AbsoluteTime[] - tFirstBatch;",
            '  Print["  first batch (", canonBatchSize, " terms): ", Round[firstBatchTime, 0.1], "s"];',
            "",
            "  If[firstBatchTime >= 5.0,",
            "    (* BYPASS: ToCanonical too expensive for this theory.                *)",
            "    (* Component decomposition will handle simplification per-component. *)",
            '    Print["ToCanonical too expensive (", Round[firstBatchTime, 0.1],',
            '      "s for ", canonBatchSize, " terms) — bypassing"];',
            '    Print["  Estimated full cost: ",',
            '      Round[firstBatchTime * nTerms / canonBatchSize, 1], "s"];',
            '    Print["  Component-level decomposition will handle simplification"];',
            "    l2Raw = Expand[l2Raw],",
            "",
            "    (* NORMAL: adaptive per-term canonicalization *)",
            "    If[firstBatchTime > 2.0 && canonBatchSize > 5,",
            "      canonBatchSize = Max[5, Floor[canonBatchSize * 2.0 / firstBatchTime]]",
            "    ];",
            "    If[firstBatchTime < 0.5 && canonBatchSize < nTerms,",
            "      canonBatchSize = Min[100, Floor[canonBatchSize * 2.0 / Max[firstBatchTime, 0.01]]]",
            "    ];",
            '    Print["  adaptive batch size: ", canonBatchSize];',
            "    batchStart = canonBatchSize + 1;",
            "    While[batchStart <= nTerms,",
            "      batchEnd = Min[batchStart + canonBatchSize - 1, nTerms];",
            f"      canonResult += ContractMetric[ToCanonical[Total[l2Terms[[batchStart ;; batchEnd]]]], {ctx.metric}];",
            "      batchStart = batchEnd + 1;",
            "    ];",
            "    l2Raw = canonResult;",
            "  ];",
            "  Clear[l2Terms, canonResult];",
            "];",
            "",
            "(* L^(2) = delta^2(sqrt|g| L) / (2 sqrt|g0|) *)",
            f"{p}Lagrangian = l2Raw / 2;",
            f'Print["L^(2) set: ", If[Head[{p}Lagrangian]===Plus, Length[{p}Lagrangian], 1], " terms"];',
            "",
            "(* Memory cleanup: free perturbation intermediates *)",
            "Clear[l2Raw, lOriginal, lDensity];",
            "Share[];",
            'Print["[TIMING] L^(2) cleanup (volume, LI, background, contortion): ", Round[AbsoluteTime[] - tSubPhase, 0.1], " seconds"];',
            "tSubPhase = AbsoluteTime[];",
            _wls_mem_print("After perturbation cleanup"),
            "",
        ]
    )

    # ------------------------------------------------------------------
    # Step 1b: Type A gauge fixing (if any) — add L_gf to L^(2)
    # Perturbation-level gauge terms (current presets) are already
    # quadratic in the perturbation field, so they are added directly
    # to L^(2) without further xPert expansion.
    # ------------------------------------------------------------------
    if ctx.gauge and any(
        _resolve_gauge_mechanism(g) == "lagrangian_term" for g in ctx.gauge
    ):
        lines.extend(_wls_gauge_fixing_type_a(ctx))

    # ------------------------------------------------------------------
    # Step 2: EOM via VarD for each dynamical field
    # ------------------------------------------------------------------
    lines.extend(
        [
            "(* ------------------------------------------------------------ *)",
            "(* EOM: VarD[field, CD][L^(2)] for each dynamical field         *)",
            "(* ------------------------------------------------------------ *)",
            "",
            "(* Expand Scalar[x]^n into products with renamed dummies.       *)",
            "(* This fixes VarD index collision from Fierz-Pauli (tr h)^2:   *)",
            "(*   Scalar[eta^ab H_ab]^2 -> Scalar[eta^ab H_ab]*Scalar[eta^cd H_cd] *)",
            "(* NOTE: Use /. (single pass) NOT //. — RenameDummies may return    *)",
            "(* the same canonical form for symmetric tensors (e.g. eta[a,b]),    *)",
            "(* causing //. to loop 65536 times before hitting $IterationLimit.   *)",
            f"l2ForVarD = {p}Lagrangian;",
            "",
        ]
    )

    # ------------------------------------------------------------------
    # For TT-gauged fields, impose tracelessness (η^{ab} h_{ab} = 0)
    # in L^(2) before VarD.  This eliminates Scalar[H[a,-a]]^2 trace
    # terms that cause ReplaceRepeated loops and memory blowup during
    # ToBasis decomposition.  Physically correct: TT gauge ⇒ tr h = 0.
    #
    # IMPORTANT: Must run BEFORE the Scalar[x]^n expansion below,
    # because that expansion wraps contents in RenameDummies[] which
    # prevents the pattern Scalar[H[a_,-a_]] from matching.
    # ------------------------------------------------------------------
    tt_fields = {entry["field"] for entry in ctx.gauge if entry["type"] == "tt"}
    for field in ctx.fields:
        if field["name"] in tt_fields and field.get("symmetry") == "symmetric":
            head = f"{p}{field['name'].capitalize()}"
            lines.extend(
                [
                    f"(* TT traceless: set tr({field['name']}) = 0 in L^(2) before VarD *)",
                    f"l2ForVarD = l2ForVarD //. Scalar[{head}[a_, -a_]] :> 0;",
                    f'Print["Imposed tr({field["name"]}) = 0: ", Short[l2ForVarD, 5]];',
                    "",
                ]
            )

    lines.extend(
        [
            "(* Evaluate Scalar[metric] → dimension (constant trace).             *)",
            f"l2ForVarD = l2ForVarD /. Scalar[{ctx.metric}[a_, b_]] :> {ctx.dim};",
            "(* NOTE: Scalar[x]^n terms from xPert's ExpandPerturbation are left   *)",
            "(* as-is. xPert already handled the Leibniz rule on R² during          *)",
            "(* Perturbation[lDensity, 2], producing 2*(δR)² + 2R·δ²R terms.     *)",
            "(* VarD handles Scalar[x]^n via its own chain rule.                   *)",
            "(* Manual expansion with RenameDummies/DummyIn is NOT needed and      *)",
            "(* causes index inhomogeneity problems with R̃² torsion theories.     *)",
            'Print["L^(2) for VarD: ", Short[l2ForVarD, 5]];',
            "",
        ]
    )

    if not matter_pert_info and not ctx.torsion:
        # --- Single-field linearization (metric only, original path) ---
        if not has_metric_pert:
            msg = "Single-field linearization requires metric perturbation"
            raise RuntimeError(msg)
        lines.extend(
            [
                "(* Vary L^(2) with respect to perturbation field H *)",
                f"eomLin = VarD[{field_head}[-a, -b], {ctx.cd}][l2ForVarD];",
                # Term-by-term ToCanonical: xperm segfaults on large VarD sums
                "If[Head[eomLin] === Plus,",
                "  eomLin = Total[ToCanonical /@ List @@ eomLin],",
                "  eomLin = ToCanonical[eomLin]",
                "];",
                f"eomLin = ContractMetric[eomLin, {ctx.metric}];",
                'Print["Linearized EOM: ", Short[eomLin, 5]];',
                "",
            ]
        )

        # Pre-decomposition TT zeroing for single-field path
        lines.extend(_wls_pre_decomposition_tt_zeroing(ctx))

        # Build SkipTuples for TT-gauged field
        skip_opt = ""
        tt_fields = {entry["field"] for entry in ctx.gauge if entry["type"] == "tt"}
        if pert_field_name in tt_fields:
            skip_tuples = ", ".join(f"{{{0},{mu}}}" for mu in range(ctx.dim))
            skip_opt = f', "SkipTuples" -> {{{skip_tuples}}}'

        lines.extend(
            [
                "(* Decompose to components *)",
                f'componentEqs = DecomposeToComponents[eomLin, {fexpr}, {ctx.chart}, {{}}, "MetricMatrix" -> {p}MetricMatrix{skip_opt}];',
                'Print["Components: ", Length[componentEqs]];',
            ]
        )
        lines.extend(_wls_vector_background_substitution(ctx, "componentEqs"))
        lines.extend(_wls_validate_backgrounds_after_decompose(ctx, "componentEqs"))

        # Build tensor component metadata for LaTeX export
        lines.append("componentMetadata = <||>;")
        lines.extend(
            _wls_component_metadata(pert_field_name, fexpr, "componentEqs", ctx.dim)
        )

        # Build fieldEquations table
        lines.extend(
            [
                "",
                "fieldEquations = Table[",
                f'  {{"{pert_field_name}_" <> ToString[componentEqs[[k, 1]]], componentEqs[[k, 2]]}},',
                "  {k, Length[componentEqs]}",
                "];",
                "",
                "(* Memory cleanup: free abstract EOM and component array *)",
                "Clear[eomLin, componentEqs, l2ForVarD];",
                "Share[];",
                "",
            ]
        )
    else:
        # --- Multi-field linearization (metric + matter perturbations) ---
        # Build list of all dynamical fields and their VarD expressions
        dyn_fields: list[dict[str, Any]] = []
        if has_metric_pert:
            dyn_fields.append(
                {
                    "name": pert_field_name,
                    "field": pert_field,
                    "head": field_head,
                    "fexpr": fexpr,
                    "vard_expr": f"{field_head}[-a, -b]",
                }
            )
        for mpi in matter_pert_info:
            mp_field_dict: dict[str, Any] = {
                "name": mpi["pert_name"],
                "type": mpi["field_type"],
            }
            if mpi["field_type"] not in {"scalar", "vector"}:
                mp_field_dict["rank"] = int(mpi["field_rank"])
            mp_fexpr = _field_expression(
                mp_field_dict,
                ctx.prefix,
                head_override=mpi["pert_head"],
            )
            dyn_fields.append(
                {
                    "name": mpi["pert_name"],
                    "field": mp_field_dict,
                    "head": mpi["pert_head"],
                    "fexpr": mp_fexpr,
                    "vard_expr": mp_fexpr,
                }
            )

        # Add torsion perturbation as a dynamical field for VarD.
        # TorsionPert[LI[1], a, -b, -c] is the perturbation of TorsionCDT,
        # auto-registered by the torsion pipeline.  VarD w.r.t. TorsionPert
        # gives the torsion field equations from the PGT Lagrangian.
        if ctx.torsion:
            # Use the perturbation field (no LI index) for VarD and decomposition.
            # SlotsOfTensor with LI = 4 (wrong), without LI = 3 (correct rank).
            torsion_pert_name = ctx.torsion["perturbation_name"]
            torsion_field_head = f"{ctx.prefix}{torsion_pert_name.capitalize()}"
            torsion_fexpr = f"{torsion_field_head}[a, -b, -c]"
            dyn_fields.append(
                {
                    "name": torsion_pert_name,
                    "field": {
                        "name": torsion_pert_name,
                        "type": "tensor",
                        "rank": 3,
                        "symmetry": "antisymmetric_23",
                    },
                    "head": torsion_field_head,
                    "fexpr": torsion_fexpr,
                    "vard_expr": torsion_fexpr,
                }
            )

        # Pre-decomposition TT zeroing: set ComponentValue = 0 for h_{0,mu}
        # before DecomposeToComponents so ToBasis substitutes zeros during
        # expansion, reducing peak memory.
        lines.extend(_wls_pre_decomposition_tt_zeroing(ctx))

        lines.extend(
            (
                'Print["[TIMING] L^(2) finalized: ", Round[AbsoluteTime[] - tSubPhase, 0.1], " seconds"];',
                "tSubPhase = AbsoluteTime[];",
            )
        )

        eom_lines = _wls_multi_field_eom(ctx, dyn_fields)
        if eom_lines is not None:
            lines.extend(eom_lines)

    # ------------------------------------------------------------------
    # Step 3: Type B gauge fixing (if any) — constraints on fieldEquations
    # ------------------------------------------------------------------
    # For the Lagrangian path (ctx.lagrangian_expr), fieldEquations is empty
    # here — Component E-L fills it later (after Phase A decomposition).
    # Type B gauge is deferred to the post-Component-E-L block in
    # _wls_metadata_and_export.  Only apply here for the legacy (non-
    # Lagrangian) linearization path where VarD already filled fieldEquations.
    if (
        not ctx.lagrangian_expr
        and ctx.gauge
        and any(_resolve_gauge_mechanism(g) == "constraint" for g in ctx.gauge)
    ):
        lines.extend(_wls_gauge_fixing_type_b(ctx))

    lines.append(
        _wls_timing_end("tLinearize", "Linearization (xPert L^(2) + EOM decomposition)")
    )
    return lines


def _resolve_gauge_mechanism(entry: dict[str, Any]) -> str:
    """Return ``'lagrangian_term'`` or ``'constraint'`` for a gauge entry."""
    if entry["type"] == "custom":
        return str(entry["mechanism"])
    return _GAUGE_PRESETS[entry["type"]]["mechanism"]


# --- WLS: Gauge fixing ---


def _wls_gauge_fixing_type_a(ctx: _WlsContext) -> list[str]:
    """Generate Type A gauge-fixing terms (added to Lagrangian before EL).

    For named presets, emits a call to the corresponding ``Build*GaugeTerm``
    function in ``GaugeFix.wl``.  For custom expressions, emits the user's
    Wolfram expression directly (after field-name substitution).
    """
    # Build perturbation head lookup for matter perturbation fields
    # to avoid name collision (e.g., perturbation "a" → "gea", not "geA")
    pert_head_map: dict[str, str] = {}
    lin = ctx.linearization or {}
    for mp in lin.get("matter_perturbations", []):
        mp_name = mp["perturbation_name"]
        mf_name = mp["field"]
        candidate = f"{ctx.prefix}{mp_name.capitalize()}"
        mf_head = f"{ctx.prefix}{mf_name.capitalize()}"
        if candidate == mf_head:
            pert_head_map[mp_name] = f"{ctx.prefix}{mp_name}"
        else:
            pert_head_map[mp_name] = candidate

    lines: list[str] = ["(* Gauge fixing: Lagrangian terms *)"]
    for entry in ctx.gauge:
        if _resolve_gauge_mechanism(entry) != "lagrangian_term":
            continue
        field_name: str = entry["field"]
        xi = entry.get("xi", 1)
        # Use perturbation head if this field is a perturbation with name collision
        pfx_field = pert_head_map.get(
            field_name, f"{ctx.prefix}{field_name.capitalize()}"
        )

        if entry["type"] == "custom":
            expr = _substitute_field_names(
                entry["expression"],
                ctx.fields,
                ctx.prefix,
                derived_fields=ctx.derived_fields,
                background_fields=ctx.background_fields,
            )
            lines.append(f"{ctx.prefix}GaugeTerm = ({expr});")
        else:
            builder = _GAUGE_PRESETS[entry["type"]]["builder"]
            lines.append(
                f"{ctx.prefix}GaugeTerm = {builder}[{pfx_field}, {ctx.metric}, {ctx.cd}, {xi}];"
            )
        lines.extend(
            (
                f"{ctx.prefix}Lagrangian = AddGaugeFixingTerm["
                f"{ctx.prefix}Lagrangian, {ctx.prefix}GaugeTerm];",
                f"{ctx.prefix}Lagrangian = ToCanonical[{ctx.prefix}Lagrangian];",
                f"{ctx.prefix}Lagrangian = ContractMetric[{ctx.prefix}Lagrangian, {ctx.metric}];",
                f'Print["Gauge-fixed Lagrangian '
                f'({entry["type"]} on {field_name}): ", {ctx.prefix}Lagrangian];',
                "",
            )
        )
    return lines


def _type_b_coulomb_constraint(
    ctx: _WlsContext,
    field_name: str,
    comp_pfx: str,
) -> list[str]:
    """Generate WLS to add a Coulomb gauge constraint (spatial divergence = 0).

    Appends an additional constraint equation ``div A_spatial = 0`` to
    ``fieldEquations`` after component decomposition.
    """
    coords = _COORDS[ctx.dim]
    coord_args = ", ".join(f"{c}[]" for c in coords)

    # Build sum of spatial first derivatives: d_1 A_1 + d_2 A_2 + ...
    deriv_terms: list[str] = []
    for i in range(1, ctx.dim):
        deriv_indices = ", ".join("1" if j == i else "0" for j in range(ctx.dim))
        comp = f"{comp_pfx}{i}"
        deriv_terms.append(f"Derivative[{deriv_indices}][{comp}][{coord_args}]")

    constraint_expr = " + ".join(deriv_terms)
    constraint_name = f"{field_name}_coulomb"

    return [
        f"(* Coulomb gauge: spatial divergence of {field_name} = 0 *)",
        f'AppendTo[fieldEquations, {{"{constraint_name}", {constraint_expr}}}];',
        f'Print["Added Coulomb constraint: div {field_name} = 0"];',
        "",
    ]


def _sym_flat_index(a: int, b: int, dim: int) -> int:
    """Flat index for symmetric pair ``(a, b)`` in ``dim``-dimensional space.

    Matches the lexicographic ordering used by ``ReplaceRank2FieldComponents``
    in ``ComponentDecompose.wl``: pairs ``(i, j)`` with ``i <= j``, flattened
    to ``0, 1, 2, ...``
    """
    if a > b:
        a, b = b, a
    return a * (2 * dim - 1 - a) // 2 + b


def _tt_transverse_constraints(
    dim: int,
    comp_pfx: str,
    field_name: str,
    coord_args: str,
    prop_spatial: int,
) -> list[str]:
    """Generate WLS for TT transverse constraints: d^i h_{i,j} = 0 per spatial j.

    For each spatial direction j, replaces ``h_{prop,j}``'s dynamical equation
    with the transverse constraint ``Σ_i ∂_i h_{ij} = 0``.  After plane-wave
    reduction (only propagation-axis derivatives survive), each constraint
    simplifies to ``∂_{prop}(h_{prop,j}) = 0``, enabling the Python constraint
    eliminator to detect the self-referencing gradient-zero and set h_{prop,j}=0.

    Parameters
    ----------
    prop_spatial : int
        1-based spatial coordinate index of the propagation axis
        (e.g. 1 for ``x``, 2 for ``y``, 3 for ``z`` in 4D).

    The constraint for ``j = dim-1`` is skipped when ``prop_spatial == dim-1``
    because that field (h_{dim-1, dim-1}) is already eliminated by the traceless
    condition and has no equation to replace.  For other propagation axes, all
    j from 1 to dim-1 (inclusive) are covered.
    """
    coords = _COORDS[dim]
    lines: list[str] = []
    for j in range(1, dim):
        # Skip when j == prop_spatial == dim-1: h_{prop,prop} = h_{last,last}
        # is already eliminated by the traceless condition — no EOM to replace.
        if j == dim - 1 and prop_spatial == dim - 1:
            continue

        deriv_terms: list[str] = []
        for i in range(1, dim):
            flat_idx = _sym_flat_index(i, j, dim)
            comp = f"{comp_pfx}{flat_idx}"
            deriv_indices = ", ".join("1" if k == i else "0" for k in range(dim))
            deriv_terms.append(f"Derivative[{deriv_indices}][{comp}][{coord_args}]")

        constraint_expr = " + ".join(deriv_terms)
        coord_label = coords[j] if j < len(coords) else str(j)

        # Replace h_{prop, j}'s dynamical equation with the transverse constraint.
        # After plane-wave reduction along prop_axis, only ∂_{prop}(h_{prop,j})
        # survives → self-referencing gradient-zero detected by Python eliminator.
        constrained_idx = _sym_flat_index(prop_spatial, j, dim)
        constrained_field = f"{field_name}_{constrained_idx}"

        lines.extend(
            [
                f"(* TT transverse: d^i h_{{i,{j}}} = 0"
                f" — replaces EOM of {constrained_field} *)",
                "fieldEquations = Table[",
                f'  If[fieldEquations[[k, 1]] === "{constrained_field}",',
                f'    {{"{constrained_field}", {constraint_expr}}},',
                "    fieldEquations[[k]]",
                "  ],",
                "  {k, Length[fieldEquations]}",
                "];",
                f'Print["TT transverse: d^i h_{{i,{coord_label}}} = 0'
                f' -> replaced {constrained_field} EOM"];',
                "",
            ]
        )
    return lines


def _tt_traceless_substitution(
    dim: int,
    comp_pfx: str,
    field_name: str,
    coord_args: str,
    metric_diagonal: list[str] | None = None,
) -> list[str]:
    """Generate WLS for TT traceless substitution: h_{d-1,d-1} → -(metric-weighted sum).

    The TT traceless condition is ``g^{ij} h_{ij} = 0`` using the background
    spatial metric.  Solving for the last spatial diagonal:

    ``h_{last} = -Σ_{i≠last} (g_{last} / g_i) * h_i``

    For flat (Minkowski/Cartesian) metrics all diagonal entries are ±1 so the
    weights are all 1 and the formula reduces to ``h_{last} → -(Σ others)``.
    For curved diagonal metrics (e.g. spherical ``diag[-1, 1, r², r²sin²θ]``)
    the weights are metric-dependent expressions in the Wolfram coordinates.

    The constraint equation written into ``fieldEquations`` is the numerator
    of ``g^{ij} h_{ij} = 0`` multiplied by ``g_{last}`` so that h_{last}
    appears with self-coefficient 1 (enabling the Python degenerate-constraint
    detector after gradient-zero fields have been zeroed first):

    ``h_last_EOM = w_0 * h_0 + w_1 * h_1 + ... + h_last``

    where ``w_i = g_{last} / g_i`` (Wolfram Simplify is applied to each weight).

    Parameters
    ----------
    metric_diagonal : list[str] | None
        Diagonal metric entries as Wolfram expression strings, including the
        time entry at index 0 (e.g. ``["-1", "1", "x[]^2", "x[]^2*Sin[y[]]^2"]``).
        Killed-coordinate values should already have been substituted by
        :func:`_apply_coord_values` before this function is called, so that
        weight expressions like ``Simplify[(x[]^2*1)/(x[]^2)]`` collapse
        immediately to constants rather than carrying ``Sin[y[]]^2`` through
        the full xPert/Christoffel computation.
        When ``None`` or empty, the flat-metric behavior (all weights = 1) is used.
    """
    spatial_diag_indices = [_sym_flat_index(i, i, dim) for i in range(1, dim)]
    last_diag_idx = spatial_diag_indices[-1]
    other_diag_indices = spatial_diag_indices[:-1]
    last_comp = f"{comp_pfx}{last_diag_idx}"

    # Compute metric weights for each spatial diagonal component.
    # For diagonal metric with entries g_i, the traceless condition is:
    #   Σ_i g^{ii} h_{ii} = Σ_i (1/g_i) h_{ii} = 0
    # Solve for h_{last}: h_{last} = -Σ_{i≠last} (g_{last}/g_i) * h_i
    if metric_diagonal and len(metric_diagonal) == dim:
        spatial_metric = metric_diagonal[1:]  # drop time component
        last_g = spatial_metric[-1]  # e.g. "x[]^2*Sin[y[]]^2"
        other_g = spatial_metric[:-1]  # e.g. ["1", "x[]^2"]
        # Wolfram Simplify collapses ratios; killed-coordinate values have already
        # been substituted in Python (via _apply_coord_values), so e.g.
        # (x[]^2*Sin[Pi/2]^2)/(x[]^2) = 1 trivially — no Sin[y[]] survives.
        weights = [f"Simplify[({last_g})/({g})]" for g in other_g]
    else:
        # Flat metric: all weights = 1
        weights = ["1"] * len(other_diag_indices)

    repl_sum = " + ".join(
        f"{w} * {comp_pfx}{idx}[args]"
        for w, idx in zip(weights, other_diag_indices, strict=False)
    )
    deriv_repl_sum = " + ".join(
        f"{w} * Derivative[d][{comp_pfx}{idx}][args]"
        for w, idx in zip(weights, other_diag_indices, strict=False)
    )
    # Constraint equation: g_{last} * Sum g^{ii} h_{ii} = 0
    # Written with h_{last} on both sides so self-coeff = 1 for Python detector.
    # The EOM for h_last is a weighted sum of the other diagonal components plus h_last.
    trace_terms = (
        " + ".join(
            f"{w} * {comp_pfx}{idx}[{coord_args}]"
            for w, idx in zip(weights, other_diag_indices, strict=False)
        )
        + f" + {last_comp}[{coord_args}]"
    )

    return [
        f"(* TT traceless: substitute {last_comp} → -(metric-weighted sum of other diags) *)",
        "(* Condition: g^{ij} h_{ij} = 0  ⟹  h_last = -Σ (g_last/g_i) h_i *)",
        f"fieldEquations = fieldEquations /. {{"
        f"{last_comp}[args___] :> -({repl_sum}), "
        f"Derivative[d__][{last_comp}][args___] :> -({deriv_repl_sum})}};",
        "",
        "(* Expand all equations to simplify kinetic terms after traceless sub *)",
        "fieldEquations = Table[{fieldEquations[[k, 1]], "
        "Expand[fieldEquations[[k, 2]]]},"
        " {k, Length[fieldEquations]}];",
        "",
        f"(* Replace {last_comp} equation with algebraic traceless constraint *)",
        "(* Form: h_last_EOM = Σ w_i h_i + h_last  (self-coeff=1 for Python detector) *)",
        f'Do[If[fieldEquations[[k, 1]] === "{field_name}_{last_diag_idx}",'
        f' fieldEquations[[k]] = {{"{field_name}_{last_diag_idx}", {trace_terms}}}],'
        f" {{k, Length[fieldEquations]}}];",
        f'Print["Applied TT traceless: {last_comp} → -(metric-weighted sum), '
        f'eq replaced with constraint"];',
        "",
    ]


def _type_b_tt_gauge(
    ctx: _WlsContext,
    field_name: str,
    comp_pfx: str,
) -> list[str]:
    """Generate WLS for TT (transverse-traceless) gauge on a symmetric rank-2 tensor.

    Applies three conditions (dimension-aware):

    1. **Temporal** (``h_{0,mu} = 0``): zero the first ``dim`` components
       (all pairs with one time index) via ``ReplaceAll``.
    2. **Transverse** (``partial^i h_{ij} = 0``): append ``dim - 1`` spatial
       divergence constraint equations.
    3. **Traceless** (``eta^{ij} h_{ij} = 0``): substitute last spatial diagonal
       ``h_{d-1,d-1} → -(sum of other spatial diags)`` throughout all equations,
       ``Expand`` to simplify, then replace h_{d-1,d-1}'s equation with the
       algebraic traceless constraint.

    Ordering matters: transverse constraints are appended *before* the traceless
    substitution so that ``h_{d-1,d-1}`` references inside the transverse
    constraints also get substituted.  Each equation is ``Expand``-ed after
    substitution to ensure the kinetic matrix simplifies correctly (analytically
    diagonal after TT substitution).
    """
    dim = ctx.dim
    coord_args = ", ".join(f"{c}[]" for c in _COORDS[dim])
    lines: list[str] = []

    # Determine propagation axis spatial index (1-based).
    # For z-propagation (default, last axis): prop_spatial = dim-1.
    # For x-propagation (spherical radial): prop_spatial = 1.
    if ctx.reduction is not None:
        prop_coord = ctx.reduction["propagation_axis"]  # e.g. "x"
        prop_spatial = _COORDS[dim].index(prop_coord)  # e.g. 1 for [t,x,y,z]
    else:
        prop_spatial = dim - 1  # default: last spatial axis

    # --- 1. Temporal: zero h_{0,mu} for mu = 0 .. dim-1 ---
    for mu in range(dim):
        idx = _sym_flat_index(0, mu, dim)
        lines.extend(wl_zero_component(f"{comp_pfx}{idx}", field_name, "TT-temporal"))

    # --- 2. Transverse: d^i h_{i,j} = 0 per spatial j ---
    # Uses propagation axis so constraints are self-referencing after reduction.
    lines.extend(
        _tt_transverse_constraints(dim, comp_pfx, field_name, coord_args, prop_spatial)
    )

    # --- 3. Traceless: h_{d-1,d-1} → -(metric-weighted sum of other spatial diags) ---
    # Uses diagonal metric entries for curved backgrounds; flat weights for Minkowski.
    lines.extend(
        _tt_traceless_substitution(
            dim, comp_pfx, field_name, coord_args, ctx.metric_diagonal
        )
    )

    return lines


def _wls_pre_decomposition_tt_zeroing(ctx: _WlsContext) -> list[str]:
    """Set ``ComponentValue = 0`` for TT-zeroed components before decomposition.

    When TT gauge is active on a symmetric rank-2 tensor, the temporal
    components ``h_{0,mu} = 0`` are known a priori.  Setting them as
    ``ComponentValue`` rules lets xCoba's ``ToBasis`` substitute zeros
    *during* expansion of ALL equations (including cross-field terms),
    reducing peak memory and intermediate expression size.

    Must be called *before* ``DecomposeToComponents`` / ``_wls_multi_field_eom``.
    The post-decomposition Type B gauge (``_type_b_tt_gauge``) remains as
    defense-in-depth for transverse + traceless constraints.
    """
    lines: list[str] = []
    for entry in ctx.gauge:
        if entry["type"] != "tt":
            continue
        field_name: str = entry["field"]
        p = ctx.prefix
        chart = ctx.chart
        head = f"{p}{field_name.capitalize()}"
        for mu in range(ctx.dim):
            lines.append(wl_component_value(head, [(0, chart), (mu, chart)], "0"))
            if mu > 0:  # symmetric: also set (mu, 0)
                lines.append(wl_component_value(head, [(mu, chart), (0, chart)], "0"))
    if lines:
        lines.insert(
            0,
            "(* Pre-decomposition TT zeroing: h_{0,mu} = 0 as ComponentValue *)",
        )
        lines.append("")
    return lines


def _wls_gauge_fixing_type_b(ctx: _WlsContext) -> list[str]:
    """Generate Type B gauge-fixing code (constraints applied after decomposition).

    Type B constraints modify ``fieldEquations`` after component decomposition:

    - **temporal**: substitute ``A_0 → 0`` (and all derivatives) everywhere
    - **axial**: substitute last spatial component → 0
    - **coulomb**: append ``div A_spatial = 0`` constraint equation
    - **tt**: temporal zeroing + traceless substitution + transverse constraints
    """
    lines: list[str] = ["(* Gauge fixing: post-decomposition constraints *)"]
    for entry in ctx.gauge:
        if _resolve_gauge_mechanism(entry) != "constraint":
            continue
        gauge_type: str = entry["type"]
        field_name: str = entry["field"]
        comp_pfx = f"{ctx.prefix}{field_name.capitalize()}"

        if gauge_type == "temporal":
            lines.extend(wl_zero_component(f"{comp_pfx}0", field_name, "temporal"))
        elif gauge_type == "axial":
            last_spatial = ctx.dim - 1
            lines.extend(
                wl_zero_component(f"{comp_pfx}{last_spatial}", field_name, "axial"),
            )
        elif gauge_type == "coulomb":
            lines.extend(_type_b_coulomb_constraint(ctx, field_name, comp_pfx))
        elif gauge_type == "tt":
            lines.extend(_type_b_tt_gauge(ctx, field_name, comp_pfx))
    return lines


# --- WLS: Plane-wave reduction ---


def _build_plane_wave_deriv_rules(ctx: _WlsContext) -> list[str]:
    """Build Wolfram substitution rules to zero transverse Derivative patterns.

    Returns a list of rule strings for killed (transverse) axes. Each rule
    zeros ``Derivative[ords__][f_][args___]`` where a transverse slot has
    order > 0.  Reused by Lagrangian reduction and VarD EOM reduction.
    """
    if ctx.reduction is None:
        return []
    prop_axis = ctx.reduction["propagation_axis"]
    coords = ctx.coords
    killed = [c for c in coords[1:] if c != prop_axis]
    rules: list[str] = []
    for c in killed:
        slot = coords.index(c) + 1
        rules.append(
            f"  Derivative[ords__][f_][args___] /; Length[{{ords}}] >= {slot}"
            f" && {{ords}}[[{slot}]] > 0 :> 0"
        )
    return rules


def _wls_plane_wave_reduction_lagrangian(ctx: _WlsContext) -> list[str]:
    """Generate Wolfram code to zero transverse derivatives in the Lagrangian.

    For plane-wave reduction, all spatial derivatives except along the
    propagation axis are set to zero.  This is applied to ``{prefix}Lagrangian``
    **before** E-L derivation / component decomposition for efficiency.

    The Wolfram substitution zeros any ``Derivative`` whose order in a
    transverse coordinate slot is > 0.
    """
    if ctx.reduction is None:
        return []

    prop_axis = ctx.reduction["propagation_axis"]
    killed = [c for c in ctx.coords[1:] if c != prop_axis]
    rules = _build_plane_wave_deriv_rules(ctx)

    p = ctx.prefix
    lines: list[str] = [
        "",
        "(* ============================================================ *)",
        "(* Plane-wave reduction: zero transverse derivatives            *)",
        f"(* Propagation axis: {prop_axis}, killed axes: {', '.join(killed):<20s}*)",
        "(* Applied BEFORE E-L derivation for efficiency                  *)",
        "(* ============================================================ *)",
        "",
        f'Print["Applying plane-wave reduction (propagation axis: {prop_axis})"];',
        f'Print["Zeroing derivatives w.r.t. transverse axes: {", ".join(killed)}"];',
        "",
        f"{p}Lagrangian = {p}Lagrangian /. {{",
        ",\n".join(rules),
        "};",
        f"{p}Lagrangian = Expand[{p}Lagrangian];",
        f'Print["Lagrangian after plane-wave reduction: ", Short[{p}Lagrangian, 5]];',
        "",
    ]
    return lines


def _wls_plane_wave_reduction_equations(ctx: _WlsContext) -> list[str]:
    """Generate Wolfram code to zero transverse derivatives in fieldEquations.

    For the linearization path, the abstract Lagrangian uses covariant
    derivatives (``CD``) that don't match the ``Derivative`` pattern.
    This function applies the reduction after component decomposition,
    when equations are in explicit ``Derivative[ords__][f_][args__]`` form.

    Also used as a fallback/defense-in-depth for the non-linearization path.
    """
    if ctx.reduction is None:
        return []

    prop_axis = ctx.reduction["propagation_axis"]
    coords = ctx.coords
    killed = [c for c in coords[1:] if c != prop_axis]

    rules: list[str] = []
    for c in killed:
        slot = coords.index(c) + 1
        # Guard: short-arity derivatives (e.g. {2,0} in 3+1D) must not
        # index beyond their length — Mathematica emits Part::partw.
        rules.append(
            f"  Derivative[ords__][f_][args___] /; Length[{{ords}}] >= {slot} && {{ords}}[[{slot}]] > 0 :> 0"
        )

    lines: list[str] = [
        "",
        "(* === Plane-wave reduction: zero transverse derivatives in fieldEquations === *)",
        f'Print["Zeroing transverse derivatives in component equations ({", ".join(killed)})"];',
        "",
        "fieldEquations = fieldEquations /. {",
        ",\n".join(rules),
        "};",
        "fieldEquations = Table[",
        "  {fieldEquations[[k, 1]], Expand[fieldEquations[[k, 2]]]},",
        "  {k, Length[fieldEquations]}",
        "];",
        "",
    ]
    return lines


def _wls_plane_wave_field_elimination(ctx: _WlsContext) -> list[str]:
    """Generate Wolfram code to iteratively eliminate zero fields.

    After plane-wave reduction and component decomposition, some fields may
    have identically zero equations (all their terms were killed).  This
    function iteratively:

    1. Finds fields whose RHS is identically zero.
    2. Substitutes those fields (and their derivatives) with zero in all
       remaining equations.
    3. Removes the zero-field equations from ``fieldEquations``.
    4. Repeats until no more fields can be eliminated.
    """
    if ctx.reduction is None:
        return []

    lines: list[str] = [
        "",
        "(* === Plane-wave reduction: iterative zero-field elimination === *)",
        'Print["Eliminating fields with identically zero equations..."];',
        "",
        "stable = False;",
        "While[!stable,",
        "  zeroFieldNames = {};",
        "  Do[",
        "    If[fieldEquations[[k, 2]] === 0,",
        "      AppendTo[zeroFieldNames, fieldEquations[[k, 1]]]",
        "    ],",
        "    {k, Length[fieldEquations]}",
        "  ];",
        "",
        "  If[Length[zeroFieldNames] == 0,",
        "    stable = True,",
        "    (* else: build substitution rules and eliminate *)",
        '    Print["Eliminating zero fields: ", zeroFieldNames];',
        "    zeroRules = {};",
        "    allHeads = Union@Cases[",
        "      fieldEquations[[All, 2]],",
        "      f_Symbol[__] :> f, {0, Infinity}",
        "    ];",
        "    Do[",
        "      Do[",
        '        If[StringContainsQ[ToString[head], StringReplace[zfn, "_" -> ""]],',
        "          AppendTo[zeroRules, head[___] :> 0];",
        "          AppendTo[zeroRules, Derivative[__][head][___] :> 0]",
        "        ],",
        "        {head, allHeads}",
        "      ],",
        "      {zfn, zeroFieldNames}",
        "    ];",
        "    fieldEquations = fieldEquations /. zeroRules;",
        "    fieldEquations = Select[fieldEquations,",
        "      !MemberQ[zeroFieldNames, #[[1]]] &];",
        "    fieldEquations = Table[",
        "      {fieldEquations[[k, 1]], Expand[fieldEquations[[k, 2]]]},",
        "      {k, Length[fieldEquations]}",
        "    ];",
        "  ]",
        "];",
        "",
        'Print["Fields after reduction: ", Length[fieldEquations]];',
        'Print["Surviving fields: ", fieldEquations[[All, 1]]];',
        "",
    ]
    return lines


def _wls_plane_wave_coordinate_evaluation(ctx: _WlsContext) -> list[str]:
    """Generate Wolfram code to evaluate killed coordinates at specified values.

    After plane-wave reduction zeros transverse derivatives, metric-dependent
    coefficients (e.g. ``Sin[y[]]^2`` from curved TT-traceless substitution,
    ``x[]^2`` from ``g_{θθ}`` in spherical coordinates) may survive as
    position-dependent prefactors in the surviving equations.

    The ``[reduction]`` TOML section supports an optional ``coordinate_values``
    map specifying Wolfram expressions for each killed coordinate:

    .. code-block:: toml

        [reduction]
        type = "plane_wave"
        propagation_axis = "x"
        coordinate_values = {y = "Pi/2"}   # equatorial plane, Sin[Pi/2]=1

    This function generates a ``ReplaceAll`` step that substitutes those values
    into all component equations and applies ``Simplify`` to reduce surviving
    factors to constants.  When ``coordinate_values`` is absent or empty, no
    code is generated (no-op).

    Must be called *after* both ``_wls_plane_wave_reduction_equations`` and
    ``_wls_plane_wave_field_elimination`` so that all derivative-zeroing and
    field-elimination steps have already been applied.
    """
    if ctx.reduction is None:
        return []

    coord_values: dict[str, str] = ctx.reduction.get("coordinate_values", {})
    if not coord_values:
        return []

    # Build Wolfram replacement rules: {y[] -> Pi/2, z[] -> 0, ...}
    rules = ", ".join(f"{coord}[] -> {val}" for coord, val in coord_values.items())

    return [
        "",
        "(* === Plane-wave: evaluate killed coordinates at specified values === *)",
        f"(* coordinate_values from [reduction] TOML: {coord_values} *)",
        "(* ReplaceAll triggers Wolfram's built-in evaluation rules immediately:",
        "   Sin[Pi/2]->1, Cos[0]->1, etc. Expand[] then collapses products.",
        "   Simplify is NOT used here: it is very slow on symbolic expressions",
        "   with free parameters (Bpeak, z0, x[]) and provides no benefit over",
        "   Expand when coordinate values are known numeric constants. *)",
        f'Print["Evaluating killed coordinates: {coord_values}"];',
        "",
        f"fieldEquations = fieldEquations /. {wl_list(rules)};",
        "fieldEquations = Table[",
        "  {fieldEquations[[k, 1]], Expand[fieldEquations[[k, 2]]]},",
        "  {k, Length[fieldEquations]}",
        "];",
        'Print["After coordinate evaluation: ", Length[fieldEquations], " equations"];',
        "",
    ]


# --- WLS: Euler-Lagrange & decomposition ---


def _wls_euler_lagrange_multi(ctx: _WlsContext) -> list[str]:
    """Generate Euler-Lagrange, decomposition, and export lines for multi-field."""
    lines: list[str] = ["(* Step 4: Euler-Lagrange equations *)"]

    for field in ctx.fields:
        fname = field["name"]
        fexpr = _field_expression(field, ctx.prefix)
        eom_var = f"eom{fname.capitalize()}"
        lines.extend(
            (
                f"{eom_var} = VarD[{fexpr}, {ctx.cd}][{ctx.prefix}Lagrangian];",
                f"{eom_var} = ToCanonical[{eom_var}];",
                f"{eom_var} = ContractMetric[{eom_var}, {ctx.metric}];",
            )
        )
        # Scalar background fields need explicit substitution (ToBasis won't touch them)
        lines.extend(_wls_scalar_background_substitution(ctx, eom_var))
        lines.extend(
            (
                f'Print["EOM {fname}: ", {eom_var}];',
                "",
            )
        )

    # Step 5: Decompose
    lines.append("(* Step 5: Decompose to components *)")
    # Build set of TT-gauged field names for SkipTuples optimization
    tt_fields = {entry["field"] for entry in ctx.gauge if entry["type"] == "tt"}
    for i, field in enumerate(ctx.fields):
        fname = field["name"]
        fexpr = _field_expression(field, ctx.prefix)
        eom_var = f"eom{fname.capitalize()}"
        comp_var = f"comp{fname.capitalize()}"

        other_exprs = [
            _field_expression(f, ctx.prefix) for j, f in enumerate(ctx.fields) if j != i
        ]
        others_str = ", ".join(other_exprs) if other_exprs else ""

        skip_opt = ""
        if fname in tt_fields:
            skip_tuples = ", ".join(f"{{{0},{mu}}}" for mu in range(ctx.dim))
            skip_opt = f', "SkipTuples" -> {{{skip_tuples}}}'

        lines.extend(
            (
                f'{comp_var} = DecomposeToComponents[{eom_var}, {fexpr}, {ctx.chart}, {{{others_str}}}, "MetricMatrix" -> {ctx.prefix}MetricMatrix{skip_opt}];',
                f'Print["{fname} components: ", Length[{comp_var}]];',
            )
        )
        # Substitute vector/tensor backgrounds explicitly (ToBasis unreliable)
        lines.extend(_wls_vector_background_substitution(ctx, comp_var))
        # Validate all backgrounds resolved
        lines.extend(_wls_validate_backgrounds_after_decompose(ctx, comp_var))
        lines.append("")

    # Build tensor component metadata for LaTeX export
    lines.append("componentMetadata = <||>;")
    for field in ctx.fields:
        fname = field["name"]
        fexpr = _field_expression(field, ctx.prefix)
        comp_var = f"comp{fname.capitalize()}"
        lines.extend(_wls_component_metadata(fname, fexpr, comp_var, ctx.dim))

    # Step 6: Export
    lines.extend(("(* Step 6: Build and export JSON *)", "fieldEquations = Flatten[{"))
    for i, field in enumerate(ctx.fields):
        fname = field["name"]
        comp_var = f"comp{fname.capitalize()}"
        comma = "," if i < len(ctx.fields) - 1 else ""
        lines.append(
            f'  Table[{{"{fname}_" <> ToString[{comp_var}[[k, 1]]], {comp_var}[[k, 2]]}}, {{k, Length[{comp_var}]}}]{comma}'
        )
    lines.extend(("}, 1];", ""))

    return lines


def _wls_euler_lagrange_single(ctx: _WlsContext) -> list[str]:
    """Generate Euler-Lagrange, decomposition, and export lines for single field."""
    field = ctx.fields[0]
    fname = field["name"]
    fexpr = _field_expression(field, ctx.prefix)

    lines = [
        "(* Step 4: Euler-Lagrange equations *)",
        f"eom = EulerLagrangeEquation[{ctx.prefix}Lagrangian, {fexpr}, {ctx.cd}];",
        'Print["EOM: ", eom];',
        "",
    ]

    # Scalar background fields need explicit substitution (ToBasis won't touch them)
    lines.extend(_wls_scalar_background_substitution(ctx, "eom"))

    lines.extend(
        [
            "(* Step 5: Decompose to components *)",
            f'componentEqs = DecomposeToComponents[eom, {fexpr}, {ctx.chart}, {{}}, "MetricMatrix" -> {ctx.prefix}MetricMatrix];',
            'Print["Components: ", Length[componentEqs]];',
        ]
    )
    lines.extend(_wls_vector_background_substitution(ctx, "componentEqs"))
    lines.extend(_wls_validate_backgrounds_after_decompose(ctx, "componentEqs"))

    # Build tensor component metadata for LaTeX export
    lines.append("componentMetadata = <||>;")
    lines.extend(_wls_component_metadata(fname, fexpr, "componentEqs", ctx.dim))

    lines.extend(
        [
            "",
            "fieldEquations = Table[",
            f'  {{"{fname}_" <> ToString[componentEqs[[k, 1]]], componentEqs[[k, 2]]}},',
            "  {k, Length[componentEqs]}",
            "];",
            "",
        ]
    )

    return lines


def _wls_linearization(ctx: _WlsContext, *, include_bg: bool = False) -> list[str]:
    """Generate xPert linearization, decomposition, and export lines.

    Raises
    ------
    ValueError
        If ``ctx.linearization`` is ``None``.
    """
    if ctx.linearization is None:
        msg = "_wls_linearization called without linearization config"
        raise ValueError(msg)
    lin = ctx.linearization
    pert_field_name = lin["perturbation_field"]
    pert_field = next(f for f in ctx.fields if f["name"] == pert_field_name)
    fexpr = _field_expression(pert_field, ctx.prefix)

    # Auto-generated internal symbols
    pert_sym = f"{ctx.prefix}hpert"
    eps_sym = f"{ctx.prefix}Epsilon"

    # Substitute CD and field names in the expression
    prefixed_expr = _substitute_field_names(
        lin["expression"].strip(),
        ctx.fields,
        ctx.prefix,
        derived_fields=ctx.derived_fields,
        background_fields=ctx.background_fields,
    )

    bg_name = f"{ctx.prefix}Bg"

    lines: list[str] = [
        "(* Step 3: xPert metric perturbation setup *)",
        f"{pert_sym}Tensor = SetupMetricPerturbation[{ctx.metric}, {pert_sym}, {eps_sym}];",
        f'Print["Perturbation: metric = background + {eps_sym} * {pert_field_name}"];',
        "",
        "(* Step 4: Linearize tensor expression *)",
        f"linExpr = LinearizeTensorExpression[{prefixed_expr}];",
        'Print["Linearized expression: ", Short[linExpr, 3]];',
        "",
        "(* Convert xPert notation to plain tensor for pipeline *)",
        f"linExprPlain = linExpr /. {pert_sym}[LI[1], idx__] :> {ctx.prefix}{pert_field_name.capitalize()}[idx];",
    ]

    # Replace bg tensor with metric (bg = background metric by construction)
    if include_bg:
        lines.extend(
            (
                "(* Replace bg with metric — bg is the background by construction *)",
                f"linExprPlain = linExprPlain /. {bg_name} -> {ctx.metric};",
            )
        )

    lines.append("linExprPlain = Simplify[linExprPlain];")

    # Scalar background fields need explicit substitution (ToBasis won't touch them)
    lines.extend(_wls_scalar_background_substitution(ctx, "linExprPlain"))

    lines.extend(
        [
            'Print["Converted to plain tensor: ", Short[linExprPlain, 3]];',
            "",
            "(* Step 5: Decompose to components *)",
            f'componentEqs = DecomposeToComponents[linExprPlain, {fexpr}, {ctx.chart}, {{}}, "MetricMatrix" -> {ctx.prefix}MetricMatrix];',
            'Print["Components: ", Length[componentEqs]];',
        ]
    )
    lines.extend(_wls_vector_background_substitution(ctx, "componentEqs"))
    lines.extend(_wls_validate_backgrounds_after_decompose(ctx, "componentEqs"))

    # Build tensor component metadata for LaTeX export
    lines.append("componentMetadata = <||>;")
    lines.extend(
        _wls_component_metadata(pert_field_name, fexpr, "componentEqs", ctx.dim)
    )

    lines.extend(
        [
            "",
            "fieldEquations = Table[",
            f'  {{"{pert_field_name}_" <> ToString[componentEqs[[k, 1]]], componentEqs[[k, 2]]}},',
            "  {k, Length[componentEqs]}",
            "];",
            "",
        ]
    )

    return lines


# --- WLS: Canonical momentum & Hamiltonian ---


def _wls_constraint_metadata(
    cs_config: dict[str, Any], spatial_coords: list[str]
) -> list[str]:
    """Generate Wolfram metadata lines for constraint solver configuration.

    Translates the ``[constraint_solver]`` TOML section into metadata keys
    that ``ConstraintSolverHints`` in ExportJSON.wl reads at export time.
    """
    method = cs_config.get("method", "auto")
    max_iter = cs_config.get("max_iterations", 30)
    tol = cs_config.get("tolerance", 1e-10)

    # Format tolerance in Wolfram scientific notation: 1e-10 → 1*^-10
    tol_str = f"{tol:.0e}"
    tol_str = tol_str.replace("e+0", "*^").replace("e-0", "*^-")
    tol_str = tol_str.replace("e+", "*^").replace("e-", "*^-")

    lines = [
        '  "solve_constraints" -> True,',
        f'  "constraint_method" -> "{method}",',
        f'  "constraint_max_iterations" -> {max_iter},',
        f'  "constraint_tolerance" -> {tol_str},',
    ]

    # Build boundary conditions Association
    bcs = cs_config.get("boundary_conditions", {})
    if bcs:
        bc_parts: list[str] = []
        for coord in spatial_coords:
            bc_info = bcs.get(coord, {})
            bc_type = bc_info.get("type", "periodic")
            if "value" in bc_info:
                bc_parts.append(
                    f'    "{coord}" -> <|"type" -> "{bc_type}", '
                    f'"value" -> {bc_info["value"]}|>'
                )
            else:
                bc_parts.append(f'    "{coord}" -> <|"type" -> "{bc_type}"|>')
        bc_str = ",\n".join(bc_parts)
        lines.append(f'  "constraint_boundary_conditions" -> <|\n{bc_str}\n  |>,')

    return lines


def _matter_pert_head_map(ctx: _WlsContext) -> dict[str, str]:
    """Build field_name → Wolfram head mapping for matter perturbation fields.

    For matter perturbation theories (Einstein-Maxwell, etc.), the L^(2)
    Lagrangian uses perturbation field heads (e.g. ``mmpa``), not the
    original field heads (e.g. ``mmpA``).  This helper replicates the
    collision-avoidance logic from ``_wls_matter_perturbation_setup``.

    Returns a dict mapping perturbation field names to their Wolfram heads,
    plus a set of original field names that should be excluded from the
    canonical pipeline (they no longer appear in L^(2)).
    """
    if not ctx.linearization:
        return {}
    mp_list = ctx.linearization.get("matter_perturbations", [])
    if not mp_list:
        return {}

    p = ctx.prefix
    head_map: dict[str, str] = {}
    for mp in mp_list:
        mf_name = mp["field"]
        mp_name = mp["perturbation_name"]
        mf_head = f"{p}{mf_name.capitalize()}"
        mp_head_candidate = f"{p}{mp_name.capitalize()}"
        # Collision avoidance: if capitalize(pert) == capitalize(field),
        # use the perturbation name as-is (e.g. "a" → "mmpa" not "mmpA").
        mp_head = f"{p}{mp_name}" if mp_head_candidate == mf_head else mp_head_candidate
        head_map[mp_name] = mp_head
    return head_map


def _matter_pert_originals(ctx: _WlsContext) -> set[str]:
    """Return the set of original field names that are replaced by perturbations.

    These fields should be excluded from the canonical pipeline because they
    no longer appear as dynamical variables in L^(2).
    """
    if not ctx.linearization:
        return set()
    return {mp["field"] for mp in ctx.linearization.get("matter_perturbations", [])}


def _canonical_field_heads(ctx: _WlsContext) -> tuple[str, str]:
    """Return (heads_str, all_heads_str) for canonical pipeline WLS generation.

    For matter perturbation theories, uses the perturbation heads (e.g.
    ``mmpa``) and excludes original fields that were replaced by perturbations.
    """
    p = ctx.prefix
    pert_heads = _matter_pert_head_map(ctx)
    originals = _matter_pert_originals(ctx)

    field_heads: list[str] = []
    for f in ctx.fields:
        fname = f["name"]
        if fname in originals:
            # Original field replaced by perturbation — skip
            continue
        if fname in pert_heads:
            field_heads.append(pert_heads[fname])
        else:
            field_heads.append(f"{p}{fname.capitalize()}")

    # Include torsion field head if present (from [torsion] TOML section).
    # Without this, the canonical pipeline can't resolve torsion tensor
    # components in R̃² Lagrangian terms, leaving 45+ abstract dummy pairs.
    if ctx.torsion is not None:
        torsion_pert_name = ctx.torsion["perturbation_name"]
        field_heads.append(f"{p}{torsion_pert_name.capitalize()}")

    all_heads = list(field_heads)
    all_heads.extend(f"{p}{bg['name'].capitalize()}" for bg in ctx.background_fields)
    return ", ".join(field_heads), ", ".join(all_heads)


def _field_component_count(field: dict[str, Any], dim: int) -> int:
    """Return the number of independent components for *field* in *dim* spacetime."""
    ftype = field["type"]
    if ftype == "scalar":
        return 1
    if ftype == "vector":
        return dim
    rank = field.get("rank", 2)
    sym = field.get("symmetry")
    if sym == "symmetric":
        # C(dim + rank - 1, rank) = (dim+rank-1)! / (rank! * (dim-1)!)
        from math import comb

        return comb(dim + rank - 1, rank)
    if sym == "antisymmetric":
        # C(dim, rank) = dim! / (rank! * (dim-rank)!)
        from math import comb

        return comb(dim, rank)
    return dim**rank


def _wls_constraint_elimination() -> list[str]:
    """Generate Wolfram code to detect and eliminate constraint fields.

    Detects two constraint types in ``fieldEquations``:

    1. **Gradient-zero**: equation where the RHS is a single spatial
       derivative of the equation's own field (no time derivatives, no
       cross-field terms).  With periodic/Neumann BCs and zero IC,
       ``gradient(f) = 0`` implies ``f = 0``.

    2. **Degenerate algebraic**: equation with ``lhsTimeOrder = 0`` where
       the identity coefficient of the self-field sums to ≈ 1.0, meaning
       the field cancels from both sides.  The residual constraint
       ``0 = Σ c_i * other_i`` is solved for the field with largest
       |coefficient| (dep_field), and both the equation's own field
       AND dep_field are eliminated.

    Substitutes detected constraints into ``lagComp`` using Mathematica's
    exact symbolic algebra — replacing the fragile Python string-algebra
    post-processing in ``reduction.py``.  Also filters ``fieldEquations``
    to exclude eliminated fields (affects ``allCompNames`` downstream).

    Emits a warning if a nonlinear constraint is detected (cannot happen
    in the linearized regime but aids debugging for pipeline extensions).

    Requires ``compToFunc``, ``fieldFuncList``, ``coordSyms``, and
    ``lagComp`` to be defined in the Wolfram session.  Sets
    ``eliminatedFromCanonical`` (list of eliminated field names).
    """
    return [
        "",
        "(* === Constraint Elimination (Wolfram-side) === *)",
        "(* Detect and eliminate constraint fields from lagComp and          *)",
        "(* fieldEquations BEFORE the Legendre transform, so the Hamiltonian *)",
        "(* is computed for surviving fields only with exact symbolic         *)",
        "(* coefficients.                                                    *)",
        'Print[""];',
        'Print["Detecting constraint fields for Hamiltonian elimination..."];',
        "",
        "(* Phase 0: Zero out gauge-eliminated fields in lagComp.            *)",
        "(* TT temporal gauge skips h_0..h_3 during EOM decomposition via    *)",
        "(* SkipTuples, so they never appear in fieldEquations.  But lagComp *)",
        "(* is decomposed from the full abstract Lagrangian and still has    *)",
        "(* these terms.  Zero them here so the Hamiltonian is clean.        *)",
        "(* SKIP when VarD is skipped (fieldEquations = {}): all fields are  *)",
        "(* alive — the Component E-L will determine which are dynamical.    *)",
        "Module[{eqNames, allNames, gaugeElim, gaugeRules},",
        "  eqNames = If[Length[fieldEquations] > 0, fieldEquations[[All, 1]], Keys[compToFunc]];",
        "  allNames = Keys[compToFunc];",
        "  gaugeElim = Complement[allNames, eqNames];",
        "  If[Length[gaugeElim] > 0,",
        "    gaugeRules = Flatten[Table[",
        "      {compToFunc[name][___] -> 0,",
        "       Derivative[__][compToFunc[name]][___] -> 0},",
        "      {name, gaugeElim}",
        "    ]];",
        "    lagComp = lagComp /. gaugeRules;",
        "    fieldFuncList = Select[fieldFuncList,",
        "      !MemberQ[compToFunc /@ gaugeElim, #] &];",
        '    Print["Zeroed gauge-eliminated fields in Lagrangian: ",',
        "      gaugeElim];",
        "  ];",
        "];",
        "",
        "(* Phase 1: Gradient-zero constraints *)  ",
        "(* Equations where RHS is spatial derivative of own field only. *)",
        "(* With periodic BCs + zero IC: gradient(f)=0 => f=0.          *)",
        "gradZeroFields = {};",
        "Do[",
        "  Module[{name, eq, ownHead},",
        "    name = fieldEquations[[k, 1]];",
        "    eq = fieldEquations[[k, 2]];",
        "    ownHead = compToFunc[name];",
        "    (* Skip if has time derivatives *)",
        "    If[!FreeQ[eq, Derivative[n_, ___][_][___] /; n >= 1], Continue[]];",
        "    (* Skip if references other fields *)",
        "    If[!FreeQ[eq,",
        "         f_Symbol[___] /; MemberQ[fieldFuncList, f] && f =!= ownHead],",
        "       Continue[]];",
        "    (* Must have spatial derivatives of own field *)",
        "    If[!FreeQ[eq, Derivative[0, __][ownHead][___]],",
        "      AppendTo[gradZeroFields, name]",
        "    ]",
        "  ],",
        "  {k, Length[fieldEquations]}",
        "];",
        "",
        "(* Phase 2: Degenerate algebraic constraints *)",
        "(* When selfCoeff ~ 1.0, the equation's own field cancels from    *)",
        "(* both sides, leaving a residual constraint among OTHER fields.   *)",
        "(* Solve for the dep field with largest |coefficient| and          *)",
        "(* eliminate both the equation field AND the dep field.            *)",
        "(* Example: h_9 = h_4 + h_7 + h_9 => 0 = h_4 + h_7 => h_4=-h_7 *)",
        "(* Warns if nonlinear constraint detected (cannot happen in       *)",
        "(* linearized regime — quadratic L => linear PDEs).               *)",
        "degenerateEqFields = {};",
        "algebraicDepFields = <||>;",
        "Do[",
        "  Module[{name, eq, ownHead, selfCoeff, otherExpr,",
        "          otherHeads, depHead, depCoeff, depName, remainExpr},",
        "    name = fieldEquations[[k, 1]];",
        "    eq = fieldEquations[[k, 2]];",
        "    ownHead = compToFunc[name];",
        "    If[MemberQ[gradZeroFields, name], Continue[]];",
        "    (* Must be purely algebraic — no derivatives at all *)",
        "    If[!FreeQ[eq, Derivative[__][_][___]], Continue[]];",
        "",
        "    (* Nonlinearity check: warn if field powers or products detected *)",
        "    If[!FreeQ[eq, Power[f_Symbol[___], n_Integer] /;",
        "          n >= 2 && MemberQ[fieldFuncList, f]],",
        '      Print["WARNING: Nonlinear constraint (field power) for ", name,',
        '        ". Constraint elimination assumes linear equations ",',
        '        "(quadratic Lagrangian). Skipping."];',
        "      Continue[]",
        "    ];",
        "    (* Also detect field products: f[...]*g[...] within a single    *)",
        "    (* multiplicative term. A linear constraint h4+h7+h9 has multiple *)",
        "    (* field symbols in the sum, but only 1 per additive term.        *)",
        "    Module[{addTerms, hasProduct = False},",
        "      addTerms = If[Head[eq] === Plus, List @@ eq, {eq}];",
        "      Do[",
        "        Module[{factors, fieldCount},",
        "          factors = If[Head[term] === Times, List @@ term, {term}];",
        "          fieldCount = Count[factors,",
        "            f_Symbol[___] /; MemberQ[fieldFuncList, f]];",
        "          If[fieldCount >= 2, hasProduct = True; Break[]]",
        "        ],",
        "        {term, addTerms}",
        "      ];",
        "      If[hasProduct,",
        '        Print["WARNING: Nonlinear constraint (field products) for ",',
        '          name, ". Skipping."];',
        "        Continue[]",
        "      ]",
        "    ];",
        "",
        "    (* Extract coefficient of identity(self) *)",
        "    selfCoeff = Coefficient[eq, ownHead[Sequence @@ coordSyms]];",
        "    If[Abs[N[selfCoeff] - 1.0] > 10^-12, Continue[]];",
        "",
        "    (* Self cancels.  Residual constraint: 0 = otherExpr *)",
        "    otherExpr = eq - selfCoeff * ownHead[Sequence @@ coordSyms];",
        "    If[otherExpr === 0, Continue[]];  (* trivial identity *)",
        "",
        "    AppendTo[degenerateEqFields, name];",
        "",
        "    (* Find field heads in the residual constraint *)",
        "    otherHeads = Select[fieldFuncList,",
        "      !FreeQ[otherExpr, #] && # =!= ownHead &];",
        "    If[Length[otherHeads] == 0, Continue[]];",
        "",
        "    (* Pick dep field: largest |coefficient| for stability *)",
        "    depHead = First[SortBy[otherHeads,",
        "      -Abs[Coefficient[otherExpr,",
        "        #[Sequence @@ coordSyms]]] &]];",
        "    depCoeff = Coefficient[otherExpr,",
        "      depHead[Sequence @@ coordSyms]];",
        "    If[Abs[N[depCoeff]] < 10^-12, Continue[]];",
        "",
        "    (* Find the component name for depHead *)",
        "    depName = First[Select[Keys[compToFunc],",
        "      compToFunc[#] === depHead &]];",
        "",
        "    (* Solve: depHead = -(remainExpr)/depCoeff *)",
        "    remainExpr = otherExpr -",
        "      depCoeff * depHead[Sequence @@ coordSyms];",
        "    (* Guard: skip if dep field already targeted by another constraint *)",
        "    If[KeyExistsQ[algebraicDepFields, depName],",
        '      Print["WARNING: dep field ", depName,',
        '        " already targeted by another constraint. Skipping."];',
        "      Continue[]",
        "    ];",
        "    algebraicDepFields[depName] = -remainExpr / depCoeff",
        "  ],",
        "  {k, Length[fieldEquations]}",
        "];",
        "",
        "(* Build the full list of eliminated fields *)",
        "eliminatedFromCanonical = Join[gradZeroFields,",
        "  degenerateEqFields, Keys[algebraicDepFields]];",
        "(* Remove duplicates (dep field might also be an eq field) *)",
        "eliminatedFromCanonical = DeleteDuplicates[eliminatedFromCanonical];",
        "",
        "If[Length[eliminatedFromCanonical] > 0,",
        "  Module[{zeroRules, algRules, allSubRules,",
        "          newLagComp, newFieldEqs},",
        "",
        "    (* Zero rules: gradient-zero + degenerate-eq fields -> 0 *)",
        "    zeroRules = Flatten[Table[",
        "      {compToFunc[name][___] -> 0,",
        "       Derivative[__][compToFunc[name]][___] -> 0},",
        "      {name, Join[gradZeroFields, degenerateEqFields]}",
        "    ]];",
        "",
        "    (* Algebraic rules: dep fields -> linear combination.     *)",
        "    (* Extract ALL derivative orders from lagComp for each     *)",
        "    (* eliminated head, then create explicit rules using D[].  *)",
        "    algRules = {};",
        "    Do[",
        "      Module[{head, subExpr, derivOrders},",
        "        head = compToFunc[depName];",
        "        subExpr = algebraicDepFields[depName];",
        "",
        "        (* Identity rule *)",
        "        AppendTo[algRules,",
        "          head[Sequence @@ coordSyms] -> subExpr];",
        "",
        "        (* Find all derivative orders for this head in lagComp *)",
        "        derivOrders = DeleteDuplicates[Cases[lagComp,",
        "          Derivative[ords__][head][___] :> {ords}, Infinity]];",
        "",
        "        (* Also scan fieldEquations for derivative patterns *)",
        "        derivOrders = DeleteDuplicates[Join[derivOrders,",
        "          Cases[fieldEquations,",
        "            Derivative[ords__][head][___] :> {ords},",
        "            Infinity]",
        "        ]];",
        "",
        "        (* Filter: only keep orders matching coordinate count *)",
        "        derivOrders = Select[derivOrders,",
        "          Length[#] == Length[coordSyms] &];",
        "",
        "        (* Explicit rule for each derivative order *)",
        "        Do[",
        "          AppendTo[algRules,",
        "            Derivative[Sequence @@ dOrd][head][",
        "              Sequence @@ coordSyms] ->",
        "            D[subExpr, Sequence @@ MapThread[",
        "              {#1, #2} &, {coordSyms, dOrd}]]",
        "          ],",
        "          {dOrd, derivOrders}",
        "        ]",
        "      ],",
        "      {depName, Keys[algebraicDepFields]}",
        "    ];",
        "",
        "    allSubRules = Join[zeroRules, algRules];",
        "",
        "    (* Substitute in lagComp — no Expand, let Legendre handle *)",
        "    newLagComp = lagComp /. allSubRules;",
        "    lagComp = newLagComp;",
        "",
        "    (* Update field function list to surviving fields only *)",
        "    fieldFuncList = Select[fieldFuncList,",
        "      !MemberQ[compToFunc /@ eliminatedFromCanonical, #] &];",
        "",
        "    (* ============================================================ *)",
        "    (* RE-DERIVE EOM from reduced Lagrangian via Euler-Lagrange.    *)",
        "    (* Simply substituting constraints into the original equations  *)",
        "    (* misses chain-rule corrections:                               *)",
        "    (*   EOM_k^reduced = EOM_k + Σ_i (∂constraint_i/∂f_k) EOM_i   *)",
        "    (* Re-deriving from lagComp is exact and extensible.            *)",
        "    (* ============================================================ *)",
        '    Print["Re-deriving EOM from reduced Lagrangian..."];',
        "",
        "    (* Build list of {name, func} for surviving fields *)",
        "    Module[{survPairs, newFieldEqs},",
        "      survPairs = Select[",
        "        Table[{fieldEquations[[k, 1]], compToFunc[fieldEquations[[k, 1]]]},",
        "          {k, Length[fieldEquations]}],",
        "        !MemberQ[eliminatedFromCanonical, #[[1]]] &",
        "      ];",
        "",
        "      (* Euler-Lagrange operator for component-level Lagrangian.    *)",
        "      (* For field f(coords): dL/df = Sum_a (-1)^|a| D^a(dL/d(D^a f))*)",
        "      (* We enumerate all derivative multi-indices present in       *)",
        "      (* lagComp for each field, then sum the E-L contributions.    *)",
        "      newFieldEqs = Table[",
        "        Module[{name, func, derivOrders, eom},",
        "          name = survPairs[[j, 1]];",
        "          func = survPairs[[j, 2]];",
        "",
        "          (* Collect all derivative orders of this field in lagComp *)",
        "          derivOrders = DeleteDuplicates[Join[",
        "            {{Sequence @@ Table[0, {Length[coordSyms]}]}},  (* identity *)",
        "            Cases[lagComp,",
        "              Derivative[ords__][func][___] :> {ords},",
        "              Infinity]",
        "          ]];",
        "",
        "          (* Filter: only keep multi-indices matching coord count *)",
        "          derivOrders = Select[derivOrders,",
        "            Length[#] == Length[coordSyms] &];",
        "",
        "          (* Compute E-L variation: Sum_a (-1)^|a| D^a[dL/d(D^a f)] *)",
        "          eom = Sum[",
        "            Module[{alpha, fDeriv, partialL, sign, diffSpec},",
        "              alpha = derivOrders[[m]];",
        "              fDeriv = Derivative[Sequence @@ alpha][func][",
        "                Sequence @@ coordSyms];",
        "              partialL = D[lagComp, fDeriv];",
        "              sign = (-1)^Total[alpha];",
        "              (* Build differentiation specification for D[] *)",
        "              diffSpec = Flatten[MapThread[",
        "                Table[{#1}, {#2}] &,",
        "                {coordSyms, alpha}",
        "              ]];",
        "              If[Length[diffSpec] > 0,",
        "                sign * D[partialL, Sequence @@ diffSpec],",
        "                sign * partialL",
        "              ]",
        "            ],",
        "            {m, Length[derivOrders]}",
        "          ];",
        "",
        "          {name, Expand[eom]}",
        "        ],",
        "        {j, Length[survPairs]}",
        "      ];",
        "",
        "      fieldEquations = newFieldEqs",
        "    ];",
        '    Print["EOM re-derived for ", Length[fieldEquations], " surviving fields"];',
        "",
        '    Print["Eliminated from canonical: ",',
        '      Length[eliminatedFromCanonical], " field(s): ",',
        "      eliminatedFromCanonical];",
        "  ],",
        "",
        '  Print["No constraint fields detected."];',
        "];",
        "",
        "(* EOM substitution for d2_t: DISABLED (issue #195).                  *)",
        "",
    ]


def _wls_canonical_phase_a(ctx: _WlsContext, all_heads_str: str) -> list[str]:  # noqa: C901, PLR0914, PLR0915
    """Generate WLS code for canonical Phase A: decompose Lagrangian + constraint elimination.

    Decomposes the abstract Lagrangian into component form (``lagComp``),
    builds the component-to-function mapping (``compToFunc``), normalizes
    derivative arities, applies plane-wave reduction, and runs constraint
    elimination on ``fieldEquations`` + ``lagComp``.

    This phase MUST run before ``BuildMultiFieldJSONStructure`` so that
    ``fieldEquations`` contains only surviving fields.

    Sets up WLS variables: ``lagComp``, ``compToFunc``, ``fieldFuncList``,
    ``velOrders``, ``coordSyms``, ``eliminatedFromCanonical``.
    """
    p = ctx.prefix

    # Build BackgroundFieldRules for DecomposeScalarExpression (same as EOM path).
    # Enables early background field evaluation during batched TraceBasisDummy,
    # critical for memory reduction in Einstein-Maxwell theories.
    bg_rules_entries: list[str] = []
    for bf in ctx.background_fields:
        if bf["type"] != "scalar" and bf.get("components"):
            bg_head = f"{p}{bf['name'].capitalize()}"
            comps_str = ", ".join(str(c) for c in bf["components"])
            contra_comps = _compute_contra_components(
                bf["components"], ctx.metric_diagonal
            )
            contra_str = ", ".join(contra_comps)
            bg_rules_entries.append(wl_bg_rule_entry(bg_head, comps_str, contra_str))
    bg_rules_opt = ""
    if bg_rules_entries:
        bg_rules_str = ", ".join(bg_rules_entries)
        bg_rules_opt = f', "BackgroundFieldRules" -> {wl_list(bg_rules_str)}'

    # NOTE: ComponentValue PD zeroing was attempted for the canonical path but
    # REMOVED.  For scalar Lagrangians, all indices are dummy — TraceBasisDummy
    # must enumerate all combinations regardless.  Extra xAct rules add
    # pattern-matching overhead that outweighs any benefit.
    # Measured: 8 rules +56% slower, 40 rules +9% slower on batch[1:50/82].
    # The per-term plane-wave reduction (after ConvertCDToDerivatives) is the
    # correct and effective approach for the canonical path.

    lines: list[str] = [
        "",
        "(* === Phase K: Canonical Momentum & Hamiltonian (component-level) === *)",
        'Print[""];',
        'Print["Computing canonical momenta and Hamiltonian..."];',
        "",
        "(* Copy Lagrangian for canonical analysis *)",
        f"lagForCanon = {p}Lagrangian;",
        "",
    ]

    # --- Self-energy sector filtering for multi-field theories with torsion ---
    # For theories where full decomposition would timeout, filter L^(2) to
    # contain only single-perturbation-field self-energy terms (h-only, a-only).
    # Background fields (Ābar) are KEPT as they are coefficients, not dynamical.
    # Cross-coupling (h*a, h*T, a*T) and torsion self-energy (T*T) are removed.
    # This produces the correct self-energy Hamiltonian for each field sector.
    if ctx.torsion is not None and len(ctx.fields) > 1:
        # Build list of perturbation field heads (NOT background fields)
        pert_heads = _matter_pert_head_map(ctx)
        originals = _matter_pert_originals(ctx)
        pert_field_heads: list[str] = []
        for f in ctx.fields:
            fname = f["name"]
            if fname in originals:
                continue
            head = pert_heads.get(fname, f"{p}{fname.capitalize()}")
            pert_field_heads.append(head)
        torsion_head = f"{p}{ctx.torsion['perturbation_name'].capitalize()}"
        pert_field_heads.append(torsion_head)

        # Generate Wolfram code to filter L^(2) by perturbation field content
        heads_wl = ", ".join(pert_field_heads)
        lines.extend(
            [
                "(* === Self-energy sector filtering ===                                 *)",
                "(* For multi-field torsion theories: keep ONLY single-perturbation-field *)",
                "(* self-energy terms. Background fields (Abar) are coefficients, kept.  *)",
                "(* Cross-coupling (h*a), torsion interactions (h*T, a*T), and torsion   *)",
                "(* self-energy (T*T, zero for non-propagating) are filtered out.        *)",
                "(* This produces correct self-energy H_hh + H_aa for conversion         *)",
                "(* measurements without the expensive torsion contraction terms.         *)",
                f"Module[{{pertHeads = {wl_list(heads_wl)}, lagTerms, selfTerms}},",
                "  lagTerms = If[Head[lagForCanon] === Plus, List @@ lagForCanon, {lagForCanon}];",
                '  Print["L^(2) before self-energy filter: ", Length[lagTerms], " terms"];',
                "",
                "  (* Keep only h-self and a-self sectors: terms with exactly 1 perturbation *)",
                "  (* field that is NOT the torsion field. Background fields (Abar) are NOT  *)",
                "  (* perturbation fields and do NOT affect this classification.              *)",
                f"  Module[{{torsionHead = {torsion_head}}},",
                "    selfTerms = Select[lagTerms, Function[{term},",
                "      Module[{presentPerts = Select[pertHeads, !FreeQ[term, #]&]},",
                "        (* Keep: terms with 0 pert fields (pure constant) or 1 non-torsion pert field *)",
                "        Length[presentPerts] == 0 ||",
                "        (Length[presentPerts] == 1 && presentPerts[[1]] =!= torsionHead)",
                "      ]",
                "    ]];",
                "  ];",
                "",
                "  lagForCanon = Total[selfTerms];",
                '  Print["L^(2) after self-energy filter: ", Length[selfTerms], " terms",',
                '    " (removed ", Length[lagTerms] - Length[selfTerms], " cross-coupling/torsion terms)"];',
                "];",
                "",
            ]
        )

    lines.extend(
        [
            "(* Apply CD shorthand rules to canonical Lagrangian.                     *)",
            "(* ComponentValues are pre-computed (supervisor's pattern) — ToValues    *)",
            "(* resolves CDN-field tensors in O(1) during StaggeredToBasis, reducing *)",
            "(* TraceBasisDummy K from ~4+ contracted pairs to ~0.                   *)",
            "(* Apply highest-order first for correct nesting:                       *)",
            "(*   CD[-a]@CD[-b]@field -> CD[-a]@CD1field -> CD2field                *)",
            "If[ListQ[$CDShorthandRules] && Length[$CDShorthandRules] > 0,",
            '  Print["Applying CD shorthand rules to canonical Lagrangian..."];',
            "  Do[lagForCanon = lagForCanon /. rule, {rule, Reverse[$CDShorthandRules]}];",
            "  (* Canonicalize with TimeConstrained — always attempt, revert if expansion. *)",
            "  (* The supervisor ALWAYS canonicalizes; old cost-bypass threshold could     *)",
            "  (* skip entirely, leaving consolidatable terms separate.                   *)",
            "  (* Guard: only keep result if term count decreases (expansion = revert).  *)",
            "  Module[{nBefore, tCanon, canonResult, nAfter},",
            "    nBefore = If[Head[lagForCanon]===Plus, Length[lagForCanon], 1];",
            "    tCanon = AbsoluteTime[];",
            f"    canonResult = Quiet[TimeConstrained[ContractMetric[ToCanonical[lagForCanon], {ctx.metric}], 30, $Failed], ToCanonical::noident];",
            "    If[canonResult =!= $Failed,",
            "      nAfter = If[Head[canonResult]===Plus, Length[canonResult], 1];",
            "      If[nAfter <= nBefore,",
            "        lagForCanon = canonResult;",
            '        Print["  ToCanonical+ContractMetric: ", nBefore, " -> ", nAfter, " terms (",',
            '          Round[AbsoluteTime[] - tCanon, 0.1], "s)"],',
            '        Print["  ToCanonical+ContractMetric: expanded ", nBefore, " -> ", nAfter,',
            '          " terms — reverted (", Round[AbsoluteTime[] - tCanon, 0.1], "s)"]',
            "      ],",
            '      Print["  ToCanonical+ContractMetric: timed out at 30s (", nBefore, " terms unchanged)"]',
            "    ];",
            "  ];",
            '  Print["After CD shorthand: ", If[Head[lagForCanon]===Plus, Length[lagForCanon], 1], " terms"];',
            "];",
            "",
            "(* NOTE: SeparateMetric is NOT applied here — it's applied inside         *)",
            "(* DecomposeScalarExpression (via SeparateFieldMetrics at line 1913 and  *)",
            "(* SeparateMetric in StaggeredToBasis at line 458). Applying it here too *)",
            "(* was a triple application that added 2-5s per term * N terms overhead. *)",
        ]
    )

    # Apply scalar BG substitutions before decomposition
    lines.extend(_wls_scalar_background_substitution(ctx, "lagForCanon"))

    lines.extend(
        [
            "",
            "(* NOTE: Together[] pre-simplification was tested (commit b2883fd) but    *)",
            "(* caused +30% regression on flat Gertsenshtein (1554s vs 1193s) even   *)",
            "(* when term count was unchanged (34→34).  Root cause: Together+Expand   *)",
            "(* subtly restructures coefficient expressions (e.g. rational vs product *)",
            "(* form), making downstream ToBasis+TraceBasisDummy slower.              *)",
            "(* Only beneficial for GW 1D (1→16 terms, 13% faster) where the single  *)",
            "(* complex EH term benefits from splitting into simpler sub-terms.        *)",
            "",
            "(* Pre-resolve Scalar[] wrappers in canonical Lagrangian.               *)",
            "(* For R̃²-decomposed torsion theories, Scalar wrappers contain         *)",
            "(* contracted tensor products that create 45+ abstract dummy pairs.     *)",
            "(* Resolve each Scalar independently via DecomposeScalarExpression       *)",
            "(* BEFORE splitting into additive terms, matching the EOM pipeline.     *)",
            "If[!FreeQ[lagForCanon, Scalar],",
            '  Print["Pre-resolving Scalar wrappers in canonical Lagrangian..."];',
            "  lagForCanon = lagForCanon /. Scalar[x_] :> Module[{resolved},",
            "    resolved = Quiet[Catch[",
            f"      DecomposeScalarExpression[x, {ctx.chart}, {wl_list(all_heads_str)},",
            f'        "MetricMatrix" -> {p}MetricMatrix]',
            "    ], {Validate::repeated, Validate::inhom}];",
            "    If[StringQ[resolved] || resolved === Null, Scalar[x], resolved]",
            "  ];",
            '  Print["Canonical Scalar resolution: ",',
            '    If[FreeQ[lagForCanon, Scalar], "all resolved", "some remain"]];',
            "];",
            "",
            "(* Re-apply CD shorthand rules AFTER Scalar pre-resolution.             *)",
            "(* Resolved Scalar[] contents may contain new CD operators               *)",
            "(* (e.g., CD[-a][CD1H[...]] from expanded Ricci scalar) that the        *)",
            "(* first pass couldn't see inside opaque Scalar wrappers.                *)",
            "(* Converting CD[CD1H] → CD2H enables ToValues O(1) resolution.         *)",
            "If[ListQ[$CDShorthandRules] && Length[$CDShorthandRules] > 0,",
            "  Do[lagForCanon = lagForCanon /. rule, {rule, Reverse[$CDShorthandRules]}];",
            "  (* Canonicalize after 2nd CD shorthand pass — only if expression still    *)",
            "  (* has abstract tensor indices (not component-form functions).            *)",
            "  (* Scalar pre-resolution converts some terms to component form;          *)",
            "  (* ToCanonical on those emits noident warnings and expands uselessly.    *)",
            "  If[!FreeQ[lagForCanon, _?AbstractIndexQ],",
            "    Module[{n2Before, t2Canon, canon2Result, n2After},",
            "      n2Before = If[Head[lagForCanon]===Plus, Length[lagForCanon], 1];",
            "      t2Canon = AbsoluteTime[];",
            f"      canon2Result = Quiet[TimeConstrained[ContractMetric[ToCanonical[lagForCanon], {ctx.metric}], 30, $Failed], ToCanonical::noident];",
            "      If[canon2Result =!= $Failed,",
            "        n2After = If[Head[canon2Result]===Plus, Length[canon2Result], 1];",
            "        If[n2After <= n2Before,",
            "          lagForCanon = canon2Result;",
            '          Print["  2nd-pass ToCanonical: ", n2Before, " -> ", n2After, " terms (",',
            '            Round[AbsoluteTime[] - t2Canon, 0.1], "s)"],',
            '          Print["  2nd-pass ToCanonical: expanded ", n2Before, " -> ", n2After,',
            '            " terms — reverted (", Round[AbsoluteTime[] - t2Canon, 0.1], "s)"]',
            "        ],",
            '        Print["  2nd-pass ToCanonical: timed out at 30s"]',
            "      ];",
            "    ],",
            '    Print["  2nd-pass ToCanonical: skipped (expression in component form)"];',
            "  ];",
            '  Print["After 2nd CD shorthand pass: ",',
            '    If[Head[lagForCanon]===Plus, Length[lagForCanon], 1], " terms"];',
            "];",
            "",
            "(* === Sector-aware Lagrangian splitting (theory-agnostic) ===          *)",
            "(* Split expanded Lagrangian by field content: single-field terms are  *)",
            "(* fast (low K), cross-field terms are the expensive part. Processing  *)",
            "(* each sector independently reduces per-term dummy pair count.        *)",
            "(* GatherBy uses actual field head presence, not hardcoded categories. *)",
            "(* For single-field theories: 1 group (no overhead).                  *)",
            "(* Ref: issue #198 — combined torsion+EM 10x slower than sum of parts.*)",
            f"Module[{{lagExp, lagList, allDynHeads = {wl_list(all_heads_str)},",
            "         termSectors, sectorGroups, sectorTerms, validationSum},",
            "  lagExp = Expand[lagForCanon];",
            "  lagList = If[Head[lagExp] === Plus, List @@ lagExp, {lagExp}];",
            '  Print["Expanded Lagrangian: ", Length[lagList], " additive terms"];',
            "",
            "  (* Build extended head list: base heads + CD shorthand heads *)",
            "  (* After CD shorthand substitution, CD[-a]@H[...] becomes CD1H[...], *)",
            "  (* so FreeQ[term, H] returns True even for H-containing terms.       *)",
            "  (* Include CD1-CD4 prefixed heads in the detection.                  *)",
            "  Module[{allExtHeads},",
            "    allExtHeads = Flatten[Map[Function[{h}, Join[{h},",
            '      Select[{Symbol["CD1" <> ToString[h]], Symbol["CD2" <> ToString[h]],',
            '              Symbol["CD3" <> ToString[h]], Symbol["CD4" <> ToString[h]]},',
            "        xTensorQ]]",
            "    ], allDynHeads]];",
            "",
            "    (* Classify each term by which BASE field heads it contains *)",
            "    (* A term matches head H if it contains H OR any of CD1H..CD4H *)",
            "    termSectors = Map[Function[{term},",
            "      Select[allDynHeads, Function[{h},",
            "        Module[{extH = Select[allExtHeads, StringMatchQ[ToString[#],",
            "          __ ~~ ToString[h]] &]},",
            "          !AllTrue[extH, FreeQ[term, #]&]",
            "        ]",
            "      ]]",
            "    ], lagList];",
            "  ];",
            "",
            "  (* Group terms by sector signature *)",
            "  sectorGroups = GatherBy[",
            "    Transpose[{lagList, termSectors}], Last];",
            "",
            "  (* Separate already-decomposed terms (empty sector = no tensor heads). *)",
            "  (* These arise from Scalar[] pre-resolution and are already in         *)",
            "  (* component form — skip DecomposeScalarExpression entirely.           *)",
            "  Module[{emptyPos, nonEmptyGroups},",
            "    emptyPos = Position[sectorGroups, g_ /; g[[1, 2]] === {}, {1}];",
            "    If[Length[emptyPos] > 0,",
            "      $tidalPreDecomposed = Total[sectorGroups[[emptyPos[[1, 1]], All, 1]]];",
            '      Print["  ", Length[sectorGroups[[emptyPos[[1, 1]]]]],'
            '        " terms already in component form (skipping decomposition)"];',
            "      nonEmptyGroups = Delete[sectorGroups, emptyPos];",
            "      sectorGroups = nonEmptyGroups,",
            "      $tidalPreDecomposed = 0",
            "    ];",
            "  ];",
            "",
            "  (* Flatten sectors into individual terms (not sector sums).          *)",
            "  (* Processing individual terms bounds per-term TraceBasisDummy cost  *)",
            "  (* to each term's own contraction depth K, rather than K of the sum.*)",
            "  (* Benchmarked: sector sums take 464s for 9-term tnrH group;       *)",
            "  (* individual terms take ~5-12s each = ~45-108s total.             *)",
            "  lagTerms = Flatten[Map[#[[All, 1]]&, sectorGroups]];",
            "  lagTerms = Select[lagTerms, # =!= 0 &];",
            "",
            "  (* Validate: sector sum must equal original *)",
            "  validationSum = Total[lagTerms] + $tidalPreDecomposed - lagExp // Expand;",
            "  If[validationSum =!= 0,",
            '    Print["WARNING: sector split validation failed — falling back to single term"];',
            "    lagTerms = {lagForCanon}; $tidalPreDecomposed = 0,",
            '    Print["Sector split: ", Length[lagTerms], " groups from ",',
            '      Length[lagList], " terms"];',
            "    Do[",
            '      Print["  sector ", k, ": ", ',
            '        sectorGroups[[k, 1, 2]], " (",',
            '        Length[sectorGroups[[k]]], " terms)"],',
            "      {k, Min[Length[sectorGroups], Length[lagTerms]]}",
            "    ]",
            "  ];",
            "];",
            "",
            'Print["Decomposing Lagrangian: ", Length[lagTerms], " sector groups"];',
            _wls_timing_start("tCanonDecomp"),
        ]
    )

    # NOTE: PD DownValues approach was tested (both ComponentValue and direct
    # DownValues with specific integers) and BOTH caused regressions:
    # - ComponentValue: +9-56% slower (xAct internal tables, not auto-evaluating)
    # - Direct DownValues: +66% slower on EM curved (adds overhead to every
    #   TraceBasisDummy evaluation step; metric DownValues work because metrics
    #   are multiplicative factors that zero products before inner evaluation,
    #   but PD wraps inner expressions that are already evaluated)
    # The correct approach is post-ToBasis ReplaceAll in DecomposeScalarExpression,
    # which zeros PD[{killedAxis,-chart}][...] terms BEFORE TraceBasisDummy.
    # This is implemented in ComponentDecompose.wl step 2.5.

    lines.extend(
        [
            "lagComp = 0;",
            "Do[",
            "  Module[{termComp, tTerm = AbsoluteTime[]},",
            "    (* Block[{Print=Null}]: suppress xAct internal prints during         *)",
            "    (* decomposition. Output buffering overhead adds ~1-2s per term.    *)",
            "    (* Supervisor wraps all expensive ops this way (SphericalEuclidean). *)",
            f"    termComp = Block[{{Print = Null}}, Quiet[Catch[DecomposeScalarExpression[lagTerms[[k]], {ctx.chart}, {wl_list(all_heads_str)}, "
            f'"MetricMatrix" -> {p}MetricMatrix{bg_rules_opt}]], {{Validate::repeated, Validate::inhom}}]];',
            "    If[termComp === Null, termComp = 0];",
        ]
    )

    # Per-term plane-wave reduction: zero transverse Derivative patterns
    # and apply coordinate_values (e.g. y→π/2) BEFORE accumulation.
    if ctx.reduction is not None:
        prop_axis = ctx.reduction["propagation_axis"]
        coords = ctx.coords
        killed = [c for c in coords[1:] if c != prop_axis]
        deriv_rules: list[str] = []
        for c in killed:
            slot = coords.index(c) + 1
            deriv_rules.append(
                f"  Derivative[ords__][f_][args___] /; Length[{{ords}}] >= {slot}"
                f" && {{ords}}[[{slot}]] > 0 :> 0"
            )
        lines.append(f"    termComp = termComp /. {wl_list(','.join(deriv_rules))};")

        coord_values: dict[str, str] = ctx.reduction.get("coordinate_values", {})
        if coord_values:
            cv_rules = ", ".join(
                f"{coord}[] -> {val}" for coord, val in coord_values.items()
            )
            lines.append(f"    termComp = termComp /. {wl_list(cv_rules)};")

        # NOTE: No Expand here — R̃² 4th-order products in 4D exceed
        # $RecursionLimit in Expand (TerminatedEvaluation crashes the script).
        # The unexpanded sum is mathematically identical; downstream IBP and
        # Legendre transform handle unexpanded Plus forms correctly.

    lines.extend(
        [
            "    lagComp += termComp;",
            "    Share[];",
            '    Print["  term ", k, "/", Length[lagTerms], ": ",',
            '      Round[AbsoluteTime[] - tTerm, 0.1], "s, ",',
            '      Round[MemoryInUse[]/1024.^2], " MB"];',
            "  ],",
            "  {k, Length[lagTerms]}",
            "];",
        ]
    )

    lines.extend(
        [
            "(* Add pre-decomposed terms (already in component form from Scalar resolution) *)",
            "If[$tidalPreDecomposed =!= 0,",
            '  Print["Adding ", If[Head[$tidalPreDecomposed]===Plus, Length[$tidalPreDecomposed], 1],',
            '    " pre-decomposed terms to lagComp"];',
            "  lagComp += $tidalPreDecomposed;",
            "];",
            "Clear[lagTerms, $tidalPreDecomposed];",
            _wls_timing_end("tCanonDecomp", "Canonical Lagrangian decomposition"),
            "",
            "(* Free abstract Lagrangian — only component form needed from here *)",
            f"Clear[lagForCanon, {p}Lagrangian]; Share[];",
            _wls_mem_print("After Lagrangian decomposition"),
        ]
    )

    # Apply vector BG substitutions after decomposition
    lines.extend(_wls_vector_background_substitution(ctx, "lagComp"))

    lines.extend(
        [
            'Print["L (components): ", Short[lagComp, 5]];',
            "",
            "(* Build velocity-order pattern: {1, 0, ...} for d_t *)",
            f"coordSyms = ScalarsOfChart[{ctx.chart}];",
            "nCoords = Length[coordSyms];",
            "velOrders = Table[If[i == 1, 1, 0], {i, nCoords}];",
            "",
            "(* Normalize all Derivative arities to full dimension.  *)",
            "(* ConvertCDToDerivatives may produce mixed arities     *)",
            "(* (e.g. Derivative[1,0] and Derivative[0,0,1] in 2+1D) *)",
            "(* which causes D[] to fail matching. Pad to nCoords.   *)",
            "lagComp = lagComp /. Derivative[orders__][g_][args__] /;",
            "  Length[{orders}] < nCoords :>",
            "  Derivative[Sequence @@ PadRight[{orders}, nCoords, 0]][g][args];",
            "",
        ]
    )

    # NOTE: Plane-wave reduction (transverse Derivative zeroing + coordinate_values)
    # is now applied per-term inside the Do loop above, keeping lagComp small for
    # faster IBP + Legendre transform.  Pre-TraceBasisDummy PD zeroing (step 2.5
    # in DecomposeScalarExpression) handles the xAct-level reduction.

    lines.extend(
        [
            "(* Map component names to Wolfram function symbols *)",
            "compToFunc = <||>;",
        ]
    )

    # Build component-to-function mapping from Python field definitions.
    # For matter perturbation fields, use the perturbation head (e.g. mmpa0)
    # instead of the original head (e.g. mmpA0).
    pert_heads = _matter_pert_head_map(ctx)
    originals = _matter_pert_originals(ctx)
    for field in ctx.fields:
        fname = field["name"]
        if fname in originals:
            # Original field replaced by perturbation — skip
            continue
        head = pert_heads[fname] if fname in pert_heads else f"{p}{fname.capitalize()}"
        n_comps = _field_component_count(field, ctx.dim)
        lines.extend(f'compToFunc["{fname}_{j}"] = {head}{j};' for j in range(n_comps))

    # Include torsion field components in compToFunc (from [torsion] TOML).
    # Without this, constraint elimination and Hamiltonian construction
    # fail on torsion component references.
    if ctx.torsion is not None:
        torsion_pert_name = ctx.torsion["perturbation_name"]
        thead = f"{p}{torsion_pert_name.capitalize()}"
        # Build a field dict for _field_component_count
        torsion_field = {
            "name": torsion_pert_name,
            "type": "tensor",
            "rank": 3,
            "symmetry": "antisymmetric_23",
        }
        t_ncomps = _field_component_count(torsion_field, ctx.dim)
        lines.extend(
            f'compToFunc["{torsion_pert_name}_{j}"] = {thead}{j};'
            for j in range(t_ncomps)
        )

    # Build field function list early — needed by constraint elimination
    # and IBP (was previously only built for IBP).
    lines.extend(
        [
            "",
            "fieldFuncList = Values[compToFunc];",
        ]
    )

    # NOTE: Constraint elimination was previously called here, but for the
    # Lagrangian path fieldEquations is empty at this point (Component E-L
    # fills it later).  Constraint elimination is now deferred to the
    # post-Component-E-L block in _wls_metadata_and_export where it runs
    # on populated fieldEquations.

    return lines


def _wls_canonical_vard_eom_path(  # noqa: C901, PLR0912, PLR0914, PLR0915
    ctx: _WlsContext,
    _all_heads_str: str,
) -> list[str]:
    """VarD-based canonical path: abstract EOM → per-term component projection.

    Alternative to ``_wls_canonical_phase_a`` that avoids the scalar Lagrangian
    decomposition bottleneck by:
    1. Computing abstract EOM via ``VarD[field, CD][L^(2)]`` (fast, O(N))
    2. Projecting each EOM term to components via ``SplinterToArray``
       (rank-N tensors → ``ComponentArray`` expands free indices)
    3. Reconstructing ``lagComp`` from bilinear identity: L = 1/2 Sum q_i * EOM_i

    This follows the supervisor's pattern (SphericalEuclidean.m lines 131-428):
    VarD at abstract level → shorthand substitution → per-term projection.

    Sets up the SAME WLS variables as ``_wls_canonical_phase_a``:
    ``lagComp``, ``compToFunc``, ``fieldFuncList``, ``velOrders``,
    ``coordSyms``, ``nCoords``, ``eliminatedFromCanonical``.
    """
    p = ctx.prefix

    # --- Build compToFunc mapping (same as Phase A lines 3841-3888) ---
    # Needed BEFORE the VarD loop for lagComp reconstruction.
    pert_heads = _matter_pert_head_map(ctx)
    originals = _matter_pert_originals(ctx)

    comp_to_func_lines: list[str] = [
        "compToFunc = <||>;",
    ]
    for field in ctx.fields:
        fname = field["name"]
        if fname in originals:
            continue
        head = pert_heads[fname] if fname in pert_heads else f"{p}{fname.capitalize()}"
        n_comps = _field_component_count(field, ctx.dim)
        comp_to_func_lines.extend(
            f'compToFunc["{fname}_{j}"] = {head}{j};' for j in range(n_comps)
        )
    if ctx.torsion is not None:
        torsion_pert_name = ctx.torsion["perturbation_name"]
        thead = f"{p}{torsion_pert_name.capitalize()}"
        torsion_field = {
            "name": torsion_pert_name,
            "type": "tensor",
            "rank": 3,
            "symmetry": "antisymmetric_23",
        }
        t_ncomps = _field_component_count(torsion_field, ctx.dim)
        comp_to_func_lines.extend(
            f'compToFunc["{torsion_pert_name}_{j}"] = {thead}{j};'
            for j in range(t_ncomps)
        )
    comp_to_func_lines.extend(
        [
            "fieldFuncList = Values[compToFunc];",
            f"coordSyms = ScalarsOfChart[{ctx.chart}];",
            "nCoords = Length[coordSyms];",
            "velOrders = Table[If[i == 1, 1, 0], {i, nCoords}];",
        ]
    )

    # --- Build the list of fields to VarD ---
    _, _field_heads_str = _canonical_field_heads(ctx)

    # --- Build VarD expressions for each dynamical field ---
    # Reuses the pattern from _wls_linearize_from_lagrangian lines 2068-2098
    dyn_field_info: list[dict[str, str]] = []

    for field in ctx.fields:
        fname = field["name"]
        if fname in originals:
            continue
        head = pert_heads[fname] if fname in pert_heads else f"{p}{fname.capitalize()}"
        ftype = field["type"]
        if ftype == "scalar":
            fexpr = f"{head}[]"
        elif ftype == "vector":
            fexpr = f"{head}[-a]"
        elif ftype == "tensor":
            rank = field.get("rank", 2)
            idx = ", ".join(["-a", "-b", "-c", "-d", "-e", "-f"][:rank])
            fexpr = f"{head}[{idx}]"
        else:
            fexpr = f"{head}[-a, -b]"

        dyn_field_info.append(
            {"name": fname, "head": head, "fexpr": fexpr, "type": ftype}
        )

    # Also include torsion field
    if ctx.torsion is not None:
        tname = ctx.torsion["perturbation_name"]
        thead = f"{p}{tname.capitalize()}"
        dyn_field_info.append(
            {
                "name": tname,
                "head": thead,
                "fexpr": f"{thead}[a, -b, -c]",
                "type": "tensor",
            }
        )

    # --- Generate the VarD EOM path ---
    lines: list[str] = [
        "",
        "(* === VarD EOM canonical path (supervisor's pattern) ===               *)",
        "(* Abstract VarD → per-term SplinterToArray projection → lagComp       *)",
        "(* reconstruction. Avoids scalar Lagrangian decomposition bottleneck.  *)",
        "(* Ref: supervisor's SphericalEuclidean.m lines 131-428.              *)",
        'Print[""];',
        'Print["VarD canonical path: abstract EOM → component projection"];',
        "",
        "(* Setup: coordinate symbols, compToFunc, fieldFuncList *)",
    ]
    lines.extend(comp_to_func_lines)
    lines.extend(
        [
            "",
            f"lagForCanon = {p}Lagrangian;",
            "(* Apply SeparateMetric before VarD.                                   *)",
            "(* Tested: removing SeparateMetric reduces term count (194→148) but   *)",
            "(* creates denser terms where individual SplinterToArray calls take    *)",
            "(* >7 minutes (vs <2s with SeparateMetric). The metric distribution   *)",
            "(* across terms keeps each term's contraction complexity bounded.     *)",
            f"lagForCanon = SeparateMetric[{ctx.metric}][lagForCanon];",
            "",
            _wls_timing_start("tVarDCanon"),
            "fieldEquations = {};",
            "",
        ]
    )

    # Build BackgroundFieldRules for DecomposeToComponents
    bg_rules_entries: list[str] = []
    for bf in ctx.background_fields:
        if bf["type"] != "scalar" and bf.get("components"):
            bg_head = f"{p}{bf['name'].capitalize()}"
            comps_str = ", ".join(str(c) for c in bf["components"])
            contra_comps = _compute_contra_components(
                bf["components"], ctx.metric_diagonal
            )
            contra_str = ", ".join(contra_comps)
            bg_rules_entries.append(wl_bg_rule_entry(bg_head, comps_str, contra_str))
    bg_rules_opt = ""
    if bg_rules_entries:
        bg_rules_str = ", ".join(bg_rules_entries)
        bg_rules_opt = f', "BackgroundFieldRules" -> {wl_list(bg_rules_str)}'

    # Apply scalar BG substitutions before VarD
    lines.extend(_wls_scalar_background_substitution(ctx, "lagForCanon"))

    # Apply plane-wave Derivative zeroing to L^(2) BEFORE VarD.
    # This zeros transverse Derivative[...] patterns that exist in
    # pre-resolved Scalar[] wrappers and CD shorthand expansions.
    # VarD on the reduced L^(2) produces fewer EOM terms.
    pw_rules = _build_plane_wave_deriv_rules(ctx)
    if pw_rules:
        lines.extend(
            [
                "",
                "(* Plane-wave reduction on L^(2) before VarD — zero transverse derivatives *)",
                "lagForCanon = lagForCanon /. {",
                ",\n".join(pw_rules),
                "};",
                "lagForCanon = Expand[lagForCanon];",
                'Print["L^(2) after plane-wave reduction: ",',
                '  If[Head[lagForCanon]===Plus, Length[lagForCanon], 1], " terms"];',
                "",
            ]
        )

    # For each dynamical field, generate VarD + projection code
    for dfi in dyn_field_info:
        head = dfi["head"]
        fexpr = dfi["fexpr"]
        fname = dfi["name"]
        ftype = dfi["type"]

        lines.extend(
            [
                f'Print["  VarD w.r.t. {fname} ({fexpr})..."];',
                "Module[{eomAbs, eomTerms, nTerms, tField = AbsoluteTime[]},",
                f"  eomAbs = VarD[{fexpr}, {ctx.cd}][lagForCanon];",
                "  (* Per-term ToCanonical: avoids xperm segfault on large sums *)",
                "  eomAbs = If[Head[eomAbs] === Plus,",
                "    Total[ToCanonical /@ List @@ eomAbs],",
                "    ToCanonical[eomAbs]];",
                f"  eomAbs = ContractMetric[eomAbs, {ctx.metric}];",
                "",
                "  (* NOTE: ScreenDollarIndices and CollectTensors were benchmarked here   *)",
                "  (* but REVERTED. ScreenDollarIndices reduces term count (194→87) but    *)",
                "  (* INCREASES per-term cost (5.5s→16.7s avg) by creating denser terms.  *)",
                "  (* CollectTensors is even worse (12-35s/term). All term-density-        *)",
                "  (* reducing approaches backfire because TraceBasisDummy is O(dim^{2K})  *)",
                "  (* in contraction depth K, not in term count.                           *)",
                "",
                "  (* Apply CD shorthand rules *)",
                "  If[ListQ[$CDShorthandRules] && Length[$CDShorthandRules] > 0,",
                "    Do[eomAbs = eomAbs /. rule, {rule, Reverse[$CDShorthandRules]}]",
                "  ];",
                "",
                '  Print["    ", If[Head[eomAbs]===Plus, Length[eomAbs], 1], " abstract EOM terms"];',
            ]
        )

        # Zero transverse Derivative patterns in abstract EOM
        if pw_rules:
            lines.extend(
                [
                    "  (* Plane-wave: zero transverse Derivative patterns in EOM *)",
                    "  eomAbs = eomAbs /. {",
                    ",\n".join(f"  {r}" for r in pw_rules),
                    "  };",
                    "  eomAbs = Expand[eomAbs];",
                    '  Print["    ", If[Head[eomAbs]===Plus, Length[eomAbs], 1], " EOM terms (after plane-wave)"];',
                ]
            )

        lines.extend(
            [
                "",
            ]
        )

        # Build additionalFields list: all OTHER dynamical fields + background fields
        other_fexprs = [d["fexpr"] for d in dyn_field_info if d["name"] != fname]
        # Add background field expressions
        for bf in ctx.background_fields:
            bf_head = f"{ctx.prefix}{bf['name'].capitalize()}"
            if bf["type"] == "scalar":
                other_fexprs.append(f"{bf_head}[]")
            elif bf["type"] == "vector":
                other_fexprs.append(f"{bf_head}[-a]")
            else:
                bf_rank = bf.get("rank", 2)
                bf_idx = ", ".join(f"-{chr(97 + i)}" for i in range(bf_rank))
                other_fexprs.append(f"{bf_head}[{bf_idx}]")
        additional_str = ", ".join(other_fexprs) if other_fexprs else ""

        # Use DecomposeToComponents with TermByTerm for per-term projection
        # This is the tested, complete pipeline (ComponentDecompose.wl:1020-1154)
        lines.extend(
            [
                "  (* Decompose EOM to components via DecomposeToComponents TermByTerm *)",
                "  (* Uses the tested pipeline: SplinterToArray per-term + post-processing *)",
                "  Module[{componentEqs},",
                f"    componentEqs = DecomposeToComponents[eomAbs, {fexpr}, {ctx.chart},",
                f"      {wl_list(additional_str)},",
                '      "TermByTerm" -> True,',
                f'      "MetricMatrix" -> {p}MetricMatrix{bg_rules_opt}];',
                '    Print["    ", Length[componentEqs], " component equations"];',
                "    Do[",
                "      Module[{compName, compEq},",
                "        compName = componentEqs[[k, 1]];",
                "        compEq = componentEqs[[k, 2]];",
                "        (* Normalize derivative arity *)",
                "        compEq = compEq /. Derivative[orders__][g_][args__] /;",
                "          Length[{orders}] < nCoords :>",
                "          Derivative[Sequence @@ PadRight[{orders}, nCoords, 0]][g][args];",
                f'        AppendTo[fieldEquations, {{"{fname}_" <> ToString[compName], compEq}}];',
                "      ];,",
                "      {k, Length[componentEqs]}",
                "    ];",
                "  ];",
            ]
        )

        # Plane-wave reduction on field equations (if applicable)
        if ctx.reduction is not None:
            prop_axis = ctx.reduction["propagation_axis"]
            coords = ctx.coords
            killed = [c for c in coords[1:] if c != prop_axis]
            deriv_rules: list[str] = []
            for c in killed:
                slot = coords.index(c) + 1
                deriv_rules.append(
                    f"  Derivative[ords__][f_][args___] /; Length[{{ords}}] >= {slot}"
                    f" && {{ords}}[[{slot}]] > 0 :> 0"
                )
            if deriv_rules:
                lines.extend(
                    (
                        "  (* Plane-wave reduction: zero transverse derivatives *)",
                        "  fieldEquations = Map[Function[{eq},",
                        f"    {{eq[[1]], eq[[2]] /. {wl_list(','.join(deriv_rules))}}}],",
                        "    fieldEquations];",
                    )
                )
            coord_values: dict[str, str] = ctx.reduction.get("coordinate_values", {})
            if coord_values:
                cv_rules = ", ".join(
                    f"{coord}[] -> {val}" for coord, val in coord_values.items()
                )
                lines.append(
                    f"  fieldEquations = fieldEquations /. {wl_list(cv_rules)};",
                )

        lines.extend(
            [
                '  Print["    done in ", Round[AbsoluteTime[] - tField, 0.1], "s"];',
                "];",
                "",
            ]
        )

    # --- lagComp reconstruction from bilinear identity ---
    lines.extend(
        [
            "(* === Reconstruct lagComp from bilinear identity ===                  *)",
            "(* For quadratic L^(2) = 1/2 q^T M q, EOM = Mq, so L = 1/2 Sum q_i*EOM_i *)",
            "(* Exact for linearized (quadratic) Lagrangians.                       *)",
            "lagComp = 0;",
            "Do[",
            "  Module[{name, eom, func},",
            "    {name, eom} = fieldEquations[[k]];",
            "    If[KeyExistsQ[compToFunc, name],",
            "      func = compToFunc[name];",
            "      lagComp += func[Sequence @@ coordSyms] * eom",
            "    ];",
            "  ],",
            "  {k, Length[fieldEquations]}",
            "];",
            "lagComp = Expand[lagComp / 2];",
            'Print["lagComp (reconstructed): ", LeafCount[lagComp], " leaves, ",',
            '  If[Head[lagComp]===Plus, Length[lagComp], 1], " terms"];',
            "(* Diagnostic: check lagComp has velocity terms for Legendre transform *)",
            "Module[{d2tHeads, d1tHeads, fieldHeads},",
            "  d2tHeads = Cases[lagComp, Derivative[n_,__][f_][__] /; n>=2 :> f, {0,Infinity}] // DeleteDuplicates;",
            "  d1tHeads = Cases[lagComp, Derivative[1,___][f_][__] :> f, {0,Infinity}] // DeleteDuplicates;",
            "  fieldHeads = Cases[lagComp, (f_Symbol)[__] /; MemberQ[fieldFuncList, f] :> f, {0,Infinity}] // DeleteDuplicates;",
            '  Print["  d^2_t field heads: ", d2tHeads];',
            '  Print["  d_t (velocity) heads: ", d1tHeads];',
            '  Print["  field heads in lagComp: ", fieldHeads];',
            '  Print["  fieldFuncList: ", fieldFuncList];',
            '  Print["  Sample term: ", Short[If[Head[lagComp]===Plus, lagComp[[1]], lagComp], 2]];',
            "];",
            "lagComp = lagComp /. Derivative[orders__][g_][args__] /;",
            "  Length[{orders}] < nCoords :>",
            "  Derivative[Sequence @@ PadRight[{orders}, nCoords, 0]][g][args];",
            "",
            f"Clear[lagForCanon, {p}Lagrangian]; Share[];",
            _wls_timing_end("tVarDCanon", "VarD canonical path (EOM + lagComp)"),
            _wls_mem_print("After VarD canonical path"),
            "",
            "eliminatedFromCanonical = {};",
            "",
        ]
    )

    # Apply vector BG substitutions after decomposition
    lines.extend(_wls_vector_background_substitution(ctx, "lagComp"))

    return lines


def _wls_canonical_phase_b(ctx: _WlsContext, _all_heads_str: str) -> list[str]:
    """Generate WLS code for canonical Phase B: IBP + Legendre transform + Hamiltonian.

    Integrates by parts on the time variable to remove second-time-derivative
    terms, performs the Legendre transform to obtain the canonical Hamiltonian,
    and parses it into structured ``hamiltonianTerms``.

    Must run AFTER ``BuildMultiFieldJSONStructure`` (so it can inject the
    canonical section into the JSON).  Reads WLS variables set by Phase A:
    ``lagComp``, ``compToFunc``, ``fieldFuncList``, ``velOrders``, ``coordSyms``.

    Sets up WLS variables: ``allCompNames``, ``piCompList``, ``canonicalH``,
    ``hamiltonianTerms``.
    """
    lines: list[str] = []

    # --- Integration by parts on time variable ---
    # The Ricci scalar R contains second derivatives of the metric (∂²g),
    # so the linearized Lagrangian L^(2) has terms like h·∂²_t h.  The
    # standard Legendre transform H = pi*v - L requires L to depend only on
    # (q, ∂_t q), not on accelerations.  We integrate by parts on the time
    # variable to convert all second-time-derivative terms to first-order
    # form: f*d^2_t g -> -(d_t f)*(d_t g).  This is the component-level
    # analogue of the Gibbons-Hawking-York boundary term in GR.
    #
    # Ref: Gibbons & Hawking (1977, Phys. Rev. D 15, 2752)
    lines.extend(
        [
            "",
            "(* Increase recursion limit for IBP on higher-derivative theories *)",
            "$RecursionLimit = 65536;",
            'Print["Phase B: starting IBP (lagComp LeafCount=", LeafCount[lagComp], ", terms=", If[Head[lagComp]===Plus, Length[lagComp], 1], ")"];',
            "",
            "(* --- Integration by parts: reduce second time derivatives ---          *)",
            "(* The Ricci scalar contains d^2 g, so L^(2) has f*d^2_t g terms.       *)",
            "(* The Legendre transform requires L = L(q, dq/dt) only.  IBP on the    *)",
            "(* time variable converts f*d^2_t g -> -(d_t f)*(d_t g), the component- *)",
            "(* level analogue of the Gibbons-Hawking-York boundary term.             *)",
            "(* Ref: Gibbons & Hawking (1977, Phys. Rev. D 15, 2752)                 *)",
            "tVar = coordSyms[[1]];",
            "",
            "(* Build O(1) lookup set for field functions — MemberQ on 34-element    *)",
            "(* list causes O(n) per-call overhead; 321 terms x factors x MemberQ    *)",
            "(* creates 10K+ evaluations that corrupt Wolfram's internal heap.       *)",
            "fieldFuncSet = Association @@ Table[f -> True, {f, fieldFuncList}];",
            "isFieldFunc[f_] := KeyExistsQ[fieldFuncSet, f];",
            "",
            "(* IBP helper: for a single additive term, find the factor with          *)",
            "(* time-derivative order >= 2 and integrate by parts once.               *)",
            "(* Uses explicit Head/Part checks instead of MatchQ with /; conditions  *)",
            "(* to avoid Wolfram's recursive pattern evaluation (segfault on >300    *)",
            "(* terms).                                                               *)",
            "(* hasHighTimeDerivative: detect Derivative[n,s1,...][f][args] where       *)",
            "(* n >= 1 AND (n >= 2 OR max(s1,...) > 0). This catches:                  *)",
            "(*   - Derivative[2,0,...] = ∂²_t f  (pure second time derivative)         *)",
            "(*   - Derivative[1,1,...] = ∂_t∂_x f (mixed time-space derivative)        *)",
            "(* but NOT Derivative[1,0,...] = ∂_t f (velocity — keep for Legendre).     *)",
            "(* IBP reduces time order by 1, converting mixed ∂_t∂_x f → ∂_x f terms  *)",
            "(* in the Legendre-ready Lagrangian L(q, ∂_t q).                          *)",
            "(* Head structure: Head[Derivative[n,...][f][args]] = Derivative[n,...][f], *)",
            "(*   Head[Derivative[n,...][f]] = Derivative[n,...],                       *)",
            "(*   Head[Derivative[n,...]] = Derivative.                                 *)",
            "hasHighTimeDerivative[expr_] := Module[{h = Head[expr], dOp, orders},",
            "  dOp = Head[h];              (* Derivative[n,...] *)",
            "  If[Head[dOp] =!= Derivative, Return[False]];",
            "  If[!isFieldFunc[h[[1]]], Return[False]];",
            "  orders = List @@ dOp;",
            "  (* Only flag PURE second+ time derivatives (∂²_t or higher).          *)",
            "  (* Mixed time-space like ∂_t∂_x are valid velocity-spatial gradients  *)",
            "  (* that the Legendre transform handles natively (L depends on ∂_t q). *)",
            "  (* Previous condition orders[[1]]>=1 && (>=2 || spatial>0) incorrectly *)",
            "  (* flagged mixed derivatives, causing unbounded recursion in the       *)",
            "  (* Power expansion branch → C-level segfault for R̃² theories (#184). *)",
            "  orders[[1]] >= 2",
            "];",
            "",
            "(* dFactor: time-derivative of a single factor WITHOUT D[] recursion.    *)",
            "(* For Derivative[n,...][f][args], directly increment time order:          *)",
            "(*   Derivative[n,s1,s2,...][f][args] -> Derivative[n+1,s1,s2,...][f][args]*)",
            "(* For bare field f[args], use D[f[args], t] (trivial, no recursion).      *)",
            "(* For numeric/symbolic constants, return 0.                               *)",
            "(* This avoids Wolfram's recursive chain-rule on mixed d^2_t d^2_x terms  *)",
            "(* that causes segfaults for R-tilde-squared theories (issue #184).        *)",
            "dFactor[fj_, tv_] := Module[{h = Head[fj]},",
            "  If[Head[h] === Symbol && MemberQ[fieldFuncList, h],",
            "    (* Bare field f[args] — simple derivative *)",
            "    Derivative[Sequence @@ ReplacePart[Table[0, {Length[coordSyms]}], 1 -> 1]][h][Sequence @@ (List @@ fj)],",
            "    If[Head[Head[h]] === Derivative,",
            "      (* Already Derivative[n,...][f][args] — increment time order *)",
            "      Module[{ord = List @@ h[[0]], fn = h[[1]], ar = List @@ fj},",
            "        Derivative[Sequence @@ ReplacePart[ord, 1 -> ord[[1]] + 1]][fn][Sequence @@ ar]",
            "      ],",
            "      (* Numeric/symbolic constant or other — derivative is 0 *)",
            "      0",
            "    ]",
            "  ]",
            "];",
            "",
            "ibpOneTerm[term_] := Module[",
            "  {factors, idx, highFactor, orders, newOrders, rest, head, args},",
            "  factors = If[Head[term] === Times, List @@ term, {term}];",
            "  (* Find first factor with time-deriv >= 2 or mixed time-space derivative *)",
            "  idx = 0;",
            "  Do[",
            "    If[hasHighTimeDerivative[factors[[i]]],",
            "      idx = i; Break[]",
            "    ],",
            "    {i, Length[factors]}",
            "  ];",
            "  If[idx == 0,",
            "    (* Also check for Power[Derivative[...][f][...], 2] *)",
            "    Do[",
            "      If[Head[factors[[i]]] === Power && factors[[i, 2]] == 2 && hasHighTimeDerivative[factors[[i, 1]]],",
            "        idx = i; Break[]",
            "      ],",
            "      {i, Length[factors]}",
            "    ]",
            "  ];",
            "  If[idx == 0, Return[term]];  (* No high time deriv — keep as-is *)",
            "",
            "  highFactor = factors[[idx]];",
            "  rest = Times @@ Delete[factors, idx];",
            "",
            "  (* Handle squared case: Power[Derivative[n,...][f][args], 2] *)",
            "  (* CANNOT IBP (∂²_t f)² to lower order — it would introduce ∂³_t.  *)",
            "  (* These are Ostrogradsky kinetic terms requiring auxiliary variables.*)",
            "  (* Skip them; they remain as-is in the Hamiltonian.                  *)",
            "  If[Head[highFactor] === Power && highFactor[[2]] == 2,",
            "    Return[term]  (* Keep as-is — squared high deriv cannot be IBP'd *)",
            "  ];",
            "",
            "  (* Standard case: rest * Derivative[n, s1, s2, ...][f][args]           *)",
            "  (* In Mathematica, Derivative[n,s1,...][f][args] has Part structure:    *)",
            "  (*   expr[[0]]     = Derivative[n,s1,...][f]     (the head)             *)",
            "  (*   expr[[0,0]]   = Derivative[n,s1,...]        (the deriv operator)   *)",
            "  (*   expr[[0,1]]   = f                           (the function symbol)  *)",
            "  (*   List @@ expr  = {args}                      (the arguments)        *)",
            "  orders = List @@ highFactor[[0, 0]];     (* {n, s1, s2, ...} *)",
            "  newOrders = ReplacePart[orders, 1 -> orders[[1]] - 1];",
            "  f = highFactor[[0, 1]];                  (* function symbol, e.g. geH5 *)",
            "  args = List @@ highFactor;               (* {t[], x[], y[], z[]} *)",
            "  (* IBP: rest * D^n_t[f] -> -D_t[rest] * D^{n-1}_t[f] *)",
            "  (* Non-recursive product rule: Wolfram's D[] recurses proportionally   *)",
            "  (* to factor count in products. For R̃² with 5+ factors, depth >4096. *)",
            "  (* Explicit product rule: D[f1*f2*...*fn, t] = Σ_i (Π_{j≠i} fj)*D[fi,t] *)",
            "  (* CRITICAL FIX (#184): dFactor directly increments time derivative    *)",
            "  (* order instead of calling D[], avoiding recursive chain-rule on       *)",
            "  (* mixed d²_t d²_x terms. Derivative[n,...][f][args] → Derivative[n+1,...]. *)",
            "  Module[{restFactors, dRest},",
            "    dRest = If[Head[rest] === Plus,",
            "      Total[Function[{term},",
            "        restFactors = If[Head[term]===Times, List@@term, {term}];",
            "        Total[Table[",
            "          Times @@ Delete[restFactors, j] * dFactor[restFactors[[j]], tVar],",
            "          {j, Length[restFactors]}",
            "        ]]",
            "      ] /@ (List @@ rest)],",
            "      (* Single term *)",
            "      restFactors = If[Head[rest]===Times, List@@rest, {rest}];",
            "      Total[Table[",
            "        Times @@ Delete[restFactors, j] * dFactor[restFactors[[j]], tVar],",
            "        {j, Length[restFactors]}",
            "      ]]",
            "    ];",
            "    -dRest * Derivative[Sequence @@ newOrders][f][Sequence @@ args]",
            "  ]",
            "];",
            "",
            "(* Apply IBP to full Lagrangian (iterate until no pure d²_t terms)    *)",
            "(* No Expand inside loop — terms stay factored. FreeQ works on any    *)",
            "(* expression structure. Post-loop Expand is safe (lower order).      *)",
            "Module[{oldLag, iter = 0, maxIter = 5},",
            "  While[iter < maxIter,",
            "    oldLag = lagComp;",
            "    Module[{ibpTerms, ibpResult, nTerms},",
            "      ibpTerms = If[Head[lagComp] === Plus, List @@ lagComp, {lagComp}];",
            "      nTerms = Length[ibpTerms];",
            '      Print["  IBP round ", iter+1, ": ", nTerms, " terms"];',
            "      (* Use Table instead of AppendTo — AppendTo on growing lists is *)",
            "      (* O(n²) in Wolfram and causes memory fragmentation/segfault.   *)",
            "      (* Use Map for batch processing — safer than Table for kernel  *)",
            "      (* Per-term timeout: if a term hangs, keep original.              *)",
            "      ibpResult = Table[",
            '        If[ii > 200, Print["  IBP term ", ii, ": ", Short[ibpTerms[[ii]], 2], " hasHigh=", AnyTrue[Cases[ibpTerms[[ii]], Derivative[__][_][__], {0,Infinity}], hasHighTimeDerivative]]];',
            "        TimeConstrained[ibpOneTerm[ibpTerms[[ii]]], 30, ibpTerms[[ii]]],",
            "        {ii, Length[ibpTerms]}",
            "      ];",
            "      lagComp = Total[ibpResult];",
            "    ];",
            "    iter++;",
            "    (* Convergence check: extract ALL Derivative subexpressions        *)",
            "    (* structurally (no /; conditions), then test each individually.   *)",
            "    (* CRITICAL: FreeQ[expr, x_ /; cond] with conditions causes        *)",
            "    (* segfault on 261-term R̃² expressions (C-level stack overflow    *)",
            "    (* from deep pattern recursion). Structural Cases[] is safe.       *)",
            "    Module[{allDerivs, hasHigh},",
            "      allDerivs = Cases[lagComp, Derivative[__][_][__], {0, Infinity}];",
            "      hasHigh = AnyTrue[allDerivs, hasHighTimeDerivative];",
            "      If[!hasHigh, Break[]]",
            "    ];",
            "  ];",
            "  If[iter > 0,",
            '    Print["[IBP] Applied ", iter, " round(s) of time-variable IBP"];',
            "  ];",
            "  (* Post-loop check: same structural extraction approach *)",
            "  Module[{allDerivs, hasHigh},",
            "    allDerivs = Cases[lagComp, Derivative[__][_][__], {0, Infinity}];",
            "    hasHigh = AnyTrue[allDerivs, hasHighTimeDerivative];",
            "    If[hasHigh,",
            '      Print["WARNING: High time-derivative terms remain after IBP (",',
            '        maxIter, " iterations)."];',
            "    ]",
            "  ];",
            "];",
            "(* Post-loop Expand: safe because IBP reduced derivative order *)",
            "(* NO post-IBP Expand — R̃² 4th-order products exceed $RecursionLimit *)",
            "(* Downstream ParseHamiltonianExpression handles unexpanded Plus forms. *)",
            'Print["L after time-IBP: ", If[Head[lagComp]===Plus, Length[lagComp], 1], " terms"];',
            'Print["Proceeding to Legendre transform..."];',
            "",
        ]
    )

    lines.extend(
        [
            "",
            "(* Compute canonical momenta: pi_i = dL/d(d_t q_i) *)",
            _wls_timing_start("tLegendre"),
            "allCompNames = fieldEquations[[All, 1]];",
            "piCompList = {};",
            "canonicalH = 0;",
            "Do[",
            "  Module[{compName, compFunc, vel, piComp},",
            "    compName = allCompNames[[k]];",
            "    compFunc = compToFunc[compName];",
            "    vel = Derivative[Sequence @@ velOrders][compFunc][Sequence @@ coordSyms];",
            "    piComp = D[lagComp, vel];",
            "    AppendTo[piCompList, {compName, piComp}];",
            "    canonicalH += piComp * vel;",
            '    Print["pi(", compName, "): ", piComp];',
            "  ],",
            "  {k, Length[allCompNames]}",
            "];",
            "",
            "(* Legendre transform: H = Sigma pi_i * vel_i - L *)",
            "(* No Expand: R̃² 4th-order products exceed $RecursionLimit in Expand *)",
            "canonicalH = canonicalH - lagComp;",
            _wls_timing_end("tLegendre", "Legendre transform (momenta + H)"),
            _wls_mem_print("After Legendre transform"),
            'Print["H (components): ", Short[canonicalH, 5]];',
            "",
            "(* --- Spatial IBP on the Hamiltonian --- *)",
            "(* Convert f * d^2_x g -> -(d_x f)*(d_x g) for ALL spatial axes.      *)",
            "(* The Hamiltonian from the Legendre transform can contain f·∇² g       *)",
            "(* terms alongside |∇f|² terms.  These must CANCEL via IBP to give      *)",
            "(* the correct coefficient on |∇f|².  Evaluating f·∇²g numerically      *)",
            "(* loses precision because the two large, nearly-cancelling terms are    *)",
            "(* computed independently.  Symbolic IBP eliminates this issue.          *)",
            "(* Ref: energy conservation requires all spatial-derivative terms in     *)",
            "(* |∇f|² form for stable numerical evaluation.                          *)",
            "Module[{spatialVars, ibpSpatialOneTerm, oldH, iter = 0, maxIter = 5},",
            "  spatialVars = Rest[coordSyms];",
            "",
            "  ibpSpatialOneTerm[term_, sVar_] := Module[",
            "    {factors, idx, highFactor, orders, newOrders, rest, f, args, posVar,",
            "     otherFieldFactors},",
            "    factors = If[Head[term] === Times, List @@ term, {term}];",
            "    posVar = Position[coordSyms, sVar][[1, 1]];",
            "    (* Find factor with spatial-deriv order >= 2 in sVar *)",
            "    idx = 0;",
            "    Do[",
            "      If[MatchQ[factors[[i]],",
            "           Derivative[__][g_][__] /; isFieldFunc[g]] &&",
            "         (List @@ factors[[i, 0, 0]])[[posVar]] >= 2,",
            "        idx = i; Break[]",
            "      ],",
            "      {i, Length[factors]}",
            "    ];",
            "    If[idx == 0, Return[term]];",
            "",
            "    (* Only IBP if the OTHER field factor has zero spatial derivative *)",
            "    (* (i.e., it is 'identity' or 'time_derivative' — no spatial deriv). *)",
            "    (* This prevents oscillation: ∂f·∂f ↔ f·∂²f.                       *)",
            "    otherFieldFactors = Select[Delete[factors, idx],",
            "      MatchQ[#, Derivative[__][g_][__] /; isFieldFunc[g]] &];",
            "    If[Length[otherFieldFactors] > 0,",
            "      Module[{otherOrds = List @@ otherFieldFactors[[1, 0, 0]]},",
            "        If[otherOrds[[posVar]] > 0, Return[term]]",
            "      ]",
            "    ];",
            "",
            "    highFactor = factors[[idx]];",
            "    rest = Times @@ Delete[factors, idx];",
            "    orders = List @@ highFactor[[0, 0]];",
            "    newOrders = ReplacePart[orders, posVar -> orders[[posVar]] - 1];",
            "    f = highFactor[[0, 1]];",
            "    args = List @@ highFactor;",
            "    (* IBP: ∫ rest · D^n f dz → -∫ D[rest,z] · D^{n-1}f dz             *)",
            "    (* (total derivative vanishes for periodic BCs)                      *)",
            "    -D[rest, sVar] * Derivative[Sequence @@ newOrders][f][Sequence @@ args]",
            "  ];",
            "",
            "  (* Iterate IBP over all spatial variables until no ∂²_x terms remain *)",
            "  While[iter < maxIter,",
            "    oldH = canonicalH;",
            "    Do[",
            "      canonicalH = Total[ibpSpatialOneTerm[#, sv] & /@ ",
            "        If[Head[canonicalH] === Plus, List @@ canonicalH, {canonicalH}]];",
            "    , {sv, spatialVars}];",
            "    iter++;",
            "    If[canonicalH === oldH, Break[]];",
            "  ];",
            "  If[iter > 0,",
            '    Print["[IBP-spatial] Applied ", iter, " round(s) of spatial IBP on H"];',
            "  ];",
            "];",
            'Print["H after spatial IBP: ", If[Head[canonicalH]===Plus, Length[canonicalH], 1], " terms"];',
            "",
            "(* Resolve abstract tensor traces in H (e.g., h^a_a = Tr(h)).          *)",
            "(* These arise from EH sector trace terms that survived decomposition.    *)",
            "(* Replace h[a_, -a_] with the explicit metric trace using MetricMatrix:  *)",
            "(*   Tr(h) = Sum_i MetricMatrix[[i+1,i+1]] * h_{ii}(coords)              *)",
            "(* Theory-agnostic: uses MetricMatrix diagonal for any signature.         *)",
            "If[!FreeQ[canonicalH, _?xTensorQ],",
            "  Module[{tensorHeads},",
            "    tensorHeads = Cases[canonicalH, (f_?xTensorQ)[_, -_] :> f, {0, Infinity}] // DeleteDuplicates;",
            '    Print["Resolving ", Length[tensorHeads], " abstract tensor trace(s) in H"];',
            "    Do[Module[{compTuples, diagPositions, traceSum, trRule},",
            "      compTuples = EnumerateComponentTuples[th, nCoords];",
            "      diagPositions = Select[Range[Length[compTuples]],",
            "        Length[compTuples[[#]]] == 2 && compTuples[[#, 1]] == compTuples[[#, 2]] &];",
            "      If[Length[diagPositions] > 0,",
            "        traceSum = Total[Table[",
            "          Module[{diagIdx = compTuples[[dp, 1]], compName, func},",
            '            compName = StringJoin[ToString[th], "_", ToString[dp - 1]];',
            "            func = compToFunc[compName];",
            "            If[func =!= compToFunc[compName],",
            f"              {ctx.prefix}MetricMatrix[[diagIdx+1, diagIdx+1]] * func[Sequence @@ coordSyms],",
            "              0",
            "            ]",
            "          ],",
            "          {dp, diagPositions}]];",
            "        trRule = th[idx_, -idx_] :> traceSum;",
            "        canonicalH = canonicalH /. trRule;",
            '        Print["  Resolved Tr(", th, ") -> ", Short[traceSum, 2]];',
            "      ];",
            "    ], {th, tensorHeads}];",
            "    canonicalH = Expand[canonicalH];",
            '    Print["H after trace resolution: ", If[Head[canonicalH]===Plus, Length[canonicalH], 1], " terms"];',
            "  ];",
            "];",
            "",
            "(* Parse H into structured quadratic terms *)",
            _wls_timing_start("tParseH"),
            "hamiltonianTerms = ParseHamiltonianExpression[canonicalH, allCompNames];",
            _wls_timing_end("tParseH", "ParseHamiltonianExpression"),
            'Print["Hamiltonian terms: ", Length[hamiltonianTerms]];',
            "",
        ]
    )
    return lines


def _wls_volume_element_code(ctx: _WlsContext) -> list[str]:
    """Generate Wolfram code to compute spatial volume element.

    When plane-wave reduction is active, computes a **factored** volume
    element: for each diagonal metric entry, drop multiplicative factors
    that depend on killed coordinates.  This gives the correct reduced
    volume element for energy integration (e.g. r² for spherical → radial).

    When no reduction, computes the standard ``√|det(g_spatial)|``.
    """
    p = ctx.prefix
    if ctx.reduction is None:
        return [
            "(* Compute spatial volume element sqrt|det(g_spatial)| for energy integration *)",
            f"sqrtDetGSpatial = Simplify[Sqrt[Abs[Det[{p}MetricMatrix[[2;;, 2;;]]]]]];",
            'Print["sqrt|g_spatial|: ", sqrtDetGSpatial];',
        ]

    # Plane-wave reduction: factored volume element
    prop_axis = ctx.reduction["propagation_axis"]
    coords = ctx.coords
    killed = [c for c in coords[1:] if c != prop_axis]
    killed_vars = ", ".join(f"{c}[]" for c in killed)

    return [
        "(* Compute REDUCED volume element for plane-wave reduction *)",
        "(* For each diagonal metric entry, drop factors depending on killed coordinates *)",
        f"spatialMetric = {p}MetricMatrix[[2;;, 2;;]];",
        f"killedVars = {wl_list(killed_vars)};",
        "reducedDet = 1;",
        "Do[",
        "  gii = spatialMetric[[k, k]];",
        "  (* Factorize and keep only factors free of killed variables *)",
        "  factors = If[Head[gii] === Times, List @@ gii, {gii}];",
        "  survivingFactors = Select[factors,",
        "    FreeQ[#, Alternatives @@ killedVars] &];",
        "  reducedDet = reducedDet * Times @@ survivingFactors,",
        "  {k, Length[spatialMetric]}",
        "];",
        "sqrtDetGSpatial = Simplify[Sqrt[Abs[reducedDet]]];",
        f'Print["Reduced volume element (killed: {", ".join(killed)}): ", sqrtDetGSpatial];',
    ]


def _wls_json_plane_wave_reduction(ctx: _WlsContext) -> list[str]:
    """Generate WLS code to remap JSON structure from ND to 1+1D.

    After ``BuildMultiFieldJSONStructure`` + canonical injection, the JSON
    uses the original ND coordinate names (e.g. ``laplacian_z``,
    ``z[]``).  This function remaps:

    1. Operator names: ``laplacian_{prop} → laplacian_x``, etc.
    2. Coordinate references in ``coefficient_symbolic``: ``{prop}[] → x[]``
    3. Bare coordinate names in ``coordinate_dependent``: ``["z"] → ["x"]``
    4. Spacetime metadata: dimension → 2, signature → [-1,1], coordinates → [t,x]
    5. Reduction provenance in metadata

    All replacement logic is self-contained in Wolfram — Python only passes
    the propagation axis name and spatial coordinate list.  Coordinate-dependent
    arrays are remapped at the Association level (before serialization).
    Operator names and symbolic coordinate references are remapped via
    ``ExportString → StringReplace`` on the JSON string.
    """
    if ctx.reduction is None:
        return []

    prop = ctx.reduction["propagation_axis"]  # e.g. "z"
    spatial = [c for c in ctx.coords if c != "t"]
    spatial_str = "{" + ", ".join(f'"{c}"' for c in spatial) + "}"

    return [
        "",
        "(* === Plane-wave reduction: remap JSON to 1+1D === *)",
        f'pwPropAxis = "{prop}";',
        f"pwSpatialAxes = {spatial_str};",
        "pwKilledAxes = DeleteCases[pwSpatialAxes, pwPropAxis];",
        "",
        'Print["Remapping JSON: " <> ToString[Length[pwSpatialAxes]+1]',
        '  <> "D → 2D (" <> pwPropAxis <> " → x)"];',
        "",
        "(* --- Update spacetime metadata --- *)",
        'jsonStructure["spacetime", "dimension"] = 2;',
        'jsonStructure["spacetime", "signature"] = {-1, 1};',
        'jsonStructure["spacetime", "coordinates"] = {"t", "x"};',
        "",
        "(* --- Store reduction provenance --- *)",
        'jsonStructure["metadata", "reduction"] = <|',
        '  "type" -> "plane_wave",',
        f'  "original_dimension" -> {ctx.dim},',
        '  "propagation_axis" -> pwPropAxis,',
        '  "eliminated_fields" -> eliminatedFromCanonical',
        "|>;",
        "",
        "(* --- Remap coordinate_dependent arrays in the Association --- *)",
        "(* Fix coordinate names BEFORE ExportString, so no fragile     *)",
        "(* string-level patching is needed.                            *)",
        "pwCoordMap = {};",
        'If[pwPropAxis =!= "x", AppendTo[pwCoordMap, pwPropAxis -> "x"]];',
        'Do[If[k =!= "x", AppendTo[pwCoordMap, k -> "x"]], {k, pwKilledAxes}];',
        "",
        "(* Remap coordinate_dependent in equations — use [[...]] Part    *)",
        '(* syntax throughout; ["key"][[i]] is a function call, not Part. *)',
        "Do[",
        "  Module[{nTerms, cd},",
        '    nTerms = Length[jsonStructure[["equations", i, "rhs", "terms"]]];',
        "    Do[",
        '      cd = jsonStructure[["equations", i, "rhs", "terms", j, "coordinate_dependent"]];',
        "      If[ListQ[cd] && Length[cd] > 0,",
        '        jsonStructure[["equations", i, "rhs", "terms", j, "coordinate_dependent"]] = cd /. pwCoordMap],',
        "      {j, nTerms}",
        "    ]",
        "  ],",
        '  {i, Length[jsonStructure[["equations"]]]}',
        "];",
        "",
        "(* Remap coordinate_dependent in hamiltonian_terms *)",
        'If[KeyExistsQ[jsonStructure, "canonical"],',
        "  Module[{nTerms, cd},",
        '    nTerms = Length[jsonStructure[["canonical", "hamiltonian_terms"]]];',
        "    Do[",
        '      cd = jsonStructure[["canonical", "hamiltonian_terms", j, "coordinate_dependent"]];',
        "      If[ListQ[cd] && Length[cd] > 0,",
        '        jsonStructure[["canonical", "hamiltonian_terms", j, "coordinate_dependent"]] = cd /. pwCoordMap],',
        "      {j, nTerms}",
        "    ]",
        "  ]",
        "];",
        "",
        "(* --- Build string replacement rules for operators and coords --- *)",
        "pwStringRules = {};",
        "",
        "(* Operator renaming: prop_axis → x *)",
        "(* Standard operators: laplacian_{axis} (underscore-separated) *)",
        "Do[",
        "  AppendTo[pwStringRules,",
        '    pfx <> "_" <> pwPropAxis -> pfx <> "_x"],',
        '  {pfx, {"laplacian", "gradient", "first_derivative", "derivative_3",',
        '         "cross_derivative"}}',
        "];",
        "(* Mixed operators: mixed_T{n}_S{m}{axis} (axis appended directly) *)",
        "Do[",
        "  AppendTo[pwStringRules,",
        '    pfx <> pwPropAxis -> pfx <> "x"],',
        '  {pfx, {"mixed_T1_S1", "mixed_T2_S1", "mixed_T2_S2",',
        '         "mixed_T3_S1"}}',
        "];",
        "",
        "(* Coordinate references in symbolic expressions: prop[] → x[] *)",
        'If[pwPropAxis =!= "x",',
        '  AppendTo[pwStringRules, pwPropAxis <> "[]" -> "x[]"]',
        "];",
        'Do[AppendTo[pwStringRules, k <> "[]" -> "x[]"], {k, pwKilledAxes}];',
        "",
        "(* Apply string rules to the JSON text (operators + symbolic coords) *)",
        'jsonStringFinal = ExportString[jsonStructure, "JSON"];',
        "jsonStringFinal = StringReplace[jsonStringFinal, pwStringRules];",
        "",
        'Print["JSON remapped to 1+1D."];',
        "",
    ]


def _wls_canonical_injection(ctx: _WlsContext) -> list[str]:
    """Generate WLS code to inject canonical structure into JSON.

    Validates Hamiltonian terms, computes volume element, and injects
    the canonical section (``hamiltonian_terms``, ``volume_element``,
    ``eliminated_from_canonical``) into ``jsonStructure``.

    Must run after Phase B has computed ``hamiltonianTerms``.
    """
    return [
        "(* === E-L Velocity Form: Inject Canonical Structure === *)",
        "(* E-L equations are preserved as-is in equations[] array. *)",
        "(* Only hamiltonian_terms are injected for energy measurement. *)",
        "",
        "(* Validate that Hamiltonian computation succeeded *)",
        "If[!ListQ[hamiltonianTerms] || Length[hamiltonianTerms] === 0,",
        '  Print["ERROR: Canonical Hamiltonian computation produced no terms."];',
        '  Print["This is required for correct energy measurement."];',
        '  Print["Check that the Lagrangian has quadratic kinetic terms."];',
        "  Exit[1]",
        "];",
        "",
        *_wls_volume_element_code(ctx),
        "",
        "(* Inject canonical structure into JSON *)",
        "canonicalSection = <|",
        '  "hamiltonian_terms" -> hamiltonianTerms',
        "|>;",
        "(* Only include volume_element when non-trivial (curved coordinates) *)",
        "If[sqrtDetGSpatial =!= 1,",
        '  canonicalSection["volume_element"] = ToString[sqrtDetGSpatial, InputForm]',
        "];",
        "(* Record Wolfram-side constraint elimination *)",
        "If[Length[eliminatedFromCanonical] > 0,",
        '  canonicalSection["wolfram_constraint_elimination"] = True;',
        '  canonicalSection["eliminated_from_canonical"] = eliminatedFromCanonical',
        "];",
        'jsonStructure["canonical"] = canonicalSection;',
        "",
        'Print["Canonical structure (hamiltonian_terms only) injected into JSON."];',
        'Print["E-L equations preserved (no Hamilton equation injection)."];',
        'Print[""];',
        "",
    ]


# --- WLS: Metadata & JSON export ---


def _wls_metadata_and_export(  # noqa: C901, PLR0912, PLR0914, PLR0915
    config: dict[str, Any],
    ctx: _WlsContext,
) -> list[str]:
    """Generate metadata and JSON export lines.

    Raises
    ------
    ValueError
        If dimension > 4 and no explicit signature provided.
    """
    sig_default = _MINKOWSKI_SIGNATURES.get(ctx.dim)
    if sig_default is None and "signature" not in config.get("spacetime", {}):
        msg = f"Dimension {ctx.dim}: must specify [spacetime].signature (no default for dim > 4)"
        raise ValueError(msg)
    sig_str = ", ".join(
        str(s) for s in config["spacetime"].get("signature", sig_default or [])
    )
    coord_str = ", ".join(f'"{c}"' for c in ctx.coords)

    is_linearization = ctx.linearization is not None

    if is_linearization and ctx.lagrangian_expr:
        # Lagrangian-first linearization: store the Lagrangian (not the EOM)
        raw_expr: str = ctx.lagrangian_expr
        linearized_str = "True"
    elif is_linearization and ctx.linearization is not None:
        # Legacy: direct EOM linearization
        raw_expr = str(ctx.linearization["expression"])
        linearized_str = "True"
    else:
        raw_expr = ctx.lagrangian_expr
        linearized_str = "False"

    escaped_expr = (
        raw_expr.replace("\\", "\\\\").replace('"', '\\"').replace("\n", " ").strip()
    )

    # Build gauge metadata string
    if ctx.gauge:
        gauge_parts: list[str] = []
        for g in ctx.gauge:
            if g["type"] == "custom":
                gauge_parts.append(f"custom({g['field']})")
            else:
                gauge_parts.append(f"{g['type']}({g['field']})")
        gauge_val = "+".join(gauge_parts)
    else:
        gauge_val = "none"

    lines: list[str] = [
        "metadata = <|",
        '  "source" -> "xAct",',
        f'  "lagrangian_expr" -> "{escaped_expr}",',
        '  "derived_from" -> "Euler-Lagrange",',
        f'  "gauge" -> "{gauge_val}",',
        f'  "linearized" -> {linearized_str},',
        f'  "dimension" -> {ctx.dim},',
        f'  "signature" -> {{{sig_str}}},',
        f'  "coordinates" -> {{{coord_str}}},',
        f'  "metric_type" -> "{ctx.metric_type}"',
    ]

    # Add constraint solver metadata if configured in TOML
    if ctx.constraint_solver is not None:
        # Spatial coordinates are all coords except "t"
        spatial_coords = [c for c in ctx.coords if c != "t"]
        cs_lines = _wls_constraint_metadata(ctx.constraint_solver, spatial_coords)
        # The last base metadata line needs a trailing comma
        lines[-1] += ","
        lines.extend(cs_lines)
        # Strip trailing comma from last line to avoid Null in Association
        lines[-1] = lines[-1].removesuffix(",")

    lines.extend(["|>;", ""])

    # Ensure eliminatedFromCanonical is always defined — plane-wave reduction
    # metadata references it even when no canonical pipeline runs.
    lines.extend(("eliminatedFromCanonical = {};", ""))

    # --- Canonical Phase A: Lagrangian decomposition + constraint elimination ---
    # Must run BEFORE BuildMultiFieldJSONStructure so that fieldEquations
    # contains only surviving fields (constraints already eliminated).
    all_heads_str = ""
    if ctx.lagrangian_expr:
        _, all_heads_str = _canonical_field_heads(ctx)
        # Memory preparation before canonical pipeline.  The EOM pass's
        # cached kernel state (Christoffel symbols, metric DownValues,
        # background field DownValues) is reused by DecomposeScalarExpression.
        # Share[] deduplicates subexpressions to reclaim memory for the
        # Lagrangian decomposition while preserving these caches.
        lines.extend(
            (
                "(* Memory cleanup before canonical pipeline — preserves cached *)",
                "(* Christoffels, metric DownValues, background field DownValues *)",
                "Share[];",
                _wls_mem_print("Before canonical pipeline"),
                "",
            )
        )
        # --- Dual-path canonical pipeline: VarD EOM (fast) with Phase A fallback ---
        # The VarD path projects rank-N EOM to components (supervisor's pattern).
        # Phase A decomposes scalar L to components (slower for complex theories).
        # Both paths produce the same WLS variables: lagComp, fieldEquations,
        # compToFunc, fieldFuncList, velOrders, coordSyms, nCoords.
        # Runtime dispatch: try VarD first, fall back to Phase A if it fails.
        has_tensor = (
            any(
                f["type"] != "scalar"
                for f in ctx.fields
                if f["name"] not in _matter_pert_originals(ctx)
            )
            or ctx.torsion is not None
        )

        if has_tensor:
            # Generate both VarD path and Phase A, with runtime LeafCount dispatch.
            # Small L^(2) → Phase A (fast, produces correct lagComp for Hamiltonian).
            # Large L^(2) → VarD path (avoids timeout on expensive decomposition).
            # Threshold: LeafCount > 50000 indicates complex theory where Phase A
            # would timeout (e.g., non-minimal R̃[μν]F coupling: LeafCount ~119K).
            lines.extend(
                [
                    "(* === Runtime dispatch: Phase A vs VarD based on L^(2) complexity === *)",
                    "(* Heuristic: LeafCount of unexpanded L^(2) + whether L^(2) has many *)",
                    "(* top-level additive terms. Large LeafCount indicates complex tensor *)",
                    "(* products that Phase A DecomposeScalarExpression handles slowly.    *)",
                    "(* VarD projects rank-N EOM (not scalar L) — fundamentally faster.    *)",
                    "(* No Expand needed for dispatch — just Length + LeafCount.           *)",
                    f"Module[{{l2LC = LeafCount[{ctx.prefix}Lagrangian],",
                    f"         l2TopTerms = If[Head[{ctx.prefix}Lagrangian] === Plus,",
                    f"           Length[{ctx.prefix}Lagrangian], 1]}},",
                    '  Print["L^(2) canonical dispatch: ", l2TopTerms, " top terms, LeafCount=", l2LC];',
                    "  If[l2TopTerms > 500 || l2LC > 50000,",
                    '    Print["Complex L^(2) — using VarD EOM path"];',
                ]
            )
            # VarD path (for large expressions)
            vard_lines = _wls_canonical_vard_eom_path(ctx, all_heads_str)
            lines.extend(f"    {line}" for line in vard_lines)
            lines.extend(
                [
                    "    $tidalVarDSuccess = True,",
                    "",
                    '    Print["Few terms — using Phase A (Lagrangian decomposition)"];',
                ]
            )
            # Phase A (for small expressions)
            phase_a_lines = _wls_canonical_phase_a(ctx, all_heads_str)
            lines.extend(f"    {line}" for line in phase_a_lines)
            lines.extend(
                [
                    "    $tidalVarDSuccess = False",
                    "  ];",
                    "];",
                    "",
                ]
            )
        else:
            # Scalar-only: use Phase A directly
            lines.extend(_wls_canonical_phase_a(ctx, all_heads_str))

        # --- Component-level E-L (only needed for Phase A fallback path) ---
        # VarD path produces fieldEquations directly; Phase A needs Component E-L.
        if ctx.lagrangian_expr:
            lines.extend(
                [
                    "",
                    "(* Skip Component E-L if VarD path already produced fieldEquations *)",
                    "If[!TrueQ[$tidalVarDSuccess],",
                    "",
                    "(* ============================================================ *)",
                    "(* Component-level Euler-Lagrange: recompute fieldEquations     *)",
                    "(* from lagComp (component Lagrangian). This replaces the       *)",
                    "(* trivial equations from abstract VarD + DecomposeToComponents *)",
                    "(* which fail for R̃² torsion theories.                        *)",
                    "(* ============================================================ *)",
                    _wls_timing_start("tCompEL"),
                    'Print["Component E-L: computing equations from component Lagrangian..."];',
                    'Print["  lagComp LeafCount: ", LeafCount[lagComp]];',
                    "",
                    "(* lagComp should now be fully resolved — ReplaceHigherRankField- *)",
                    "(* Components handles all symmetry permutations and Derivative    *)",
                    "(* forms since the symmetry-aware rewrite.                        *)",
                    "",
                    "(* Build field function list from lagComp.                     *)",
                    "(* Match BOTH bare f[args] AND Derivative[...][f][args] forms.  *)",
                    "(* Fields that appear only in derivatives (e.g. photon in F²)   *)",
                    "(* would be missed by the bare pattern alone.                    *)",
                    "fieldFuncsEL = Join[",
                    "  Cases[lagComp,",
                    "    (f_Symbol)[__] /; MemberQ[fieldFuncList, f] :> f,",
                    "    {0, Infinity}],",
                    "  Cases[lagComp,",
                    "    Derivative[__][f_Symbol][__] /; MemberQ[fieldFuncList, f] :> f,",
                    "    {0, Infinity}]",
                    "] // DeleteDuplicates // Sort;",
                    'Print["  Field functions: ", Length[fieldFuncsEL], " — ", fieldFuncsEL];',
                    "",
                    "(* Compute Euler-Lagrange equations *)",
                    "eomListEL = ComponentEulerLagrange[lagComp, fieldFuncsEL, coordSyms];",
                    "",
                    "(* Replace fieldEquations with component E-L results *)",
                    "fieldEquations = {};",
                    "Do[",
                    "  Module[{fname, eqExpr, nameRules},",
                    "    fname = ToString[fieldFuncsEL[[i]]];",
                    "    (* Map function name to field label *)",
                    "    nameRules = {};",
                    "    Scan[Function[{key},",
                    "      Module[{head = ToString[compToFunc[key]], fld = key},",
                    "        If[fname === head, AppendTo[nameRules, fname -> fld]]",
                    "      ]],",
                    "      Keys[compToFunc]",
                    "    ];",
                    "    If[Length[nameRules] > 0, fname = fname /. nameRules];",
                    "    eqExpr = eomListEL[[i]];",
                    "    (* Resolve any remaining tensor refs created by D[] *)",
                    f"    Do[eqExpr = ReplaceTensorFieldComponents[eqExpr, afh, {ctx.chart}, coordSyms, nCoords],",
                    f"      {wl_list(f'afh, {wl_list(all_heads_str)}')}];",
                    "    AppendTo[fieldEquations, {fname, eqExpr}];",
                    "  ],",
                    "  {i, Length[fieldFuncsEL]}",
                    "];",
                    _wls_timing_end("tCompEL", "Component E-L (differentiate lagComp)"),
                    'Print["Component E-L: ", Length[fieldEquations], " equations"];',
                    "",
                    "]; (* End If[!TrueQ[$tidalVarDSuccess]] *)",
                    "",
                ]
            )

    # --- Post-Component-E-L pipeline for linearization+Lagrangian path ---
    # Type B gauge, plane-wave reduction, and constraint elimination MUST
    # run AFTER Component E-L fills fieldEquations.  Previously these ran
    # inside _wls_linearize_from_lagrangian (on empty fieldEquations) or
    # inside _wls_canonical_phase_a (also empty).  Now deferred to here.
    #
    # Ordering:
    #   1. Type B gauge — creates TT constraint eqs for constraint elimination
    #   2. Plane-wave reduction — zeros transverse derivs, eliminates zero fields
    #   3. Constraint elimination — detects/removes constraint fields from both
    #      fieldEquations and lagComp
    if ctx.lagrangian_expr and ctx.linearization is not None:
        has_type_b_post = bool(ctx.gauge) and any(
            _resolve_gauge_mechanism(g) == "constraint" for g in ctx.gauge
        )
        if has_type_b_post:
            lines.extend(_wls_gauge_fixing_type_b(ctx))
        if ctx.reduction is not None:
            lines.extend(_wls_plane_wave_reduction_equations(ctx))
            lines.extend(_wls_plane_wave_field_elimination(ctx))
            lines.extend(_wls_plane_wave_coordinate_evaluation(ctx))
        lines.extend(_wls_constraint_elimination())

    # Build JSON — always use multi-field builder since fieldEquations
    # Inject tensor component metadata into metadata for LaTeX export
    lines.extend(
        (
            "(* Inject tensor component metadata for LaTeX export *)",
            "If[Length[componentMetadata] > 0,",
            '  metadata["component_metadata"] = componentMetadata;',
            "];",
            "",
        )
    )

    # is constructed with proper labels by both single and multi-field paths.
    # fieldEquations now contains only surviving fields (constraints
    # eliminated by Phase A above), so the JSON is born correct.
    #
    # Substitute parameter values before export: R̃² 4th-order equations
    # have symbolic coefficients (e.g. 5*b5/kappa^2) that cause Expand
    # inside BuildMultiFieldJSONStructure to exceed $RecursionLimit.
    # Numeric substitution collapses these to floats, avoiding deep recursion.
    param_rules = []
    if ctx.parameters:
        for pname, pval in ctx.parameters.items():
            # Constants are DefConstantSymbol (no prefix), not prefixed fields
            param_rules.append(f"{pname} -> {pval}")
    if param_rules:
        lines.extend(
            (
                "(* NOTE: Parameter substitution deferred to after Phase B.           *)",
                "(* Setting DefConstantSymbol values before IBP/Legendre crashes the   *)",
                "(* Wolfram kernel (segfault) — xAct's evaluation engine conflicts     *)",
                "(* with numeric constant symbols during pattern matching.              *)",
                "",
            )
        )
    lines.extend(
        (
            "jsonStructure = BuildMultiFieldJSONStructure[fieldEquations, metadata];",
            "",
        )
    )

    # Inject runtime parameter defaults into JSON metadata
    if ctx.parameters:
        param_entries = ", ".join(f'"{k}" -> {v}' for k, v in ctx.parameters.items())
        lines.extend(
            (
                f'jsonStructure["metadata", "parameters"] = <|{param_entries}|>;',
                "",
            )
        )

    # --- Canonical Phase B: IBP + Legendre transform + Hamiltonian ---
    # Must run AFTER BuildMultiFieldJSONStructure (injects canonical section
    # into jsonStructure).  Reads allCompNames from fieldEquations.
    # Skip for R̃² torsion: IBP fails on mixed d²_t d²_x terms (D[] recurses
    # >4096 deep on nested R̃² expressions). Phase A's ParseHamiltonian-
    # Expression already provides the quadratic Hamiltonian (14 terms) which
    # is sufficient for energy measurements. Full canonical Hamiltonian
    # (Ostrogradsky-reduced) is tracked in issue #158.
    if ctx.lagrangian_expr:
        lines.extend(_wls_canonical_phase_b(ctx, all_heads_str))
        lines.extend(_wls_canonical_injection(ctx))

    # Free fieldEquations now that both JSON structure and canonical
    # pipeline have finished using it.
    lines.extend(("Clear[fieldEquations]; Share[];", ""))

    # --- Deferred parameter substitution (AFTER Phase B) ---
    # Setting DefConstantSymbol values BEFORE IBP/Legendre crashes the
    # Wolfram kernel — xAct's evaluation conflicts with numeric constants.
    # Substitute parameters in the JSON structure AFTER canonical processing.
    if ctx.parameters:
        param_entries_wl = ", ".join(f"{k} -> {v}" for k, v in ctx.parameters.items())
        lines.extend(
            (
                "(* Substitute parameter values in JSON structure *)",
                f"paramRules = {wl_list(param_entries_wl)};",
                "jsonStructure = jsonStructure /. paramRules;",
                'Print["Parameter substitution applied to JSON structure."];',
                "",
            )
        )

    # --- Plane-wave reduction: remap JSON to 1+1D ---
    # After BuildMultiFieldJSONStructure + canonical injection, the JSON
    # still uses the original ND coordinate names (e.g. laplacian_z).
    # Remap to 1+1D: surviving axis → "x", killed axes removed.
    if ctx.reduction:
        lines.extend(_wls_json_plane_wave_reduction(ctx))

    # Export
    escaped_output = str(ctx.output_path).replace("\\", "\\\\").replace('"', '\\"')
    lines.extend(
        (
            f'outputPath = "{escaped_output}";',
            "outputDir = DirectoryName[outputPath];",
            'If[outputDir =!= "" && !DirectoryQ[outputDir], CreateDirectory[outputDir]];',
            "",
            'Print["JSON Output:"];',
            # When reduction is active, jsonStringFinal holds the remapped
            # JSON string (never re-imported to avoid Association issues).
            "Print[jsonStringFinal];"
            if ctx.reduction
            else 'Print[ExportString[jsonStructure, "JSON"]];',
            "",
            # Write directly from string when reduction was applied
            "WriteString[outputPath, jsonStringFinal]; Close[outputPath];"
            if ctx.reduction
            else 'Export[outputPath, jsonStructure, "JSON"];',
            'Print[""];',
            'Print["Exported to: ", outputPath];',
            "",
            f'Print["*** {ctx.theory_name} derivation complete! ***"];',
        )
    )

    return lines


# --- WLS assembly & execution ---


def generate_wls(
    config: dict[str, Any],
    output_override: str | None = None,
    *,
    config_dir: Path | None = None,
) -> str:
    """Generate a complete .wls script from a TOML config.

    Parameters
    ----------
    config : dict
        Parsed TOML configuration.
    output_override : str | None
        Override output JSON path.
    config_dir : Path | None
        Directory of the TOML config file.  Relative output paths are resolved
        against this directory (falls back to CWD if *None*).

    Returns
    -------
    str
        Complete Wolfram Language script.

    """
    _validate_config(config)

    prefix = _make_prefix(config)
    dim = config["spacetime"]["dimension"]

    # Resolve pipeline path to absolute so the WLS script works from any location
    wolfram_dir = Path(__file__).resolve().parent.parent / "wolfram"

    # Resolve output path to absolute — relative paths resolve against config_dir
    raw_output = output_override or config.get("output", {}).get("path", "output.json")
    resolved_output = Path(raw_output)
    if not resolved_output.is_absolute():
        base = config_dir if config_dir is not None else Path.cwd()
        resolved_output = (base / resolved_output).resolve()

    ctx = _WlsContext(
        prefix=prefix,
        dim=dim,
        fields=config["fields"],
        constants=config.get("constants", {}).get("names", []),
        coords=_COORDS[dim],
        manifold=f"{prefix}M{dim}",
        metric=f"{prefix}Eta",
        cd=f"{prefix}CD",
        chart=f"{prefix}Cart",
        torsion=config.get("torsion"),  # [torsion] section (None if absent)
        cdt=f"{prefix}CDT",
        theory_name=config.get("theory", {}).get("name", "Custom Theory"),
        output_path=str(resolved_output),
        lagrangian_expr=config.get("lagrangian", {}).get("expression", "").strip(),
        is_multi=len(config["fields"]) > 1,
        pipeline_path=str(wolfram_dir),
        parameters={k: float(v) for k, v in config.get("parameters", {}).items()},
        derived_fields=config.get("derived_fields", []),
        background_fields=config.get("background_fields", []),
        linearization=config.get("linearization"),
        constraint_solver=config.get("constraint_solver"),
        gauge=config.get("gauge", []),
        reduction=config.get("reduction"),
        metric_diagonal=_apply_coord_values(
            [str(e) for e in config["spacetime"].get("diagonal", [])],
            config.get("reduction", {}).get("coordinate_values", {}),
        ),
        metric_type=config["spacetime"].get("metric", "minkowski"),
    )

    is_linearization = ctx.linearization is not None
    has_gauge = bool(ctx.gauge)
    has_type_a = has_gauge and any(
        _resolve_gauge_mechanism(g) == "lagrangian_term" for g in ctx.gauge
    )
    has_type_b = has_gauge and any(
        _resolve_gauge_mechanism(g) == "constraint" for g in ctx.gauge
    )

    lines: list[str] = []
    lines.extend(_wls_header(ctx))
    lines.extend(
        _wls_packages(
            ctx.pipeline_path,
            load_xpert=is_linearization,
            load_gauge=has_type_a,  # Only load GaugeFix.wl for Type A
        )
    )
    lines.extend(_wls_spacetime(config, ctx))
    lines.extend(_wls_fields(ctx, include_bg=_needs_bg_tensor(config)))
    lines.extend(_wls_derived_fields(ctx))

    if is_linearization and ctx.lagrangian_expr:
        # Lagrangian-first linearization: single-path via L^(2)
        # L -> L^(2) = Perturbation[L,2]/2 -> VarD[H,CD][L^(2)] -> EOM
        # Same L^(2) feeds canonical pipeline (pi, H)
        lines.extend(_wls_lagrangian(ctx))
        lines.extend(
            _wls_linearize_from_lagrangian(
                ctx,
                include_bg=_needs_bg_tensor(config),
            )
        )
        # EOM computed inside _wls_linearize_from_lagrangian
        # Reduction applied inside _wls_linearize_from_lagrangian (after L^(2))
    elif is_linearization:
        # Legacy: direct EOM linearization (deprecated — no [lagrangian])
        lines.extend(_wls_linearization(ctx, include_bg=_needs_bg_tensor(config)))
    else:
        lines.extend(_wls_lagrangian(ctx))
        # Plane-wave reduction: zero transverse derivatives in Lagrangian
        # BEFORE gauge fixing and EL derivation for efficiency
        if ctx.reduction is not None:
            lines.extend(_wls_plane_wave_reduction_lagrangian(ctx))
        # Type A gauge fixing: modify Lagrangian before EL derivation
        if has_type_a:
            lines.extend(_wls_gauge_fixing_type_a(ctx))
        lines.extend(_wls_pre_decomposition_tt_zeroing(ctx))
        if ctx.is_multi:
            lines.extend(_wls_euler_lagrange_multi(ctx))
        else:
            lines.extend(_wls_euler_lagrange_single(ctx))
        # Type B gauge fixing: constraints applied after decomposition
        if has_type_b:
            lines.extend(_wls_gauge_fixing_type_b(ctx))

    # Plane-wave reduction: zero transverse derivatives in fieldEquations.
    # For the linearization+Lagrangian path, fieldEquations is empty here
    # (Component E-L fills it later inside _wls_metadata_and_export).
    # Defer plane-wave reduction to the post-Component-E-L block.
    # For the non-linearization path, fieldEquations is already populated
    # by VarD + DecomposeToComponents, so run it here.
    if ctx.reduction is not None and not (is_linearization and ctx.lagrangian_expr):
        lines.extend(_wls_plane_wave_reduction_equations(ctx))
        lines.extend(_wls_plane_wave_field_elimination(ctx))
        lines.extend(_wls_plane_wave_coordinate_evaluation(ctx))

    lines.extend(_wls_metadata_and_export(config, ctx))

    script = "\n".join(lines) + "\n"

    # Pre-execution bracket validation (#191): catch codegen errors before
    # the 30-min wolframscript execution.  Unbalanced brackets from f-string
    # quoting bugs are the most common codegen failure mode.
    bracket_errors = validate_wls_brackets(script)
    if bracket_errors:
        import logging

        logger = logging.getLogger(__name__)
        logger.warning(
            "Generated .wls script has %d bracket error(s):\n  %s",
            len(bracket_errors),
            "\n  ".join(bracket_errors[:10]),
        )

    return script


# --- Execution ---


def _run_wolframscript(script_path: Path, *, timeout: int = 0) -> int:
    """Run wolframscript on a .wls file.

    Parameters
    ----------
    script_path
        Path to the .wls script.
    timeout
        Maximum seconds to wait. 0 means no timeout.

    Returns
    -------
    int
        Exit code from wolframscript.
    """
    if shutil.which("wolframscript") is None:
        from tidal.cli._console import error_with_hint as _cerror_hint

        _cerror_hint(
            "'wolframscript' not found on PATH.",
            [
                "Install Wolfram Engine (free for development): https://www.wolfram.com/engine/",
                "Or use --dry-run to see the generated script without execution.",
            ],
        )
        return 1

    # Kill any stale WolframKernel processes to avoid license seat exhaustion.
    # The Wolfram Engine license limits concurrent sessions; orphaned kernels
    # from previous crashed/interrupted runs consume seats and cause
    # "license error" segfaults on subsequent invocations.
    subprocess.run(["pkill", "-f", "WolframKernel"], capture_output=True, check=False)

    # Warn if no swap — large derivations may crash without sufficient memory
    try:
        with Path("/proc/swaps").open(encoding="ascii") as f:
            if len(f.readlines()) <= 1:  # header only
                from tidal.cli._console import warn as _cwarn

                _cwarn(
                    "No swap space available. Large derivations may crash. "
                    "Docker Desktop: Settings → Resources → Memory → 16 GB+"
                )
    except OSError:
        pass

    timeout_effective = timeout if timeout > 0 else None
    timeout_str = f" (timeout: {timeout}s)" if timeout_effective else ""
    print(f"Running: wolframscript -file {script_path}{timeout_str}")
    print()

    try:
        result = subprocess.run(
            ["wolframscript", "-file", str(script_path)],
            capture_output=False,
            check=False,
            timeout=timeout_effective,
        )
    except subprocess.TimeoutExpired:
        print(
            f"\nError: wolframscript timed out after {timeout}s.",
            file=sys.stderr,
        )
        print(
            "Investigate which pipeline stage is the bottleneck and optimize "
            "the algorithm. Use --timeout 0 to disable if truly needed.",
            file=sys.stderr,
        )
        # Clean up orphaned WolframKernel after timeout
        subprocess.run(
            ["pkill", "-f", "WolframKernel"], capture_output=True, check=False
        )
        return 1

    return result.returncode


def _derive_from_toml(config_path: Path, args: Namespace) -> int:  # noqa: C901, PLR0912, PLR0915
    """Run derivation from a TOML config file.

    Parameters
    ----------
    config_path : Path
        Path to the TOML config file.
    args : Namespace
        Parsed CLI arguments.

    Returns
    -------
    int
        Exit code.
    """
    with config_path.open("rb") as f:
        config = tomllib.load(f)

    script_content = generate_wls(
        config, output_override=args.output, config_dir=config_path.parent.resolve()
    )

    if args.dry_run:
        print(script_content)
        return 0

    # --- Derivation caching ---
    # Hash the generated WLS script (which captures all TOML config, Wolfram
    # pipeline code paths, and field/parameter definitions).  If the output
    # JSON already exists with a matching hash, skip the expensive Wolfram run.
    script_hash = hashlib.sha256(script_content.encode("utf-8")).hexdigest()
    raw_output = args.output or config.get("output", {}).get("path", "output.json")
    resolved_out = Path(raw_output)
    if not resolved_out.is_absolute():
        resolved_out = (config_path.parent.resolve() / resolved_out).resolve()

    force = getattr(args, "force_derive", False)
    if not force and resolved_out.exists():
        try:
            existing = _json_mod.loads(resolved_out.read_text(encoding="utf-8"))
            existing_hash = existing.get("metadata", {}).get("derivation_hash", "")
            if existing_hash == script_hash:
                print(f"Derivation cache hit: {resolved_out.name}")
                print(
                    "Generated script unchanged — skipping wolframscript. "
                    "Use --force-derive to re-run."
                )
                return 0
        except Exception:  # noqa: BLE001, S110
            pass  # Corrupted JSON or missing fields — re-derive

    if args.save_script:
        save_path = Path(args.save_script)
        save_path.write_text(script_content, encoding="utf-8")
        print(f"Saved script to: {save_path.resolve()}")
        if shutil.which("wolframscript") is None:
            print(
                "Note: wolframscript not found. Run the script manually when available."
            )
            return 0
        return _run_wolframscript(save_path, timeout=args.timeout)

    # Use temp file
    with tempfile.NamedTemporaryFile(
        encoding="utf-8", mode="w", suffix=".wls", delete=False, prefix="tidal_derive_"
    ) as tmp:
        tmp.write(script_content)
        tmp_path = Path(tmp.name)

    try:
        ret = _run_wolframscript(tmp_path, timeout=args.timeout)
    finally:
        tmp_path.unlink(missing_ok=True)

    # Use resolved output path from cache check above
    resolved = resolved_out

    # Wolfram Engine in Docker sometimes exits with non-zero code due to
    # license cleanup errors AFTER successfully writing the JSON.  If the
    # output file exists with valid structure, treat as success for
    # post-processing (reduction, validation, hash injection).
    if ret != 0 and resolved.exists():
        try:
            probe = _json_mod.loads(resolved.read_text(encoding="utf-8"))
            if probe.get("equations") and len(probe["equations"]) > 0:
                print(
                    f"\nNote: wolframscript exited with code {ret} but "
                    f"JSON was exported successfully — proceeding with "
                    f"post-processing.",
                    file=sys.stderr,
                )
                ret = 0
        except Exception:  # noqa: BLE001, S110
            pass  # JSON missing or corrupt — honour the non-zero exit code

    # NOTE: Plane-wave reduction (coordinate remapping, operator renaming,
    # dimension change) is now handled entirely in Wolfram via
    # _wls_json_plane_wave_reduction(). No Python post-processing needed.

    # Post-validate output JSON if wolframscript succeeded
    if ret == 0 and resolved.exists():
        try:
            from tidal.symbolic.json_loader import (
                load_equation_system,
            )

            spec = load_equation_system(resolved)
            print()
            print(
                f"Validation: JSON loaded successfully ({spec.n_components} components)"
            )
        except Exception as exc:  # noqa: BLE001
            from tidal.cli._console import warn as _cwarn

            _cwarn(f"JSON validation failed: {exc}")
            ret = 1

    # Inject derivation hash into JSON metadata for future cache checks
    if ret == 0 and resolved.exists():
        try:
            spec_data = _json_mod.loads(resolved.read_text(encoding="utf-8"))
            spec_data.setdefault("metadata", {})["derivation_hash"] = script_hash
            resolved.write_text(
                _json_mod.dumps(spec_data, indent="\t"), encoding="utf-8"
            )
        except Exception:  # noqa: BLE001, S110
            pass  # Non-critical — derivation succeeded, hash injection is optional

    return ret


def derive_command(args: Namespace) -> int:
    """Execute the derive command.

    Parameters
    ----------
    args : Namespace
        Parsed CLI arguments.

    Returns
    -------
    int
        Exit code.
    """
    from tidal.cli._console import error as _cerror

    config_path = Path(args.config)
    if not config_path.exists():
        _cerror(f"file not found: {config_path}")
        return 1

    ext = config_path.suffix.lower()

    if ext == ".wls":
        return _run_wolframscript(config_path, timeout=args.timeout)

    if ext in {".toml", ".tml"}:
        return _derive_from_toml(config_path, args)

    _cerror(
        f"unsupported file extension '{ext}'. Use .toml for config or .wls for script."
    )
    return 1
