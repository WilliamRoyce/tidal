"""``tidal measure`` — Extract physics measurements from simulation output.

Loads an NPZ file produced by ``tidal simulate --output`` and runs
measurement analyses: energy decomposition, conversion probability,
mixing length, spectral analysis, and energy conservation diagnostics.

The JSON equation spec can be auto-discovered from the NPZ metadata
(stored by ``tidal simulate``) or provided explicitly via ``--spec``.

Physics notes
-------------
- **Energy**: Canonical Hamiltonian energy from the Lagrangian structure.
  Decomposed into per-field (kinetic + gradient + mass) and interaction terms.
- **Conservation**: Relative drift ``|dE/E_0|`` over the simulation window.
  Systems with no explicit time-dependence should conserve energy.
- **Conversion**: ``P(t) = E_target(t) / E_source(0)`` — fraction of initial
  source energy transferred to the target field(s).
- **Mixing length**: ``L_mix = pi / omega_dom`` — half-period of the dominant
  oscillation frequency in ``P(t)``.  Uncertainty from HWHM of the spectral
  peak: ``dL = (pi / omega^2) * FWHM/2``.
- **Spectrum**: Spatial power spectrum ``|hat{phi}(k)|^2`` at initial and
  final snapshots.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

import numpy as np

if TYPE_CHECKING:
    from argparse import Namespace

    from tidal.measurement._conversion import ConversionResult
    from tidal.measurement._io import SimulationData

_VALID_MEASUREMENTS = frozenset({
    "summary",
    "energy",
    "conversion",
    "mixing",
    "spectrum",
    "spectral_conversion",
    "dispersion",
    "conservation",
})


# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------


def _resolve_spec_path(data_path: Path, spec_arg: str | None) -> Path:  # noqa: C901
    """Resolve the JSON spec path from CLI flag, directory metadata, or NPZ.

    Resolution order:
    1. Explicit ``--spec`` flag (highest priority)
    2. ``spec_path`` from ``metadata.json`` (snapshot directory)
    3. ``_spec_path`` stored in the NPZ by ``tidal simulate``
    4. Error — no spec found

    Raises
    ------
    FileNotFoundError
        If the resolved path does not exist.
    ValueError
        If no spec path can be determined.
    """
    if spec_arg is not None:
        p = Path(spec_arg)
        if not p.exists():
            msg = f"Spec file not found: {p}"
            raise FileNotFoundError(msg)
        return p

    # Try snapshot directory metadata
    if data_path.is_dir():
        metadata_file = data_path / "metadata.json"
        if metadata_file.exists():
            import json

            metadata = json.loads(metadata_file.read_text())
            if "spec_path" in metadata:
                p = Path(str(metadata["spec_path"]))
                if p.exists():
                    return p
                # Try relative to directory
                relative = data_path.parent / p.name
                if relative.exists():
                    return relative

        msg = (
            f"Cannot determine JSON spec for {data_path.name} — "
            f"use --spec to provide the path explicitly"
        )
        raise ValueError(msg)

    # Try NPZ metadata
    data = np.load(str(data_path), allow_pickle=False)
    if "_spec_path" in data:
        stored = str(data["_spec_path"])
        p = Path(stored)
        if p.exists():
            return p
        # Stored path doesn't exist — try relative to NPZ directory
        relative = data_path.parent / p.name
        if relative.exists():
            return relative

    msg = (
        f"Cannot determine JSON spec for {data_path.name} — "
        f"use --spec to provide the path explicitly"
    )
    raise ValueError(msg)


def _parse_measurements(raw: str | None) -> set[str]:
    """Parse ``--what=energy,mixing`` into a set of measurement names.

    Returns ``{"summary"}`` when *raw* is None (default behavior).

    Raises
    ------
    ValueError
        If any requested measurement name is unknown.
    """
    if raw is None:
        return {"summary"}

    names = {s.strip() for s in raw.split(",")}
    unknown = names - _VALID_MEASUREMENTS
    if unknown:
        msg = (
            f"Unknown measurement(s): {', '.join(sorted(unknown))}. "
            f"Valid: {', '.join(sorted(_VALID_MEASUREMENTS))}"
        )
        raise ValueError(msg)
    return names


def _parse_field_list(raw: str | None) -> tuple[str, ...] | None:
    """Parse comma-separated field names into a tuple, or None if absent."""
    if raw is None:
        return None
    return tuple(s.strip() for s in raw.split(",") if s.strip())


def _load_data(
    data_path: Path,
    spec_path: Path,
    param_overrides: list[str],
) -> SimulationData:
    """Load simulation data from NPZ or snapshot directory.

    Auto-detects format: directories use memory-mapped lazy loading (O(1)
    RAM), NPZ files are loaded eagerly.

    Merges ``--param`` overrides with parameters stored in the data.
    """
    from tidal.cli._simulate import _parse_params  # pyright: ignore[reportPrivateUsage]
    from tidal.measurement import SimulationData
    from tidal.symbolic import load_equation_system

    spec = load_equation_system(spec_path)
    data = SimulationData.load(data_path, spec)

    # Merge CLI param overrides
    if param_overrides:
        merged = dict(data.parameters)
        cli_params = _parse_params(param_overrides, spec)
        merged.update(cli_params)
        from dataclasses import replace

        data = replace(data, parameters=merged)

    return data


# ------------------------------------------------------------------
# Individual measurement runners
# ------------------------------------------------------------------


def _run_energy(data: SimulationData) -> dict[str, Any]:
    """Compute per-field and total energy timeseries."""
    from tidal.measurement import compute_energy_timeseries

    times, per_field, interaction, total = compute_energy_timeseries(data)
    return {
        "times": times.tolist(),
        "per_field": {k: v.tolist() for k, v in per_field.items()},
        "interaction": interaction.tolist(),
        "total": total.tolist(),
    }


def _run_conservation(data: SimulationData, threshold: float) -> dict[str, Any]:
    """Check energy conservation."""
    from tidal.measurement import check_energy_conservation

    diag = check_energy_conservation(data, threshold=threshold)
    return {
        "max_relative_error": diag.max_relative_error,
        "is_conserved": diag.is_conserved,
        "threshold": threshold,
    }


def _run_conversion(
    data: SimulationData,
    source: tuple[str, ...] | None,
    target: tuple[str, ...] | None,
) -> dict[str, Any]:
    """Compute conversion probability.

    Auto-detects source/target when not specified: first dynamical field
    as source, all remaining dynamical fields as target.

    The internal ``_result_obj`` key stores the ConversionResult for
    downstream mixing computation; it is stripped before serialization.

    Raises
    ------
    ValueError
        If no dynamical fields or auto-detection fails.
    """
    from tidal.measurement import (
        compute_conversion_probability,
        compute_group_conversion,
    )

    dyn = data.dynamical_fields
    if not dyn:
        msg = "No dynamical fields found — conversion requires at least 2 fields"
        raise ValueError(msg)

    # Auto-detect source/target
    if source is None:
        source = (dyn[0],)
    if target is None:
        remaining = tuple(f for f in dyn if f not in source)
        if not remaining:
            msg = (
                f"Cannot auto-detect target: all dynamical fields ({', '.join(dyn)}) "
                f"are in the source set. Use --target explicitly."
            )
            raise ValueError(msg)
        target = remaining

    # Single-field or group conversion
    if len(source) == 1 and len(target) == 1:
        result = compute_conversion_probability(data, source[0], target[0])
    else:
        result = compute_group_conversion(data, list(source), list(target))

    peak_idx = int(np.argmax(result.probability))
    return {
        "source": list(source),
        "target": list(target),
        "peak_probability": float(result.probability[peak_idx]),
        "peak_time": float(result.times[peak_idx]),
        "_result_obj": result,
    }


def _run_mixing(conversion_result: ConversionResult) -> dict[str, Any]:
    """Extract mixing length from conversion probability timeseries."""
    from tidal.measurement import compute_mixing_length, compute_mixing_spectrum

    mixing = compute_mixing_length(conversion_result)
    spectrum = compute_mixing_spectrum(conversion_result)

    return {
        "mixing_length": mixing.mixing_length,
        "mixing_length_uncertainty": mixing.mixing_length_uncertainty,
        "dominant_frequency": mixing.dominant_frequency,
        "frequency_fwhm": mixing.frequency_fwhm,
        "max_conversion": mixing.max_conversion,
        "n_peaks": len(mixing.peaks),
        "rayleigh_resolution": spectrum.rayleigh_resolution,
        "_spectrum_obj": spectrum,
    }


def _run_spectrum(data: SimulationData) -> dict[str, Any]:
    """Compute spatial power spectrum at initial and final snapshots."""
    from tidal.measurement import compute_spectrum

    result: dict[str, Any] = {}
    for name in data.fields:
        snap_initial = compute_spectrum(
            data.fields[name][0], data.grid_spacing, data.periodic
        )
        snap_final = compute_spectrum(
            data.fields[name][-1], data.grid_spacing, data.periodic
        )
        result[name] = {
            "initial": {
                "wavenumbers": snap_initial.wavenumbers.tolist(),
                "power": snap_initial.power_spectrum.tolist(),
            },
            "final": {
                "wavenumbers": snap_final.wavenumbers.tolist(),
                "power": snap_final.power_spectrum.tolist(),
            },
        }
    return result


def _run_spectral_conversion(
    data: SimulationData,
    source: tuple[str, ...] | None,
    target: tuple[str, ...] | None,
) -> dict[str, Any]:
    """Compute per-mode spectral conversion P(k,t).

    Raises
    ------
    ValueError
        If no dynamical fields, or auto-detection fails.
    """
    from tidal.measurement import (
        compute_group_spectral_conversion,
        compute_spectral_conversion,
    )

    dyn = data.dynamical_fields
    if not dyn:
        msg = "No dynamical fields found — spectral conversion requires at least 2 fields"
        raise ValueError(msg)
    if source is None:
        source = (dyn[0],)
    if target is None:
        remaining = tuple(f for f in dyn if f not in source)
        if not remaining:
            msg = (
                f"Cannot auto-detect target: all dynamical fields ({', '.join(dyn)}) "
                f"are in the source set. Use --target explicitly."
            )
            raise ValueError(msg)
        target = remaining

    if len(source) == 1 and len(target) == 1:
        result = compute_spectral_conversion(data, source[0], target[0])
    else:
        result = compute_group_spectral_conversion(data, list(source), list(target))

    n_active = int(result.active_modes.sum())
    peak_k_idx = int(np.argmax(result.probability[-1])) if n_active > 0 else 0

    return {
        "source": list(source),
        "target": list(target),
        "n_modes": len(result.wavenumbers),
        "n_active_modes": n_active,
        "peak_mode_k": float(result.wavenumbers[peak_k_idx]) if n_active > 0 else None,
        "peak_mode_P": float(result.probability[-1, peak_k_idx]) if n_active > 0 else None,
        "_result_obj": result,
    }


def _run_dispersion(
    data: SimulationData,
    field_name: str,
) -> dict[str, Any]:
    """Compute dispersion relation omega(k) for a single field."""
    from tidal.measurement import compute_dispersion

    result = compute_dispersion(data, field_name)
    n_active = int(np.count_nonzero(result.peak_frequencies > 0.0))

    return {
        "field": field_name,
        "n_modes": len(result.wavenumbers),
        "n_active_modes": n_active,
        "rayleigh_resolution": result.rayleigh_resolution,
        "_result_obj": result,
    }


def _run_summary(
    data: SimulationData,
    threshold: float,
) -> dict[str, Any]:
    """Run full summary: energy + conservation + auto-detect conversion + mixing.

    Auto-detects conversion fields from the dynamical field list (first
    field as source, rest as target).  Mixing is attempted but failures
    are reported as errors, not crashes.
    """
    results: dict[str, Any] = {}

    results["energy"] = _run_energy(data)
    results["conservation"] = _run_conservation(data, threshold)

    # Auto-detect conversion (needs >= 2 dynamical fields)
    dyn = data.dynamical_fields
    if len(dyn) >= 2:  # noqa: PLR2004
        try:
            results["conversion"] = _run_conversion(data, None, None)
            # Chain: mixing from conversion
            conv_result: ConversionResult = results["conversion"]["_result_obj"]
            try:
                results["mixing"] = _run_mixing(conv_result)
            except ValueError as e:
                results["mixing"] = {"error": str(e)}
        except ValueError as e:
            results["conversion"] = {"error": str(e)}

    return results


# ------------------------------------------------------------------
# Output formatting
# ------------------------------------------------------------------


def _strip_internal_keys(obj: dict[str, Any]) -> dict[str, Any]:
    """Remove keys starting with ``_`` from a dict, recursively."""
    result: dict[str, Any] = {}
    for k, v in obj.items():
        if k.startswith("_"):
            continue
        if isinstance(v, dict):
            result[k] = _strip_internal_keys(cast("dict[str, Any]", v))
        else:
            result[k] = v
    return result


def _format_json(results: dict[str, Any], data: SimulationData) -> str:
    """Produce JSON output with simulation metadata."""
    output: dict[str, Any] = {
        "simulation": {
            "time_range": [float(data.times[0]), float(data.times[-1])],
            "n_snapshots": data.n_snapshots,
            "fields": list(data.fields.keys()),
            "parameters": data.parameters,
        },
    }
    output.update(_strip_internal_keys(results))
    return json.dumps(output, indent=2, default=float)


def _format_text_section_conservation(lines: list[str], cons: dict[str, Any]) -> None:
    """Append conservation section to *lines*."""
    if "error" in cons:
        lines.append(f"Energy Conservation: ERROR ({cons['error']})")
    else:
        status = "PASS" if cons["is_conserved"] else "FAIL"
        lines.extend([
            f"Energy Conservation: {status}",
            f"  max |dE/E| = {cons['max_relative_error']:.2e}",
            f"  threshold  = {cons['threshold']:.0e}",
        ])
    lines.append("")


def _format_text_section_energy(lines: list[str], eng: dict[str, Any]) -> None:
    """Append energy section to *lines*."""
    if "error" in eng:
        lines.append(f"Energy: ERROR ({eng['error']})")
    else:
        lines.append("Per-Field Energy (final):")
        for name, series in eng["per_field"].items():
            lines.append(f"  {name}: {series[-1]:.6f}")
        lines.extend([
            f"  interaction: {eng['interaction'][-1]:.6f}",
            f"  total:       {eng['total'][-1]:.6f}",
        ])
    lines.append("")


def _format_text_section_conversion(lines: list[str], conv: dict[str, Any]) -> None:
    """Append conversion section to *lines*."""
    if "error" in conv:
        lines.append(f"Conversion: ERROR ({conv['error']})")
    else:
        src = ", ".join(conv["source"])
        tgt = ", ".join(conv["target"])
        lines.extend([
            f"Conversion ({src} -> {tgt}):",
            f"  Peak P(t) = {conv['peak_probability']:.6f}",
            f"  at t = {conv['peak_time']:.2f}",
        ])
    lines.append("")


def _format_text_section_mixing(lines: list[str], mix: dict[str, Any]) -> None:
    """Append mixing section to *lines*."""
    if "error" in mix:
        lines.append(f"Mixing Length: not extracted ({mix['error']})")
    else:
        lines.extend([
            "Mixing Length:",
            f"  L_mix     = {mix['mixing_length']:.4f}"
            f" +/- {mix['mixing_length_uncertainty']:.4f}",
            f"  omega_dom = {mix['dominant_frequency']:.4f}"
            f"  (FWHM = {mix['frequency_fwhm']:.4f})",
            f"  Rayleigh  = {mix['rayleigh_resolution']:.4f}",
        ])
    lines.append("")


def _format_text_section_spectrum(lines: list[str], spec: dict[str, Any]) -> None:
    """Append spectrum section to *lines*."""
    if "error" in spec:
        lines.append(f"Spectrum: ERROR ({spec['error']})")
    else:
        lines.append("Power Spectrum:")
        for name, snap in spec.items():
            n_bins = len(snap["initial"]["wavenumbers"])
            lines.append(f"  {name}: {n_bins} frequency bins")
    lines.append("")


def _format_text_section_spectral_conversion(
    lines: list[str], sc: dict[str, Any],
) -> None:
    """Append spectral conversion section to *lines*."""
    if "error" in sc:
        lines.append(f"Spectral Conversion: ERROR ({sc['error']})")
    else:
        src = ", ".join(sc["source"])
        tgt = ", ".join(sc["target"])
        lines.extend([
            f"Spectral Conversion ({src} -> {tgt}):",
            f"  Active k-modes: {sc['n_active_modes']} / {sc['n_modes']}",
        ])
        if sc["peak_mode_k"] is not None:
            lines.extend([
                f"  Peak P(k, t_final) at |k| = {sc['peak_mode_k']:.4f}",
                f"    P(k) = {sc['peak_mode_P']:.6f}",
            ])
    lines.append("")


def _format_text_section_dispersion(
    lines: list[str], disp: dict[str, Any],
) -> None:
    """Append dispersion section to *lines*."""
    if "error" in disp:
        lines.append(f"Dispersion: ERROR ({disp['error']})")
    else:
        lines.extend([
            f"Dispersion ({disp['field']}):",
            f"  Active k-modes: {disp['n_active_modes']} / {disp['n_modes']}",
            f"  Rayleigh resolution: {disp['rayleigh_resolution']:.4f} rad/time",
        ])
    lines.append("")


def _format_text(results: dict[str, Any], data: SimulationData) -> str:
    """Produce human-readable aligned text output."""
    lines: list[str] = []
    sep = "=" * 64

    lines.extend([
        sep,
        f"Measurement: {', '.join(data.fields.keys())}",
        sep,
        "",
        "Simulation:",
        f"  Time range: {float(data.times[0]):.1f} -> {float(data.times[-1]):.1f}"
        f"  ({data.n_snapshots} snapshots)",
        f"  Fields: {', '.join(data.fields.keys())}",
    ])
    if data.parameters:
        param_str = ", ".join(f"{k}={v}" for k, v in data.parameters.items())
        lines.append(f"  Parameters: {param_str}")
    lines.append("")

    if "conservation" in results:
        _format_text_section_conservation(lines, results["conservation"])
    if "energy" in results:
        _format_text_section_energy(lines, results["energy"])
    if "conversion" in results:
        _format_text_section_conversion(lines, results["conversion"])
    if "mixing" in results:
        _format_text_section_mixing(lines, results["mixing"])
    if "spectrum" in results:
        _format_text_section_spectrum(lines, results["spectrum"])
    if "spectral_conversion" in results:
        _format_text_section_spectral_conversion(lines, results["spectral_conversion"])
    if "dispersion" in results:
        _format_text_section_dispersion(lines, results["dispersion"])

    lines.append(sep)
    return "\n".join(lines)


# ------------------------------------------------------------------
# Main entry point
# ------------------------------------------------------------------


def _run_individual_measurements(  # noqa: C901, PLR0912, PLR0915
    measurements: set[str],
    data: SimulationData,
    args: Namespace,
    threshold: float,
) -> dict[str, Any] | int:
    """Run individual measurements (non-summary mode).

    Returns a results dict on success, or 1 on early error.
    """
    results: dict[str, Any] = {}

    if "energy" in measurements:
        try:
            results["energy"] = _run_energy(data)
        except ValueError as e:
            results["energy"] = {"error": str(e)}

    if "conservation" in measurements:
        try:
            results["conservation"] = _run_conservation(data, threshold)
        except ValueError as e:
            results["conservation"] = {"error": str(e)}

    if "conversion" in measurements or "mixing" in measurements:
        source = _parse_field_list(getattr(args, "source", None))
        target = _parse_field_list(getattr(args, "target", None))

        if "conversion" in measurements and source is None:
            print(
                "Error: --source required for --what=conversion "
                "(or use --what=summary for auto-detection)",
                file=sys.stderr,
            )
            return 1

        try:
            conv = _run_conversion(data, source, target)
            results["conversion"] = conv

            if "mixing" in measurements:
                conv_result: ConversionResult = conv["_result_obj"]
                try:
                    results["mixing"] = _run_mixing(conv_result)
                except ValueError as e:
                    results["mixing"] = {"error": str(e)}
        except ValueError as e:
            results["conversion"] = {"error": str(e)}
            if "mixing" in measurements:
                results["mixing"] = {"error": f"conversion failed: {e}"}

    if "spectrum" in measurements:
        try:
            results["spectrum"] = _run_spectrum(data)
        except ValueError as e:
            results["spectrum"] = {"error": str(e)}

    if "spectral_conversion" in measurements:
        source = _parse_field_list(getattr(args, "source", None))
        target = _parse_field_list(getattr(args, "target", None))

        if source is None:
            print(
                "Error: --source required for --what=spectral_conversion "
                "(or use --what=summary for auto-detection)",
                file=sys.stderr,
            )
            return 1

        try:
            results["spectral_conversion"] = _run_spectral_conversion(
                data, source, target,
            )
        except ValueError as e:
            results["spectral_conversion"] = {"error": str(e)}

    if "dispersion" in measurements:
        source = _parse_field_list(getattr(args, "source", None))
        # Default to first dynamical field if no --source
        dyn = data.dynamical_fields
        if source is not None:
            field_name = source[0]
        elif dyn:
            field_name = dyn[0]
        else:
            print("Error: no dynamical fields for dispersion", file=sys.stderr)
            return 1

        try:
            results["dispersion"] = _run_dispersion(data, field_name)
        except ValueError as e:
            results["dispersion"] = {"error": str(e)}

    return results


def measure_command(args: Namespace) -> int:
    """Run ``tidal measure`` subcommand.

    Flow:
    1. Validate NPZ exists
    2. Resolve JSON spec path (--spec or NPZ metadata)
    3. Parse --what measurement set
    4. Load SimulationData
    5. Run requested measurements (dependency-ordered)
    6. Format and output (text, JSON, or plot)

    Returns 0 on success, 1 on error.
    """
    data_path = Path(args.npz_path)
    if not data_path.exists():
        print(f"Error: data path not found: {data_path}", file=sys.stderr)
        return 1

    spec_path = _resolve_spec_path(data_path, getattr(args, "spec", None))
    measurements = _parse_measurements(getattr(args, "what", None))
    quiet: bool = getattr(args, "quiet", False)

    if not quiet:
        fmt_label = "directory" if data_path.is_dir() else "NPZ"
        print(f"Loading: {data_path.name} ({fmt_label})")
        print(f"Spec:    {spec_path.name}")

    data = _load_data(data_path, spec_path, getattr(args, "param", None) or [])

    if not quiet:
        print(
            f"  {data.n_snapshots} snapshots, "
            f"{len(data.fields)} fields, "
            f"t=[{float(data.times[0]):.1f}, {float(data.times[-1]):.1f}]"
        )

    threshold: float = getattr(args, "energy_threshold", 1e-3)

    # Run measurements
    if "summary" in measurements:
        results = _run_summary(data, threshold)
    else:
        outcome = _run_individual_measurements(measurements, data, args, threshold)
        if isinstance(outcome, int):
            return outcome
        results = outcome

    # Output
    output_path: str | None = getattr(args, "output", None)
    json_mode: bool = getattr(args, "json_output", False)

    if output_path is not None and output_path.endswith((".png", ".pdf")):
        from tidal.cli._measure_plot import save_measurement_plot

        save_measurement_plot(Path(output_path), data, results)
        if not quiet:
            print(f"  Saved plot to: {output_path}")
    elif json_mode:
        print(_format_json(results, data))
    else:
        print(_format_text(results, data))

    return 0
