"""``tidal measure`` — Extract physics measurements from simulation output.

Loads simulation output from a snapshot directory produced by
``tidal simulate --output`` and runs measurement analyses: energy
decomposition, conversion probability, mixing length, spectral analysis,
and energy conservation diagnostics.

The JSON equation spec can be auto-discovered from ``metadata.json``
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

_VALID_MEASUREMENTS = frozenset(
    {
        "summary",
        "energy",
        "conversion",
        "mixing",
        "spectrum",
        "dispersion",
        "conservation",
        "effective_mass",
        "asymptotic",
        "peak_conversion",
        "velocity",
        "resonance",
    },
)


# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------


def _resolve_spec_path(data_path: Path, spec_arg: str | None) -> Path:
    """Resolve the JSON spec path from CLI flag or directory metadata.

    Resolution order:
    1. Explicit ``--spec`` flag (highest priority)
    2. ``spec_path`` from ``metadata.json`` (snapshot directory)
    3. Error — no spec found

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


def _resolve_source_target(
    data: SimulationData,
    source: tuple[str, ...] | None,
    target: tuple[str, ...] | None,
    *,
    require_both: bool = False,
    measurement_name: str = "measurement",
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Resolve source/target fields with optional auto-detection.

    When *require_both* is ``False`` (default), auto-detects missing fields
    from dynamical field list: source defaults to the first dynamical field,
    target defaults to the remaining dynamical fields.

    When *require_both* is ``True``, raises ``ValueError`` if either is None
    (used for velocity mismatch and resonance where auto-detection is
    ambiguous).

    Raises
    ------
    ValueError
        If auto-detection fails (no dynamical fields, or all fields in source
        set) or if *require_both* and a field is missing.
    """
    if require_both:
        if source is None:
            msg = f"--source required for {measurement_name}"
            raise ValueError(msg)
        if target is None:
            msg = f"--target required for {measurement_name}"
            raise ValueError(msg)
        return source, target

    dyn = data.dynamical_fields
    if not dyn:
        msg = (
            f"No dynamical fields found — {measurement_name} requires at least 2 fields"
        )
        raise ValueError(msg)

    if source is None:
        source = (dyn[0],)
    if target is None:
        remaining = tuple(f for f in dyn if f not in source)
        if not remaining:
            msg = (
                f"Cannot auto-detect target: all dynamical fields "
                f"({', '.join(dyn)}) are in the source set. "
                f"Use --target explicitly."
            )
            raise ValueError(msg)
        target = remaining

    return source, target


def _run_measurement_safe(
    func: Any,  # noqa: ANN401
    *args: Any,  # noqa: ANN401
    **kwargs: Any,  # noqa: ANN401
) -> dict[str, Any]:
    """Run a measurement function with ValueError error handling.

    Returns the function result on success, or ``{"error": str(exc)}``
    if a ``ValueError`` is raised.
    """
    try:
        return func(*args, **kwargs)  # type: ignore[no-any-return]
    except ValueError as e:
        return {"error": str(e)}


def _load_data(
    data_path: Path,
    spec_path: Path,
    param_overrides: list[str],
) -> SimulationData:
    """Load simulation data from a snapshot directory (memory-mapped, O(1) RAM).

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
    """
    from tidal.measurement import (
        compute_conversion_probability,
        compute_group_conversion,
    )

    source, target = _resolve_source_target(
        data,
        source,
        target,
        measurement_name="conversion",
    )

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
            data.fields[name][0],
            data.grid_spacing,
            data.periodic,
        )
        snap_final = compute_spectrum(
            data.fields[name][-1],
            data.grid_spacing,
            data.periodic,
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


def _run_dispersion(
    data: SimulationData,
    field_names: list[str],
) -> dict[str, Any]:
    """Compute dispersion relation omega(k) for a field group.

    Spectral power is summed over all fields in *field_names*, making the
    measurement rotationally covariant within the group.
    """
    from tidal.measurement import compute_dispersion

    result = compute_dispersion(data, field_names)
    n_active = int(np.count_nonzero(result.peak_frequencies > 0.0))

    return {
        "field": result.field_name,
        "n_modes": len(result.wavenumbers),
        "n_active_modes": n_active,
        "rayleigh_resolution": result.rayleigh_resolution,
        "_result_obj": result,
    }


def _run_effective_mass(
    data: SimulationData,
    field_names: list[str],
) -> dict[str, Any]:
    """Compute effective mass from dispersion relation.

    The effective mass is the Lorentz-invariant 4-momentum norm:
    ``m²_eff = ω² - k²`` at the dominant frequency per k-bin.
    """
    from tidal.measurement import compute_effective_mass

    result = compute_effective_mass(data, field_names)
    return {
        "field": result.field_name,
        "m2_eff": result.m2_eff,
        "m2_eff_std": result.m2_eff_std,
        "n_active_modes": result.n_active_modes,
        "_result_obj": result,
    }


def _run_asymptotic(
    data: SimulationData,
    source: tuple[str, ...] | None,
    target: tuple[str, ...] | None,
) -> dict[str, Any]:
    """Compute asymptotic scattering observables.

    Forward/reflected split is defined by the source field's initial
    propagation direction (spectral centroid), making it independent
    of coordinate axis choice.
    """
    from tidal.measurement import compute_asymptotic_conversion

    source, target = _resolve_source_target(
        data,
        source,
        target,
        measurement_name="asymptotic",
    )

    result = compute_asymptotic_conversion(data, list(source), list(target))
    return {
        "source": list(source),
        "target": list(target),
        "P_final": result.P_final,
        "P_transmitted": result.P_transmitted,
        "P_reflected": result.P_reflected,
        "source_wavevector": list(result.source_wavevector),
        "_result_obj": result,
    }


def _run_peak_conversion(
    data: SimulationData,
    source: tuple[str, ...] | None,
    target: tuple[str, ...] | None,
) -> dict[str, Any]:
    """Extract scalar conversion summary: P_max, t_peak, P_final.

    Reuses :func:`_run_conversion` and distills to sweep-friendly scalars.
    """
    conv = _run_conversion(data, source, target)
    result_obj: ConversionResult = conv["_result_obj"]
    peak_idx = int(np.argmax(result_obj.probability))
    return {
        "source": conv["source"],
        "target": conv["target"],
        "P_max": float(result_obj.probability[peak_idx]),
        "P_max_time": float(result_obj.times[peak_idx]),
        "P_final": float(result_obj.probability[-1]),
    }


def _run_velocity(
    data: SimulationData,
    field_names: list[str],
) -> dict[str, Any]:
    """Compute group and phase velocities from dispersion relation."""
    from tidal.measurement import compute_velocities

    result = compute_velocities(data, field_names)
    return {
        "field": result.field_name,
        "n_active_modes": result.n_active_modes,
        "v_group_mean": result.group_velocity_mean,
        "v_phase_mean": result.phase_velocity_mean,
        "_result_obj": result,
    }


def _run_velocity_mismatch(
    data: SimulationData,
    source: tuple[str, ...] | None,
    target: tuple[str, ...] | None,
) -> dict[str, Any]:
    """Compute velocity mismatch between source and target field groups."""
    from tidal.measurement import compute_velocity_mismatch

    source, target = _resolve_source_target(
        data,
        source,
        target,
        require_both=True,
        measurement_name="velocity mismatch",
    )

    result = compute_velocity_mismatch(data, list(source), list(target))
    return {
        "source": list(source),
        "target": list(target),
        "v_group_mean_source": result.source_velocity.group_velocity_mean,
        "v_group_mean_target": result.target_velocity.group_velocity_mean,
        "v_mismatch_max": result.max_mismatch,
        "v_mismatch_mean": result.mean_mismatch,
        "_result_obj": result,
    }


def _run_resonance(
    data: SimulationData,
    source: tuple[str, ...] | None,
    target: tuple[str, ...] | None,
) -> dict[str, Any]:
    """Compute resonance analysis between source and target fields."""
    from tidal.measurement import compute_resonance_analysis

    source, target = _resolve_source_target(
        data,
        source,
        target,
        require_both=True,
        measurement_name="resonance",
    )

    result = compute_resonance_analysis(data, list(source), list(target))
    return {
        "source": result.source_field,
        "target": result.target_field,
        "n_resonant_modes": result.n_resonant_modes,
        "conversion_bandwidth": result.conversion_bandwidth,
        "peak_conversion_k": result.peak_conversion_k,
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
    from tidal.cli._console import key_value, pass_fail

    if "error" in cons:
        lines.append(f"Energy Conservation: ERROR ({cons['error']})")
    else:
        lines.extend(
            [
                pass_fail("Energy Conservation", passed=cons["is_conserved"]),
                key_value("max |dE/E|", f"{cons['max_relative_error']:.2e}"),
                key_value("threshold", f"{cons['threshold']:.0e}"),
            ],
        )
    lines.append("")


def _format_text_section_energy(lines: list[str], eng: dict[str, Any]) -> None:
    """Append energy section to *lines*."""
    from tidal.cli._console import key_value

    if "error" in eng:
        lines.append(f"Energy: ERROR ({eng['error']})")
    else:
        lines.append("Per-Field Energy (final):")
        for name, series in eng["per_field"].items():
            lines.append(key_value(name, f"{series[-1]:.6f}"))
        lines.extend(
            (
                key_value("interaction", f"{eng['interaction'][-1]:.6f}"),
                key_value("total", f"{eng['total'][-1]:.6f}"),
            ),
        )
    lines.append("")


def _format_text_section_conversion(lines: list[str], conv: dict[str, Any]) -> None:
    """Append conversion section to *lines*."""
    from tidal.cli._console import key_value

    if "error" in conv:
        lines.append(f"Conversion: ERROR ({conv['error']})")
    else:
        src = ", ".join(conv["source"])
        tgt = ", ".join(conv["target"])
        lines.extend(
            (
                f"Conversion ({src} -> {tgt}):",
                key_value("Peak P(t)", f"{conv['peak_probability']:.6f}"),
                key_value("at t", f"{conv['peak_time']:.2f}"),
            ),
        )
    lines.append("")


def _format_text_section_mixing(lines: list[str], mix: dict[str, Any]) -> None:
    """Append mixing section to *lines*."""
    if "error" in mix:
        lines.append(f"Mixing Length: not extracted ({mix['error']})")
    else:
        lines.extend(
            [
                "Mixing Length:",
                f"  L_mix     = {mix['mixing_length']:.4f}"
                f" +/- {mix['mixing_length_uncertainty']:.4f}",
                f"  omega_dom = {mix['dominant_frequency']:.4f}"
                f"  (FWHM = {mix['frequency_fwhm']:.4f})",
                f"  Rayleigh  = {mix['rayleigh_resolution']:.4f}",
            ],
        )
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


def _format_text_section_dispersion(
    lines: list[str],
    disp: dict[str, Any],
) -> None:
    """Append dispersion section to *lines*."""
    if "error" in disp:
        lines.append(f"Dispersion: ERROR ({disp['error']})")
    else:
        lines.extend(
            [
                f"Dispersion ({disp['field']}):",
                f"  Active k-modes: {disp['n_active_modes']} / {disp['n_modes']}",
                f"  Rayleigh resolution: {disp['rayleigh_resolution']:.4f} rad/time",
            ],
        )
    lines.append("")


def _format_text_section_effective_mass(
    lines: list[str],
    em: dict[str, Any],
) -> None:
    """Append effective mass section to *lines*."""
    if "error" in em:
        lines.append(f"Effective Mass: ERROR ({em['error']})")
    else:
        lines.extend(
            [
                f"Effective Mass ({em['field']}):",
                f"  m²_eff  = {em['m2_eff']:.6f} +/- {em['m2_eff_std']:.6f}",
                f"  Active modes: {em['n_active_modes']}",
            ],
        )
    lines.append("")


def _format_text_section_asymptotic(
    lines: list[str],
    asym: dict[str, Any],
) -> None:
    """Append asymptotic conversion section to *lines*."""
    if "error" in asym:
        lines.append(f"Asymptotic Conversion: ERROR ({asym['error']})")
    else:
        src = ", ".join(asym["source"])
        tgt = ", ".join(asym["target"])
        lines.extend(
            [
                f"Asymptotic Conversion ({src} -> {tgt}):",
                f"  P_final     = {asym['P_final']:.6f}",
                f"  P_transmitted = {asym['P_transmitted']:.6f}",
                f"  P_reflected = {asym['P_reflected']:.6f}",
                f"  Source k    = ({', '.join(f'{k:.3f}' for k in asym['source_wavevector'])})",
            ],
        )
    lines.append("")


def _format_text_section_peak_conversion(
    lines: list[str],
    pc: dict[str, Any],
) -> None:
    """Append peak conversion section to *lines*."""
    if "error" in pc:
        lines.append(f"Peak Conversion: ERROR ({pc['error']})")
    else:
        src = ", ".join(pc["source"])
        tgt = ", ".join(pc["target"])
        lines.extend(
            [
                f"Peak Conversion ({src} -> {tgt}):",
                f"  P_max   = {pc['P_max']:.6f}  at t = {pc['P_max_time']:.2f}",
                f"  P_final = {pc['P_final']:.6f}",
            ],
        )
    lines.append("")


def _format_text(results: dict[str, Any], data: SimulationData) -> str:  # noqa: C901
    """Produce human-readable aligned text output."""
    from tidal.cli._console import header as _header

    lines: list[str] = []

    lines.append(_header(f"Measurement: {', '.join(data.fields.keys())}"))
    lines.extend(
        [
            "",
            "Simulation:",
            f"  Time range: {float(data.times[0]):.1f} -> {float(data.times[-1]):.1f}"
            f"  ({data.n_snapshots} snapshots)",
            f"  Fields: {', '.join(data.fields.keys())}",
        ],
    )
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
    if "dispersion" in results:
        _format_text_section_dispersion(lines, results["dispersion"])
    if "effective_mass" in results:
        _format_text_section_effective_mass(lines, results["effective_mass"])
    if "asymptotic" in results:
        _format_text_section_asymptotic(lines, results["asymptotic"])
    if "peak_conversion" in results:
        _format_text_section_peak_conversion(lines, results["peak_conversion"])

    lines.append("=" * 64)
    return "\n".join(lines)


# ------------------------------------------------------------------
# Main entry point
# ------------------------------------------------------------------


def _filter_to_dynamical(
    source: tuple[str, ...] | None,
    data: SimulationData,
    measurement_name: str,
) -> list[str]:
    """Filter source fields to dynamical-only, falling back to all dynamical."""
    dyn = list(data.dynamical_fields)
    if source is None:
        return dyn
    filtered = [f for f in source if f in set(dyn)]
    if filtered:
        return filtered
    if dyn:
        print(
            f"Note: {measurement_name}: no dynamical fields in --source list; "
            f"using all dynamical fields: {dyn}",
            file=sys.stderr,
        )
    return dyn


def _run_individual_measurements(  # noqa: C901, PLR0912
    measurements: set[str],
    data: SimulationData,
    args: Namespace,
    threshold: float,
) -> dict[str, Any] | int:
    """Run individual measurements (non-summary mode).

    Returns a results dict on success, or 1 on early error.
    """
    results: dict[str, Any] = {}

    # Parse source/target once for all measurements that use them.
    source = _parse_field_list(getattr(args, "source", None))
    target = _parse_field_list(getattr(args, "target", None))

    # Measurements that require --source (error if missing).
    require_source = {
        "conversion",
        "asymptotic",
        "peak_conversion",
        "resonance",
    }
    needs_source = measurements & require_source
    if needs_source and source is None:
        from tidal.cli._console import error_with_hint

        names = ", ".join(sorted(needs_source))
        error_with_hint(
            f"--source required for --what={names} "
            f"(or use --what=summary for auto-detection)",
            ["Example: `--what conversion --source phi --target psi`"],
        )
        return 1

    if "energy" in measurements:
        results["energy"] = _run_measurement_safe(_run_energy, data)

    if "conservation" in measurements:
        results["conservation"] = _run_measurement_safe(
            _run_conservation,
            data,
            threshold,
        )

    if "conversion" in measurements or "mixing" in measurements:
        try:
            conv = _run_conversion(data, source, target)
            results["conversion"] = conv

            if "mixing" in measurements:
                conv_result: ConversionResult = conv["_result_obj"]
                results["mixing"] = _run_measurement_safe(_run_mixing, conv_result)
        except ValueError as e:
            results["conversion"] = {"error": str(e)}
            if "mixing" in measurements:
                results["mixing"] = {"error": f"conversion failed: {e}"}

    if "spectrum" in measurements:
        results["spectrum"] = _run_measurement_safe(_run_spectrum, data)

    if "dispersion" in measurements:
        dyn_in_source = _filter_to_dynamical(source, data, "dispersion")
        if not dyn_in_source:
            from tidal.cli._console import error_with_hint

            error_with_hint(
                "no dynamical fields for dispersion",
                ["Check spec with `tidal inspect <json>`"],
            )
            return 1
        results["dispersion"] = _run_measurement_safe(
            _run_dispersion,
            data,
            dyn_in_source,
        )

    if "effective_mass" in measurements:
        dyn_in_source = _filter_to_dynamical(source, data, "effective_mass")
        results["effective_mass"] = _run_measurement_safe(
            _run_effective_mass,
            data,
            dyn_in_source,
        )

    if "asymptotic" in measurements:
        results["asymptotic"] = _run_measurement_safe(
            _run_asymptotic,
            data,
            source,
            target,
        )

    if "peak_conversion" in measurements:
        results["peak_conversion"] = _run_measurement_safe(
            _run_peak_conversion,
            data,
            source,
            target,
        )

    if "velocity" in measurements:
        dyn_in_source = _filter_to_dynamical(source, data, "velocity")
        results["velocity"] = _run_measurement_safe(_run_velocity, data, dyn_in_source)
        if "error" not in results["velocity"] and target is not None:
            results["velocity_mismatch"] = _run_measurement_safe(
                _run_velocity_mismatch,
                data,
                tuple(dyn_in_source),
                target,
            )

    if "resonance" in measurements:
        results["resonance"] = _run_measurement_safe(
            _run_resonance,
            data,
            source,
            target,
        )

    return results


def measure_command(args: Namespace) -> int:  # noqa: C901, PLR0911, PLR0912, PLR0915
    """Run ``tidal measure`` subcommand.

    Flow:
    1. Validate data path exists
    2. Resolve JSON spec path (--spec or metadata.json)
    3. Parse --what measurement set
    4. Load SimulationData
    5. Run requested measurements (dependency-ordered)
    6. Format and output (text, JSON, or plot)

    Returns 0 on success, 1 on error.
    """
    if getattr(args, "list_types", False):
        print("Available measurement types:")
        print("  summary              Overview of all measurements")
        print("  energy               Per-field and total energy")
        print("  conservation         Energy conservation check (PASS/FAIL)")
        print("  conversion           Field conversion probability")
        print("  mixing               Mixing length and dominant frequency")
        print("  spectrum             Power spectrum analysis")
        print("  dispersion           Dispersion relation extraction")
        print("  effective_mass       Effective mass from dispersion")
        print("  asymptotic           Late-time asymptotic behavior")
        print("  peak_conversion      Peak conversion probability")
        print("  velocity             Group/phase velocity measurement")
        print("  resonance            Resonance detection")
        return 0

    from tidal.cli._console import error_with_hint

    if args.data_path is None:
        from tidal.cli._console import error as _error

        _error("data_path is required")
        return 1

    data_path = Path(args.data_path)
    if not data_path.exists():
        error_with_hint(
            f"data path not found: {data_path}",
            ["Run `tidal simulate ... --output DIR` first, then `tidal measure DIR`"],
        )
        return 1

    try:
        spec_path = _resolve_spec_path(data_path, getattr(args, "spec", None))
    except FileNotFoundError:
        error_with_hint(
            f"spec file not found for {data_path.name}",
            ["Use `tidal list` to find specs"],
        )
        return 1
    except ValueError:
        error_with_hint(
            f"cannot determine spec for {data_path.name}",
            ["Provide explicitly: `tidal measure DIR --spec spec.json`"],
        )
        return 1

    try:
        measurements = _parse_measurements(getattr(args, "what", None))
    except ValueError:
        error_with_hint(
            f"unknown measurement type in --what={getattr(args, 'what', '')}",
            ["Available: energy, conversion, mixing, spectrum, conservation, etc."],
        )
        return 1
    quiet: bool = getattr(args, "quiet", False)

    if not quiet:
        print(f"Loading: {data_path.name}")
        print(f"Spec:    {spec_path.name}")

    data = _load_data(data_path, spec_path, getattr(args, "param", None) or [])

    if not quiet:
        print(
            f"  {data.n_snapshots} snapshots, "
            f"{len(data.fields)} fields, "
            f"t=[{float(data.times[0]):.1f}, {float(data.times[-1]):.1f}]",
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
            print(f"  Saved plot to: {Path(output_path).resolve()}")
    elif json_mode:
        print(_format_json(results, data))
    else:
        print(_format_text(results, data))

    return 0
