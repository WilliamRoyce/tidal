"""``tidal sweep --config`` — TOML configuration for parameter sweeps.

Loads a TOML file that declaratively specifies the full sweep setup:
parameters, sweep ranges, simulation settings, measurements, output,
and execution options.  CLI flags override TOML values.
"""

from __future__ import annotations

import sys
import tomllib
from dataclasses import dataclass, field
from pathlib import Path  # noqa: TC003 — used at runtime
from typing import TYPE_CHECKING, Any, cast

import numpy as np

if TYPE_CHECKING:
    from argparse import Namespace


# ---------------------------------------------------------------------------
# Configuration data model
# ---------------------------------------------------------------------------


@dataclass
class SweepConfig:
    """Parsed sweep configuration from TOML."""

    spec_path: Path
    swept_params: dict[str, list[float]]
    fixed_params: dict[str, float]
    sim_settings: dict[str, Any]
    measurements: list[str]
    source: str | None
    target: str | None
    output: Path
    parallel: int | None
    resume: bool
    energy_threshold: float
    converge_sizes: list[int] | None
    force_large_sweep: bool
    dry_run: bool
    adaptive_config: dict[str, dict[str, Any]] = field(
        default_factory=lambda: cast("dict[str, dict[str, Any]]", {}),
    )
    sweep_strategy: str | None = None
    n_samples: int | None = None
    n_replicates: int = 1
    base_seed: int = 42
    ic_perturbation: float | None = None
    param_noise: dict[str, float] | None = None


_KNOWN_SECTIONS = frozenset(
    {
        "spec",
        "parameters",
        "sweep",
        "simulation",
        "measurement",
        "output",
        "execution",
        "convergence",
        "ensemble",
    },
)

_SIM_KEY_MAP: dict[str, str] = {
    "grid_shape": "grid_shape",
    "bounds": "bounds",
    "periodic": "periodic",
    "bc": "bc",
    "ic": "ic",
    "ic_formula": "ic_formula",
    "ic_center": "ic_center",
    "ic_width": "ic_width",
    "ic_amplitude": "ic_amplitude",
    "ic_component": "ic_component",
    "ic_wavevector": "ic_wavevector",
    "ic_formula_velocity": "ic_formula_velocity",
    "ic_field": "ic_field",
    "ic_file": "ic_file",
    "ic_noise_seed": "ic_noise_seed",
    "t_end": "t_end",
    "dt": "dt",
    "scheme": "scheme",
    "rtol": "rtol",
    "atol": "atol",
    "method": "method",
    "max_step": "max_step",
    "snapshots": "snapshots",
    "mode": "mode",
    "allow_inconsistent_ic": "allow_inconsistent_ic",
}

# ---------------------------------------------------------------------------
# TOML section parsers
# ---------------------------------------------------------------------------

# Keys that should be converted from TOML list to CLI-style string
_LIST_TO_COLON_KEYS = frozenset({"bounds"})
_LIST_TO_COMMA_KEYS = frozenset({"grid_shape", "ic_center", "ic_wavevector"})


def _parse_sweep_section(
    name: str,
    section: dict[str, Any],
    *,
    strategy: str | None = None,
) -> tuple[list[float], dict[str, Any]]:
    """Parse a ``[sweep.NAME]`` TOML section into values + adaptive config.

    Returns
    -------
    tuple[list[float], dict[str, Any]]
        (parameter values, adaptive_config dict — empty if not adaptive)

    Raises
    ------
    ValueError
        If the section is malformed.
    """
    if "values" in section:
        vals = cast("list[Any]", section["values"])
        if not isinstance(vals, list) or len(vals) < 2:  # noqa: PLR2004  # type: ignore[reportUnnecessaryIsInstance]
            msg = f"[sweep.{name}].values must be a list with >= 2 entries"
            raise ValueError(msg)
        return [float(v) for v in vals], {}

    for key in ("start", "stop"):
        if key not in section:
            msg = f"[sweep.{name}] missing required key '{key}'"
            raise ValueError(msg)

    start = float(section["start"])
    stop = float(section["stop"])
    scale = str(section.get("scale", "linear")).lower()

    if scale == "adaptive":
        return _parse_adaptive(name, section, start, stop)

    # Space-filling strategies (LHS, Sobol) only need bounds — count
    # comes from n_samples at the sweep level.
    if strategy in {"latin_hypercube", "sobol"} and "count" not in section:
        return [start, stop], {}

    return _parse_range(name, section, start, stop, scale)


def _parse_adaptive(
    name: str,
    section: dict[str, Any],
    start: float,
    stop: float,
) -> tuple[list[float], dict[str, Any]]:
    """Parse adaptive sweep section.

    Raises
    ------
    ValueError
        If initial_count < 2.
    """
    initial_count = int(section.get("initial_count", 5))
    if initial_count < 2:  # noqa: PLR2004
        msg = f"[sweep.{name}].initial_count must be >= 2"
        raise ValueError(msg)
    values = np.linspace(start, stop, initial_count).tolist()
    adaptive = {
        "bounds": (start, stop),
        "initial_count": initial_count,
        "max_count": int(section.get("max_count", 20)),
        "metric": section.get("metric"),
        "threshold": float(section.get("threshold", 0.01)),
    }
    return values, adaptive


def _parse_range(
    name: str,
    section: dict[str, Any],
    start: float,
    stop: float,
    scale: str,
) -> tuple[list[float], dict[str, Any]]:
    """Parse linear/log range sweep section.

    Raises
    ------
    ValueError
        If count is missing, too small, or scale is unknown.
    """
    if "count" not in section:
        msg = f"[sweep.{name}] missing required key 'count' (or use 'values' for explicit list)"
        raise ValueError(msg)
    count = int(section["count"])
    if count < 2:  # noqa: PLR2004
        msg = f"[sweep.{name}].count must be >= 2, got {count}"
        raise ValueError(msg)
    if scale == "linear":
        return np.linspace(start, stop, count).tolist(), {}
    if scale == "log":
        if start <= 0 or stop <= 0:
            msg = (
                f"[sweep.{name}] log scale requires positive bounds, got {start}:{stop}"
            )
            raise ValueError(msg)
        return cast(
            "list[float]",
            np.logspace(np.log10(start), np.log10(stop), count).tolist(),
        ), {}
    msg = f"[sweep.{name}] unknown scale '{scale}' (expected: linear, log, adaptive)"
    raise ValueError(msg)


_SWEEP_META_KEYS = frozenset({"strategy", "n_samples"})


def _parse_sweeps(
    data: dict[str, Any],
) -> tuple[dict[str, list[float]], dict[str, dict[str, Any]], str | None, int | None]:
    """Parse all ``[sweep.*]`` sections plus sweep-level metadata.

    Returns
    -------
    tuple
        (swept_params, adaptive_config, strategy, n_samples)
    """
    swept: dict[str, list[float]] = {}
    adaptive: dict[str, dict[str, Any]] = {}
    strategy: str | None = None
    n_samples: int | None = None
    raw = data.get("sweep")
    if raw and isinstance(raw, dict):
        raw_dict = cast("dict[str, Any]", raw)
        strategy = cast("str | None", raw_dict.get("strategy"))
        n_samples_raw = raw_dict.get("n_samples")
        if n_samples_raw is not None:
            n_samples = int(n_samples_raw)
        for name, section in raw_dict.items():
            if name in _SWEEP_META_KEYS:
                continue
            if not isinstance(section, dict):
                continue
            values, ac = _parse_sweep_section(
                name,
                cast("dict[str, Any]", section),
                strategy=strategy,
            )
            swept[name] = values
            if ac:
                adaptive[name] = ac
    return swept, adaptive, strategy, n_samples


def _parse_convergence(
    data: dict[str, Any],
    *,
    has_sweeps: bool,
) -> list[int] | None:
    """Parse [convergence] section.

    Raises
    ------
    ValueError
        If grid_sizes has < 2 entries or both sweep and convergence are set.
    """
    if "convergence" not in data:
        return None
    conv = data["convergence"]
    if "grid_sizes" not in conv:
        return None
    sizes = sorted(int(s) for s in conv["grid_sizes"])
    if len(sizes) < 2:  # noqa: PLR2004
        msg = "convergence.grid_sizes must have >= 2 entries"
        raise ValueError(msg)
    if has_sweeps:
        msg = "[sweep.*] and [convergence] are mutually exclusive"
        raise ValueError(msg)
    return sizes


def _parse_measurement(
    data: dict[str, Any],
) -> tuple[list[str], str | None, str | None, float]:
    """Parse [measurement] section."""
    if "measurement" not in data:
        return [], None, None, 1e-3
    meas = data["measurement"]
    return (
        list(meas.get("types", [])),
        meas.get("source"),
        meas.get("target"),
        float(meas.get("energy_threshold", 1e-3)),
    )


def _parse_output(data: dict[str, Any], toml_dir: Path) -> Path:
    """Parse [output] section."""
    if "output" not in data:
        return toml_dir / "sweep_output"
    out = data["output"]
    if isinstance(out, str):
        return toml_dir / out
    if isinstance(out, dict) and "path" in out:
        return toml_dir / str(cast("dict[str, Any]", out)["path"])
    return toml_dir / "sweep_output"


def _parse_execution(
    data: dict[str, Any],
) -> tuple[int | None, bool, bool, bool, int, int]:
    """Parse [execution] section.

    Returns (parallel, resume, force_large, dry_run, n_replicates, base_seed).
    """
    if "execution" not in data:
        return None, False, False, False, 1, 42
    exe = data["execution"]
    return (
        exe.get("parallel"),
        bool(exe.get("resume", False)),
        bool(exe.get("force_large_sweep", False)),
        bool(exe.get("dry_run", False)),
        int(exe.get("n_replicates", 1)),
        int(exe.get("base_seed", 42)),
    )


def _parse_ensemble(
    data: dict[str, Any],
) -> tuple[float | None, dict[str, float] | None]:
    """Parse [ensemble] section for IC perturbation and parameter noise.

    Returns (ic_perturbation, param_noise).
    """
    if "ensemble" not in data:
        return None, None
    ens = data["ensemble"]
    ic_pert = ens.get("ic_perturbation")
    if ic_pert is not None:
        ic_pert = float(ic_pert)
    param_noise: dict[str, float] | None = None
    raw_pn = ens.get("param_noise")
    if raw_pn and isinstance(raw_pn, dict):
        param_noise = {str(k): float(v) for k, v in raw_pn.items()}  # type: ignore[reportUnknownVariableType, reportUnknownArgumentType]
    return ic_pert, param_noise


def _parse_sim_settings(data: dict[str, Any], filename: str) -> dict[str, Any]:
    """Parse [simulation] section, warning on unknown keys."""
    sim_settings: dict[str, Any] = {}
    for k, v in data.get("simulation", {}).items():
        if k in _SIM_KEY_MAP:
            sim_settings[k] = v
        else:
            print(
                f"  Warning: unknown simulation key '{k}' in {filename}",
                file=sys.stderr,
            )
    return sim_settings


# ---------------------------------------------------------------------------
# Loader and CLI integration
# ---------------------------------------------------------------------------


def load_sweep_config(path: Path) -> SweepConfig:
    """Load and validate a sweep TOML configuration file.

    Parameters
    ----------
    path : Path
        Path to the TOML configuration file.

    Returns
    -------
    SweepConfig
        Parsed configuration.

    Raises
    ------
    FileNotFoundError
        If the TOML file doesn't exist.
    ValueError
        If required keys are missing or values are invalid.
    """
    if not path.exists():
        msg = f"Config file not found: {path}"
        raise FileNotFoundError(msg)

    with path.open("rb") as f:
        data = tomllib.load(f)

    toml_dir = path.parent

    for key in data:
        if key not in _KNOWN_SECTIONS:
            print(
                f"  Warning: unknown TOML section '{key}' in {path.name}",
                file=sys.stderr,
            )

    if "spec" not in data:
        msg = f"Missing required key 'spec' in {path.name}"
        raise ValueError(msg)

    swept_params, adaptive_config, strategy, n_samples = _parse_sweeps(data)
    converge_sizes = _parse_convergence(data, has_sweeps=bool(swept_params))

    if not swept_params and converge_sizes is None:
        msg = f"No [sweep.*] sections or [convergence] found in {path.name}"
        raise ValueError(msg)

    meas, source, target, energy_threshold = _parse_measurement(data)
    parallel, resume, force_large, dry_run, n_replicates, base_seed = _parse_execution(
        data,
    )
    ic_perturbation, param_noise = _parse_ensemble(data)

    return SweepConfig(
        spec_path=toml_dir / data["spec"],
        swept_params=swept_params,
        fixed_params={k: float(v) for k, v in data.get("parameters", {}).items()},
        sim_settings=_parse_sim_settings(data, path.name),
        measurements=meas,
        source=source,
        target=target,
        output=_parse_output(data, toml_dir),
        parallel=parallel,
        resume=resume,
        energy_threshold=energy_threshold,
        converge_sizes=converge_sizes,
        force_large_sweep=force_large,
        dry_run=dry_run,
        adaptive_config=adaptive_config,
        sweep_strategy=strategy,
        n_samples=n_samples,
        n_replicates=n_replicates,
        base_seed=base_seed,
        ic_perturbation=ic_perturbation,
        param_noise=param_noise,
    )


def _convert_toml_val(toml_key: str, toml_val: Any) -> Any:  # noqa: ANN401
    """Convert TOML types to CLI-compatible types."""
    if toml_key in _LIST_TO_COLON_KEYS and isinstance(toml_val, list):
        return ":".join(str(x) for x in cast("list[Any]", toml_val))
    if toml_key in _LIST_TO_COMMA_KEYS:
        if isinstance(toml_val, (list, tuple)):
            return ",".join(str(x) for x in cast("list[Any]", toml_val))
        return str(toml_val)
    return toml_val


def _apply_sim_settings(config: SweepConfig, args: Namespace) -> None:
    """Fill simulation settings from config into args (CLI wins)."""
    for toml_key, attr_name in _SIM_KEY_MAP.items():
        if toml_key not in config.sim_settings:
            continue
        toml_val = config.sim_settings[toml_key]

        if toml_key == "ic_field" and isinstance(toml_val, list):
            current: list[Any] = getattr(args, "ic_field", None) or []
            if not current:
                args.ic_field = toml_val
            continue

        toml_val = _convert_toml_val(toml_key, toml_val)
        cli_attr = attr_name.replace("-", "_")
        if getattr(args, cli_attr, None) is None:
            setattr(args, cli_attr, toml_val)


def _apply_execution_settings(config: SweepConfig, args: Namespace) -> None:  # noqa: C901
    """Fill execution settings from config into args (CLI wins)."""
    if getattr(args, "parallel", None) is None and config.parallel is not None:
        args.parallel = config.parallel
    if not getattr(args, "resume", False) and config.resume:
        args.resume = True
    if not getattr(args, "force_large_sweep", False) and config.force_large_sweep:
        args.force_large_sweep = True
    if not getattr(args, "dry_run", False) and config.dry_run:
        args.dry_run = True
    if config.adaptive_config:
        args._adaptive_config = config.adaptive_config  # noqa: SLF001
    if getattr(args, "sweep_strategy", None) is None and config.sweep_strategy:
        args.sweep_strategy = config.sweep_strategy
    if getattr(args, "n_samples", None) is None and config.n_samples is not None:
        args.n_samples = config.n_samples
    # Ensemble / replicate settings (CLI wins)
    if getattr(args, "n_replicates", None) is None and config.n_replicates > 1:
        args.n_replicates = config.n_replicates
    if getattr(args, "base_seed", None) is None and config.base_seed != 42:  # noqa: PLR2004
        args.base_seed = config.base_seed
    if (
        getattr(args, "ic_perturbation", None) is None
        and config.ic_perturbation is not None
    ):
        args.ic_perturbation = config.ic_perturbation
    if getattr(args, "param_noise", None) is None and config.param_noise:
        # Convert dict to CLI-style list: ["B0=0.01", "m2=0.05"]
        args.param_noise = [f"{k}={v}" for k, v in config.param_noise.items()]


def apply_config_to_args(
    config: SweepConfig,
    args: Namespace,
) -> tuple[dict[str, list[float]], list[int] | None]:
    """Merge SweepConfig into an argparse Namespace, respecting CLI overrides.

    CLI flags that were explicitly set (non-default) take priority over TOML.
    TOML values fill in anything not set on the CLI.

    Parameters
    ----------
    config : SweepConfig
        Parsed TOML configuration.
    args : Namespace
        Parsed CLI arguments (mutated in place).

    Returns
    -------
    tuple[dict[str, list[float]], list[int] | None]
        (swept_params, converge_sizes) from TOML.
    """
    if not getattr(args, "json_path", None):
        args.json_path = str(config.spec_path)
    if not getattr(args, "output", None):
        args.output = str(config.output)

    # Fixed parameters: merge (CLI --param overrides TOML)
    existing_params: list[str] = getattr(args, "param", None) or []
    cli_param_names: set[str] = {
        p.split("=", 1)[0].strip() for p in existing_params if "=" in p
    }
    for name, val in config.fixed_params.items():
        if name not in cli_param_names:
            existing_params.append(f"{name}={val}")
    args.param = existing_params

    _apply_sim_settings(config, args)

    if not getattr(args, "measure", None) and config.measurements:
        args.measure = ",".join(config.measurements)
    if not getattr(args, "source", None) and config.source:
        args.source = config.source
    if not getattr(args, "target", None) and config.target:
        args.target = config.target
    if getattr(args, "energy_threshold", None) is None:
        args.energy_threshold = config.energy_threshold

    _apply_execution_settings(config, args)

    return config.swept_params, config.converge_sizes
