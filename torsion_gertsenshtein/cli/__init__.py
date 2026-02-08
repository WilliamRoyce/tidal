"""Command-line interface for the Lagrangian-to-PDE pipeline.

Entry point: ``tg`` command with subcommands:

- ``tg derive``   — Derive equations from Lagrangian via Wolfram/xAct
- ``tg inspect``  — Display equation system information from JSON
- ``tg simulate`` — Run PDE simulation from JSON specification
- ``tg list``     — List available JSON specifications
- ``tg validate`` — Validate a JSON equation specification
"""

from __future__ import annotations

import argparse
import sys

__all__ = ["main"]


def _build_parser() -> argparse.ArgumentParser:
    """Build the top-level argument parser with subcommands."""
    parser = argparse.ArgumentParser(
        prog="tg",
        description="Lagrangian-to-PDE pipeline: derive, inspect, simulate, list.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {_get_version()}",
    )
    sub = parser.add_subparsers(dest="command", help="Available commands")

    # --- derive ---
    derive_parser = sub.add_parser(
        "derive",
        help="Derive equations from Lagrangian via Wolfram/xAct",
        description=(
            "Generate equations of motion from a Lagrangian. "
            "Accepts a TOML config file (.toml) or a Wolfram script (.wls)."
        ),
        epilog=(
            "Examples:\n"
            "  tg derive theory.toml                     # run derivation\n"
            "  tg derive theory.toml --dry-run            # preview .wls without running\n"
            "  tg derive theory.toml --save-script eq.wls # save generated script\n"
            "  tg derive script.wls                       # run existing .wls directly"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    derive_parser.add_argument(
        "config",
        help="Path to .toml config or .wls script (auto-detected by extension)",
    )
    derive_parser.add_argument(
        "--save-script",
        metavar="PATH",
        default=None,
        help="Save generated .wls script to this path (TOML mode only)",
    )
    derive_parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print generated .wls without running wolframscript",
    )
    derive_parser.add_argument(
        "--output",
        metavar="PATH",
        default=None,
        help="Override output JSON path from config",
    )

    # --- inspect ---
    inspect_parser = sub.add_parser(
        "inspect",
        help="Display equation system information from JSON",
        description="Load a JSON specification and display its contents.",
        epilog=(
            "Examples:\n"
            "  tg inspect examples/data/klein_gordon_1d.json\n"
            "  tg inspect spec.json --params    # show default parameter values"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    inspect_parser.add_argument(
        "json_path",
        help="Path to the JSON equation specification",
    )
    inspect_parser.add_argument(
        "--params",
        action="store_true",
        help="Show default parameter values from metadata",
    )
    inspect_parser.add_argument(
        "--json",
        action="store_true",
        dest="json_output",
        help="Output machine-readable JSON instead of human-readable text",
    )

    # --- simulate ---
    sim_parser = sub.add_parser(
        "simulate",
        help="Run PDE simulation from JSON specification",
        description="Build and run a PDE simulation from a JSON equation specification.",
        epilog=(
            "Examples:\n"
            "  tg simulate spec.json --param m2=1.0 --t-end 10\n"
            "  tg simulate spec.json --ic gaussian --ic-width 2.0 --output result.png\n"
            "  tg simulate spec.json --mode constraint --bc neumann\n"
            "  tg simulate spec.json --ic formula --ic-formula 'exp(-(x-5)**2)'"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sim_parser.add_argument(
        "json_path",
        help="Path to the JSON equation specification",
    )
    # Parameters
    sim_parser.add_argument(
        "--param",
        action="append",
        default=[],
        metavar="KEY=VAL",
        help="Set symbolic parameter (repeatable, e.g. --param m2=1.0)",
    )
    # Grid
    sim_parser.add_argument(
        "--grid-shape",
        default=None,
        metavar="N[,N,N]",
        help="Grid points per axis (e.g. 32 or 32,32,32). Default: auto from dimension",
    )
    sim_parser.add_argument(
        "--bounds",
        default=None,
        metavar="LO:HI[,...]",
        help="Domain bounds per axis (e.g. 0:20 or 0:20,0:10). Default: 0:10",
    )
    sim_parser.add_argument(
        "--periodic",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Use periodic boundary conditions (default: True)",
    )
    sim_parser.add_argument(
        "--bc",
        default=None,
        metavar="BC[,BC,BC]",
        help="Per-axis boundary conditions: periodic, neumann (comma-separated). Overrides --periodic.",
    )
    # Initial conditions
    sim_parser.add_argument(
        "--ic",
        choices=["gaussian", "plane-wave", "zero", "formula"],
        default="gaussian",
        help="Initial condition type (default: gaussian)",
    )
    sim_parser.add_argument(
        "--ic-formula",
        default=None,
        metavar="EXPR",
        help="Math expression for --ic=formula. Variables: coordinate names (e.g. x,y,z), np (numpy), pi.",
    )
    sim_parser.add_argument(
        "--ic-center", default=None, metavar="X[,X,X]",
        help="Gaussian center position (default: domain midpoint)",
    )
    sim_parser.add_argument(
        "--ic-width", type=float, default=None, metavar="W",
        help="Gaussian width (default: domain_size/10)",
    )
    sim_parser.add_argument(
        "--ic-amplitude", type=float, default=1.0, metavar="A",
        help="IC peak amplitude (default: 1.0)",
    )
    sim_parser.add_argument(
        "--ic-component", default=None, metavar="NAME",
        help="Field component for IC (default: first field)",
    )
    sim_parser.add_argument(
        "--ic-wavevector", default=None, metavar="K[,K,K]",
        help="Wavevector for plane-wave IC (e.g. 0.1,0.0,0.0)",
    )
    # Mode
    sim_parser.add_argument(
        "--mode",
        choices=["evolve", "constraint"],
        default="evolve",
        help="'evolve' = time evolution (default), 'constraint' = single constraint solve",
    )
    # Solver
    sim_parser.add_argument(
        "--t-end", type=float, default=10.0,
        help="Simulation duration (default: 10.0)",
    )
    sim_parser.add_argument(
        "--dt", type=float, default=None,
        help="Time step (default: auto from CFL)",
    )
    sim_parser.add_argument(
        "--scheme", choices=["runge-kutta", "scipy"], default="runge-kutta",
        help="Solver scheme (default: runge-kutta)",
    )
    sim_parser.add_argument(
        "--snapshots", type=float, default=None, metavar="DT",
        help="Snapshot interval (default: t_end/100)",
    )
    # Output
    sim_parser.add_argument(
        "--output", default=None, metavar="PATH",
        help="Output file path (default: {spec_dir}/{stem}_output.png)",
    )
    sim_parser.add_argument(
        "--format",
        choices=["png", "npz", "summary"],
        default=None,
        dest="output_format",
        help="Output format (default: inferred from --output extension)",
    )
    sim_parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Print summary only, skip visualization",
    )
    sim_parser.add_argument(
        "--quiet", "-q",
        action="store_true",
        help="Suppress progress messages (results and errors still shown)",
    )

    # --- list ---
    list_parser = sub.add_parser(
        "list",
        help="List available JSON specifications",
        description="Scan a directory for JSON equation specifications and display summaries.",
        epilog=(
            "Examples:\n"
            "  tg list                          # scan default examples/data/\n"
            "  tg list --dir /path/to/specs      # scan custom directory"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    list_parser.add_argument(
        "--dir",
        default=None,
        metavar="PATH",
        help="Directory to scan (default: examples/data/)",
    )

    # --- validate ---
    validate_parser = sub.add_parser(
        "validate",
        help="Validate a JSON equation specification",
        description="Check a JSON specification for errors (unknown operators, bad references, etc.).",
        epilog=(
            "Examples:\n"
            "  tg validate examples/data/klein_gordon_1d.json\n"
            "  tg validate spec.json"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    validate_parser.add_argument(
        "json_path",
        help="Path to the JSON equation specification to validate",
    )

    return parser


def _get_version() -> str:
    """Return the package version string from installed metadata."""
    from importlib.metadata import PackageNotFoundError, version

    try:
        return version("torsion_gertsenshtein")
    except PackageNotFoundError:
        return "unknown"


def _dispatch(args: argparse.Namespace) -> int:
    """Lazily import and run the appropriate command handler.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed CLI arguments with ``command`` attribute.

    Returns
    -------
    int
        Exit code from the command handler.

    Raises
    ------
    ValueError
        If ``args.command`` is not a recognized subcommand.
    """
    if args.command == "derive":
        from torsion_gertsenshtein.cli._derive import derive_command

        return derive_command(args)
    if args.command == "inspect":
        from torsion_gertsenshtein.cli._inspect import inspect_command

        return inspect_command(args)
    if args.command == "simulate":
        from torsion_gertsenshtein.cli._simulate import simulate_command

        return simulate_command(args)
    if args.command == "list":
        from torsion_gertsenshtein.cli._list import list_command

        return list_command(args)
    if args.command == "validate":
        from torsion_gertsenshtein.cli._validate import validate_command

        return validate_command(args)
    msg = f"Unknown command: {args.command}"
    raise ValueError(msg)


def main(argv: list[str] | None = None) -> int:
    """CLI entry point.

    Parameters
    ----------
    argv : list[str] | None
        Command-line arguments. If None, uses sys.argv.

    Returns
    -------
    int
        Exit code (0 for success, non-zero for errors).
    """
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        return 0

    try:
        return _dispatch(args)
    except KeyboardInterrupt:
        print("\nInterrupted.")
        return 130
    except Exception as exc:  # noqa: BLE001
        print(f"Error: {exc}", file=sys.stderr)
        return 1
