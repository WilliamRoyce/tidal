"""Command-line interface for the Lagrangian-to-PDE pipeline.

Entry point: ``tg`` command with subcommands:

- ``tg derive``   — Derive equations from Lagrangian via Wolfram/xAct
- ``tg inspect``  — Display equation system information from JSON
- ``tg simulate`` — Run PDE simulation from JSON specification
- ``tg list``     — List available JSON specifications
"""

from __future__ import annotations

import argparse
import sys


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

    # --- simulate ---
    sim_parser = sub.add_parser(
        "simulate",
        help="Run PDE simulation from JSON specification",
        description="Build and run a PDE simulation from a JSON equation specification.",
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
    # Initial conditions
    sim_parser.add_argument(
        "--ic",
        choices=["gaussian", "plane-wave", "zero"],
        default="gaussian",
        help="Initial condition type (default: gaussian)",
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
        help="Snapshot interval (default: t_end/20)",
    )
    # Output
    sim_parser.add_argument(
        "--output", default=None, metavar="PATH",
        help="Output file path (default: outputs/{stem}_output.png)",
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

    # --- list ---
    list_parser = sub.add_parser(
        "list",
        help="List available JSON specifications",
        description="Scan a directory for JSON equation specifications and display summaries.",
    )
    list_parser.add_argument(
        "--dir",
        default=None,
        metavar="PATH",
        help="Directory to scan (default: examples/data/)",
    )

    return parser


def _get_version() -> str:
    """Return the package version string."""
    try:
        from torsion_gertsenshtein import __version__

        return str(__version__)
    except (ImportError, AttributeError):
        return "unknown"


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
    except KeyboardInterrupt:
        print("\nInterrupted.")
        return 130
    except Exception as exc:  # noqa: BLE001
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    return 0
