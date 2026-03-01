"""Command-line interface for the Lagrangian-to-PDE pipeline.

Entry point: ``tidal`` command with subcommands:

- ``tidal derive``   — Derive equations from Lagrangian via Wolfram/xAct
- ``tidal inspect``  — Display equation system information from JSON
- ``tidal simulate`` — Run PDE simulation from JSON specification
- ``tidal measure``  — Extract physics measurements from simulation output
- ``tidal plot``     — Generate individual plots from simulation output
- ``tidal list``     — List available JSON specifications
- ``tidal validate`` — Validate a JSON equation specification
"""

from __future__ import annotations

import argparse
import sys

from tidal.solver._defaults import DEFAULT_ATOL, DEFAULT_RTOL

__all__ = ["main"]


def _build_parser() -> argparse.ArgumentParser:  # noqa: PLR0915
    """Build the top-level argument parser with subcommands."""
    parser = argparse.ArgumentParser(
        prog="tidal",
        description="Lagrangian-to-PDE pipeline: derive, inspect, simulate, measure, plot, list.",
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
            "  tidal derive theory.toml                     # run derivation\n"
            "  tidal derive theory.toml --dry-run            # preview .wls without running\n"
            "  tidal derive theory.toml --save-script eq.wls # save generated script\n"
            "  tidal derive script.wls                       # run existing .wls directly"
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
            "  tidal inspect examples/data/klein_gordon_1d.json\n"
            "  tidal inspect spec.json --params    # show default parameter values"
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
            "  tidal simulate spec.json --param m2=1.0 --t-end 10\n"
            "  tidal simulate spec.json --ic gaussian --ic-width 2.0 --output result.png\n"
            "  tidal simulate spec.json --mode constraint --bc neumann\n"
            "  tidal simulate spec.json --ic formula --ic-formula 'exp(-(x-5)**2)'"
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
        choices=["gaussian", "plane-wave", "zero", "formula", "file", "noise"],
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
        "--ic-center",
        default=None,
        metavar="X[,X,X]",
        help="Gaussian center position (default: domain midpoint)",
    )
    sim_parser.add_argument(
        "--ic-width",
        type=float,
        default=None,
        metavar="W",
        help="Gaussian width (default: domain_size/10)",
    )
    sim_parser.add_argument(
        "--ic-amplitude",
        type=float,
        default=1.0,
        metavar="A",
        help="IC peak amplitude (default: 1.0)",
    )
    sim_parser.add_argument(
        "--ic-component",
        default=None,
        metavar="NAME",
        help="Field component for IC (default: first field)",
    )
    sim_parser.add_argument(
        "--ic-wavevector",
        default=None,
        metavar="K[,K,K]",
        help="Wavevector for plane-wave or gaussian IC (e.g. 3 or 0.1,0.0,0.0). "
        "With gaussian: creates a travelling wave packet (positive k = right-mover)",
    )
    sim_parser.add_argument(
        "--ic-formula-velocity",
        default=None,
        metavar="EXPR",
        help="Velocity (time derivative) expression for --ic=formula. "
        "Same namespace as --ic-formula (x, y, z, sin, cos, exp, ...).",
    )
    sim_parser.add_argument(
        "--ic-field",
        action="append",
        default=[],
        metavar="FIELD:EXPR",
        help="Per-field IC formula override (repeatable). "
        "Format: FIELD:EXPR or FIELD:velocity:EXPR. "
        "Applied after --ic. Example: --ic-field 'chi_0:0.1*sin(x)'",
    )
    sim_parser.add_argument(
        "--ic-file",
        default=None,
        metavar="PATH",
        help="Path to .npy file or simulation output directory for --ic=file.",
    )
    sim_parser.add_argument(
        "--ic-noise-seed",
        type=int,
        default=None,
        metavar="N",
        help="Random seed for --ic=noise (reproducible noise).",
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
        "--t-end",
        type=float,
        default=10.0,
        help="Simulation duration (default: 10.0)",
    )
    sim_parser.add_argument(
        "--dt",
        type=float,
        default=None,
        help="Time step (default: auto from CFL)",
    )
    sim_parser.add_argument(
        "--scheme",
        choices=["ida", "leapfrog", "cvode", "scipy", "auto"],
        default="auto",
        help=(
            "Solver scheme (default: auto). "
            "'auto' selects the best adaptive solver: cvode for wave systems, "
            "ida for DAE/dissipative systems. "
            "'cvode' uses SUNDIALS/CVODE — adaptive BDF/Adams for wave systems. "
            "'ida' uses SUNDIALS/IDA — handles all equation types including constraints. "
            "'scipy' uses scipy.integrate.solve_ivp — DOP853/RK45/Radau/BDF. "
            "'leapfrog' uses symplectic Störmer-Verlet (fixed dt, zero energy drift)."
        ),
    )
    sim_parser.add_argument(
        "--rtol",
        type=float,
        default=DEFAULT_RTOL,
        help=f"Relative tolerance for adaptive solvers (default: {DEFAULT_RTOL})",
    )
    sim_parser.add_argument(
        "--atol",
        type=float,
        default=DEFAULT_ATOL,
        help=f"Absolute tolerance for adaptive solvers (default: {DEFAULT_ATOL})",
    )
    sim_parser.add_argument(
        "--method",
        type=str,
        default=None,
        help=(
            "Integration method (default: auto). "
            "cvode: 'BDF' (default) or 'Adams'. "
            "scipy: 'DOP853' (default), 'RK45', 'Radau', 'BDF', 'RK23', 'LSODA'."
        ),
    )
    sim_parser.add_argument(
        "--max-step",
        type=float,
        default=None,
        help=(
            "Maximum step size for adaptive solvers. "
            "Default: unbounded for cvode/ida, CFL dt for scipy."
        ),
    )
    sim_parser.add_argument(
        "--snapshots",
        type=float,
        default=None,
        metavar="DT",
        help="Snapshot interval (default: t_end/100)",
    )
    # Output
    sim_parser.add_argument(
        "--output",
        default=None,
        metavar="PATH",
        help="Output file path (default: {spec_dir}/{stem}_output.png)",
    )
    sim_parser.add_argument(
        "--format",
        choices=["png", "summary"],
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
        "--quiet",
        "-q",
        action="store_true",
        help="Suppress progress messages (results and errors still shown)",
    )
    sim_parser.add_argument(
        "--require-stable",
        action="store_true",
        default=False,
        help=(
            "Abort if the pre-simulation stability check detects an unstable "
            "mass matrix (negative eigenvalue). Default: warn only."
        ),
    )
    sim_parser.add_argument(
        "--allow-inconsistent-ic",
        action="store_true",
        default=False,
        help=(
            "Allow inconsistent initial conditions for constraint equations. "
            "When set, constraint violations produce warnings instead of errors. "
            "Default: error if constraint ICs cannot be made consistent."
        ),
    )

    # --- list ---
    list_parser = sub.add_parser(
        "list",
        help="List available JSON specifications",
        description="Scan a directory for JSON equation specifications and display summaries.",
        epilog=(
            "Examples:\n"
            "  tidal list                          # scan default examples/data/\n"
            "  tidal list --dir /path/to/specs      # scan custom directory"
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
            "  tidal validate examples/data/klein_gordon_1d.json\n"
            "  tidal validate spec.json"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    validate_parser.add_argument(
        "json_path",
        help="Path to the JSON equation specification to validate",
    )

    # --- measure ---
    measure_parser = sub.add_parser(
        "measure",
        help="Extract physics measurements from simulation output",
        description=(
            "Load simulation output from 'tidal simulate --output' and run "
            "measurement analyses (energy, conversion, mixing length, etc.)."
        ),
        epilog=(
            "Examples:\n"
            "  tidal measure result_dir/ --spec spec.json\n"
            "  tidal measure result_dir/ --json\n"
            "  tidal measure result_dir/ --what conversion --source phi_0 --target chi_0\n"
            "  tidal measure result_dir/ --what energy,conservation\n"
            "  tidal measure result_dir/ --output measurement.png"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    measure_parser.add_argument(
        "data_path",
        help="Path to the simulation output directory",
    )
    measure_parser.add_argument(
        "--spec",
        default=None,
        metavar="PATH",
        help="Path to JSON equation spec (auto-discovered from metadata.json if omitted)",
    )
    measure_parser.add_argument(
        "--what",
        default=None,
        metavar="TYPE[,TYPE,...]",
        help=(
            "Measurements to run (comma-separated). "
            "Options: summary, energy, conversion, mixing, spectrum, spectral_conversion, dispersion, conservation. "
            "Default: summary"
        ),
    )
    measure_parser.add_argument(
        "--source",
        default=None,
        metavar="FIELD[,FIELD,...]",
        help="Source field(s) for conversion measurement (comma-separated)",
    )
    measure_parser.add_argument(
        "--target",
        default=None,
        metavar="FIELD[,FIELD,...]",
        help="Target field(s) for conversion measurement (comma-separated)",
    )
    measure_parser.add_argument(
        "--param",
        action="append",
        default=[],
        metavar="KEY=VAL",
        help="Override parameter value (repeatable, e.g. --param m2=1.0)",
    )
    measure_parser.add_argument(
        "--energy-threshold",
        type=float,
        default=1e-3,
        metavar="T",
        help="Energy conservation threshold (default: 1e-3)",
    )
    measure_parser.add_argument(
        "--output",
        default=None,
        metavar="PATH",
        help="Save measurement plot (.png or .pdf)",
    )
    measure_parser.add_argument(
        "--json",
        action="store_true",
        dest="json_output",
        help="Output measurements as JSON instead of text",
    )
    measure_parser.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        help="Suppress progress messages",
    )

    # --- plot ---
    plot_parser = sub.add_parser(
        "plot",
        help="Generate individual plots from simulation output",
        description=(
            "Read disk-backed simulation output and produce a single focused plot. "
            "Users compose visualizations via multiple calls in shell scripts."
        ),
        epilog=(
            "Examples:\n"
            "  tidal plot output_dir/ --type heatmap --field phi_0\n"
            "  tidal plot output_dir/ --type snapshot --time-index -1\n"
            "  tidal plot output_dir/ --type amplitude --overlay 'exp(-0.1*t)'\n"
            "  tidal plot output_dir/ --type energy --fields phi_0,chi_0\n"
            "  tidal plot output_dir/ --type profile --cross-section y=25.0\n"
            "  tidal plot output_dir/ --type compare --output compare.png"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    plot_parser.add_argument(
        "data_dir",
        help="Path to simulation output directory (from 'tidal simulate --output')",
    )
    plot_parser.add_argument(
        "--type",
        required=True,
        choices=[
            "heatmap",
            "snapshot",
            "amplitude",
            "energy",
            "profile",
            "compare",
            "hamiltonian",
            "conservation",
        ],
        help="Plot type to generate",
    )
    # Field selection
    plot_parser.add_argument(
        "--field",
        default=None,
        metavar="NAME",
        help="Single field for heatmap/snapshot/profile (default: first field)",
    )
    plot_parser.add_argument(
        "--fields",
        default=None,
        metavar="NAME[,NAME,...]",
        help="Multiple fields for amplitude/energy/compare (default: all)",
    )
    # Time selection
    plot_parser.add_argument(
        "--time-index",
        type=int,
        default=None,
        metavar="N",
        help="Snapshot index for snapshot type (supports negative: -1 = last)",
    )
    plot_parser.add_argument(
        "--time-indices",
        default=None,
        metavar="N,N,...",
        help="Comma-separated time indices for profile (default: 5 evenly spaced)",
    )
    # Spatial slicing
    plot_parser.add_argument(
        "--cross-section",
        default=None,
        metavar="AXIS=VAL",
        help="Slice 2D data along axis (e.g. y=25.0) for profile/compare",
    )
    # Reference overlay
    plot_parser.add_argument(
        "--overlay",
        default=None,
        metavar="EXPR",
        help="Analytic formula vs time for amplitude plot (e.g. 'exp(-0.1*t)')",
    )
    # Conservation threshold
    plot_parser.add_argument(
        "--threshold",
        type=float,
        default=None,
        metavar="T",
        help="Conservation threshold for --type=conservation (default: 1e-3)",
    )
    # Output options
    plot_parser.add_argument(
        "--output",
        default=None,
        metavar="PATH",
        help="Output file path (default: DATA_DIR/{type}.png)",
    )
    plot_parser.add_argument(
        "--dpi",
        type=int,
        default=None,
        metavar="N",
        help="Output resolution in DPI (default: 150)",
    )
    plot_parser.add_argument(
        "--figsize",
        default=None,
        metavar="W,H",
        help="Figure size in inches (e.g. 8,6)",
    )
    plot_parser.add_argument(
        "--cmap",
        default=None,
        metavar="NAME",
        help="Colormap for heatmap/snapshot (default: RdBu_r)",
    )
    plot_parser.add_argument(
        "--title",
        default=None,
        help="Custom figure title",
    )
    # Metadata
    plot_parser.add_argument(
        "--spec",
        default=None,
        metavar="PATH",
        help="Override JSON spec path (auto-discovered from metadata.json)",
    )
    plot_parser.add_argument(
        "--param",
        action="append",
        default=[],
        metavar="KEY=VAL",
        help="Override parameter value (repeatable, e.g. --param m2=1.0)",
    )
    plot_parser.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        help="Suppress progress messages",
    )

    return parser


def _get_version() -> str:
    """Return the package version string from installed metadata."""
    from importlib.metadata import PackageNotFoundError, version

    try:
        return version("tidal")
    except PackageNotFoundError:
        return "unknown"


def _dispatch(args: argparse.Namespace) -> int:  # noqa: PLR0911
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
        from tidal.cli._derive import derive_command

        return derive_command(args)
    if args.command == "inspect":
        from tidal.cli._inspect import inspect_command

        return inspect_command(args)
    if args.command == "simulate":
        from tidal.cli._simulate import simulate_command

        return simulate_command(args)
    if args.command == "list":
        from tidal.cli._list import list_command

        return list_command(args)
    if args.command == "validate":
        from tidal.cli._validate import validate_command

        return validate_command(args)
    if args.command == "measure":
        from tidal.cli._measure import measure_command

        return measure_command(args)
    if args.command == "plot":
        from tidal.cli._plot_command import plot_command

        return plot_command(args)
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
    except (ValueError, FileNotFoundError, TypeError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    except RuntimeError:
        import traceback

        traceback.print_exc()
        return 1
