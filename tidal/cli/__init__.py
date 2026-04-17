"""Command-line interface for the Lagrangian-to-PDE pipeline.

Entry point: ``tidal`` command with subcommands:

- ``tidal derive``   — Derive equations from Lagrangian via Wolfram/xAct
- ``tidal inspect``  — Display equation system information from JSON
- ``tidal simulate`` — Run PDE simulation from JSON specification
- ``tidal measure``  — Extract physics measurements from simulation output
- ``tidal plot``     — Generate individual plots from simulation output
- ``tidal sweep``    — Run parameter sweeps and convergence studies
- ``tidal list``     — List available JSON specifications
- ``tidal validate`` — Validate a JSON equation specification
"""

from __future__ import annotations

import argparse
import os
import sys

from tidal.solver._defaults import DEFAULT_ATOL, DEFAULT_RTOL

__all__ = ["main"]


def _build_parser() -> argparse.ArgumentParser:  # noqa: PLR0915
    """Build the top-level argument parser with subcommands."""
    # Use rich-argparse for colored help if available, fallback to default
    try:
        from rich_argparse import RichHelpFormatter

        formatter = RichHelpFormatter
    except ImportError:
        formatter = argparse.HelpFormatter

    parser = argparse.ArgumentParser(
        prog="tidal",
        description="Lagrangian-to-PDE pipeline: derive, inspect, simulate, measure, plot, sweep, sample, list.",
        formatter_class=formatter,
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {_get_version()}",
    )
    parser.add_argument(
        "--no-banner",
        action="store_true",
        default=False,
        help="Suppress the startup banner",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        default=False,
        help="Show detailed diagnostic output (solver selection, CFL, Jacobian tier)",
    )
    parser.add_argument(
        "--cite",
        action="store_true",
        default=False,
        help="Print citation information for TIDAL and its dependencies",
    )
    parser.add_argument(
        "--time",
        action="store_true",
        default=False,
        help="Show wall-clock elapsed time after command completes",
    )
    sub = parser.add_subparsers(dest="command", help="Available commands")

    # --- derive ---
    derive_parser = sub.add_parser(
        "derive",
        aliases=["der"],
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
        "--force-derive",
        action="store_true",
        help="Force re-derivation even if the generated script is unchanged",
    )
    derive_parser.add_argument(
        "--output",
        metavar="PATH",
        default=None,
        help="Override output JSON path from config",
    )
    derive_parser.add_argument(
        "--timeout",
        type=int,
        default=600,
        metavar="SECONDS",
        help=(
            "Maximum time (seconds) for wolframscript execution. "
            "Default: 600 (10 min). If exceeded, investigate which pipeline "
            "stage is the bottleneck and optimize. Use 0 for no timeout."
        ),
    )

    # --- inspect ---
    inspect_parser = sub.add_parser(
        "inspect",
        aliases=["insp"],
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
    inspect_parser.add_argument(
        "--latex",
        action="store_true",
        help="Output equations as LaTeX math",
    )
    inspect_parser.add_argument(
        "--latex-format",
        choices=["align", "document", "raw"],
        default="align",
        dest="latex_format",
        help="LaTeX output format (default: align). 'document' wraps in standalone .tex",
    )

    # --- simulate ---
    sim_parser = sub.add_parser(
        "simulate",
        aliases=["sim"],
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
        nargs="?",
        default=None,
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
        choices=["ida", "leapfrog", "cvode", "scipy", "modal", "auto"],
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
    sim_parser.add_argument(
        "--fd-order",
        type=int,
        choices=[2, 4, 6],
        default=4,
        help=(
            "Finite-difference accuracy order for spatial operators (default: 4). "
            "Higher orders use wider stencils: 2->3pt, 4->5pt, 6->7pt. "
            "Order 4 is recommended: O(dx^4) convergence at ~2x per-point cost "
            "vs order 2 — typically enables halving N for the same accuracy. "
            "Ref: Fornberg (1988), Mathematics of Computation 51(184)."
        ),
    )
    sim_parser.add_argument(
        "--leapfrog-order",
        type=int,
        choices=[2, 4],
        default=None,
        help=(
            "Leapfrog integrator order. Default: auto-detect (4 for time-independent "
            "non-dissipative systems, 2 otherwise). "
            "2 = Störmer-Verlet (1 force eval/step). "
            "4 = Yoshida triple-composition (3 force evals/step, O(dt⁴) accuracy). "
            "Yoshida allows ~2x larger dt for the same accuracy. "
            "Only applies when --scheme leapfrog is used. "
            "Ref: Yoshida (1990), Physics Letters A 150(5-7), pp. 262-268."
        ),
    )
    sim_parser.add_argument(
        "--spectral",
        action=argparse.BooleanOptionalAction,
        default=None,
        help=(
            "Use FFT-based spectral operators instead of finite-difference stencils. "
            "Default: auto-detect (enabled when all BCs are periodic and solver is "
            "not IDA). Use --no-spectral to force finite-difference stencils. "
            "Achieves exponential convergence for smooth problems -- typically "
            "N=64-128 instead of N=512-1024. "
            "Incompatible with --scheme ida (use leapfrog, cvode, or scipy). "
            "Ref: Burns et al. (2020), Phys. Rev. Research 2:023068."
        ),
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
        "--dry-run",
        action="store_true",
        default=False,
        help="Preview simulation setup without running (grid, solver, IC, memory estimate)",
    )
    sim_parser.add_argument(
        "--list-schemes",
        action="store_true",
        default=False,
        help="List available solver schemes and exit",
    )
    sim_parser.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Overwrite existing --output directory without prompting",
    )
    sim_parser.add_argument(
        "--report",
        type=str,
        default=None,
        metavar="PATH",
        help="Generate self-contained HTML summary report at PATH",
    )
    sim_parser.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        help="Suppress progress messages (results and errors still shown)",
    )
    # Note: --require-stable removed. The modal solver's eigenvalue pre-check
    # in _evolve_per_mode provides accurate instability detection using the
    # full evolution matrix. The old mass-only check produced false positives
    # for theories with non-unit kinetic coefficients (e.g., PGT torsion).
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
    # Resume from checkpoint
    sim_parser.add_argument(
        "--resume",
        default=None,
        metavar="DIR",
        help=(
            "Resume simulation from a snapshot directory. "
            "Inherits grid, parameters, and BC from saved metadata. "
            "Loads final snapshot by default (use --snapshot N to pick)."
        ),
    )
    sim_parser.add_argument(
        "--snapshot",
        type=int,
        default=None,
        metavar="N",
        help="Snapshot index to resume from (default: last). Requires --resume.",
    )
    sim_parser.add_argument(
        "--t-additional",
        type=float,
        default=None,
        metavar="T",
        help=(
            "Additional simulation time beyond the checkpoint "
            "(alternative to --t-end when using --resume)."
        ),
    )

    # --- list ---
    list_parser = sub.add_parser(
        "list",
        aliases=["ls"],
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
        aliases=["val"],
        help="Validate a JSON equation specification",
        description="Check a JSON specification for errors (unknown operators, bad references, etc.).",
        epilog=(
            "Examples:\n"
            "  tidal validate examples/data/klein_gordon_1d.json\n"
            "  tidal validate spec.json --stability\n"
            "  tidal validate spec.json --stability --param m2=1.0"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    validate_parser.add_argument(
        "json_path",
        help="Path to the JSON equation specification to validate",
    )
    validate_parser.add_argument(
        "--stability",
        action="store_true",
        default=False,
        help=(
            "Run stability checks: tachyon detection (negative mass-matrix eigenvalues) "
            "and ghost detection (wrong-sign kinetic terms from Hamiltonian). "
            "Uses a 1-point grid; exact for constant-coefficient systems."
        ),
    )
    validate_parser.add_argument(
        "--param",
        action="append",
        default=[],
        metavar="KEY=VAL",
        help="Set symbolic parameter for stability check (repeatable, e.g. --param m2=1.0)",
    )

    # --- measure ---
    measure_parser = sub.add_parser(
        "measure",
        aliases=["meas"],
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
        nargs="?",
        default=None,
        help="Path to the simulation output directory",
    )
    measure_parser.add_argument(
        "--list-types",
        action="store_true",
        default=False,
        help="List available measurement types and exit",
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
            "Options: summary, energy, conversion, mixing, spectrum, dispersion, conservation. "
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
            "  tidal plot output_dir/ --type compare --output compare.png\n"
            "  tidal plot sweep_dir/ --type sweep --metric P_final "
            "--overlay 'sin(kappa * B0 * t_end / 2)**2'\n"
            "  tidal plot sweep_dir/ --type sweep --metric P_final "
            "--overlay 'sin(kappa * Bpeak * R * sqrt(pi/2))**2'"
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
            "sweep",
            "sweep-compare",
            "sweep-parallel",
            "sweep-tornado",
            "sweep-scatter",
            "sweep-grouped",
            "sweep-histogram",
            "sweep-divergence",
            "campaign",
            "convergence",
            "replicate-convergence",
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
        help=(
            "Analytic reference formula for overlay. "
            "For --type=amplitude: variable 't' (time), e.g. 'exp(-0.1*t)'. "
            "For --type=sweep: swept parameter names, fixed parameters, and numeric "
            "sim settings (t_end, grid_shape) are available as variables. "
            "1D example: 'sin(kappa * B0 * t_end / 2)**2'. "
            "2D example: 'sin(kappa * Bpeak * R * sqrt(pi/2))**2' (3-panel figure)."
        ),
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
        "--log-scale",
        dest="log_scale",
        action="store_true",
        help="Use logarithmic colorbar for heatmaps (useful for inv_B_min)",
    )
    plot_parser.add_argument(
        "--log-y",
        dest="log_y",
        action="store_true",
        help="Use logarithmic y-axis for 1D sweep plots",
    )
    plot_parser.add_argument(
        "--log-x",
        dest="log_x",
        action="store_true",
        help="Use logarithmic x-axis for 1D sweep plots",
    )
    plot_parser.add_argument(
        "--x-param",
        dest="x_param",
        default=None,
        metavar="PARAM",
        help="Swept parameter for the x-axis (required for --type sweep-grouped)",
    )
    plot_parser.add_argument(
        "--group-by",
        dest="group_by",
        default=None,
        metavar="PARAM",
        help=(
            "Swept parameter to group 1D lines by (required for --type sweep-grouped). "
            "Produces one connected line per unique value of the grouping parameter."
        ),
    )
    plot_parser.add_argument(
        "--hline",
        action="append",
        default=[],
        metavar="VALUE[:LABEL[:COLOR]]",
        help=(
            "Add horizontal reference line to 1D sweep plots (repeatable). "
            "E.g. --hline 1.0:P=1:red --hline 0.1::orange"
        ),
    )
    plot_parser.add_argument(
        "--divergent-center",
        type=float,
        default=None,
        metavar="VALUE",
        help=(
            "Center colormap at VALUE using TwoSlopeNorm "
            "(e.g. P_EM baseline for amplification plots). "
            "Implies --cmap RdBu_r unless --cmap is set."
        ),
    )
    # Sweep-specific options
    plot_parser.add_argument(
        "--metric",
        default=None,
        metavar="NAME[,NAME,...]",
        help=(
            "Metric column(s) for sweep plots (e.g. P_max, L_mix, max_energy_error). "
            "Derived metrics A, log10_A, C0, P_EM, P_valid are auto-computed "
            "from --baseline-formula when requested."
        ),
    )
    plot_parser.add_argument(
        "--baseline-formula",
        dest="baseline_formula",
        default=None,
        metavar="EXPR",
        help=(
            "Analytical baseline formula for derived metrics (A, C0, log10_A). "
            "Default: 'sin(kappa * B0 * t_end / 2)**2' (Gertsenshtein). "
            "Uses swept/fixed param names and sim settings as variables."
        ),
    )
    plot_parser.add_argument(
        "--baseline-sweep",
        dest="baseline_sweep",
        default=None,
        metavar="DIR",
        help=(
            "Path to baseline sweep directory for paired amplification. "
            "Computes A_paired = P(signal)/P(baseline) by matching "
            "parameter values. Use with --metric A_paired."
        ),
    )
    plot_parser.add_argument(
        "--title",
        default=None,
        help="Custom figure title",
    )
    # Scatter matrix options
    plot_parser.add_argument(
        "--params",
        default=None,
        metavar="P1,P2,...",
        help="Subset of swept parameters for sweep-scatter (comma-separated)",
    )
    plot_parser.add_argument(
        "--log-diagonal",
        dest="log_diagonal",
        action="store_true",
        help="Use log y-axis on scatter matrix diagonal panels",
    )
    # 1D sweep options
    plot_parser.add_argument(
        "--scatter",
        action="store_true",
        help="Use scatter plot instead of line plot for 1D sweep",
    )
    plot_parser.add_argument(
        "--annotate-extremes",
        dest="annotate_extremes",
        action="store_true",
        help="Annotate max/min values with arrows on 1D sweep",
    )
    # Axis transform options
    plot_parser.add_argument(
        "--arctan-axes",
        dest="arctan_axes",
        action="store_true",
        help=(
            "Apply arctan transform to parameter axes (compresses ±∞ to bounded range). "
            "Tick labels show actual parameter values."
        ),
    )
    # 2D heatmap options
    plot_parser.add_argument(
        "--clamp-color",
        dest="clamp_color",
        default=None,
        metavar="MIN:MAX",
        help="Clamp colorbar range (e.g. --clamp-color=-8:5)",
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
    # Campaign-specific options
    plot_parser.add_argument(
        "--layout",
        default=None,
        metavar="RxC",
        help="Grid layout for campaign plots (e.g. 3x3). Auto-detected if omitted.",
    )
    plot_parser.add_argument(
        "--classify",
        default=None,
        metavar="RULE",
        help=(
            "Classify campaign subplots (e.g. 'active:A>1.01'). "
            "Active panels get warm background, null panels get grey."
        ),
    )

    # --- sweep ---
    sweep_parser = sub.add_parser(
        "sweep",
        aliases=["sw"],
        help="Run parameter sweeps and convergence studies",
        description=(
            "Run simulate + measure across a parameter grid and aggregate "
            "scalar metrics into portable CSV/JSON output."
        ),
        epilog=(
            "Examples:\n"
            "  tidal sweep spec.json --sweep 'g0=0.1:1.0:10' --measure conversion\n"
            "  tidal sweep spec.json --sweep 'g0=0.1:1.0:5' --sweep 'mChi2=0.5:4.0:8'\n"
            "  tidal sweep spec.json --converge '32,64,128,256' --measure conservation\n"
            "  tidal sweep spec.json --sweep 'g0=0.1,0.5,1.0' --resume --output sweep_out/"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sweep_parser.add_argument(
        "json_path",
        nargs="?",
        default=None,
        help="Path to JSON equation specification (optional if --config provides spec)",
    )
    # Sweep specification
    sweep_parser.add_argument(
        "--sweep",
        action="append",
        default=[],
        metavar="PARAM=SPEC",
        help=(
            "Parameter sweep spec (repeatable). "
            "Formats: NAME=START:STOP:N (linear), NAME=START:STOP:N:log, "
            "NAME=V1,V2,V3,..."
        ),
    )
    sweep_parser.add_argument(
        "--converge",
        default=None,
        metavar="N1,N2,...",
        help="Grid convergence study (comma-separated grid sizes). Mutually exclusive with --sweep.",
    )
    # Measurement config
    sweep_parser.add_argument(
        "--measure",
        default=None,
        metavar="TYPE[,TYPE,...]",
        help=(
            "Measurements to run per simulation (comma-separated). "
            "Options: summary, energy, conversion, mixing, dispersion, "
            "conservation, effective_mass, asymptotic, peak_conversion. "
            "Default: summary"
        ),
    )
    sweep_parser.add_argument(
        "--source",
        default=None,
        metavar="FIELD[,FIELD,...]",
        help="Source field(s) for conversion measurement (comma-separated)",
    )
    sweep_parser.add_argument(
        "--target",
        default=None,
        metavar="FIELD[,FIELD,...]",
        help="Target field(s) for conversion measurement (comma-separated)",
    )
    sweep_parser.add_argument(
        "--energy-threshold",
        type=float,
        default=1e-3,
        metavar="T",
        help="Energy conservation threshold (default: 1e-3)",
    )
    # Simulation passthrough flags (same as simulate)
    sweep_parser.add_argument(
        "--param",
        action="append",
        default=[],
        metavar="KEY=VAL",
        help="Set fixed parameter (repeatable, e.g. --param m2=1.0)",
    )
    sweep_parser.add_argument(
        "--grid-shape",
        default=None,
        metavar="N[,N,N]",
        help="Grid points per axis (default: auto)",
    )
    sweep_parser.add_argument(
        "--bounds",
        default=None,
        metavar="LO:HI[,...]",
        help="Domain bounds per axis (e.g. 0:20 or -50:50)",
    )
    sweep_parser.add_argument(
        "--periodic",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Use periodic boundary conditions (default: True)",
    )
    sweep_parser.add_argument(
        "--bc",
        default=None,
        metavar="BC[,BC,BC]",
        help="Per-axis boundary conditions (periodic, neumann). Overrides --periodic.",
    )
    sweep_parser.add_argument(
        "--ic",
        choices=["gaussian", "plane-wave", "zero", "formula", "file", "noise"],
        default="gaussian",
        help="Initial condition type (default: gaussian)",
    )
    sweep_parser.add_argument(
        "--ic-formula",
        default=None,
        metavar="EXPR",
        help="IC formula expression (for --ic=formula)",
    )
    sweep_parser.add_argument(
        "--ic-center",
        default=None,
        metavar="X[,X,X]",
        help="Gaussian center position (default: domain midpoint)",
    )
    sweep_parser.add_argument(
        "--ic-width",
        type=float,
        default=None,
        metavar="W",
        help="Gaussian width (default: domain_size/10)",
    )
    sweep_parser.add_argument(
        "--ic-amplitude",
        type=float,
        default=1.0,
        metavar="A",
        help="IC amplitude (default: 1.0)",
    )
    sweep_parser.add_argument(
        "--ic-component",
        default=None,
        metavar="NAME",
        help="Field component for IC (default: first field)",
    )
    sweep_parser.add_argument(
        "--ic-wavevector",
        default=None,
        metavar="K[,K,K]",
        help="Wavevector for plane-wave or gaussian IC",
    )
    sweep_parser.add_argument(
        "--ic-formula-velocity",
        default=None,
        metavar="EXPR",
        help="Velocity formula for --ic=formula",
    )
    sweep_parser.add_argument(
        "--ic-field",
        action="append",
        default=[],
        metavar="FIELD:EXPR",
        help="Per-field IC formula override (repeatable)",
    )
    sweep_parser.add_argument(
        "--ic-file",
        default=None,
        metavar="PATH",
        help="IC file path (for --ic=file)",
    )
    sweep_parser.add_argument(
        "--ic-noise-seed",
        type=int,
        default=None,
        metavar="N",
        help="Random seed for --ic=noise",
    )
    sweep_parser.add_argument(
        "--t-end",
        type=float,
        default=10.0,
        help="Simulation duration (default: 10.0)",
    )
    sweep_parser.add_argument(
        "--dt",
        type=float,
        default=None,
        help="Time step (default: auto from CFL)",
    )
    sweep_parser.add_argument(
        "--scheme",
        choices=["ida", "leapfrog", "cvode", "scipy", "modal", "auto"],
        default="auto",
        help="Solver scheme (default: auto)",
    )
    sweep_parser.add_argument(
        "--rtol",
        type=float,
        default=DEFAULT_RTOL,
        help=f"Relative tolerance (default: {DEFAULT_RTOL})",
    )
    sweep_parser.add_argument(
        "--atol",
        type=float,
        default=DEFAULT_ATOL,
        help=f"Absolute tolerance (default: {DEFAULT_ATOL})",
    )
    sweep_parser.add_argument(
        "--method",
        type=str,
        default=None,
        help="Integration method (default: auto)",
    )
    sweep_parser.add_argument(
        "--max-step",
        type=float,
        default=None,
        help="Maximum step size for adaptive solvers",
    )
    sweep_parser.add_argument(
        "--snapshots",
        type=float,
        default=None,
        metavar="DT",
        help="Snapshot interval (default: t_end/100)",
    )
    sweep_parser.add_argument(
        "--fd-order",
        type=int,
        choices=[2, 4, 6],
        default=4,
        help=(
            "Finite-difference accuracy order for spatial operators (default: 4). "
            "Higher orders use wider stencils: 2->3pt, 4->5pt, 6->7pt. "
            "Order 4 is recommended: O(dx^4) convergence at ~2x per-point cost "
            "vs order 2 — typically enables halving N for the same accuracy. "
            "Ref: Fornberg (1988), Mathematics of Computation 51(184)."
        ),
    )
    sweep_parser.add_argument(
        "--leapfrog-order",
        type=int,
        choices=[2, 4],
        default=None,
        help=(
            "Leapfrog integrator order (default: auto-detect). "
            "2 = Störmer-Verlet, 4 = Yoshida. "
            "Only applies when --scheme leapfrog. "
            "Ref: Yoshida (1990), Physics Letters A 150(5-7)."
        ),
    )
    sweep_parser.add_argument(
        "--spectral",
        action=argparse.BooleanOptionalAction,
        default=None,
        help=(
            "Use FFT spectral operators (default: auto-detect). "
            "Auto-enabled when all BCs are periodic and solver is not IDA. "
            "Use --no-spectral to force finite-difference. "
            "Ref: Burns et al. (2020), Phys. Rev. Research 2:023068."
        ),
    )
    sweep_parser.add_argument(
        "--mode",
        choices=["evolve", "constraint"],
        default="evolve",
        help="Simulation mode (default: evolve)",
    )
    sweep_parser.add_argument(
        "--allow-inconsistent-ic",
        action="store_true",
        default=False,
        help="Allow inconsistent ICs for constraint equations",
    )
    # Sweep execution
    sweep_parser.add_argument(
        "--output",
        default=None,
        metavar="DIR",
        help="Output directory for sweep results (required)",
    )
    sweep_parser.add_argument(
        "--resume",
        action="store_true",
        help="Skip runs with existing metadata.json (resume interrupted sweep)",
    )
    sweep_parser.add_argument(
        "--parallel",
        type=int,
        default=None,
        metavar="N",
        help=(
            "Number of parallel worker processes (default: sequential). "
            "Always set this on HPC — sequential default costs 20-50x wall time. "
            "Rule of thumb on a 112-core sapphire node: "
            "< 5 s/point -> min(N,32) (pool startup ~40%% of run time above P=32); "
            "5-30 s/point -> min(N,56); "
            "> 30 s/point -> min(N,112) (startup <1%%, use the full node). "
            "Super-linear speedup at P in {8,16} from BLAS cache effects on short runs. "
            "Run scripts/hpc_scaling_sweep.sh to profile a new workload."
        ),
    )
    sweep_parser.add_argument(
        "--force-large-sweep",
        action="store_true",
        help="Allow sweeps with more than 10,000 runs (overrides safety limit)",
    )
    sweep_parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print sweep plan (grid size, parameter combos) without running",
    )
    sweep_parser.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Overwrite existing --output directory without prompting",
    )
    sweep_parser.add_argument(
        "--keep-runs",
        action="store_true",
        help=(
            "Keep per-run snapshot directories after measurement "
            "(default: deleted to save disk). Required for "
            "`tidal plot --type sweep-compare` overlays. "
            "Warning: ~11 MB per run."
        ),
    )
    sweep_parser.add_argument(
        "--config",
        default=None,
        metavar="TOML",
        help="Load sweep configuration from TOML file (CLI flags override TOML values)",
    )
    # Adaptive sampling (F2a)
    sweep_parser.add_argument(
        "--adaptive-metric",
        default=None,
        metavar="METRIC",
        help="Metric to drive adaptive refinement (e.g. P_max)",
    )
    sweep_parser.add_argument(
        "--adaptive-budget",
        type=int,
        default=None,
        metavar="N",
        help="Maximum total points for adaptive refinement (default: 20)",
    )
    sweep_parser.add_argument(
        "--adaptive-threshold",
        type=float,
        default=None,
        metavar="T",
        help="Stop adaptive refinement when max interval score < T (default: 0.01)",
    )

    # Ensemble / replicate runs
    sweep_parser.add_argument(
        "--n-replicates",
        type=int,
        default=None,
        metavar="N",
        help="Number of replicate runs per parameter point (default: 1, no replicates)",
    )
    sweep_parser.add_argument(
        "--base-seed",
        type=int,
        default=None,
        metavar="S",
        help="Base RNG seed for ensemble; replicate i gets seed S+i (default: 42)",
    )
    sweep_parser.add_argument(
        "--ic-perturbation",
        type=float,
        default=None,
        metavar="SCALE",
        help=(
            "Additive Gaussian noise (relative to ic-amplitude) on top of "
            "deterministic ICs. Each replicate uses a different seed."
        ),
    )
    sweep_parser.add_argument(
        "--param-noise",
        type=str,
        action="append",
        default=None,
        metavar="PARAM=SCALE",
        help=(
            "Per-parameter Gaussian noise for uncertainty quantification. "
            "Each replicate draws param += N(0, SCALE). Repeatable."
        ),
    )

    # Multi-D sampling
    sweep_parser.add_argument(
        "--sweep-strategy",
        type=str,
        default=None,
        choices=["grid", "latin_hypercube", "sobol"],
        help="Sampling strategy for multi-parameter sweeps (default: grid = Cartesian product)",
    )
    sweep_parser.add_argument(
        "--n-samples",
        type=int,
        default=None,
        metavar="N",
        help="Number of samples for latin_hypercube/sobol strategies",
    )

    # Critical field extraction (post-sweep)
    sweep_parser.add_argument(
        "--critical-field",
        dest="critical_field",
        metavar="PARAM",
        help=(
            "After sweep, extract minimum field strength for full conversion. "
            "Collapses the named swept parameter by finding the first value "
            "where --metric >= --cf-threshold (default: P_final >= 0.99). "
            "Results saved to OUTPUT/critical_field/"
        ),
    )
    sweep_parser.add_argument(
        "--cf-threshold",
        dest="cf_threshold",
        type=float,
        default=0.99,
        metavar="T",
        help="Conversion threshold for --critical-field (default: 0.99)",
    )
    sweep_parser.add_argument(
        "--no-interpolate",
        dest="no_interpolate",
        action="store_true",
        help="Disable interpolation between sweep grid points for --critical-field",
    )

    # --- analyze ---
    analyze_parser = sub.add_parser(
        "analyze",
        help="Post-hoc analysis of sweep or inference results",
        epilog="""Examples:
  tidal analyze sweep_results/ --sensitivity sobol --metric P_max
  tidal analyze sweep_results/ --critical-field B0 --metric P_final --threshold 0.99
  tidal analyze ns_output/ --inference --importance
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    analyze_parser.add_argument(
        "data_path",
        help="Path to sweep output directory",
    )
    analyze_parser.add_argument(
        "--sensitivity",
        choices=["sobol", "morris"],
        default="sobol",
        help="Sensitivity analysis method (default: sobol)",
    )
    analyze_parser.add_argument(
        "--metric",
        help="Metric to analyze (e.g. P_final, P_max, max_energy_error)",
    )
    analyze_parser.add_argument(
        "--bootstrap",
        type=int,
        default=100,
        help="Bootstrap resamples for Sobol confidence intervals (default: 100)",
    )
    # Critical field analysis
    analyze_parser.add_argument(
        "--critical-field",
        dest="critical_field",
        metavar="PARAM",
        help=(
            "Extract minimum field strength for full conversion. "
            "Collapses the named swept parameter by finding the first "
            "value where --metric >= --threshold (default: P_final >= 0.99)"
        ),
    )
    analyze_parser.add_argument(
        "--threshold",
        type=float,
        default=0.99,
        help="Metric threshold for full conversion (default: 0.99)",
    )
    analyze_parser.add_argument(
        "--output",
        help="Output directory for critical field results (required with --critical-field)",
    )
    analyze_parser.add_argument(
        "--no-interpolate",
        dest="no_interpolate",
        action="store_true",
        help="Disable interpolation between sweep grid points",
    )
    analyze_parser.add_argument(
        "--reference-formula",
        dest="reference_formula",
        choices=["boccaletti", "uniform"],
        help=(
            "Compute threshold from analytical E-M formula. "
            "'boccaletti': localized Gaussian B-field, "
            "'uniform': uniform B periodic domain. "
            "Overrides --threshold."
        ),
    )
    analyze_parser.add_argument(
        "--reference-B",
        dest="reference_B",
        type=float,
        help="Reference B value for analytical formula (required with --reference-formula)",
    )
    # Inference analysis
    analyze_parser.add_argument(
        "--inference",
        action="store_true",
        default=False,
        help="Treat data as nested sampling output (loads InferenceResult)",
    )
    analyze_parser.add_argument(
        "--importance",
        action="store_true",
        default=False,
        help="Compute parameter importance via KL divergence (requires --inference)",
    )
    analyze_parser.add_argument(
        "--n-bootstrap",
        type=int,
        default=100,
        dest="n_bootstrap_importance",
        help="Bootstrap samples for importance uncertainties (default: 100)",
    )

    # --- sample (Bayesian inference) ---
    sample_parser = sub.add_parser(
        "sample",
        aliases=["samp"],
        help="Bayesian inference via Monte Carlo and nested sampling",
        description=(
            "Run Bayesian parameter estimation using simulation-based "
            "likelihood evaluation.  Supports simple Monte Carlo and "
            "nested sampling (dynesty/PolyChord) with anesthetic visualization."
        ),
        epilog=(
            "Examples:\n"
            "  tidal sample spec.json --prior 'g0=uniform:0.01:0.5' "
            "--likelihood 'P_max:maximize' --method mc --n-samples 50 --output mc_out/\n"
            "  tidal sample spec.json --prior 'alpha=uniform:0:10' "
            "--prior 'xi=log_uniform:0.01:10' --constraint 'xi > 0' "
            "--likelihood 'P_max:gaussian:0.3:0.1' --method nested --output ns_out/"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sample_parser.add_argument(
        "json_path",
        help="Path to JSON equation specification",
    )
    # Priors
    sample_parser.add_argument(
        "--prior",
        action="append",
        default=[],
        metavar="SPEC",
        help=(
            "Parameter prior (repeatable). "
            "Format: NAME=DIST:ARG1:ARG2. "
            "Distributions: uniform, log_uniform, normal, arctan_uniform. "
            "Example: --prior 'alpha=uniform:0.01:10'"
        ),
    )
    # Constraints
    sample_parser.add_argument(
        "--constraint",
        action="append",
        default=[],
        metavar="EXPR",
        help=(
            "Hard prior constraint (repeatable). "
            "Rejects parameter combinations before simulation. "
            "Example: --constraint 'xi > 0'"
        ),
    )
    # Likelihood
    sample_parser.add_argument(
        "--likelihood",
        default=None,
        metavar="SPEC",
        help=(
            "Likelihood specification: METRIC:TYPE[:ARGS]. "
            "Types: maximize, minimize, extremize, gaussian:TARGET:SIGMA, "
            "threshold:MIN_VALUE. "
            "Example: --likelihood 'P_max:maximize'"
        ),
    )
    sample_parser.add_argument(
        "--baseline-formula",
        default=None,
        metavar="FORMULA",
        dest="baseline_formula",
        help=(
            "Baseline formula for extremize likelihood. "
            "Evaluated per-point with current parameter values. "
            'Example: --baseline-formula "sin(kappa * B0 * t_end / 2)**2"'
        ),
    )
    # Sampling method
    sample_parser.add_argument(
        "--method",
        choices=["mc", "nested"],
        default="mc",
        help="Sampling method (default: mc)",
    )
    sample_parser.add_argument(
        "--sampler",
        choices=["polychord", "dynesty"],
        default="polychord",
        help="Nested sampling backend (default: polychord)",
    )
    sample_parser.add_argument(
        "--n-samples",
        type=int,
        default=100,
        metavar="N",
        help="Number of Monte Carlo samples (default: 100)",
    )
    sample_parser.add_argument(
        "--nlive",
        type=int,
        default=100,
        metavar="N",
        help="Number of live points for nested sampling (default: 100)",
    )
    sample_parser.add_argument(
        "--dlogz",
        type=float,
        default=0.01,
        help="Evidence tolerance for nested sampling (default: 0.01)",
    )
    sample_parser.add_argument(
        "--dynamic",
        action="store_true",
        default=False,
        help="Use dynamic nested sampling (dynesty only)",
    )
    sample_parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed (default: 42)",
    )
    # Measurement
    sample_parser.add_argument(
        "--measure",
        default=None,
        metavar="TYPE[,TYPE,...]",
        help="Measurements to run (auto-detected from --likelihood if omitted)",
    )
    sample_parser.add_argument(
        "--source",
        default=None,
        metavar="FIELD[,FIELD,...]",
        help="Source field(s) for conversion measurement",
    )
    sample_parser.add_argument(
        "--target",
        default=None,
        metavar="FIELD[,FIELD,...]",
        help="Target field(s) for conversion measurement",
    )
    sample_parser.add_argument(
        "--energy-threshold",
        type=float,
        default=1e-3,
        metavar="T",
        help="Energy conservation threshold (default: 1e-3)",
    )
    # Simulation passthrough (same as sweep)
    sample_parser.add_argument(
        "--param",
        action="append",
        default=[],
        metavar="KEY=VAL",
        help="Set fixed parameter (repeatable)",
    )
    sample_parser.add_argument(
        "--grid-shape",
        default=None,
        metavar="N[,N,N]",
        help="Grid points per axis",
    )
    sample_parser.add_argument(
        "--bounds",
        default=None,
        metavar="LO:HI[,...]",
        help="Domain bounds per axis",
    )
    sample_parser.add_argument(
        "--periodic",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Use periodic boundary conditions (default: True)",
    )
    sample_parser.add_argument("--bc", default=None, metavar="BC[,BC,BC]")
    sample_parser.add_argument(
        "--ic",
        choices=["gaussian", "plane-wave", "zero", "formula", "file", "noise"],
        default="gaussian",
        help="Initial condition type (default: gaussian)",
    )
    sample_parser.add_argument("--ic-formula", default=None)
    sample_parser.add_argument("--ic-center", default=None)
    sample_parser.add_argument("--ic-width", type=float, default=None)
    sample_parser.add_argument("--ic-amplitude", type=float, default=1.0)
    sample_parser.add_argument("--ic-component", default=None)
    sample_parser.add_argument("--ic-wavevector", default=None)
    sample_parser.add_argument("--ic-formula-velocity", default=None)
    sample_parser.add_argument("--ic-field", action="append", default=[])
    sample_parser.add_argument("--ic-file", default=None)
    sample_parser.add_argument("--ic-noise-seed", type=int, default=None)
    sample_parser.add_argument("--t-end", type=float, default=10.0)
    sample_parser.add_argument("--dt", type=float, default=None)
    sample_parser.add_argument(
        "--scheme",
        choices=["ida", "leapfrog", "cvode", "scipy", "modal", "auto"],
        default="auto",
    )
    sample_parser.add_argument("--rtol", type=float, default=DEFAULT_RTOL)
    sample_parser.add_argument("--atol", type=float, default=DEFAULT_ATOL)
    sample_parser.add_argument(
        "--method-solver", type=str, default=None, dest="method_solver"
    )
    sample_parser.add_argument("--max-step", type=float, default=None)
    sample_parser.add_argument("--snapshots", type=float, default=None)
    sample_parser.add_argument("--fd-order", type=int, choices=[2, 4, 6], default=4)
    sample_parser.add_argument(
        "--leapfrog-order", type=int, choices=[2, 4], default=None
    )
    sample_parser.add_argument(
        "--spectral", action=argparse.BooleanOptionalAction, default=None
    )
    sample_parser.add_argument(
        "--mode", choices=["evolve", "constraint"], default="evolve"
    )
    sample_parser.add_argument(
        "--allow-inconsistent-ic", action="store_true", default=False
    )
    # Execution
    sample_parser.add_argument(
        "--output",
        default=None,
        metavar="DIR",
        help="Output directory for inference results (required)",
    )
    sample_parser.add_argument(
        "--parallel",
        type=int,
        default=None,
        metavar="N",
        help="Number of parallel workers",
    )
    # Plots
    sample_parser.add_argument(
        "--corner",
        action="store_true",
        default=False,
        help="Generate corner plot",
    )
    sample_parser.add_argument(
        "--trace",
        action="store_true",
        default=False,
        help="Generate trace plot",
    )
    sample_parser.add_argument(
        "--importance",
        action="store_true",
        default=False,
        help="Generate parameter importance bar chart",
    )
    # Analysis
    sample_parser.add_argument(
        "--analyze",
        action="store_true",
        default=False,
        help="Run parameter importance analysis after sampling (nested only)",
    )
    sample_parser.add_argument(
        "--nlive-auto",
        nargs="?",
        const="standard",
        default=None,
        choices=["fast", "standard", "production"],
        metavar="PRECISION",
        help=(
            "Auto-scale nlive based on number of parameters. "
            "Overrides --nlive. Levels: fast, standard (default), production"
        ),
    )
    # Resume fields needed by _build_sim_args
    sample_parser.add_argument("--resume", default=None)
    sample_parser.add_argument("--snapshot", default=None)
    sample_parser.add_argument("--t-additional", default=None)

    # --- doctor ---
    sub.add_parser(
        "doctor",
        aliases=["doc"],
        help="Check environment health (Python, dependencies, Wolfram, xAct)",
        epilog="""Examples:
  tidal doctor              # Full environment check
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Shell completion (optional dependency: shtab)
    try:
        import shtab

        shtab.add_argument_to(parser)
    except ImportError:
        pass

    return parser


# ---------------------------------------------------------------------------
# Dispatch and entry point
# ---------------------------------------------------------------------------


def _get_version() -> str:
    """Return the package version string from installed metadata."""
    from importlib.metadata import PackageNotFoundError, version

    try:
        return version("tidal")
    except PackageNotFoundError:
        return "unknown"


_COMMAND_ALIASES: dict[str, str] = {
    "der": "derive",
    "insp": "inspect",
    "sim": "simulate",
    "ls": "list",
    "val": "validate",
    "meas": "measure",
    "sw": "sweep",
    "samp": "sample",
    "doc": "doctor",
}


def _dispatch(args: argparse.Namespace) -> int:  # noqa: PLR0911, C901
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
    cmd = _COMMAND_ALIASES.get(args.command, args.command)
    if cmd == "derive":
        from tidal.cli._derive import derive_command

        return derive_command(args)
    if cmd == "inspect":
        from tidal.cli._inspect import inspect_command

        return inspect_command(args)
    if cmd == "simulate":
        from tidal.cli._simulate import simulate_command

        return simulate_command(args)
    if cmd == "list":
        from tidal.cli._list import list_command

        return list_command(args)
    if cmd == "validate":
        from tidal.cli._validate import validate_command

        return validate_command(args)
    if cmd == "measure":
        from tidal.cli._measure import measure_command

        return measure_command(args)
    if cmd == "plot":
        from tidal.cli._plot_command import plot_command

        return plot_command(args)
    if cmd == "sweep":
        from tidal.cli._sweep import sweep_command

        return sweep_command(args)
    if cmd == "analyze":
        from tidal.cli._analyze import analyze_command

        return analyze_command(args)
    if cmd == "sample":
        from tidal.cli._sample import sample_command

        return sample_command(args)
    if cmd == "doctor":
        from tidal.cli._doctor import doctor_command

        return doctor_command(args)
    msg = f"Unknown command: {cmd}"
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

    # Show banner for long-running commands and bare `tidal` (no subcommand).
    # Quick commands (list, validate, inspect, measure, plot, analyze) skip it.
    # Skip for --cite (info-only flag).
    if (
        not args.no_banner
        and not os.environ.get("TIDAL_NO_BANNER")
        and not args.cite
        and _COMMAND_ALIASES.get(args.command, args.command)
        in {None, "simulate", "derive", "sweep"}
    ):
        from tidal.banner import print_banner

        banner_theme = os.environ.get("TIDAL_BANNER_THEME", "ocean")
        banner_layout = os.environ.get("TIDAL_BANNER_LAYOUT", "auto")
        print_banner(
            theme=banner_theme,
            layout=banner_layout,
            version=_get_version(),
        )

    if args.cite:
        from tidal.cli._cite import print_citation

        print_citation()
        return 0

    if args.command is None:
        parser.print_help()
        return 0

    from tidal.cli._console import configure as _configure_console
    from tidal.cli._console import error as _console_error

    _configure_console(
        verbose=getattr(args, "verbose", False),
        quiet=getattr(args, "quiet", False),
    )

    import time

    t0 = time.perf_counter()
    try:
        rc = _dispatch(args)
    except KeyboardInterrupt:
        print("\nInterrupted.", file=sys.stderr)
        rc = 130
    except FileNotFoundError as exc:
        from tidal.cli._console import error_with_hint

        error_with_hint(
            str(exc),
            [
                "Run 'tidal list' to see available specifications.",
                "Run 'tidal derive <theory.toml>' to generate from a Lagrangian.",
            ],
        )
        rc = 1
    except (ValueError, TypeError) as exc:
        _console_error(str(exc))
        rc = 1
    except RuntimeError:
        import traceback

        traceback.print_exc()
        rc = 1

    if args.time:
        elapsed = time.perf_counter() - t0
        if elapsed >= 60:  # noqa: PLR2004
            mins = int(elapsed // 60)
            secs = elapsed % 60
            print(
                f"[TIME] {args.command} completed in {mins}m {secs:.1f}s",
                file=sys.stderr,
            )
        else:
            print(f"[TIME] {args.command} completed in {elapsed:.2f}s", file=sys.stderr)

    return rc
