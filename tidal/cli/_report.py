"""``tidal simulate --report`` — Self-contained HTML summary report.

Generates a single HTML file with embedded plots (base64 PNG), parameter
tables, solver info, energy conservation check, and a reproducible command.
No external dependencies beyond matplotlib (already required).
"""

from __future__ import annotations

import base64
import io
import sys
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from tidal.measurement._io import SimulationData
    from tidal.solver.grid import GridInfo
    from tidal.symbolic.json_loader import EquationSystem


def _fig_to_base64() -> str:
    """Render the current matplotlib figure as a base64-encoded PNG."""
    import matplotlib.pyplot as plt

    buf = io.BytesIO()
    plt.savefig(buf, format="png", dpi=120, bbox_inches="tight")
    plt.close()
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("ascii")


def _make_overview_plot(sim_data: SimulationData, grid_info: GridInfo) -> str:
    """Generate the overview heatmap/line plot and return as base64 PNG."""
    import matplotlib.pyplot as plt

    n_fields = len(sim_data.fields)
    fig, axes = plt.subplots(
        1, max(n_fields, 1), figsize=(5 * n_fields, 4), squeeze=False,
    )

    for i, (name, data) in enumerate(sim_data.fields.items()):
        ax = axes[0][i]
        if data.ndim == 2:  # noqa: PLR2004
            # 1+1D: heatmap (time x space)
            extent = [
                float(grid_info.bounds[0][0]),
                float(grid_info.bounds[0][1]),
                float(sim_data.times[0]),
                float(sim_data.times[-1]),
            ]
            ax.imshow(
                data,
                aspect="auto",
                origin="lower",
                extent=extent,
                cmap="RdBu_r",
            )
            ax.set_xlabel("x")
            ax.set_ylabel("t")
        else:
            # Multi-D: just plot time evolution of max amplitude
            peaks = [float(np.max(np.abs(data[t]))) for t in range(len(sim_data.times))]
            ax.plot(sim_data.times, peaks)
            ax.set_xlabel("t")
            ax.set_ylabel("|peak|")
        ax.set_title(name)

    fig.tight_layout()
    return _fig_to_base64()


def _energy_conservation(
    sim_data: SimulationData,
) -> dict[str, Any]:
    """Compute energy conservation metric if possible."""
    try:
        from tidal.measurement._diagnostics import check_energy_conservation

        result = check_energy_conservation(sim_data)
    except Exception as exc:  # noqa: BLE001
        return {"error": str(exc)}
    else:
        return {
            "is_conserved": result.is_conserved,
            "max_relative_error": result.max_relative_error,
            "threshold": 1e-3,
        }


def _escape(text: str) -> str:
    """HTML-escape a string."""
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def generate_report(  # noqa: PLR0913, PLR0914, PLR0917
    sim_data: SimulationData,
    spec: EquationSystem,
    params: dict[str, float],
    grid_info: GridInfo,
    scheme: str,
    report_path: str,
    argv: list[str] | None = None,
) -> None:
    """Generate a self-contained HTML report for a simulation run.

    Parameters
    ----------
    sim_data : SimulationData
        Completed simulation data.
    spec : EquationSystem
        Equation specification.
    params : dict
        Resolved parameter values.
    grid_info : GridInfo
        Grid configuration.
    scheme : str
        Solver scheme used.
    report_path : str
        Output HTML file path.
    argv : list[str] | None
        Original command-line arguments for reproducibility.
    """
    from importlib.metadata import PackageNotFoundError, version

    try:
        tidal_version = version("tidal")
    except PackageNotFoundError:
        tidal_version = "unknown"

    # Overview plot
    overview_b64 = _make_overview_plot(sim_data, grid_info)

    # Energy conservation
    cons = _energy_conservation(sim_data)

    # Build parameter table rows
    param_rows = "".join(
        f"<tr><td><code>{_escape(k)}</code></td><td>{v}</td></tr>\n"
        for k, v in sorted(params.items())
    )
    if not param_rows:
        param_rows = "<tr><td colspan='2'>No parameters</td></tr>"

    # Field summary rows
    field_rows = ""
    for name in spec.component_names:
        if name not in sim_data.fields:
            continue
        init_peak = float(np.max(np.abs(sim_data.fields[name][0])))
        final_peak = float(np.max(np.abs(sim_data.fields[name][-1])))
        ratio = final_peak / init_peak if init_peak > 0 else 0.0
        field_rows += (
            f"<tr><td><code>{_escape(name)}</code></td>"
            f"<td>{init_peak:.4f}</td><td>{final_peak:.4f}</td>"
            f"<td>{ratio:.4f}</td></tr>\n"
        )

    # Energy conservation section
    if "error" in cons:
        cons_html = f"<p class='warn'>Energy check failed: {_escape(cons['error'])}</p>"
    else:
        status = "PASS" if cons["is_conserved"] else "FAIL"
        status_class = "pass" if cons["is_conserved"] else "fail"
        cons_html = (
            f"<p class='{status_class}'><strong>{status}</strong> — "
            f"max |dE/E| = {cons['max_relative_error']:.2e} "
            f"(threshold: {cons['threshold']:.0e})</p>"
        )

    # Reproducible command
    cmd = "tidal " + " ".join(argv) if argv else " ".join(sys.argv)

    # Grid info
    grid_shape = "x".join(str(s) for s in grid_info.shape)
    bc_str = ", ".join("periodic" if p else "non-periodic" for p in grid_info.periodic)
    bounds_str = ", ".join(f"[{b[0]}, {b[1]}]" for b in grid_info.bounds)

    html = f"""\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>TIDAL Simulation Report</title>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
         max-width: 900px; margin: 2em auto; padding: 0 1em; color: #333; }}
  h1 {{ color: #0066cc; border-bottom: 2px solid #0066cc; padding-bottom: 0.3em; }}
  h2 {{ color: #444; margin-top: 1.5em; }}
  table {{ border-collapse: collapse; width: 100%; margin: 1em 0; }}
  th, td {{ border: 1px solid #ddd; padding: 8px 12px; text-align: left; }}
  th {{ background: #f5f5f5; font-weight: 600; }}
  tr:nth-child(even) {{ background: #fafafa; }}
  code {{ background: #f0f0f0; padding: 2px 6px; border-radius: 3px; font-size: 0.9em; }}
  pre {{ background: #1e1e1e; color: #d4d4d4; padding: 1em; border-radius: 6px;
         overflow-x: auto; font-size: 0.85em; }}
  .pass {{ color: #2e7d32; font-weight: bold; }}
  .fail {{ color: #c62828; font-weight: bold; }}
  .warn {{ color: #e65100; }}
  .plot {{ text-align: center; margin: 1.5em 0; }}
  .plot img {{ max-width: 100%; border: 1px solid #eee; border-radius: 4px; }}
  footer {{ margin-top: 3em; padding-top: 1em; border-top: 1px solid #eee;
            color: #999; font-size: 0.85em; }}
</style>
</head>
<body>

<h1>TIDAL Simulation Report</h1>

<h2>Specification</h2>
<table>
  <tr><th>Lagrangian</th><td>{_escape(spec.metadata.get("lagrangian_expr", "N/A"))}</td></tr>
  <tr><th>Dimension</th><td>{spec.dimension}D ({spec.spatial_dimension}+1D)</td></tr>
  <tr><th>Fields</th><td>{spec.n_components} ({", ".join(spec.component_names)})</td></tr>
  <tr><th>Equations</th><td>{len(spec.equations)}</td></tr>
</table>

<h2>Parameters</h2>
<table>
  <tr><th>Parameter</th><th>Value</th></tr>
  {param_rows}
</table>

<h2>Grid & Solver</h2>
<table>
  <tr><th>Grid</th><td>{grid_shape} points</td></tr>
  <tr><th>Bounds</th><td>{bounds_str}</td></tr>
  <tr><th>BCs</th><td>{bc_str}</td></tr>
  <tr><th>Solver</th><td>{_escape(scheme)}</td></tr>
  <tr><th>Time range</th><td>{float(sim_data.times[0]):.2f} → {float(sim_data.times[-1]):.2f} ({sim_data.n_snapshots} snapshots)</td></tr>
</table>

<h2>Field Summary</h2>
<table>
  <tr><th>Field</th><th>Initial peak</th><th>Final peak</th><th>Ratio</th></tr>
  {field_rows}
</table>

<h2>Energy Conservation</h2>
{cons_html}

<h2>Overview</h2>
<div class="plot">
  <img src="data:image/png;base64,{overview_b64}" alt="Simulation overview">
</div>

<h2>Reproduce</h2>
<pre>{_escape(cmd)}</pre>

<footer>Generated by TIDAL v{_escape(tidal_version)}</footer>

</body>
</html>
"""

    Path(report_path).write_text(html, encoding="utf-8")
