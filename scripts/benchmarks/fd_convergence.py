"""FD-stencil and spectral convergence on the smooth gradient test.

Measures L_2 error of the TIDAL gradient operator (Fornberg stencils,
orders 2/4/6, plus FFT-based pseudo-spectral) against the exact derivative
of f(x) = sin(x) on a periodic grid over [0, 2pi], sweeping N.

Serves:  manuscript/sections/appendices/numerical.tex:236 (Figure FDConvergence)
Consumes: scripts/figures/figC2_fd_convergence.py
Writes:  benchmark_results/canonical/fd_convergence.json by default
"""

from __future__ import annotations

import argparse
import datetime
import json
import platform
import socket
import subprocess  # noqa: S404
import sys
from pathlib import Path

import numpy as np
import scipy

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from tidal.solver.grid import GridInfo  # noqa: E402
from tidal.solver.operators import (  # noqa: E402
    get_fd_order,
    get_spectral,
    gradient,
    set_fd_order,
    set_spectral,
)

# Sweep parameters — N values cover the regime where round-off floor is
# reached for the spectral path (~64) and where high-order FD curves
# have asymptoted to their stencil convergence rate (>= 32).
N_VALUES = [16, 32, 64, 128, 256, 512]
FD_ORDERS = [2, 4, 6]


def _l2_error(numeric: np.ndarray, exact: np.ndarray) -> float:
    return float(np.sqrt(np.mean((numeric - exact) ** 2)))


def _git_sha() -> str:
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=REPO_ROOT,
            stderr=subprocess.DEVNULL,
        )
        return out.decode().strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown"


def _metadata() -> dict:
    return {
        "timestamp": datetime.datetime.now(datetime.UTC).isoformat(),
        "host": socket.gethostname(),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "git_sha": _git_sha(),
        "parameters": {
            "n_values": N_VALUES,
            "fd_orders": FD_ORDERS,
            "domain": [0.0, 2.0 * float(np.pi)],
            "test_function": "sin(x)",
            "exact_derivative": "cos(x)",
        },
    }


def _measure(n: int, order: int | None) -> dict:
    """One (N, scheme) measurement. `order=None` selects spectral."""
    bounds = ((0.0, 2.0 * np.pi),)
    grid = GridInfo(bounds=bounds, shape=(n,), periodic=(True,))
    x = grid.axes_coords(0)
    f = np.sin(x)
    exact = np.cos(x)

    prev_order = get_fd_order()
    prev_spectral = get_spectral()
    try:
        if order is None:
            set_spectral(True)
            scheme = "spectral"
        else:
            set_spectral(False)
            set_fd_order(order)
            scheme = f"fd_o{order}"

        numeric = gradient(f, axis=0, grid=grid)
        err = _l2_error(numeric, exact)
    finally:
        set_spectral(prev_spectral)
        set_fd_order(prev_order)

    return {"n": n, "scheme": scheme, "l2_error": err}


def run() -> dict:
    rows: list[dict] = []
    for n in N_VALUES:
        rows.extend(_measure(n, order) for order in FD_ORDERS)
        rows.append(_measure(n, None))
    return {"metadata": _metadata(), "results": rows}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out",
        type=Path,
        default=REPO_ROOT / "benchmark_results" / "canonical" / "fd_convergence.json",
        help="output JSON path",
    )
    args = parser.parse_args()

    data = run()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        json.dump(data, fh, indent=2)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
