"""Fourier modal solver convergence on Einstein–Maxwell vs Boccaletti formula.

Serves:   manuscript/sections/appendices/numerical.tex:547 (figC_solver_convergence)
Consumes: scripts/figures/figC1_solver_convergence.py
Writes:   benchmark_results/canonical/solver_convergence.json by default

Status: stub. Two-panel benchmark:

    Panel (a) — relative L_2 error of the Fourier modal solver against the
                analytic Boccaletti P_conv(t) = sin^2(kappa B_0 D / 2) for
                Einstein–Maxwell on a uniform B_0 background, as a function
                of grid resolution N in {16, 32, 64, 128, ..., 2048}.
                Expected: machine-precision plateau within a few N (modal
                solver is exact for constant-coefficient periodic systems).

    Panel (b) — modal vs CVODE vs IDA at fixed N over an rtol/atol sweep.
                Demonstrates that the modal solver's error floor is lower
                than the adaptive ODE backends' achievable accuracy.

Implementation requirements:

    1. Use ``examples/einstein_maxwell_baseline/`` as the regression
       fixture; its theory.toml is already vetted and produces the
       Boccaletti result.
    2. ``tidal simulate --scheme modal`` and ``--scheme cvode`` /
       ``--scheme ida`` with the same IC.
    3. Compare snapshot-by-snapshot against the analytic formula.

Expected output shape:

    {
      "metadata": {...},
      "panel_a": [{"n": 16, "l2_error": ...}, ...],
      "panel_b": [{"scheme": "modal", "rtol": 1e-6, "l2_error": ...}, ...]
    }
"""

from __future__ import annotations

import sys


def main() -> None:
    msg = (
        "solver_convergence: not yet implemented. See module docstring for "
        "the two-panel design and the modal/CVODE/IDA driver glue."
    )
    print(msg, file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
    main()
