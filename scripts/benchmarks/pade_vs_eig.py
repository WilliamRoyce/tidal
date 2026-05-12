"""Padé scaling-and-squaring vs eigendecomposition wall-time across theories.

Serves:   manuscript/sections/appendices/numerical.tex:348 (tab:PadeBenchmark)
Consumes: scripts/tables/tab_pade_benchmark.py
Writes:   benchmark_results/canonical/pade_vs_eig.json by default

Status: stub. The implementation needs to drive the modal solver per
theory in a way that captures the **total wall-time across all per-mode
matrix-exponential evaluations**, separately for the Padé path
(``scipy.linalg.expm``, production default) and the eigendecomposition
reference path (``V @ diag(exp(D·t)) @ V_inv``, retired but reconstructed
for the benchmark). The existing 3-row table came from an uncommitted
optimisation study; reproducing it faithfully requires:

    1. Loading each theory via ``tidal.symbolic.json_loader.load_equation_system``.
    2. Choosing a representative ``(grid, parameters, t_end)`` triple per
       theory matching the appendix's column (N=256 for coupled_scalars,
       N=512 for gertsenshtein, 16**3 for coupled_proca_3d).
    3. Building per-mode matrices via the modal solver's internal
       ``_build_per_mode_matrices`` / ``_build_evolution_matrices``
       (constraint-bearing theories need the latter).
    4. Timing both paths over the full per-mode loop, with repeats to
       reduce noise.

For 5–8 rows the extra theories from the TODO are
``gertsenshtein_localised``, ``torsion_gertsenshtein``, ``plasma_dispersion``.

Each follow-up session implementing this should produce a JSON of the form

    {
      "metadata": {"timestamp": ..., "host": ..., "git_sha": ..., ...},
      "results": [
        {"theory": "coupled_scalars", "n": 256, "block": "1d",
         "pade_s": 0.0030, "eig_s": 0.0028, "ratio": 1.07, "n_modes": 129, ...},
        ...
      ]
    }

so ``scripts/tables/tab_pade_benchmark.py`` can render rows directly.
"""

from __future__ import annotations

import sys


def main() -> None:
    msg = (
        "pade_vs_eig: not yet implemented. See module docstring for the "
        "design and the missing per-theory pipeline glue."
    )
    print(msg, file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
    main()
