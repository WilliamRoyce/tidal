"""Sparse-CSC vs dense ``expm_multiply`` wall-time on localised Gertsenshtein.

Serves:   manuscript/sections/appendices/numerical.tex:490 (tab:SparseCSCBenchmark)
Consumes: scripts/tables/tab_sparse_csc.py
Writes:   benchmark_results/canonical/sparse_csc_vs_dense.json by default

Status: stub. Reproduces the "factor of 23–27x speedup" headline by
constructing the mode-coupled position-dependent evolution matrix for
``gertsenshtein_localised`` and timing ``scipy.sparse.linalg.expm_multiply``
in dense vs CSC formats. Implementation requirements:

    1. Load ``examples/data/gertsenshtein_localised.json`` (or whichever
       JSON exposes the localised B-field background that drives the
       sparse vs dense gap).
    2. For each ``N`` in {128, 256, 512}: build the mode-coupled evolution
       matrix (see ``tidal.solver.modal._build_convolution_matrix``).
    3. Time ``expm_multiply(A_dense, y0)`` and ``expm_multiply(A_csc, y0)``.
    4. Record density (fraction of non-zero entries) per N.

Expected output shape:

    {
      "metadata": {...},
      "results": [
        {"n": 128, "density_pct": 0.78, "dense_s": ..., "sparse_s": ...,
         "speedup": ...},
        ...
      ]
    }
"""

from __future__ import annotations

import sys


def main() -> None:
    msg = (
        "sparse_csc_vs_dense: not yet implemented. See module docstring for "
        "the design and the missing modal-matrix construction glue."
    )
    print(msg, file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
    main()
