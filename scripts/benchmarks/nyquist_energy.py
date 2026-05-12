"""Relative energy drift |dE/E| across grid resolutions, pre/post-fix.

Serves:   manuscript/sections/appendices/numerical.tex:525 (tab:NyquistEnergy)
Consumes: scripts/tables/tab_nyquist_energy.py
Writes:   benchmark_results/canonical/nyquist_energy.json by default

Status: stub. Sweeps the per-snapshot relative energy drift at
N ∈ {128, 256, 512, 1024} for two theories — ``coupled_scalars`` and
``gertsenshtein`` — under the modal solver. The headline post-fix numbers
are ~2e-14 (coupled_scalars all N) and ~1.3e-14 (gertsenshtein), confirming
that the volume-element / Parseval-spectral fix lands |dE/E| at the
floating-point floor across resolutions.

The pre-fix comparison row (~1e-5) cannot be regenerated without checking
out the broken commit; pragmatically the pre-fix value is loaded from an
archived measurement at
``benchmark_results/canonical/nyquist_energy_prefix.json`` (committed once
and treated as illustrative).

Implementation requirements:

    1. Load each theory, build IC (Gaussian, plane-wave — match what the
       existing one-off measurements used).
    2. Run ``tidal simulate`` (or the modal solver directly) at each N.
    3. Call ``tidal measure --what conservation`` and pull |dE/E|.
    4. Aggregate.

Expected output shape:

    {
      "metadata": {...},
      "results": [
        {"theory": "coupled_scalars", "n": 128, "abs_dE_over_E": 2.1e-14},
        ...
      ],
      "prefix_archive": {"comment": "...", "value": 1.0e-5}
    }
"""

from __future__ import annotations

import sys


def main() -> None:
    msg = (
        "nyquist_energy: not yet implemented. See module docstring for the "
        "design and the simulate/measure pipeline glue needed."
    )
    print(msg, file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
    main()
