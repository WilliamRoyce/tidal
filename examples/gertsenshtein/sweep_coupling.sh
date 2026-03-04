#!/usr/bin/env bash
# Gertsenshtein Coupling Sweep — B0 parameter scan
#
# Sweeps B0 (background magnetic field strength) and measures conversion
# probability at each value. Compares against the analytical formula:
#   P(graviton -> photon) = sin^2(kappa * B0 * D / (4*sqrt(pi)))
#
# Ref: Palessandro & Rothman (2023), arXiv:2301.02072, Eq. 26
# Gauge: TT (graviton) + Lorenz (photon), as in the reference above
#
# Running:
#   cd examples/gertsenshtein && bash sweep_coupling.sh

set -euo pipefail
cd "$(dirname "$0")"

echo "=== Gertsenshtein B0 Coupling Sweep ==="
echo ""

# Sweep B0 from weak to strong coupling
# In natural units (kappa=1), P = sin^2(B0 * D / (4*sqrt(pi)))
# With D ~ 50 (half domain for round-trip), B0 from 0.01 to 1.0
tidal sweep ../data/gertsenshtein.json \
  --sweep "B0=0.01:0.5:10" \
  --measure conversion \
  --grid-shape 256 \
  --bounds 0:100 \
  --periodic \
  --ic gaussian \
  --ic-amplitude 0.1 \
  --ic-width 3.0 \
  --ic-center 20.0 \
  --t-end 40.0 \
  --param kappa=1.0 \
  --output ../data/gertsenshtein_sweep_B0

echo ""
echo "=== Sweep complete ==="
echo "Results in: examples/data/gertsenshtein_sweep_B0/"
