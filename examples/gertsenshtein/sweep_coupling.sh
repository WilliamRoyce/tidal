#!/usr/bin/env bash
# Gertsenshtein Coupling Sweep — B0 parameter scan (quick version)
#
# Quick sweep with fewer points and shorter integration time.
# For the full validation sweep, use sweep_B0.sh instead.
#
# Sweeps B0 and measures conversion probability P(graviton -> photon).
# Uses plane-wave IC for direct comparison with analytical formula:
#   P = sin^2(kappa * B0 * D / 2)
#
# Ref: Gertsenshtein (1962), JETP 14:84
#
# Running:
#   cd examples/gertsenshtein && bash sweep_coupling.sh

set -euo pipefail
cd "$(dirname "$0")"

echo "=== Gertsenshtein B0 Coupling Sweep (quick) ==="
echo ""

# Quick sweep: 10 points, t_end=100 (D=100, ~1 oscillation cycle)
tidal sweep ../data/gertsenshtein.json \
  --sweep "B0=0.01:0.25:10" \
  --measure conversion \
  --source h_7 --target a_2 \
  --grid-shape 512 \
  --bounds 0:100 \
  --periodic \
  --ic plane-wave \
  --ic-wavevector 2.0106 \
  --ic-amplitude 0.1 \
  --ic-component h_7 \
  --t-end 100.0 \
  --param kappa=1.0 \
  --output ../data/gertsenshtein_sweep_quick

echo ""
echo "=== Sweep complete ==="
echo "Results in: examples/data/gertsenshtein_sweep_quick/"
