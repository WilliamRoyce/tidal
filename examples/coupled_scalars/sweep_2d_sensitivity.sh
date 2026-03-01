#!/usr/bin/env bash
# 2D Sensitivity Analysis — Coupling × Mass Ratio
#
# Runs a Latin Hypercube sweep over gCpl × mChi2 and performs Sobol
# sensitivity analysis to determine which parameter controls conversion.
# Generates advanced visualizations: parallel coordinates, tornado,
# and scatter matrix.
#
# Demonstrates features F2b (LHS), F4 (Sobol), F5 (query), F8 (viz).
#
# Running this script:
#   cd examples/coupled_scalars && bash sweep_2d_sensitivity.sh

set -euo pipefail
cd "$(dirname "$0")"

OUT=../data/coupled_scalars_sensitivity

echo "=== 2D Latin Hypercube Sweep: gCpl × mChi2 (30 samples) ==="

tidal sweep --config sweep_2d_sensitivity.toml

echo ""
echo "=== Sensitivity Analysis ==="

# Sobol sensitivity (requires SALib: pip install SALib)
if python -c "import SALib" 2>/dev/null; then
  tidal analyze "$OUT" --sensitivity sobol --metric P_max
else
  echo "  Skipping Sobol (SALib not installed). Running Morris instead:"
  tidal analyze "$OUT" --sensitivity morris --metric P_max
fi

echo ""
echo "=== Generating Visualizations ==="

tidal plot "$OUT" --type sweep --metric P_max \
  --output "$OUT/sweep_P_max.png" --quiet

tidal plot "$OUT" --type sweep-parallel --metric P_max \
  --output "$OUT/parallel_P_max.png" --quiet

tidal plot "$OUT" --type sweep-tornado --metric P_max \
  --output "$OUT/tornado_P_max.png" --quiet

tidal plot "$OUT" --type sweep-scatter --metric P_max \
  --output "$OUT/scatter_P_max.png" --quiet

echo ""
echo "Plots saved to: $OUT/"
echo "  sweep_P_max.png      — 2D heatmap"
echo "  parallel_P_max.png   — parallel coordinates"
echo "  tornado_P_max.png    — parameter importance"
echo "  scatter_P_max.png    — pairwise scatter matrix"
