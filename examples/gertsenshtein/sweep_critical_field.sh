#!/usr/bin/env bash
# Gertsenshtein Critical Field Extraction
#
# 1D B0 sweep with --critical-field to extract B_min: minimum field strength
# for full conversion (P_final >= 0.99).
#
# For uniform B in periodic domain with kappa=1, t_end=50:
#   P_final = sin^2(kappa * B0 * t_end / 2) = sin^2(25 * B0)
#   Full conversion when 25 * B0 = pi/2 (first maximum)
#   B_min = pi / 50 ~ 0.0628
#
# The --critical-field flag runs the extraction automatically after the sweep,
# saving reduced results to OUTPUT/critical_field/.
#
# When torsion theories gain a background B field, the same workflow extends
# to multi-parameter sweeps:
#   tidal sweep torsion.json \
#     --sweep "B0=0.005:0.3:25" \
#     --sweep "alpha1=0.01:1.0:10:log" \
#     --critical-field B0 --cf-threshold 0.99
# producing a heatmap of 1/B_min(alpha1).
#
# IC: monochromatic plane wave cos(k*z) on h_7 (graviton h_+)
# k = 2*pi*32/100 ~ 2.011 (32 wavelengths in domain [0,100])
#
# Ref: Gertsenshtein (1962), JETP 14:84
#
# Running:
#   cd examples/gertsenshtein && bash sweep_critical_field.sh

set -euo pipefail
cd "$(dirname "$0")"

echo "=== Gertsenshtein Critical Field Extraction ==="
echo ""

tidal sweep ../data/gertsenshtein.json \
  --sweep "B0=0.005:0.15:30" \
  --measure conversion \
  --source h_7 --target a_2 \
  --grid-shape 512 \
  --bounds 0:100 \
  --periodic \
  --ic plane-wave \
  --ic-wavevector 2.0106 \
  --ic-amplitude 0.1 \
  --ic-component h_7 \
  --t-end 50.0 \
  --param kappa=1.0 \
  --critical-field B0 --cf-threshold 0.99 \
  --output ../data/gertsenshtein_sweep_critical

echo ""
echo "=== Sweep complete ==="
echo ""
echo "Results:        examples/data/gertsenshtein_sweep_critical/"
echo "Critical field: examples/data/gertsenshtein_sweep_critical/critical_field/"
echo ""
echo "Expected B_min ~ pi/50 ~ 0.0628 (analytical, uniform B, kappa=1)"
