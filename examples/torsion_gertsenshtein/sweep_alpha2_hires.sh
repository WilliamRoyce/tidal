#!/usr/bin/env bash
# High-resolution alpha2 sweep at delta1=1.0 (amplification boundary zoom)
#
# Goal: Resolve the amplification profile A(alpha2) near the instability boundary.
# The boundary at delta1=1.0 is at alpha2 ~ -0.76 (upper, tachyonic photon)
# and alpha2 ~ -1.14 (lower, Block A denominator pole).
#
# Sweeps alpha2 from -1.2 to -0.5 with 200 points for fine resolution.
# B0=0.0001 to stay in linear regime.

set -euo pipefail

SPEC="examples/data/torsion_gertsenshtein_nonminimal.json"
OUTPUT="examples/data/nonminimal_alpha2_hires"

uv run tidal sweep "$SPEC" \
  --sweep "alpha2=-1.2:-0.5:200" \
  --measure peak_conversion --source h_5 --target a_1 \
  --grid-shape 256 --bounds "0:100" --periodic --t-end 50 \
  --param delta1=1.0 --param alpha1=0 --param alpha3=1.0 --param kappa=1 --param B0=0.0001 \
  --output "$OUTPUT" --parallel 4 --resume

echo ""
echo "Postprocessing..."
uv run python examples/torsion_gertsenshtein/postprocess_sweep.py "$OUTPUT/results.csv" 0.0001
