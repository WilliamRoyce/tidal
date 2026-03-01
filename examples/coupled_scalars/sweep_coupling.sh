#!/usr/bin/env bash
# Coupled Scalars — Sweep coupling strength gCpl
#
# Sweeps the coupling constant from weak to strong coupling and measures
# conversion probability P(t) = E_chi(t) / E_phi(0) at each point.
#
# Produces:
#   sweep_coupling_output/results.csv   — tabular results
#   sweep_coupling_output/sweep.json    — provenance metadata
#   sweep_coupling_output/results.json  — JSON results
#
# Usage:
#   cd examples/coupled_scalars && bash sweep_coupling.sh
#
# Requires: pre-derived JSON spec at examples/data/coupled_scalars.json
#   (run `tidal derive theory.toml` first if missing)

set -euo pipefail
cd "$(dirname "$0")"

OUTPUT=../data/coupled_scalars_sweep

# 1D sweep: coupling strength gCpl from 0.01 to 1.0 (8 points, log scale)
tidal sweep ../data/coupled_scalars.json \
  --sweep "gCpl=0.01:1.0:8:log" \
  --param mPhi2=1.0 --param mChi2=4.0 \
  --measure conversion,conservation \
  --source phi_0 --target chi_0 \
  --grid-shape 64 \
  --bounds 0:50 \
  --periodic \
  --ic gaussian \
  --ic-component phi_0 \
  --ic-center 15.0 \
  --ic-width 3.0 \
  --t-end 10.0 \
  --output "$OUTPUT"

echo ""
echo "=== Sweep complete ==="
echo "Results: $OUTPUT/results.csv"

# Plot: P_max vs coupling strength
# tidal plot "$OUTPUT" --type sweep --metric P_max --output "$OUTPUT/sweep_P_max.png"
