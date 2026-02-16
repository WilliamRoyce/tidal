#!/usr/bin/env bash
# Coupled Proca 2+1D — Full derive → inspect → simulate → measure pipeline
#
# Physics: Two massive vector fields (A, B) in 2+1D with periodic BCs.
# Each field has 3 components (A_0, A_1, A_2) and (B_0, B_1, B_2).
# The A_0 and B_0 components are constraints (solved via coupled FFT);
# A_1, A_2, B_1, B_2 are dynamical (time-evolved).
#
# This exercises the coupled constraint solver:
#   - Coupled FFT solve (periodic BCs, operator-matrix method)
#   - Two Helmholtz scales (mA2 vs mB2)
#   - Cross-field identity coupling in constraints (gcoup)
#
# See also:
#   coupled_proca.wls          — manual Wolfram derivation
#   measure_conversion.py      — detailed measurement analysis
#
# Running this script:
#   cd examples/coupled_proca && uv run bash run.sh
#
# Or run each step manually to learn the tidal CLI:
#
#   # Step 1: Derive equations from the Lagrangian (requires wolframscript)
#   uv run tidal derive theory.toml
#
#   # Step 2: Inspect the generated equation system (6 fields, 6x6 matrices)
#   uv run tidal inspect ../data/coupled_proca_3d.json
#
#   # Step 3: Simulate (Gaussian IC, periodic BCs, 16x16 grid)
#   uv run tidal simulate ../data/coupled_proca_3d.json \
#     --param mA2=1.0 --param mB2=2.0 --param gcoup=0.5 \
#     --ic gaussian --grid-shape 16 --t-end 2.0 --dt 0.05 \
#     --bc periodic,periodic --scheme runge-kutta \
#     --output ../data/coupled_proca_output
#
#   # Step 4: Measure conversion from A-field group to B-field group
#   uv run tidal measure ../data/coupled_proca_output \
#     --what conversion,mixing \
#     --source A_0,A_1,A_2 --target B_0,B_1,B_2
#
#   # Step 5: Spectral conversion P(k,t) — per-mode A->B transfer
#   uv run tidal measure ../data/coupled_proca_output \
#     --what spectral_conversion --source A_1 --target B_1,B_2
#
#   # Step 6: Dispersion relation omega(k) for A_1
#   uv run tidal measure ../data/coupled_proca_output \
#     --what dispersion --source A_1
#
#   # Step 7: Combined plot with all measurements
#   uv run tidal measure ../data/coupled_proca_output \
#     --what conversion,mixing,spectral_conversion,dispersion \
#     --source A_0,A_1,A_2 --target B_0,B_1,B_2 \
#     --output ../data/coupled_proca_measurement.png

set -euo pipefail
cd "$(dirname "$0")"

# Step 1: Derive equations from Lagrangian (requires wolframscript)
tidal derive theory.toml

# Step 2: Inspect the equation system (A_0..A_2 + B_0..B_2)
tidal inspect ../data/coupled_proca_3d.json

# Step 3: Run simulation
# Gaussian IC with periodic BCs; constraint solver auto-detects A_0, B_0
tidal simulate ../data/coupled_proca_3d.json \
  --param mA2=1.0 --param mB2=2.0 --param gcoup=0.5 \
  --ic gaussian --grid-shape 16 --t-end 2.0 --dt 0.05 \
  --bc periodic,periodic --scheme runge-kutta \
  --output ../data/coupled_proca_output

# Step 4: Measure conversion between vector field groups and mixing length
# Source = all A components, Target = all B components (group conversion)
# P(t) = E_B(t) / E_A(0) tracks energy transfer from A to B
tidal measure ../data/coupled_proca_output \
  --what conversion,mixing \
  --source A_0,A_1,A_2 --target B_0,B_1,B_2

# Step 5: Spectral conversion P(k,t) — per-mode A->B energy transfer
tidal measure ../data/coupled_proca_output \
  --what spectral_conversion \
  --source A_1 --target B_1,B_2

# Step 6: Dispersion relation omega(k) for A_1
tidal measure ../data/coupled_proca_output \
  --what dispersion \
  --source A_1

# Step 7: Combined measurement plot (all panels in one figure)
tidal measure ../data/coupled_proca_output \
  --what conversion,mixing,spectral_conversion,dispersion \
  --source A_0,A_1,A_2 --target B_0,B_1,B_2 \
  --output ../data/coupled_proca_measurement.png
