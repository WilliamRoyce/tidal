#!/usr/bin/env bash
# Scalar-Vector Coupling 2+1D — Full derive → inspect → simulate → measure pipeline
#
# Physics: Mixed-rank cross-field coupling between a scalar phi and a vector A
# in 2+1D. This stress-tests the pipeline with:
#   - Mixed-rank cross-field coupling (scalar phi + vector A_0, A_1, A_2)
#   - Cross-field first_derivative_t and gradient operators
#   - Epsilon tensor (Chern-Simons) + cross-field coupling in same system
#   - 4 symbolic constants (phim2, Am2, kCS, gSV)
#   - 4x4 mass/coupling matrices
#
# A_0 is a constraint; phi_0, A_1, A_2 are dynamical.
#
# See also:
#   measure_conversion.py        — detailed measurement analysis
#
# Running this script:
#   cd examples/scalar_vector_coupling && uv run bash run.sh
#
# Or run each step manually to learn the tidal CLI:
#
#   # Step 1: Derive equations from the Lagrangian (requires wolframscript)
#   uv run tidal derive theory.toml
#
#   # Step 2: Inspect the generated equation system (4 fields, 4x4 matrices)
#   uv run tidal inspect ../data/scalar_vector_coupling.json
#
#   # Step 3: Simulate with Gaussian IC for phi (periodic BCs, 48x48 grid)
#   uv run tidal simulate ../data/scalar_vector_coupling.json \
#     --param phim2=1.0 --param Am2=0.5 --param kCS=0.3 --param gSV=0.2 \
#     --grid-shape 48 --bounds 0:10,0:10 --bc periodic,periodic \
#     --ic gaussian --t-end 5.0 --dt 0.005 \
#     --output ../data/scalar_vector_coupling_output
#
#   # Step 4: Measure scalar-to-vector conversion and mixing length
#   uv run tidal measure ../data/scalar_vector_coupling_output \
#     --what conversion,mixing --source phi_0 --target A_0,A_1,A_2
#
#   # Step 5: Spectral conversion P(k,t) — per-mode energy transfer
#   uv run tidal measure ../data/scalar_vector_coupling_output \
#     --what spectral_conversion --source phi_0 --target A_1,A_2
#
#   # Step 6: Dispersion relation omega(k) for phi
#   uv run tidal measure ../data/scalar_vector_coupling_output \
#     --what dispersion --source phi_0
#
#   # Step 7: Combined plot with all measurements
#   uv run tidal measure ../data/scalar_vector_coupling_output \
#     --what conversion,mixing,spectral_conversion,dispersion \
#     --source phi_0 --target A_0,A_1,A_2 \
#     --output ../data/scalar_vector_coupling_measurement.png

set -euo pipefail
cd "$(dirname "$0")"

# Step 1: Derive equations from Lagrangian (requires wolframscript)
tidal derive theory.toml

# Step 2: Inspect the equation system (phi_0 + A_0, A_1, A_2)
tidal inspect ../data/scalar_vector_coupling.json

# Step 3: Run simulation
# Gaussian IC for phi; periodic BCs ensure energy conservation
tidal simulate ../data/scalar_vector_coupling.json \
  --param phim2=1.0 --param Am2=0.5 --param kCS=0.3 --param gSV=0.2 \
  --grid-shape 48 \
  --bounds 0:10,0:10 \
  --bc periodic,periodic \
  --ic gaussian \
  --t-end 5.0 \
  --dt 0.005 \
  --output ../data/scalar_vector_coupling_output

# Step 4: Measure scalar-to-vector conversion and characteristic mixing length
# Source = scalar phi_0, Target = all vector components (group conversion)
# P(t) = E_A(t) / E_phi(0) tracks energy transfer from scalar to vector
tidal measure ../data/scalar_vector_coupling_output \
  --what conversion,mixing \
  --source phi_0 --target A_0,A_1,A_2

# Step 5: Spectral conversion P(k,t)
# Shows which Fourier modes participate in scalar-to-vector energy transfer
tidal measure ../data/scalar_vector_coupling_output \
  --what spectral_conversion \
  --source phi_0 --target A_1,A_2

# Step 6: Dispersion relation omega(k) for the scalar field
tidal measure ../data/scalar_vector_coupling_output \
  --what dispersion \
  --source phi_0

# Step 7: Combined measurement plot (all panels in one figure)
tidal measure ../data/scalar_vector_coupling_output \
  --what conversion,mixing,spectral_conversion,dispersion \
  --source phi_0 --target A_0,A_1,A_2 \
  --output ../data/scalar_vector_coupling_measurement.png
