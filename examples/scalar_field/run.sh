#!/usr/bin/env bash
# Klein-Gordon 1+1D — Full derive → inspect → simulate → plot pipeline
#
# Physics: Massive scalar field φ in 1+1D Minkowski spacetime with mass m².
# Gaussian pulse propagation with dispersion from the mass term.
#
# Running this script:
#   cd examples/scalar_field && bash run.sh
#
# Or run each step manually to learn the tidal CLI:
#
#   # Step 1: Derive equations from the Lagrangian (requires wolframscript)
#   tidal derive theory.toml
#
#   # Step 2: Inspect the generated equation system
#   tidal inspect ../data/klein_gordon_1d.json
#
#   # Step 3: Simulate — Gaussian pulse with mass-induced dispersion
#   tidal simulate ../data/klein_gordon_1d.json \
#     --param m2=1.0 --grid-shape 256 --bounds 0:100 --periodic \
#     --ic gaussian --ic-width 5.0 --t-end 30.0 --output ../data/scalar_field_output
#
#   # Step 4: Generate plots from simulation output
#   tidal plot ../data/scalar_field_output --type heatmap --quiet
#   tidal plot ../data/scalar_field_output --type amplitude --quiet

set -euo pipefail
cd "$(dirname "$0")"

# Derive equations from Lagrangian (requires wolframscript)
tidal derive theory.toml

# Or preview the generated Wolfram script without running:
# tidal derive theory.toml --dry-run

# Inspect the equation system
tidal inspect ../data/klein_gordon_1d.json

# Run simulation (Gaussian pulse, disk-backed output)
tidal simulate ../data/klein_gordon_1d.json \
  --param m2=1.0 \
  --grid-shape 256 \
  --bounds 0:100 \
  --periodic \
  --ic gaussian \
  --ic-width 5.0 \
  --t-end 30.0 \
  --output ../data/scalar_field_output

# Visualize results (plots saved into the simulation output directory)
tidal plot ../data/scalar_field_output --type heatmap --quiet
tidal plot ../data/scalar_field_output --type amplitude --quiet
tidal plot ../data/scalar_field_output --type snapshot --time-index 0 --quiet
tidal plot ../data/scalar_field_output --type snapshot --time-index -1 --quiet
tidal plot ../data/scalar_field_output --type profile --quiet
tidal plot ../data/scalar_field_output --type hamiltonian --quiet
tidal plot ../data/scalar_field_output --type conservation --quiet

# ============================================================================
# Run 2: Travelling wave packet (right-mover with carrier frequency)
# ============================================================================
# The static Gaussian above disperses symmetrically in both directions.
# Using --ic-formula-velocity, we lock the field and velocity together so
# the packet travels cleanly rightward — the bread-and-butter IC for
# wave physics simulations.
#
# For a right-mover cos(kx - wt): field = envelope*cos(kx),
# velocity = df/dt|_0 = +w*envelope*sin(kx) where w = |k| (massless limit).

tidal simulate ../data/klein_gordon_1d.json \
  --param m2=1.0 \
  --grid-shape 256 \
  --bounds 0:100 \
  --periodic \
  --ic formula \
  --ic-formula 'exp(-(x - 25)**2 / 50) * cos(2 * x)' \
  --ic-formula-velocity '2 * exp(-(x - 25)**2 / 50) * sin(2 * x)' \
  --t-end 30.0 \
  --output ../data/scalar_field_travelling

# Compare: static (symmetric spreading) vs travelling (rightward propagation)
tidal plot ../data/scalar_field_travelling --type heatmap \
  --title "Travelling wave packet (right-mover)" --quiet
tidal plot ../data/scalar_field_travelling --type profile --quiet
tidal plot ../data/scalar_field_travelling --type conservation --quiet
