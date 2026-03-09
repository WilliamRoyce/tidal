#!/usr/bin/env bash
# Maxwell Matter Perturbation (CI Test) — Derive → Inspect → Simulate → Measure
#
# Physics: Proca vector field A perturbed around position-dependent background
# Abar = (0, -B0*x) in 1+1D Minkowski. After linearization with Lorenz gauge,
# the perturbation satisfies free Proca equations (mass m2). The background B0
# validates the matter perturbation pipeline infrastructure but doesn't appear
# in the linearized equations (Proca is quadratic → linear perturbation equations
# are background-independent).
#
# This is a fast CI-level test (~seconds to derive, ~seconds to simulate).
#
# Running:
#   cd examples/maxwell_perturbation && bash run.sh

set -euo pipefail
cd "$(dirname "$0")"

echo "=== Maxwell Matter Perturbation (CI Test) ==="
echo ""

# Step 1: Derive Proca perturbation equations
echo "--- Step 1: Derive ---"
tidal derive theory.toml
echo ""

# Step 2: Inspect
echo "--- Step 2: Inspect ---"
tidal inspect ../data/maxwell_perturbation.json
echo ""

# Step 3: Simulate — Gaussian wave packet on spatial component a_1
# IC on a_1 (spatial, stable mass term) not a_0 (temporal, unstable in Lorenz gauge)
echo "--- Step 3: Simulate ---"
tidal simulate ../data/maxwell_perturbation.json \
  --grid-shape 128 \
  --bounds 0:50 \
  --periodic \
  --ic gaussian \
  --ic-amplitude 0.1 \
  --ic-width 3.0 \
  --ic-center 25.0 \
  --ic-component a_1 \
  --t-end 10.0 \
  --param m2=1.0 --param B0=0.5 \
  --output ../data/maxwell_perturbation_output

echo ""

# Step 4: Measure energy conservation
echo "--- Step 4: Measure ---"
tidal measure ../data/maxwell_perturbation_output --what conservation --param m2=1.0 --param B0=0.5
tidal measure ../data/maxwell_perturbation_output --what energy --param m2=1.0 --param B0=0.5
echo ""

echo "=== Done ==="
