#!/usr/bin/env bash
# CLI equivalents for the Massive Gravity 2+1D example (Fierz-Pauli)
# See also: simulation.py (Python simulation)
#
# Fierz-Pauli mass term: G^(1)_ab - m^2 (h_ab - eta_ab h) = 0
# The trace h = eta^cd h_cd couples diagonal metric components (h_0, h_3, h_5).
# This is the unique ghost-free linear mass term for spin-2.
#
# To run manually:  cd examples/massive_gravity

set -euo pipefail
cd "$(dirname "$0")"

# Derive equations from TOML config (xPert linearization with FP mass term)
tidal derive theory.toml

# Inspect the equation system (6 components: h_0..h_5)
tidal inspect ../data/massive_gravity_3d.json

# Simulate with Gaussian perturbation in h_3 (h_xx spatial component)
# The mass term creates dispersive propagation: omega^2 = k^2 + m^2
tidal simulate ../data/massive_gravity_3d.json \
  --param m2=1.0 \
  --grid-shape 64 \
  --bounds 0:50 \
  --periodic \
  --ic gaussian \
  --ic-component h_3 \
  --t-end 5.0 \
  --scheme scipy

# For parameter sweep (vary mass):
# for m2 in 0.5 1.0 2.0; do
#   tidal simulate ../data/massive_gravity_3d.json --param m2=$m2
# done

# For detailed simulation with physics validation, use:
# python simulation.py
