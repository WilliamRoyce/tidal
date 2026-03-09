#!/usr/bin/env bash
# Massive Gravity 2+1D (Fierz-Pauli) — Full derive → inspect → simulate → plot pipeline
#
# Physics: Fierz-Pauli mass term: G^(1)_ab - m² (h_ab - η_ab h) = 0
# The trace h = η^cd h_cd couples diagonal metric components (h_0, h_3, h_5).
# This is the unique ghost-free linear mass term for spin-2. Uses xPert
# linearization with L^(2) preprocessing.
#
# Running this script:
#   cd examples/massive_gravity && bash run.sh
#
# Or run each step manually:
#   tidal derive theory.toml
#   tidal inspect ../data/massive_gravity_3d.json
#   tidal simulate ../data/massive_gravity_3d.json --param m2=1.0 \
#     --grid-shape 32 --bounds 0:50 --periodic --ic formula \
#     --ic-component h_3 \
#     --ic-formula "np.exp(-((x-25)**2+(y-25)**2)/50) * np.sin(2*np.pi*x/50)" \
#     --t-end 2.0 --output ../data/massive_gravity_output
#   tidal plot ../data/massive_gravity_output --type amplitude --quiet

set -euo pipefail
cd "$(dirname "$0")"

# Derive equations from TOML config (xPert linearization with FP mass term)
tidal derive theory.toml

# Inspect the equation system (6 components: h_0..h_5)
tidal inspect ../data/massive_gravity_3d.json

# Simulate with sine-modulated Gaussian in h_3 (h_xx component).
# h_3 is dynamical with dispersion omega^2 = k^2 + 2*m^2.
# The h_0 (lapse) equation doesn't involve h_0 — it's a compatibility
# constraint on spatial components h_3/h_4/h_5 (the Hamiltonian constraint
# G^(1)_00 = m²(h_3+h_5) doesn't involve h_00).  h_0 stays at zero
# (temporal gauge).  The --allow-inconsistent-ic flag is used because
# the Gaussian × sin IC doesn't exactly satisfy the compatibility constraint.
tidal simulate ../data/massive_gravity_3d.json \
  --atol 1e-4 \
  --rtol 1e-6 \
  --param m2=1.0 \
  --grid-shape 32 \
  --bounds 0:50 \
  --periodic \
  --ic formula \
  --ic-component h_3 \
  --ic-formula "np.exp(-((x - 25)**2 + (y - 25)**2) / 50) * np.sin(2 * np.pi * x / 50)" \
  --t-end 2.0 \
  --allow-inconsistent-ic \
  --output ../data/massive_gravity_output

# Visualize results (plots saved into the simulation output directory)
tidal plot ../data/massive_gravity_output --type snapshot --field h_3 --time-index 0 --quiet
tidal plot ../data/massive_gravity_output --type snapshot --field h_3 --time-index -1 --quiet
tidal plot ../data/massive_gravity_output --type amplitude --overlay 'cos(sqrt(2)*t)*0.5' --quiet
tidal plot ../data/massive_gravity_output --type profile --field h_3 --cross-section y=25.0 --quiet

# For parameter sweep (vary mass):
# for m2 in 0.5 1.0 2.0; do
#   tidal simulate ../data/massive_gravity_3d.json --param m2=$m2
# done
