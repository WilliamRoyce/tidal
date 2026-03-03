#!/usr/bin/env bash
# Gravitational Waves 3+1D — Full derive → inspect → simulate pipeline
#
# Physics: xPert linearization of the Einstein equations with TT
# (transverse-traceless) gauge, reducing to physical polarization modes
# h_+ and h_×. Gaussian-modulated cosine wave packet propagating along z.
#
# TT gauge IC requirements:
#   - Traceless: h_xx + h_yy + h_zz = 0  →  h_7 = -h_4 (for h_9 = h_zz = 0)
#   - Transverse: ∂_i h_{ij} = 0  →  satisfied when h_4/h_7 depend only on z
#   For + polarization (z-propagating): set h_4 (h_xx), h_7 = -h_4 (h_yy),
#   all other spatial components zero.
#
#
# Running this script:
#   cd examples/gravitational_waves && bash run.sh
#
# Or run each step manually:
#   tidal derive theory.toml
#   tidal inspect ../data/linearized_gravity.json
#   tidal simulate ../data/linearized_gravity.json \
#     --grid-shape 4,4,32 --bounds 0:4,0:4,0:20 --periodic \
#     --ic formula --ic-component h_4 \
#     --ic-formula "np.exp(-(z - 10.0)**2 / 4.5) * np.cos(0.6283 * (z - 10.0))" \
#     --ic-field "h_7:-(np.exp(-(z - 10.0)**2 / 4.5) * np.cos(0.6283 * (z - 10.0)))" \
#     --t-end 5.0 --output ../data/linearized_gravity_output

set -euo pipefail
cd "$(dirname "$0")"

# Derive TT gauge-fixed linearized Einstein equations
tidal derive theory.toml
tidal inspect ../data/linearized_gravity.json

# Simulate with h_plus (h_4 = h_xx) initial condition
# Gaussian-modulated cosine wave packet propagating along z.
# TT gauge requires h_7 (h_yy) = -h_4 (traceless) for the constraint
# solver to cascade correctly (h_9 = -(h_4 + h_7) = 0).
# Gauge-unfixed constraints (h_1..h_3, h_transverse_*) have no self-terms
# — IDA automatically freezes these at zero (temporal gauge).
tidal simulate ../data/linearized_gravity.json \
  --rtol 1e-4 \
  --atol 1e-6 \
  --grid-shape 4,4,32 \
  --bounds 0:4,0:4,0:20 \
  --periodic \
  --ic formula \
  --ic-formula "np.exp(-(z - 10.0)**2 / 4.5) * np.cos(0.6283 * (z - 10.0))" \
  --ic-component h_4 \
  --ic-field "h_7:-(np.exp(-(z - 10.0)**2 / 4.5) * np.cos(0.6283 * (z - 10.0)))" \
  --t-end 5.0 \
  --output ../data/linearized_gravity_output

# Visualize results (plots saved into the simulation output directory)
tidal plot ../data/linearized_gravity_output --type snapshot --field h_4 --time-index 0 --quiet
tidal plot ../data/linearized_gravity_output --type snapshot --field h_4 --time-index -1 --quiet
