#!/usr/bin/env bash
# Gravitational Waves 3+1D — Full derive → inspect → simulate pipeline
#
# Physics: xPert linearization of the Einstein equations with TT
# (transverse-traceless) gauge, reducing to physical polarisation modes
# h_+ and h_×. Gaussian-modulated cosine wave packet propagating along z.
#
# NOTE: 3D data — tidal plot is for 1D/2D only.
#
# Running this script:
#   cd examples/gravitational_waves && bash run.sh
#
# Or run each step manually:
#   tidal derive theory.toml
#   tidal inspect ../data/linearized_gravity.json
#   tidal simulate ../data/linearized_gravity.json \
#     --grid-shape 4,4,64 --bounds 0:4,0:4,0:40 --periodic \
#     --ic formula --ic-component h_4 \
#     --ic-formula "np.exp(-(z - 20.0)**2 / 18.0) * np.cos(0.6283 * (z - 20.0))" \
#     --t-end 15.0 --scheme scipy

set -euo pipefail
cd "$(dirname "$0")"

# Derive TT gauge-fixed linearized Einstein equations
tidal derive theory.toml
tidal inspect ../data/linearized_gravity.json

# Simulate with h_plus (h_4 = h_xx) initial condition
# Gaussian-modulated cosine wave packet propagating along z
tidal simulate ../data/linearized_gravity.json \
  --grid-shape 4,4,64 \
  --bounds 0:4,0:4,0:40 \
  --periodic \
  --ic formula \
  --ic-formula "np.exp(-(z - 20.0)**2 / 18.0) * np.cos(0.6283 * (z - 20.0))" \
  --ic-component h_4 \
  --t-end 15.0 \
  --scheme scipy
