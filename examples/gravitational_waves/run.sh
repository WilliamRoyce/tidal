#!/usr/bin/env bash
# Gravitational Waves 3+1D — Full derive → inspect → simulate pipeline
#
# Physics: xPert linearization of the Einstein equations with TT
# (transverse-traceless) gauge, reducing to physical polarization modes
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
#     --grid-shape 4,4,32 --bounds 0:4,0:4,0:20 --periodic \
#     --ic formula --ic-component h_4 \
#     --ic-formula "np.exp(-(z - 10.0)**2 / 4.5) * np.cos(0.6283 * (z - 10.0))" \
#     --t-end 5.0

set -euo pipefail
cd "$(dirname "$0")"

# Derive TT gauge-fixed linearized Einstein equations
tidal derive theory.toml
tidal inspect ../data/linearized_gravity.json

# Simulate with h_plus (h_4 = h_xx) initial condition
# Gaussian-modulated cosine wave packet propagating along z
# Quick test: coarser grid + shorter time. For production resolution,
# use --grid-shape 4,4,64 --bounds 0:4,0:4,0:40 --t-end 15.0
#
# NOTE: This system has gauge-unfixed constraints (h_0..h_3, h_9,
# h_transverse_*) with no self-terms — IDA's IDACalcIC may fail.
# Phase B (gauge fixing) will resolve this; for now, use leapfrog:
tidal simulate ../data/linearized_gravity.json \
  --grid-shape 4,4,32 \
  --bounds 0:4,0:4,0:20 \
  --periodic \
  --scheme leapfrog \
  --ic formula \
  --ic-formula "np.exp(-(z - 10.0)**2 / 4.5) * np.cos(0.6283 * (z - 10.0))" \
  --ic-component h_4 \
  --t-end 5.0
