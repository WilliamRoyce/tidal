#!/usr/bin/env bash
# Gravitational Waves 3+1D — Full derive → inspect → simulate pipeline
#
# This example uses xPert linearization of the Einstein equations.
# The TOML config produces gauge-unfixed linearized Einstein equations (10 components).
# The simulation uses --ic formula for TT-gauge initial conditions.
#
# To run manually:  cd examples/gravitational_waves

set -euo pipefail
cd "$(dirname "$0")"

# Derive equations from TOML config (xPert linearization, gauge-unfixed)
tidal derive theory.toml

# Inspect the equation system (10 components: h_00..h_33)
tidal inspect ../data/linearized_gravity.json

# Simulate with TT-gauge h_plus (h_4 = h_xx) initial condition
# This applies a Gaussian-modulated cosine wave packet along z
tidal simulate ../data/linearized_gravity.json \
  --grid-shape 4,4,64 \
  --bounds 0:4,0:4,0:40 \
  --periodic \
  --ic formula \
  --ic-formula "np.exp(-(z - 20.0)**2 / 18.0) * np.cos(0.6283 * (z - 20.0))" \
  --ic-component h_4 \
  --t-end 15.0 \
  --scheme scipy
