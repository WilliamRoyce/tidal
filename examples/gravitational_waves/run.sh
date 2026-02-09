#!/usr/bin/env bash
# CLI equivalents for the Gravitational Waves 3+1D example
# See also: linearized_gravity.wls (manual derivation), gw_simulation.py (Python simulation)
#
# NOTE: This example uses xPert linearization of the Einstein equations,
# which requires a separate workflow beyond what `tidal derive` TOML supports.
# The derive step requires the manual .wls script.
# The simulation uses --ic formula for TT-gauge initial conditions.
# For the full multi-component TT-gauge setup, see gw_simulation.py.
#
# To run manually:  cd examples/gravitational_waves

set -euo pipefail
cd "$(dirname "$0")"

# Derive equations (xPert linearization, .wls pass-through to wolframscript)
tidal derive linearized_gravity.wls

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
  --dt 0.01

# For full TT-gauge IC with h_plus, h_cross, and h_yy = -h_plus,
# use the Python script: python gw_simulation.py
