#!/usr/bin/env bash
# CLI equivalents for the Gravitational Waves 3+1D example
# See also: linearized_gravity.wls (manual derivation), gw_simulation.py (Python simulation)
#
# NOTE: This example uses xPert linearization of the Einstein equations,
# which requires a separate workflow beyond what `tg derive` TOML supports.
# The derive step requires the manual .wls script.
# The simulation uses custom TT-gauge initial conditions (plane-wave packet)
# that are not yet expressible via CLI IC presets.
# The Python script is recommended for full-fidelity simulation.

set -euo pipefail

# Derive equations (manual .wls only — xPert linearization)
# tg derive linearized_gravity.wls    # pass-through to wolframscript

# Inspect the equation system (10 components: h_00..h_33)
tg inspect ../data/linearized_gravity.json

# Approximate simulation with Gaussian IC (not TT-gauge)
# For proper TT-gauge initial conditions, use gw_simulation.py
tg simulate ../data/linearized_gravity.json \
  --grid-shape 4,4,64 \
  --bounds 0:4,0:4,0:40 \
  --periodic \
  --ic gaussian \
  --ic-component h_5 \
  --ic-width 3.0 \
  --t-end 15.0 \
  --dt 0.01
