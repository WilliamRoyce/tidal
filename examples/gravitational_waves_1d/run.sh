#!/usr/bin/env bash
# Gravitational Waves (TT gauge + z-propagation) — Reduced 1+1D pipeline
#
# Physics: Linearized gravity in TT gauge with plane-wave reduction along z.
# The full 3+1D linearized Einstein equations (10 symmetric h_ab components,
# reduced to 5 after temporal+traceless+transverse gauge conditions) further
# reduce to just 2 uncoupled 1D wave equations for the plus and cross
# polarizations:
#   d2t(h_+) = d2z(h_+)     (plus polarization, h_xx = -h_yy)
#   d2t(h_x) = d2z(h_x)     (cross polarization, h_xy)
#
# This is the primary use case for plane-wave reduction: testing GW
# propagation along a single spatial axis at dramatically lower cost
# than a full 3D simulation.
#
# Running this script:
#   cd examples/gravitational_waves_1d && bash run.sh

set -euo pipefail
cd "$(dirname "$0")"

# Derive equations from Lagrangian (requires wolframscript)
tidal derive theory.toml

# Inspect the equation system (should show 1+1D, 2 fields)
tidal inspect ../data/gw_plane_wave_1d.json

# Run 1D simulation (Gaussian pulse for both polarizations)
tidal simulate ../data/gw_plane_wave_1d.json \
  --grid-shape 256 \
  --bounds -10:10 \
  --bc neumann \
  --ic gaussian \
  --ic-width 1.0 \
  --t-end 8.0

# Check energy conservation
tidal measure output/ --what conservation
