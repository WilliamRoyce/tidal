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
# Note: The derived JSON contains 12 fields (5 dynamical + 7 constraints).
# IDA correctly handles the TT gauge constraints. The physical DOFs are
# h_4 (h_+ = h_xx) and h_5 (h_x = h_xy), each satisfying d2t = d2z.
#
# Note: --ic-component h_4 initializes only h_4 without its traceless
# partner h_7 = -h_4. This violates the transverse subsidiary constraint
# but does not affect the wave equation dynamics (h_4 is self-contained).
# For fully TT-gauge-consistent IC, one would need h_4 = -h_7 = Gaussian.
#
# Running this script:
#   cd examples/gravitational_waves_1d && bash run.sh

set -euo pipefail
cd "$(dirname "$0")"

# Derive equations from Lagrangian (requires wolframscript)
tidal derive theory.toml

# Inspect the equation system (should show 1+1D, 12 fields with constraints)
tidal inspect ../data/gw_plane_wave_1d.json

# Run 1D simulation (Gaussian pulse on h_+ polarization)
tidal simulate ../data/gw_plane_wave_1d.json \
  --grid-shape 256 \
  --bounds=-10:10 \
  --bc neumann \
  --ic gaussian \
  --ic-width 1.0 \
  --ic-component h_4 \
  --t-end 8.0

# Energy conservation measurement (currently blocked by unsupported
# mixed_* Hamiltonian operators — see .github-issues-pending.md)
# tidal measure output/ --what conservation
