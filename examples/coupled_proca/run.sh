#!/bin/bash
# Coupled Proca: two massive vector fields in 2+1D (periodic BCs)
#
# This exercises the coupled constraint solver:
#   - Coupled FFT solve (periodic BCs)
#   - Two Helmholtz scales (mA2 vs mB2)
#   - Cross-field identity coupling in constraints
set -euo pipefail

tidal simulate examples/data/coupled_proca_3d.json \
  --param mA2=1.0 --param mB2=2.0 --param gcoup=0.5 \
  --ic gaussian --grid-size 16 --t-end 2.0 --dt 0.05 \
  --bc periodic,periodic --scheme runge-kutta \
  --output outputs/coupled_proca_output.npz

# Measure conversion between vector field groups and mixing length
tidal measure outputs/coupled_proca_output.npz \
  --spec examples/data/coupled_proca_3d.json \
  --what conversion,mixing \
  --source A_0,A_1,A_2 --target B_0,B_1,B_2
