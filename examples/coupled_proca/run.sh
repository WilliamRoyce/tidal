#!/bin/bash
# Coupled Proca: two massive vector fields in a 2+1D cavity (Dirichlet BCs)
#
# This exercises the non-periodic constraint solver:
#   - Sparse matrix solver (not FFT)
#   - Gauss-Seidel iteration for coupled A_0-B_0 constraints
#   - Dirichlet boundary conditions
set -euo pipefail

tidal simulate examples/data/coupled_proca_3d.json \
  --param mA2=1.0 --param mB2=2.0 --param gcoup=0.5 \
  --ic gaussian --grid-size 16 --t-end 2.0 --dt 0.05 \
  --bc dirichlet,dirichlet --scheme runge-kutta \
  --output outputs/coupled_proca_output.png
