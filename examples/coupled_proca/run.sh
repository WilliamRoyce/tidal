#!/usr/bin/env bash
# Coupled Proca 2+1D — Full simulate → measure pipeline
#
# Physics: Two massive vector fields (A, B) in 2+1D with periodic BCs.
# Each field has 3 components (A_0, A_1, A_2) and (B_0, B_1, B_2).
# The A_0 and B_0 components are constraints (solved via coupled FFT);
# A_1, A_2, B_1, B_2 are dynamical (time-evolved).
#
# This exercises the coupled constraint solver:
#   - Coupled FFT solve (periodic BCs, operator-matrix method)
#   - Two Helmholtz scales (mA2 vs mB2)
#   - Cross-field identity coupling in constraints (gcoup)
#
# See also:
#   coupled_proca.wls          — manual Wolfram derivation
#   measure_conversion.py      — detailed measurement analysis
#
# Running this script (from the repository root):
#   uv run bash examples/coupled_proca/run.sh
#
# Or run each step manually to learn the tidal CLI:
#
#   # Simulate (Gaussian IC, periodic BCs, 16x16 grid)
#   uv run tidal simulate examples/data/coupled_proca_3d.json \
#     --param mA2=1.0 --param mB2=2.0 --param gcoup=0.5 \
#     --ic gaussian --grid-shape 16 --t-end 2.0 --dt 0.05 \
#     --bc periodic,periodic --scheme runge-kutta \
#     --output outputs/coupled_proca_output
#
#   # Measure conversion from A-field group to B-field group
#   uv run tidal measure outputs/coupled_proca_output \
#     --spec examples/data/coupled_proca_3d.json \
#     --what conversion,mixing \
#     --source A_0,A_1,A_2 --target B_0,B_1,B_2
#
#   # Spectral conversion P(k,t) — per-mode A->B transfer
#   uv run tidal measure outputs/coupled_proca_output \
#     --spec examples/data/coupled_proca_3d.json \
#     --what spectral_conversion --source A_1 --target B_1,B_2
#
#   # Dispersion relation omega(k) for A_1
#   uv run tidal measure outputs/coupled_proca_output \
#     --spec examples/data/coupled_proca_3d.json \
#     --what dispersion --source A_1

set -euo pipefail

# Simulate (from repo root — paths are relative to repo root)
# Gaussian IC with periodic BCs; constraint solver auto-detects A_0, B_0
tidal simulate examples/data/coupled_proca_3d.json \
  --param mA2=1.0 --param mB2=2.0 --param gcoup=0.5 \
  --ic gaussian --grid-shape 16 --t-end 2.0 --dt 0.05 \
  --bc periodic,periodic --scheme runge-kutta \
  --output outputs/coupled_proca_output

# Measure conversion between vector field groups and mixing length
# Source = all A components, Target = all B components (group conversion)
# P(t) = E_B(t) / E_A(0) tracks energy transfer from A to B
tidal measure outputs/coupled_proca_output \
  --spec examples/data/coupled_proca_3d.json \
  --what conversion,mixing \
  --source A_0,A_1,A_2 --target B_0,B_1,B_2

# Spectral conversion P(k,t) — per-mode A->B energy transfer
tidal measure outputs/coupled_proca_output \
  --spec examples/data/coupled_proca_3d.json \
  --what spectral_conversion \
  --source A_1 --target B_1,B_2

# Dispersion relation omega(k) for A_1
tidal measure outputs/coupled_proca_output \
  --spec examples/data/coupled_proca_3d.json \
  --what dispersion \
  --source A_1
