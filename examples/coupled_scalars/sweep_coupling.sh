#!/usr/bin/env bash
# Gertsenshtein Effective Theory — B0 coupling strength sweep
#
# Sweeps the background magnetic field strength B0, measuring the
# graviton-photon conversion probability P(t) = sin²(κ*B0*t/2).
#
# Expected physics (massless vacuum, κ=1):
#   P(t) = sin²(B0 * t_end / 2) → Rabi oscillation in B0
#   First maximum at B0 = π/t_end ≈ 0.063 (for t_end=50)
#   Period: ΔB0 = 2π/t_end ≈ 0.126
#
# This sweeps B0 from weak to strong coupling, covering ~2 Rabi cycles.
#
# Ref: Gertsenshtein (1962), JETP 14:84
#      Raffelt & Stodolsky (1988), Phys. Rev. D 37:1237
#
# Running:
#   cd examples/coupled_scalars && bash sweep_coupling.sh

set -euo pipefail
cd "$(dirname "$0")"

OUT=../data/coupled_scalars_sweep

echo "=== Gertsenshtein Effective: B0 Coupling Sweep ==="
echo "    B0 = 0.005 to 0.25 (20 points), kappa=1, omegaP2=0, mg2=0"
echo "    Plane-wave k=2.0 on h_0, domain [0,100], N=256, t_end=50"
echo ""

tidal sweep ../data/coupled_scalars.json \
  --sweep "B0=0.005:0.25:20" \
  --param kappa=1.0 --param omegaP2=0.0 --param mg2=0.0 \
  --measure conversion,conservation \
  --source h_0 --target a_0 \
  --grid-shape 256 \
  --bounds 0:100 \
  --periodic \
  --ic plane-wave --ic-component h_0 --ic-wavevector 2.0 --ic-amplitude 0.1 \
  --t-end 50.0 \
  --output "$OUT"

echo ""

# Plot with Gertsenshtein analytical overlay
tidal plot "$OUT" --type sweep \
  --metric P_final \
  --title "Gertsenshtein P(B0) — effective 1+1D theory" \
  --overlay 'sin(kappa * B0 * t_end / 2)**2' \
  --output "$OUT/sweep.png" --quiet

echo "Results: $OUT/"
echo "Plot:    $OUT/sweep.png"
