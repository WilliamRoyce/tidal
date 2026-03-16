#!/usr/bin/env bash
# Gertsenshtein Effective — Plasma Detuning Sweep
#
# Sweeps the photon effective mass omegaP2 (= plasma frequency squared)
# from 0 to 2.0, demonstrating detuning suppression of conversion.
#
# Expected physics (Raffelt-Stodolsky master formula):
#   P = sin²(2θ) × sin²(Δ_osc * D / 2)
#   tan(2θ) = κ*B0 / |Δ|,  Δ = -omegaP2/(2ω),  Δ_osc = √(Δ² + (κB0)²)
#
# At omegaP2=0 (massless): full Rabi oscillation P = sin²(κB0D/2)
# At omegaP2>>0: suppressed as P ~ (κB0/omegaP2)² — detuned
#
# Uses --ic-wavevector 2 for a propagating graviton wavepacket.
# The velocity and resonance measurements show how the conversion
# bandwidth narrows as the mass mismatch increases.
#
# Ref: Raffelt & Stodolsky (1988), Phys. Rev. D 37:1237
#
# Running:
#   cd examples/coupled_scalars && bash sweep_mass_ratio.sh

set -euo pipefail
cd "$(dirname "$0")"

OUT=../data/sweep_mass_ratio

echo "=== Plasma Detuning Sweep: omegaP2 = 0.0 to 2.0 (15 points) ==="

tidal sweep ../data/coupled_scalars.json \
  --sweep "omegaP2=0.0:2.0:15" \
  --param kappa=1.0 --param B0=0.1 --param mg2=0.0 \
  --measure conversion,conservation \
  --source h_0 --target a_0 \
  --grid-shape 256 \
  --bounds 0:100 \
  --periodic \
  --ic plane-wave --ic-component h_0 --ic-wavevector 2.0 --ic-amplitude 0.1 \
  --t-end 50.0 \
  --output "$OUT"

echo ""
echo "Results: $OUT/"
echo ""
echo "Visualize:"
echo "  tidal plot $OUT/ --type sweep --metric P_max"
