#!/usr/bin/env bash
# Gertsenshtein with Massive Photon — mA² family oscillation traces
#
# Produces P(t) traces at fixed B₀ = 0.10 for a family of photon masses
# spanning on-resonance through off-resonance. Used to visualize the
# transition from coherent Rabi oscillation to Raffelt-Stodolsky
# suppression.
#
# Resonance for κ=1, B₀=0.10:  mA² = κ²B₀²/2 = 0.005
# Coupling scale:              2·κ·B₀·k = 0.402  (k = 2.0106)
#
# mA² values chosen to show:
#   0.000  — off by resonance/2, near-resonant Rabi
#   0.005  — exactly on resonance, peak mixing
#   0.050  — modestly detuned (~0.1·coupling²)
#   0.200  — detuned by ~half coupling, amplitude reduced
#   0.500  — detuned by ~coupling, P_max ~ 0.5
#   1.000  — strongly detuned, P_max ~ 0.15
#
# Ref: Raffelt & Stodolsky 1988, Phys. Rev. D 37:1237
#
# Running:
#   cd examples/gertsenshtein_proca && bash sweep_mA2_family.sh

set -euo pipefail
cd "$(dirname "$0")"

JSON=../data/gertsenshtein_proca.json
OUT=../data/sweep_proca_mA2_family

echo "=== Massive Photon Gertsenshtein: mA² Family (fixed B₀=0.10) ==="
echo "    6 simulations across on- and off-resonance mA² values"
echo ""

# Use explicit list of mA² values — more pedagogically useful than a
# uniform sweep because the points are chosen to span the Lorentzian.
tidal sweep "$JSON" \
  --sweep "mA2=0.0,0.005,0.05,0.2,0.5,1.0" \
  --param kappa=1.0 --param B0=0.10 \
  --measure conversion,conservation \
  --source h_7 --target a_2 \
  --grid-shape 512 \
  --bounds 0:100 \
  --periodic \
  --ic plane-wave --ic-component h_7 --ic-wavevector 2.0106 --ic-amplitude 0.1 \
  --t-end 50.0 \
  --output "$OUT" --resume

echo ""
echo "Results: $OUT/"
echo "Use reproduce_figures.sh to generate the overlay plot of P(t) vs mA²."
