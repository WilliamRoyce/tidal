#!/usr/bin/env bash
# Sweep B — 2D Heatmap: δ vs ξ at fixed on-resonance (mA²=mT²=0.5)
#
# Produces a 50×50 Cartesian grid exploring the kinetic-mixing /
# ghost-freedom plane. The ghost boundary is the analytical curve
# |δ| = √ξ/2, equivalently δ² = ξ/4. Beyond this boundary one
# kinetic eigenvalue of K_kin = [[1, -2δ], [-2δ, ξ]] flips sign
# → ghost DOF → simulation diverges.
#
# Physics:
#   At resonance mA²=mT², the photon and dark-photon dispersions match,
#   giving maximal energy transfer. The kinetic mixing δ and the kinetic
#   coefficient ξ jointly control the coupling strength:
#     E[t,a] = -2δ·mA² / (4δ²-ξ)
#   For ghost-free parameters (4δ²<ξ), E[t,a] > 0 and increases with δ.
#
# Ghost-freedom (verified symbolically):
#   K_kin = [[1, -2δ], [-2δ, ξ]], det = ξ-4δ² > 0 iff |δ| < √ξ/2.
#   At ξ=1: eigenvalues are 1±2δ, both positive iff |δ|<1/2.
#   The boundary δ²=ξ/4 should appear as the divergence frontier.
#
# Metric: A_total = P_max / sin²(κB₀t_end/2)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SPEC="${SCRIPT_DIR}/../data/dark_photon_plasma.json"
OUTPUT="${SCRIPT_DIR}/../data/dp_plasma_heatmap_dxi"

# Fixed physics (on resonance)
KAPPA=1.0
B0=0.01
MA2=0.5
MT2=0.5
K0=2.0

# Grid and simulation settings
GRID=256
BOUNDS="0:100"
T_END=50.0

echo "=== Sweep B: δ vs ξ heatmap (50×50) ==="
echo "Fixed: κ=${KAPPA}, B₀=${B0}, mA²=${MA2}, mT²=${MT2} (on resonance)"
echo "Ghost boundary: |δ| = √ξ/2"
echo "Output: ${OUTPUT}"
echo ""

tidal sweep "${SPEC}" \
  --sweep "deltam=-2.0:2.0:50" \
  --sweep "xi=0.1:5.0:50" \
  --measure peak_conversion \
  --source h_5 --target a_1 \
  --grid-shape "${GRID}" --bounds "${BOUNDS}" --periodic \
  --ic plane-wave --ic-wavevector "${K0}" --ic-amplitude 0.1 --ic-component h_5 \
  --t-end "${T_END}" \
  --param "kappa=${KAPPA}" --param "B0=${B0}" \
  --param "mA2=${MA2}" --param "mT2=${MT2}" \
  --parallel "${TIDAL_PARALLEL:-4}" --resume \
  --output "${OUTPUT}"

echo ""
echo "--- Generating plots ---"

# 2D heatmap coloured by amplification A
tidal plot "${OUTPUT}" --type sweep \
  --metric A \
  --baseline-formula "sin(kappa * B0 * ${T_END} / 2)**2" \
  --title "Plasma dark photon: A(δ, ξ) at resonance mA²=mT²=${MA2}" \
  --output "${OUTPUT}/heatmap_A_dxi.png" --quiet

# Stability filter
python "${SCRIPT_DIR}/analyze_sweep.py" "${OUTPUT}"

echo ""
echo "=== Sweep B complete ==="
echo "Plot: ${OUTPUT}/heatmap_A_dxi.png  [amplification vs (δ, ξ)]"
echo "Ghost boundary: δ²=ξ/4 (analytical)"
