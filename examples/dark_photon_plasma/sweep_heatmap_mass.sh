#!/usr/bin/env bash
# Sweep A — 2D Heatmap: mA² vs mT² at fixed (δ=0.3, ξ=1)
#
# Produces a 50×50 Cartesian grid for clean heatmap visualisation of
# the amplification factor A = P(h→a₁)/P_GR where P_GR = sin²(κB₀t/2).
#
# Physics:
#   The photon effective mass mA² (plasma proxy) detunes the Gertsenshtein
#   h↔a channel. The dark-photon mass mT² controls the Raffelt-Stodolsky
#   oscillation frequency for a↔t mixing. At resonance mA² ≈ mT², the
#   photon and dark-photon dispersions match, giving maximal energy transfer.
#
#   Ghost-freedom condition: |δ| < √ξ/2. At δ=0.3, ξ=1 → |δ|=0.3 < 0.5 ✓
#
# Metric: A_total = P_max / sin²(κB₀t_end/2) — amplification relative
#   to pure vacuum Gertsenshtein baseline.
#
# Refs:
#   An, Pospelov, Pradler (2013), arXiv:1302.3884 — conversion formulas
#   Holdom (1986), Phys. Lett. B 166, 196 — kinetic mixing triviality

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SPEC="${SCRIPT_DIR}/../data/dark_photon_plasma.json"
OUTPUT="${SCRIPT_DIR}/../data/dp_plasma_heatmap_mass"

# Fixed physics
KAPPA=1.0
B0=0.01
XI=1.0
DELTAM=0.3
ALPHA3=0.15    # alpha3 = mT2/2; alpha3=0.15 → mass^2=0.30, near mA2 resonance
K0=2.0

# Grid and simulation settings
GRID=256
BOUNDS="0:100"
T_END=50.0

echo "=== Sweep A: mA² vs alpha3 heatmap (50×50) ==="
echo "Fixed: κ=${KAPPA}, B₀=${B0}, ξ=${XI}, δ=${DELTAM}, alpha3=${ALPHA3}"
echo "Note: mass^2 = 2*alpha3; alpha3 ∈ [0.005, 2.5] ≡ mT2 ∈ [0.01, 5] in FV model"
echo "Output: ${OUTPUT}"
echo ""

tidal sweep "${SPEC}" \
  --sweep "mA2=0.01:5.0:50" \
  --sweep "alpha3=0.005:2.5:50" \
  --measure peak_conversion \
  --source h_5 --target a_1 \
  --grid-shape "${GRID}" --bounds "${BOUNDS}" --periodic \
  --ic plane-wave --ic-wavevector "${K0}" --ic-amplitude 0.1 --ic-component h_5 \
  --t-end "${T_END}" \
  --param "kappa=${KAPPA}" --param "B0=${B0}" \
  --param "xi=${XI}" --param "deltam=${DELTAM}" \
  --parallel "${TIDAL_PARALLEL:-4}" --resume \
  --output "${OUTPUT}"

echo ""
echo "--- Generating plots ---"

# 2D heatmap coloured by amplification A
tidal plot "${OUTPUT}" --type sweep \
  --metric A \
  --baseline-formula "sin(kappa * B0 * ${T_END} / 2)**2" \
  --title "Plasma dark photon: amplification A(mA², alpha3) at δ=${DELTAM}, ξ=${XI}" \
  --output "${OUTPUT}/heatmap_A.png" --quiet

# Also plot raw P_max for comparison
tidal plot "${OUTPUT}" --type sweep \
  --metric P_max --log-scale \
  --title "Plasma dark photon: P_max(mA², alpha3)" \
  --output "${OUTPUT}/heatmap_Pmax.png" --quiet

# Stability filter
python "${SCRIPT_DIR}/analyze_sweep.py" "${OUTPUT}"

echo ""
echo "=== Sweep A complete ==="
echo "Plots:"
echo "  ${OUTPUT}/heatmap_A.png     [amplification heatmap]"
echo "  ${OUTPUT}/heatmap_Pmax.png  [raw P_max heatmap]"
