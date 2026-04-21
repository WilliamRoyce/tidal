#!/usr/bin/env bash
# B₀ Sensitivity Sweep: Plasma Dark Photon Model
#
# Scans B₀ at an on-resonance mixing configuration (mA² = mT²)
# to test whether the dark-photon amplification A_dark is B₀-independent.
#
# Physics:
#   At resonance mA² ≈ mT², the photon and dark-photon dispersions match,
#   maximizing the Raffelt-Stodolsky oscillation. The dark-photon coupling
#   E[t,a] = -2δ·mA²/(4δ²-ξ) has NO B₀ dependence, so A_dark should be
#   constant across B₀ values in the perturbative (small-P) regime.
#
#   Deviation from constancy would indicate:
#   - Higher-order B₀ effects (nonlinear backreaction)
#   - Breakdown of the linearized regime (P_GR approaching 1)
#   - Numerical artefacts from the modal solver
#
# Refs:
#   An, Pospelov, Pradler (2013), arXiv:1302.3884 — conversion formulas
#   Hwang & Noh (2310.04150) — linearized regime validity

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SPEC="${SCRIPT_DIR}/../data/dark_photon_plasma.json"
OUTPUT="${SCRIPT_DIR}/../data/dark_photon_plasma_B0"

# Fixed physics at resonance
KAPPA=1.0
MA2=0.5        # photon effective mass (plasma)
MT2=0.5        # dark-photon mass (= mA2 for resonance)
XI=1.0         # canonical kinetic
DELTAM=0.3     # kinetic mixing
K0=2.0         # GW frequency

# Grid and simulation settings
GRID=256
BOUNDS="0:100"
T_END=50.0

echo "=== Plasma Dark Photon: B₀ Sensitivity Sweep ==="
echo "Fixed: kappa=${KAPPA}, mA2=${MA2}, mT2=${MT2}, xi=${XI}, deltam=${DELTAM}"
echo "Sweep: B0 = 0.001 .. 0.015 (20 points)"
echo "At resonance: mA2 = mT2 = ${MA2} (matched dispersions)"
echo "Output: ${OUTPUT}"
echo ""

tidal sweep "${SPEC}" \
  --sweep "B0=0.001:0.015:20" \
  --measure peak_conversion \
  --source h_5 --target a_1 \
  --grid-shape "${GRID}" --bounds "${BOUNDS}" --periodic \
  --ic plane-wave --ic-wavevector "${K0}" --ic-amplitude 0.1 --ic-component h_5 \
  --t-end "${T_END}" \
  --param "kappa=${KAPPA}" --param "mA2=${MA2}" --param "mT2=${MT2}" \
  --param "xi=${XI}" --param "deltam=${DELTAM}" \
  --parallel "${TIDAL_PARALLEL:-4}" --resume \
  --output "${OUTPUT}"

echo ""
echo "--- Generating plots ---"

tidal plot "${OUTPUT}" --type sweep \
  --metric P_max \
  --title "Plasma dark photon: P_max vs B₀ (on resonance mA²=mT²=${MA2})" \
  --overlay "sin(kappa * B0 * ${T_END} / 2)**2" \
  --output "${OUTPUT}/sweep_B0.png" --quiet

echo ""
echo "=== Sweep complete ==="
echo "Results: ${OUTPUT}"
echo "Plot: ${OUTPUT}/sweep_B0.png"
echo ""
echo "Expected: A_dark ≈ const (B₀-independent dark-photon amplification)"
echo "  P_GR baseline: sin²(κ·B₀·t/2)"
echo "  Perturbative regime: B₀ ≤ 0.015 gives P_GR ≤ 0.14"
