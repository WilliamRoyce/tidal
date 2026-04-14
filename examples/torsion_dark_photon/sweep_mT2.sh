#!/usr/bin/env bash
# 1D mT2 sweep: dark-photon Raffelt-Stodolsky resonance scan
#
# The dark-photon-torsion theory has two independent Gertsenshtein channels:
#   1. h ↔ a (standard photon)            coupling kappa*B0
#   2. h ↔ t (dark photon, Proca mass mT2) coupling 2*kappa*B0*deltam
# These do not cross-mix at the modal eigendecomposition level: pure h IC
# at fixed k excites only channel 1 (photon Gertsenshtein), pure t IC
# excites channel 2 (dark Gertsenshtein). We probe channel 2 directly.
#
# This sweep uses a t_1 plane-wave IC and measures P(t_1 -> h_5), the
# dark-Gertsenshtein conversion. Varying the dark-photon mass mT2
# detunes the t-photon dispersion from omega = k, producing a
# Raffelt-Stodolsky-type resonance:
#
#   P_max = sin^2(2*theta) * sin^2(Delta_osc * D / 2)
#   tan(2*theta) = 2*kappa*B0*deltam / |Delta|
#   Delta = -mT2 / (2*omega)        (omega = k for plane wave)
#   Delta_osc = sqrt(Delta^2 + (2*kappa*B0*deltam)^2)
#
# At mT2 = 0 the channel is on resonance, P_max = sin^2(kappa*B0*deltam*t)
# Detuning |Delta| > coupling shrinks the conversion towards zero.
#
# Refs:
#   Holdom (1986), Phys. Lett. B 166, 196 — kinetic mixing
#   Raffelt & Stodolsky (1988), PRD 37, 1237 — two-state resonance
#   Fabbrichesi, Gabrielli, Lanfranchi (2020), arXiv:2005.01515 — review
#   Berlin et al. (2024), arXiv:2405.08865 — dark-photon mixing methods

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SPEC="${SCRIPT_DIR}/../data/torsion_dark_photon.json"
OUTPUT="${SCRIPT_DIR}/../data/torsion_dark_photon_sweep_mT2"

KAPPA=1.0
B0=0.01
K0=2.0
XI=1.0
DELTAM=0.3   # below critical kinetic mixing |deltam| < sqrt(xi)/2 = 0.5

GRID=256
BOUNDS="0:100"
T_END=50.0

echo "=== Dark Photon: mT2 Raffelt-Stodolsky resonance sweep ==="
echo "Fixed: kappa=${KAPPA}, B0=${B0}, xi=${XI}, deltam=${DELTAM}, t_end=${T_END}"
echo "Channel: t_1 -> h_5 (dark Gertsenshtein, coupling 2*kappa*B0*deltam)"
echo "Sweep: mT2 = 1e-6 .. 1.0 (60 log-spaced points)"
COUPLING=$(python3 -c "print(2*${KAPPA}*${B0}*${DELTAM})")
echo "Effective coupling: 2*kappa*B0*deltam = ${COUPLING}"
echo "On-resonance (mT2=0): P_max = sin^2(coupling * t_end / 2) ~ sin^2($(python3 -c "print(${COUPLING}*${T_END}/2)"))"
echo "Output: ${OUTPUT}"
echo ""

uv run tidal sweep "${SPEC}" \
  --sweep "mT2=1e-6:1.0:60:log" \
  --measure peak_conversion \
  --source t_1 --target h_5 \
  --grid-shape "${GRID}" --bounds "${BOUNDS}" --periodic \
  --ic plane-wave --ic-wavevector "${K0}" --ic-amplitude 0.1 --ic-component t_1 \
  --t-end "${T_END}" \
  --param "kappa=${KAPPA}" --param "B0=${B0}" \
  --param "xi=${XI}" --param "deltam=${DELTAM}" \
  --resume \
  --parallel 4 \
  --output "${OUTPUT}"

echo ""
echo "--- Generating plots ---"
uv run tidal plot "${OUTPUT}" --type sweep \
  --metric P_max \
  --title "Dark Gertsenshtein resonance (deltam=${DELTAM}, B0=${B0})" \
  --output "${OUTPUT}/sweep_mT2.png" --quiet

echo ""
echo "=== Sweep complete ==="
echo "Results: ${OUTPUT}"
echo "Plot:    ${OUTPUT}/sweep_mT2.png"
