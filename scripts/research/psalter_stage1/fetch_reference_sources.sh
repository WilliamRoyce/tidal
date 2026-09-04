#!/usr/bin/env bash
# ==============================================================================
# Fetch the Barker reference sources that the H8 Stage-1 study cites
# ==============================================================================
# Downloads the specific upstream files that `docs/cosmology/stage1_engineering_plan.md`
# analyses, into a working directory that is NOT tracked by git.
#
# Why a fetch script rather than vendored copies: PSALTer and its supplemental
# materials are GPL-3.0-or-later, TIDAL is MIT. Committing those sources verbatim
# would be distribution -- the exact trigger recorded as a release blocker on
# GH #495. Reading and adapting them is authorized (decision D6, author permission
# explicit); redistributing them from this repository is not. So we commit the
# route, not the payload.
#
# Usage:
#   bash scripts/research/psalter_stage1/fetch_reference_sources.sh [OUTPUT_DIR]
#
# Default OUTPUT_DIR: <repo root>/third_party/psalter_reference (gitignored)
# ==============================================================================

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
OUT_DIR="${1:-${REPO_ROOT}/third_party/psalter_reference}"

# Pinned revisions -- verified 2026-09-03. Override via environment if a later
# revision is wanted, but re-verify the study's claims before trusting them.
PSALTER_COMMIT="${PSALTER_COMMIT:-bb45adb0fa21e467dbd88d4dc36ef21b84abbe6d}"  # v2.0.2
SM2506B_REF="${SM2506B_REF:-HEAD}"   # TorC / CTEG companion
SM2607_REF="${SM2607_REF:-HEAD}"     # numerical polology companion

log() { echo "[fetch-psalter-ref] $*"; }

mkdir -p "${OUT_DIR}"

# ---- PSALTer itself (the package we analyse; never installed by this script) --
if [[ ! -d "${OUT_DIR}/PSALTer/.git" ]]; then
    log "Cloning PSALTer at ${PSALTER_COMMIT}"
    git clone --quiet https://github.com/wevbarker/PSALTer "${OUT_DIR}/PSALTer"
    git -C "${OUT_DIR}/PSALTer" checkout --quiet "${PSALTER_COMMIT}"
else
    log "PSALTer checkout already present -- skipping"
fi

# ---- Individual reference files (raw download; the repos carry large blobs) ---
# SupplementalMaterials-2607 in particular is slow to clone in full.
fetch() {
    local url="$1" dest="$2"
    mkdir -p "$(dirname "${dest}")"
    if curl -sfL "${url}" -o "${dest}"; then
        log "  ok  ${dest#"${OUT_DIR}"/}"
    else
        log "  FAILED ${url}"
        return 1
    fi
}

SM6="https://raw.githubusercontent.com/wevbarker/SupplementalMaterials-2506b/${SM2506B_REF}"
SM7="https://raw.githubusercontent.com/wevbarker/SupplementalMaterials-2607/${SM2607_REF}"

log "Fetching SupplementalMaterials-2506b (TorC / CTEG)"
fetch "${SM6}/ParticleSpectroscopy/system-tests-Qtorsion/ParticleSpectrographCTEG.m" \
      "${OUT_DIR}/sm2506b/ParticleSpectrographCTEG.m"
fetch "${SM6}/ParticleSpectroscopy/ParticleSpectroscopy/PoincareGaugeTheory.m" \
      "${OUT_DIR}/sm2506b/PoincareGaugeTheory.m"
fetch "${SM6}/ParticleSpectroscopy/ParticleSpectroscopy/PoincareGaugeTheory/Linearise.m" \
      "${OUT_DIR}/sm2506b/Linearise.m"
fetch "${SM6}/ParticleSpectroscopy/ParticleSpectroscopy/PoincareGaugeTheory/LagrangianKarananasCouplings.m" \
      "${OUT_DIR}/sm2506b/LagrangianKarananasCouplings.m"
fetch "${SM6}/ParticleSpectroscopy/ParticleSpectroscopy/PoincareGaugeTheory/Models.m" \
      "${OUT_DIR}/sm2506b/Models.m"
# The committed result association -- a direct oracle for our own CTEG run.
fetch "${SM6}/ParticleSpectroscopy/ParticleSpectrographCTEG.mx" \
      "${OUT_DIR}/sm2506b/ParticleSpectrographCTEG.mx"

log "Fetching SupplementalMaterials-2607 (numerical polology)"
fetch "${SM7}/WolframLanguage/ParticleSpectroscopy.m" \
      "${OUT_DIR}/sm2607/ParticleSpectroscopy.m"
fetch "${SM7}/WolframLanguage/ParticleSpectroscopy/FieldKinematics.m" \
      "${OUT_DIR}/sm2607/FieldKinematics.m"
fetch "${SM7}/WolframLanguage/ParticleSpectroscopy/JuliaExport.m" \
      "${OUT_DIR}/sm2607/JuliaExport.m"
fetch "${SM7}/WolframLanguage/ParticleSpectroscopy/Models/VectorTheory.m" \
      "${OUT_DIR}/sm2607/VectorTheory.m"
fetch "${SM7}/WolframLanguage/ParticleSpectroscopy/Models/A23Theory.m" \
      "${OUT_DIR}/sm2607/A23Theory.m"
# The two small spectrograph exports -- candidate unit-test fixtures.
fetch "${SM7}/WolframLanguage/ParticleSpectrographVectorTheory.wxf" \
      "${OUT_DIR}/sm2607/ParticleSpectrographVectorTheory.wxf"
fetch "${SM7}/WolframLanguage/ParticleSpectrographA23Theory.wxf" \
      "${OUT_DIR}/sm2607/ParticleSpectrographA23Theory.wxf"

log "Done. Sources in ${OUT_DIR}"
log "These are GPL-3.0-or-later: read and adapt with provenance, do not commit verbatim."
