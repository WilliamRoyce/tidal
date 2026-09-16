#!/usr/bin/env bash
# ==============================================================================
# manifest.sh -- two-layer sha256 manifest of the Wolfram userbase (R-C, #567)
# ==============================================================================
# Proves that the lane added exactly the intended package directories and
# modified nothing else. verify-wolfram-setup.sh fingerprints only four
# packages' version lines, so it cannot make that claim; this can.
#
# Layer (i)  -- everything that must be byte-identical before/after:
#               Applications/** minus the new directories, Kernel/**,
#               Licensing/**, Autoload/**, SystemFiles/**, .activation_backup/**
# Layer (ii) -- directories allowed to change (kernel-maintained caches and
#               cloud-session state), reported file by file:
#               Paclets/**, FrontEnd/**, SearchIndices/**,
#               ApplicationData/CloudObject/**, ApplicationData/Credentials/**
# New        -- the six directories the lane may add.
# Any path outside all three classes that differs is a FINDING.
#
# Usage:
#   manifest.sh before [OUTDIR]     snapshot -> OUTDIR/before.*
#   manifest.sh after  [OUTDIR]     snapshot -> OUTDIR/after.*
#   manifest.sh diff   [OUTDIR]     compare before/after, print verdict lines
# OUTDIR defaults to <repo>/third_party/perturbations_runs/manifest (gitignored).
# No kernel is started. Paths are printed relative to the userbase.
# ==============================================================================
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
UB="${WOLFRAM_USERBASE:-$HOME/.local/wolfram/userbase}"
MODE="${1:-}"
OUT="${2:-$REPO_ROOT/third_party/perturbations_runs/manifest}"
NEW_DIRS=(Applications/xAct/xPand Applications/SymmetricFunctions Applications/BrauerAlgebra
          Applications/xAct/xBrauer Applications/xAct/TraceFree Applications/xAct/xMAG)
ALLOW_DIRS=(Paclets FrontEnd SearchIndices ApplicationData/CloudObject ApplicationData/Credentials)

[[ -d "$UB/Applications/xAct/xTensor" ]] || { echo "MANIFEST_ERROR=userbase not found (set WOLFRAM_USERBASE)"; exit 1; }
mkdir -p "$OUT"

classify() {  # stdin: "hash  ./path" lines -> three files
    local tag=$1
    : > "$OUT/$tag.layer1.sha256"; : > "$OUT/$tag.layer2.sha256"; : > "$OUT/$tag.new.sha256"
    while IFS= read -r line; do
        local path="${line#*  ./}" cls=layer1 d
        for d in "${NEW_DIRS[@]}";   do [[ "$path" == "$d/"* ]] && cls=new && break; done
        if [[ $cls == layer1 ]]; then
            for d in "${ALLOW_DIRS[@]}"; do [[ "$path" == "$d/"* ]] && cls=layer2 && break; done
        fi
        printf '%s\n' "$line" >> "$OUT/$tag.$cls.sha256"
    done
}

snapshot() {
    local tag=$1
    ( cd "$UB" && find . -type f -print0 | LC_ALL=C sort -z | xargs -0 sha256sum ) > "$OUT/$tag.all.sha256"
    classify "$tag" < "$OUT/$tag.all.sha256"
    echo "MANIFEST_${tag^^}_FILES=$(wc -l < "$OUT/$tag.all.sha256")"
    echo "MANIFEST_${tag^^}_LAYER1_FILES=$(wc -l < "$OUT/$tag.layer1.sha256")"
    echo "MANIFEST_${tag^^}_LAYER2_FILES=$(wc -l < "$OUT/$tag.layer2.sha256")"
    echo "MANIFEST_${tag^^}_NEW_FILES=$(wc -l < "$OUT/$tag.new.sha256")"
    echo "MANIFEST_${tag^^}_LAYER1_SHA256=$(sha256sum "$OUT/$tag.layer1.sha256" | cut -d' ' -f1)"
    echo "MANIFEST_${tag^^}_UTC=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
}

case "$MODE" in
    before|after) snapshot "$MODE" ;;
    diff)
        for f in before after; do [[ -f "$OUT/$f.all.sha256" ]] || { echo "MANIFEST_ERROR=missing $f snapshot"; exit 1; }; done
        if cmp -s "$OUT/before.layer1.sha256" "$OUT/after.layer1.sha256"; then
            echo "MANIFEST_LAYER1=identical"
        else
            echo "MANIFEST_LAYER1=DIFFERENT (FINDING)"; diff "$OUT/before.layer1.sha256" "$OUT/after.layer1.sha256" | head -50
        fi
        echo "MANIFEST_LAYER2_CHANGED_FILES=$(diff "$OUT/before.layer2.sha256" "$OUT/after.layer2.sha256" | grep -c '^[<>]' || true)"
        # diff exits 1 when the files differ; under pipefail that would abort the script here
        { diff "$OUT/before.layer2.sha256" "$OUT/after.layer2.sha256" || true; } | grep '^[<>]' | sed 's/^\([<>]\) [0-9a-f]*  \.\//\1 /' | sort -u -k2 | head -80 || true
        echo "MANIFEST_NEW_BEFORE_FILES=$(wc -l < "$OUT/before.new.sha256")"
        echo "MANIFEST_NEW_AFTER_FILES=$(wc -l < "$OUT/after.new.sha256")"
        for d in "${NEW_DIRS[@]}"; do
            n=$(grep -c "  ./$d/" "$OUT/after.new.sha256" || true); echo "MANIFEST_NEW_DIR $d files=$n"
        done
        ;;
    *) echo "usage: $0 before|after|diff [OUTDIR]"; exit 2 ;;
esac
