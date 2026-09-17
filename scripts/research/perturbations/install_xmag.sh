#!/usr/bin/env bash
# ==============================================================================
# install_xmag.sh -- additive install of xMAG and its dependency chain
# ==============================================================================
# R-C (#567), widened by the user at planning (Q4). Adds FIVE new directories,
# nothing else, each stamped with INSTALLED_COMMIT:
#   Applications/SymmetricFunctions/   } plain Mathematica packages,
#   Applications/BrauerAlgebra/        } from THelpin/xBrauer_Bundle
#   Applications/xAct/xBrauer/           from THelpin/xBrauer_Bundle
#   Applications/xAct/TraceFree/         from xAct-contrib/TraceFree (repo root)
#   Applications/xAct/xMAG/              from THelpin/xMAG (repo root)
# Sources are fetched at pinned commits (planning, 2026-09-15). Refuses if any
# target exists. No patch to any of the chain (Q4). No kernel is started unless
# the userbase cannot be found.
#
# Usage:
#   install_xmag.sh [--dest-root DIR] [--print-removal] [--help]
#   --dest-root DIR   the Applications directory (default: resolved userbase)
#   --print-removal   print (never run) the removal commands
# Exit: 0 installed and additivity proven; 1 any failure; 3 kernel live.
# ==============================================================================
set -euo pipefail
XBRAUER_REPO="${XBRAUER_REPO:-https://github.com/THelpin/xBrauer_Bundle}"
XBRAUER_COMMIT="${XBRAUER_COMMIT:-48be67e1c9037650f661c9d8ae538e807e407ba1}"
TRACEFREE_REPO="${TRACEFREE_REPO:-https://github.com/xAct-contrib/TraceFree}"
TRACEFREE_COMMIT="${TRACEFREE_COMMIT:-4e53ab3996f3be5da312576f1a24fb7be2ab5ddf}"
XMAG_REPO="${XMAG_REPO:-https://github.com/THelpin/xMAG}"
XMAG_COMMIT="${XMAG_COMMIT:-88026e47af1f8651f999775339e5c080d82d9d04}"
DEST_ROOT=""; PRINT_REMOVAL="false"
say() { printf '%s\n' "${1//$HOME/\~}"; }
die() { say "[ERROR] $1"; exit 1; }
show_help() { awk 'NR==1 && /^#!/ {next} /^#/ {sub(/^# ?/, ""); print; next} {exit}' "${BASH_SOURCE[0]}"; }
while [[ $# -gt 0 ]]; do
    case $1 in
        --dest-root) DEST_ROOT="$2"; shift 2 ;;
        --print-removal) PRINT_REMOVAL="true"; shift ;;
        --help|-h) show_help; exit 0 ;;
        *) die "unknown option: $1" ;;
    esac
done
guard() {
    if pgrep -x wolframscript >/dev/null || pgrep -x WolframKernel >/dev/null; then
        say "[ERROR] a Wolfram kernel is live; one kernel at a time (single license)"; exit 3
    fi
}
if [[ -z "$DEST_ROOT" ]]; then
    ub="${WOLFRAM_USERBASE:-$HOME/.local/wolfram/userbase}"
    if [[ ! -f "$ub/Applications/xAct/xTensor/xTensor.m" ]]; then
        guard; ub="$(wolframscript -code '$UserBaseDirectory' 2>/dev/null | tr -d '\r\n')"
        [[ -f "$ub/Applications/xAct/xTensor/xTensor.m" ]] || die "cannot locate the xAct userbase"
    fi
    DEST_ROOT="$ub/Applications"
fi
XACT="$DEST_ROOT/xAct"
say "XMAG_DEST_ROOT=$DEST_ROOT"
TARGETS=("$DEST_ROOT/SymmetricFunctions" "$DEST_ROOT/BrauerAlgebra" "$XACT/xBrauer" "$XACT/TraceFree" "$XACT/xMAG")
if [[ "$PRINT_REMOVAL" == "true" ]]; then
    for t in "${TARGETS[@]}"; do say "removal (NOT run): rm -rf '$t'   # only if $t/INSTALLED_COMMIT names this installer"; done
    exit 0
fi
[[ -d "$XACT" ]] || die "destination does not exist: $XACT"
for t in "${TARGETS[@]}"; do [[ -e "$t" ]] && die "refusing: $t already exists (additive only)"; done
command -v git >/dev/null || die "git is required"

WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT INT TERM
fetch_pinned() {  # name repo sha -> $WORK/name
    local name=$1 repo=$2 sha=$3 dir="$WORK/$1"
    git init --quiet "$dir"
    git -C "$dir" remote add origin "$repo"
    if git -C "$dir" fetch --quiet --depth 1 origin "$sha" 2>/dev/null; then
        git -C "$dir" checkout --quiet FETCH_HEAD
    else
        rm -rf "$dir"; git clone --quiet "$repo" "$dir" || die "clone failed: $repo"
        git -C "$dir" checkout --quiet "$sha" || die "commit not found: $sha in $repo"
    fi
    local got; got="$(git -C "$dir" rev-parse HEAD)"
    [[ "$got" == "$sha" ]] || die "$name resolved $got, expected $sha"
    say "XMAG_FETCHED $name@$got"
}
fetch_pinned bundle "$XBRAUER_REPO" "$XBRAUER_COMMIT"
fetch_pinned tracefree "$TRACEFREE_REPO" "$TRACEFREE_COMMIT"
fetch_pinned xmag "$XMAG_REPO" "$XMAG_COMMIT"

for f in bundle/SymmetricFunctions/SymmetricFunctions.m bundle/BrauerAlgebra/BrauerAlgebra.m bundle/xBrauer/xBrauer.m tracefree/TraceFree.m xmag/xMAG.m; do
    [[ -f "$WORK/$f" ]] || die "payload shape unexpected: $f missing"
done
version_of() { grep -m1 -oE '\$Version *= *\{"[0-9.]+", *\{[0-9, ]+\}\}' "$1" || echo "unknown"; }

before_root="$(find "$DEST_ROOT" -mindepth 1 -maxdepth 1 -printf '%f\n' | LC_ALL=C sort)"
before_xact="$(find "$XACT" -mindepth 1 -maxdepth 1 -printf '%f\n' | LC_ALL=C sort)"

copy_pkg() {  # src dest repo sha mainfile
    local src=$1 dest=$2 repo=$3 sha=$4 main=$5
    mkdir -p "$dest"
    ( cd "$src" && tar --exclude=.git -cf - . ) | ( cd "$dest" && tar -xf - )
    {
        echo "$sha"; echo "source_url=$repo"; echo "package_version=$(version_of "$dest/$main")"
        echo "installed_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"; echo "installed_by=scripts/research/perturbations/install_xmag.sh"; echo "patch=none"
    } > "$dest/INSTALLED_COMMIT"
    say "XMAG_INSTALLED $dest ($(version_of "$dest/$main"); $(find "$dest" -type f | wc -l) files)"
}
copy_pkg "$WORK/bundle/SymmetricFunctions" "$DEST_ROOT/SymmetricFunctions" "$XBRAUER_REPO" "$XBRAUER_COMMIT" SymmetricFunctions.m
copy_pkg "$WORK/bundle/BrauerAlgebra"      "$DEST_ROOT/BrauerAlgebra"      "$XBRAUER_REPO" "$XBRAUER_COMMIT" BrauerAlgebra.m
copy_pkg "$WORK/bundle/xBrauer"            "$XACT/xBrauer"                 "$XBRAUER_REPO" "$XBRAUER_COMMIT" xBrauer.m
copy_pkg "$WORK/tracefree"                 "$XACT/TraceFree"               "$TRACEFREE_REPO" "$TRACEFREE_COMMIT" TraceFree.m
copy_pkg "$WORK/xmag"                      "$XACT/xMAG"                    "$XMAG_REPO" "$XMAG_COMMIT" xMAG.m

after_root="$(find "$DEST_ROOT" -mindepth 1 -maxdepth 1 -printf '%f\n' | LC_ALL=C sort)"
after_xact="$(find "$XACT" -mindepth 1 -maxdepth 1 -printf '%f\n' | LC_ALL=C sort)"
added_root="$(comm -13 <(printf '%s\n' "$before_root") <(printf '%s\n' "$after_root") | tr '\n' ' ')"
added_xact="$(comm -13 <(printf '%s\n' "$before_xact") <(printf '%s\n' "$after_xact") | tr '\n' ' ')"
removed="$(comm -23 <(printf '%s\n' "$before_root") <(printf '%s\n' "$after_root"); comm -23 <(printf '%s\n' "$before_xact") <(printf '%s\n' "$after_xact"))"
if [[ "$added_root" == "BrauerAlgebra SymmetricFunctions " && "$added_xact" == "TraceFree xBrauer xMAG " && -z "$removed" ]]; then
    say "XMAG_ADDITIVITY_OK=Applications/{SymmetricFunctions,BrauerAlgebra} and Applications/xAct/{xBrauer,TraceFree,xMAG} added, nothing removed"
else
    die "additivity check failed: added_root=[$added_root] added_xact=[$added_xact] removed=[$removed]"
fi
say "XMAG_INSTALL_OK"
