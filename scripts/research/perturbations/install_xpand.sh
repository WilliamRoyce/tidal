#!/usr/bin/env bash
# cspell:words HDRS tolower
# ==============================================================================
# install_xpand.sh -- additive install of xPand 0.4.4 into the Wolfram userbase
# ==============================================================================
# R-C (#567). Places ONLY Applications/xAct/xPand/ under the userbase and writes
# an INSTALLED_VERSION stamp inside it. Refuses if the directory exists; never
# touches any other package (the certified xAct tree is read-only for this lane).
# Precedent: scripts/install-psalter.sh (stamp), scripts/research/interfaces/
# (route not payload). No kernel is started unless the userbase cannot be found.
#
# Source: Cyril Pitrou's page, tarball pinned by sha256 (planning, 2026-09-15).
# The archive is flat (xPand/ at the top) and carries macOS ._* / .DS_Store junk,
# which is stripped on extraction.
#
# Usage:
#   install_xpand.sh [--tarball FILE] [--dest DIR] [--wayback] [--print-removal]
#   --tarball FILE    use a local copy instead of downloading
#   --dest DIR        the Applications/xAct directory to install into
#                     (default: resolved userbase; used for mock rehearsals)
#   --wayback         request an archive.org snapshot of the tarball URL and
#                     record it (outward-facing; asked for in review round 1)
#   --print-removal   print (never run) the command that would remove the install
#   --help
# Exit: 0 installed and additivity proven; 1 any failure; 3 kernel live.
# ==============================================================================
set -euo pipefail
XPAND_VERSION="0.4.4"
XPAND_URL="${XPAND_URL:-http://www2.iap.fr/users/pitrou/xPand_0.4.4.tar.gz}"
XPAND_SHA256="${XPAND_SHA256:-26e7abcac7bb655235ec39b73850729cf4465748d0d5ea2ad03c0607aef5ceab}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
TARBALL=""; DEST=""; WAYBACK="false"; PRINT_REMOVAL="false"
say() { printf '%s\n' "${1//$HOME/\~}"; }
die() { say "[ERROR] $1"; exit 1; }
show_help() { awk 'NR==1 && /^#!/ {next} /^#/ {sub(/^# ?/, ""); print; next} {exit}' "${BASH_SOURCE[0]}"; }

while [[ $# -gt 0 ]]; do
    case $1 in
        --tarball) TARBALL="$2"; shift 2 ;;
        --dest) DEST="$2"; shift 2 ;;
        --wayback) WAYBACK="true"; shift ;;
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

resolve_dest() {
    if [[ -n "$DEST" ]]; then return; fi
    local ub="${WOLFRAM_USERBASE:-$HOME/.local/wolfram/userbase}"
    if [[ ! -f "$ub/Applications/xAct/xTensor/xTensor.m" ]]; then
        guard
        ub="$(wolframscript -code '$UserBaseDirectory' 2>/dev/null | tr -d '\r\n')"
        [[ -f "$ub/Applications/xAct/xTensor/xTensor.m" ]] || die "cannot locate the xAct userbase"
    fi
    DEST="$ub/Applications/xAct"
}

resolve_dest
say "XPAND_DEST=$DEST"
if [[ "$PRINT_REMOVAL" == "true" ]]; then
    say "removal (NOT run): rm -rf '$DEST/xPand'   # only if $DEST/xPand/INSTALLED_VERSION names this installer"
    exit 0
fi
[[ -d "$DEST" ]] || die "destination does not exist: $DEST"
[[ -e "$DEST/xPand" ]] && die "refusing: $DEST/xPand already exists (this installer is additive only)"

WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT INT TERM
if [[ -n "$TARBALL" ]]; then
    cp "$TARBALL" "$WORK/xPand.tar.gz"; say "XPAND_SOURCE=local copy"
else
    curl -fsSL --retry 3 -o "$WORK/xPand.tar.gz" "$XPAND_URL" || die "download failed: $XPAND_URL"
    say "XPAND_SOURCE=$XPAND_URL"
fi
HASH="$(sha256sum "$WORK/xPand.tar.gz" | cut -d' ' -f1)"
say "XPAND_TARBALL_SHA256=$HASH"
[[ "$HASH" == "$XPAND_SHA256" ]] || die "sha256 mismatch: expected $XPAND_SHA256"

mkdir -p "$WORK/x"
tar -xzf "$WORK/xPand.tar.gz" -C "$WORK/x" --exclude='._*' --exclude='.DS_Store' 2>/dev/null \
    || tar -xzf "$WORK/xPand.tar.gz" -C "$WORK/x" --exclude='._*' --exclude='.DS_Store' 2>&1 | grep -v 'Ignoring unknown extended header' || true
SRC="$WORK/x/xPand"
[[ -f "$SRC/xPand.m" && -f "$SRC/Kernel/init.m" ]] || die "payload shape unexpected: xPand/xPand.m or xPand/Kernel/init.m missing"
VERSION_LINE="$(grep -m1 'xAct`xPand`\$Version' "$SRC/xPand.m" || true)"
say "XPAND_VERSION_LINE=$VERSION_LINE"
[[ "$VERSION_LINE" == *"\"$XPAND_VERSION\""* ]] || die "version line does not say $XPAND_VERSION"

BEFORE="$(find "$DEST" -mindepth 1 -maxdepth 1 -printf '%f\n' | LC_ALL=C sort)"
cp -a "$SRC" "$DEST/xPand"
AFTER="$(find "$DEST" -mindepth 1 -maxdepth 1 -printf '%f\n' | LC_ALL=C sort)"
ADDED="$(comm -13 <(printf '%s\n' "$BEFORE") <(printf '%s\n' "$AFTER"))"
REMOVED="$(comm -23 <(printf '%s\n' "$BEFORE") <(printf '%s\n' "$AFTER"))"
if [[ "$ADDED" == "xPand" && -z "$REMOVED" ]]; then say "XPAND_ADDITIVITY_OK=only xPand/ added"; else die "additivity check failed: added=[$ADDED] removed=[$REMOVED]"; fi

SNAPSHOT="not requested"
if [[ "$WAYBACK" == "true" ]]; then
    # Save Page Now answers 302 with the snapshot in Location; anything else is a failure and
    # is recorded as such (never the save URL). From this container archive.org refuses the
    # TLS handshake ("alert access denied", 2026-09-16), so the request must be made elsewhere.
    HDRS="$(curl -sS --max-time 90 -D - -o /dev/null "https://web.archive.org/save/$XPAND_URL" 2>&1 || true)"
    LOC="$(printf '%s' "$HDRS" | tr -d '\r' | awk 'tolower($1)=="location:" {print $2; exit}')"
    if [[ "$LOC" == https://web.archive.org/web/* ]]; then SNAPSHOT="$LOC";
    else SNAPSHOT="request failed: $(printf '%s' "$HDRS" | head -1 | tr -d '\r' | cut -c1-120)"; fi
    say "XPAND_WAYBACK=$SNAPSHOT"
fi

{
    echo "$XPAND_VERSION"
    echo "source_url=$XPAND_URL"
    echo "tarball_sha256=$HASH"
    echo "version_line=$VERSION_LINE"
    echo "wayback=$SNAPSHOT"
    echo "installed_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "installed_by=scripts/research/perturbations/install_xpand.sh"
    echo "patch=none"
} > "$DEST/xPand/INSTALLED_VERSION"
say "XPAND_STAMP=$DEST/xPand/INSTALLED_VERSION"
say "XPAND_FILES=$(find "$DEST/xPand" -type f | wc -l)"
say "XPAND_INSTALL_OK"
