#!/usr/bin/env bash
# cspell:words notfname
# ==============================================================================
# Verification Script for Wolfram Engine, xAct/xCoba and PSALTer
# ==============================================================================
# This script checks that all Wolfram-related components are properly installed
# and working. It should be run after container creation or to diagnose issues.
#
# Checks 7-10 cover PSALTer. They WARN when PSALTer is absent, so a container
# that only needs xAct still passes, and FAIL when it is present but broken.
# Pass --require-psalter to turn absence into a failure.
#
# Usage:
#   ./scripts/verify-wolfram-setup.sh [--require-psalter]
#
# Options:
#   --require-psalter  Treat a missing PSALTer install as a failure
#   --help, -h         Show this help
#
# Exit Codes:
#   0 - All checks passed
#   2 - Install works but a capability is degraded (not certifiable)
#   1 - One or more checks failed
#
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

ERRORS=0
WARNINGS=0
DEGRADED=0

# PSALTer checks are advisory unless --require-psalter is passed.
REQUIRE_PSALTER="false"
PSALTER_PRESENT="false"
PSALTER_LOAD_TIMEOUT="${PSALTER_LOAD_TIMEOUT:-900}"
PSALTER_SMOKE_TIMEOUT="${PSALTER_SMOKE_TIMEOUT:-900}"
PSALTER_COMMIT="${PSALTER_COMMIT:-bb45adb0fa21e467dbd88d4dc36ef21b84abbe6d}"
# The engine the Tier-1 gate was certified on (GH #543). A different engine is not a
# broken install -- the gate can still be run on it to reproduce or cross-check -- but
# the result is not certifiable, so a mismatch is reported as DEGRADED (exit 2).
# Override with EXPECTED_WOLFRAM_VERSION=<major.minor.release> when certifying anew.
EXPECTED_WOLFRAM_VERSION="${EXPECTED_WOLFRAM_VERSION:-14.3.0}"
# The certified engine lives on the devcontainer mount; a kernel reporting any other
# $InstallationDirectory is not the local engine (wolframscript falls back to CLOUD
# evaluation when the pinned kernel path is absent, and '1+1' passes there too).
EXPECTED_ENGINE_DIR="${EXPECTED_ENGINE_DIR:-${HOME}/.local/wolfram/engine/${EXPECTED_WOLFRAM_VERSION%.*}}"
# The certified xAct is the CODE, not the tarball label: the four package $Version strings
# measured on the certified install (xAct 1.3.0 bundle). The "1.2.1" that earlier documents
# carried was install-xact-xcoba.sh's default, never a measurement.
EXPECTED_XACT_FINGERPRINT="${EXPECTED_XACT_FINGERPRINT:-xCore=0.6.10 xPerm=1.2.4 xTensor=1.3.0 xCoba=0.8.6}"
# The two Function Repository resources PSALTer needs, as registered by
# scripts/psalter/register_resources.wl with FIXED identities (evidence/tier1-20260911-pass).
EXPECTED_RESOURCE_UUIDS="${EXPECTED_RESOURCE_UUIDS:-d40a8dd6-658c-47d2-8719-2f5fc8e1f83d 2f89f2e6-7bc8-4491-84bf-d8e13c69addb}"
# xPand carries the FRW derivation (D-C = O1', 2026-09-17), so it is part of the certified
# configuration: absent or different means runs are possible but not certifiable (DEGRADED).
# It is NOT in EXPECTED_XACT_FINGERPRINT, which names the four-package xAct bundle that many
# documents quote; this is its own check with its own expected version.
EXPECTED_XPAND_VERSION="${EXPECTED_XPAND_VERSION:-0.4.4}"

log_pass() {
    echo -e "${GREEN}[PASS]${NC} $1"
}

log_fail() {
    echo -e "${RED}[FAIL]${NC} $1"
    ((ERRORS++)) || true
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
    ((WARNINGS++)) || true
}

log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

# A capability that is genuinely absent but does NOT mean the install is broken.
# Distinct from log_soft_fail on purpose (found 2026-09-09): PSALTer's two Wolfram
# Function Repository dependencies cannot be fetched in this container at all, and
# routing them through log_soft_fail made them hard failures under --require-psalter
# -- which is what the Tier-1 gate calls, so the gate refused to start and became
# permanently impossible to run one minute after it last ran. A degradation is reported
# loudly and carried in the exit code (2), but it does not claim the install cannot
# be exercised. Anything demanding a certifiable install still refuses on 2.
log_degraded() {
    echo -e "${YELLOW}[DEGRADED]${NC} $1"
    ((DEGRADED++)) || true
}

# A missing PSALTer install is a warning by default and a failure under
# --require-psalter: the script runs at container creation, before PSALTer
# exists, but the install gate needs it to be hard.
log_soft_fail() {
    if [[ "$REQUIRE_PSALTER" == "true" ]]; then
        log_fail "$1"
    else
        log_warn "$1"
    fi
}

# Check 1: Wolfram Engine binary available
check_wolfram_binary() {
    log_info "Checking Wolfram Engine installation..."
    
    # The mounted kernel first, and on its own line: `command -v wolframscript`
    # is satisfied by the image's own client (/usr/bin/wolframscript ->
    # /opt/Wolfram/...) with no engine installed, after which every later check
    # runs against a CLOUD evaluation -- slow, and confusing for exactly the
    # fresh user this script is supposed to orient (#559).
    local kernel="${EXPECTED_ENGINE_DIR}/Executables/WolframKernel"
    if [[ ! -x "$kernel" ]]; then
        log_fail "No Wolfram kernel at ${kernel}"
        log_fail "  The engine is not installed in its mount. Setup path:"
        log_fail "  .devcontainer/docs/WOLFRAM_GUIDE.md"
        return 1
    fi
    log_pass "Engine kernel present: ${kernel}"

    if command -v wolframscript &> /dev/null; then
        local path
        path=$(command -v wolframscript)
        log_pass "wolframscript found at: ${path}"
        return 0
    else
        log_fail "wolframscript not found in PATH"
        return 1
    fi
}

# Check 2: Wolfram Engine activation
check_wolfram_activation() {
    log_info "Checking Wolfram Engine activation..."
    
    # Compare against stdout only, and match a LINE rather than the whole
    # capture: Wolfram Engine intermittently segfaults while shutting down, after
    # producing correct results, and folding that stderr text into the compared
    # value made this check fail at random on a perfectly good install.
    local result
    if result=$(wolframscript -code '1+1' 2>/dev/null); then
        if echo "$result" | grep -qx "2"; then
            log_pass "Wolfram Engine is activated and responding correctly"
            
            # Get version info
            local version
            version=$(wolframscript -code 'ToString[$VersionNumber] <> "." <> ToString[$ReleaseNumber]' 2>/dev/null | grep -E '^[0-9]+\.[0-9]+\.[0-9]+$' | head -1 || true)
            version="${version:-unknown}"
            log_info "  Version: ${version}"
            if [[ "$version" == "$EXPECTED_WOLFRAM_VERSION" ]]; then
                log_pass "  Engine matches the certified version (${EXPECTED_WOLFRAM_VERSION})"
            else
                log_degraded "  Engine ${version} is not the certified ${EXPECTED_WOLFRAM_VERSION} (GH #543): runs are possible, not certifiable"
            fi
            
            local license
            license=$(wolframscript -code '$LicenseType' 2>/dev/null || echo "unknown")
            log_info "  License: ${license}"

            # Provenance of the kernel itself. Without this, an absent local engine is
            # invisible: wolframscript evaluates '1+1' in the cloud and the check passes.
            local instdir
            instdir=$(wolframscript -code '$InstallationDirectory' 2>/dev/null | tr -d '\r' | tail -1)
            if [[ "$instdir" == "$EXPECTED_ENGINE_DIR"* ]]; then
                log_pass "  Kernel is the local engine at ${instdir}"
            else
                log_fail "  Kernel is NOT the local engine: \$InstallationDirectory=${instdir:-<none>} (expected under ${EXPECTED_ENGINE_DIR})"
                log_info "  A cloud-evaluated or foreign kernel cannot certify anything. Check WOLFRAMSCRIPT_KERNELPATH and the engine mount."
                return 1
            fi
            return 0
        else
            log_fail "Wolfram Engine returned unexpected result: ${result}"
            return 1
        fi
    else
        log_fail "Wolfram Engine not activated or error occurred"
        log_info "  Run: ./scripts/activate-wolfram.sh"
        return 1
    fi
}

# Check 3: User Applications directory exists
check_applications_dir() {
    log_info "Checking Wolfram Applications directory..."
    
    local user_dir
    user_dir=$(wolframscript -code '$UserBaseDirectory' 2>/dev/null | tr -d '\n\r')
    local apps_dir="${user_dir}/Applications"
    
    if [[ -d "$apps_dir" ]]; then
        log_pass "Applications directory exists: ${apps_dir}"
        return 0
    else
        log_fail "Applications directory not found: ${apps_dir}"
        return 1
    fi
}

# Check 4: xAct installation
check_xact_installed() {
    log_info "Checking xAct installation..."
    
    local user_dir
    user_dir=$(wolframscript -code '$UserBaseDirectory' 2>/dev/null | tr -d '\n\r')
    local xact_dir="${user_dir}/Applications/xAct"
    
    if [[ -d "$xact_dir" ]]; then
        log_pass "xAct directory found: ${xact_dir}"
        
        # Check for key packages
        local packages=("xCore" "xPerm" "xTensor" "xCoba")
        for pkg in "${packages[@]}"; do
            if [[ -d "${xact_dir}/${pkg}" ]]; then
                log_pass "  Package ${pkg} found"
            else
                log_fail "  Package ${pkg} missing"
            fi
        done

        # The certified xAct, asserted from the code: each package's own $Version line,
        # read from its .m header (no kernel needed). Degraded, not failed -- the install
        # runs, but a different bundle is not the one the Tier-1 gate was passed on.
        local fp="" pkg ver
        for pkg in "${packages[@]}"; do
            ver=$(grep -m1 -oE '\$Version *= *\{"[0-9.]+"' "${xact_dir}/${pkg}/${pkg}.m" 2>/dev/null | grep -oE '[0-9.]+' | head -1)
            fp="${fp}${fp:+ }${pkg}=${ver:-missing}"
        done
        if [[ "$fp" == "$EXPECTED_XACT_FINGERPRINT" ]]; then
            log_pass "  xAct packages match the certified bundle (${fp})"
        else
            log_degraded "  xAct packages (${fp}) are not the certified set (${EXPECTED_XACT_FINGERPRINT}): runs are possible, not certifiable"
        fi
        return 0
    else
        log_fail "xAct directory not found: ${xact_dir}"
        log_info "  Run: ./scripts/install-xact-xcoba.sh"
        return 1
    fi
}

# Check 5: xPerm binary compatibility
check_xperm_binary() {
    log_info "Checking xPerm binary compatibility..."
    
    local user_dir
    user_dir=$(wolframscript -code '$UserBaseDirectory' 2>/dev/null | tr -d '\n\r')
    local xperm_binary="${user_dir}/Applications/xAct/xPerm/mathlink/xperm.linux.64-bit"
    
    if [[ ! -f "$xperm_binary" ]]; then
        log_fail "xPerm binary not found: ${xperm_binary}"
        return 1
    fi
    
    # Check for GLIBC compatibility
    local ldd_output
    if ldd_output=$(ldd "$xperm_binary" 2>&1); then
        if echo "$ldd_output" | grep -q "not found"; then
            log_fail "xPerm binary has missing dependencies"
            log_info "  Run: ./scripts/install-xact-xcoba.sh to recompile"
            return 1
        else
            log_pass "xPerm binary dependencies satisfied"
            return 0
        fi
    else
        log_warn "Could not check xPerm binary dependencies"
        return 0
    fi
}

# Check 6: Load xAct packages
check_xact_loads() {
    log_info "Checking xAct package loading..."
    
    local test_code='
    Quiet[
      Needs["xAct`xCore`"];
      Needs["xAct`xPerm`"];
      Needs["xAct`xTensor`"];
      Needs["xAct`xCoba`"];
    , {General::stop, UpSetDelayed::write, SetDelayed::write}];
    Print["PACKAGES_LOADED"];
    '
    
    local result
    if result=$(timeout 60 wolframscript -code "$test_code" 2>&1); then
        if echo "$result" | grep -q "PACKAGES_LOADED"; then
            log_pass "All xAct packages load successfully"
            
            # Check for connection to xPerm external executable
            if echo "$result" | grep -q "Connection established"; then
                log_pass "  xPerm external executable connected"
            elif echo "$result" | grep -q "GLIBC"; then
                log_fail "  xPerm has GLIBC compatibility issues"
                return 1
            fi
            return 0
        else
            log_fail "xAct packages failed to load"
            log_info "  Output: ${result}"
            return 1
        fi
    else
        log_fail "Timeout or error loading xAct packages"
        return 1
    fi
}

# Check 7: PSALTer installed
check_xpand_installed() {
    log_info "Checking xPand installation (the FRW derivation engine, D-C = O1')..."

    local user_dir xpand_m ver
    user_dir=$(wolframscript -code '$UserBaseDirectory' 2>/dev/null | tr -d '\n\r')
    xpand_m="${user_dir}/Applications/xAct/xPand/xPand.m"

    if [[ ! -f "$xpand_m" ]]; then
        log_degraded "xPand not installed at ${user_dir}/Applications/xAct/xPand: FRW derivations cannot run"
        log_info "  Run: ./scripts/install-xpand.sh"
        return 0
    fi
    # Read the version from the package header; no kernel needed (the same tactic as xAct's
    # fingerprint), so this check costs nothing and cannot be confused by a loaded session.
    ver=$(grep -m1 -oE 'xAct`xPand`\$Version *= *\{"[0-9.]+"' "$xpand_m" | grep -oE '[0-9.]+' | head -1)
    if [[ "$ver" == "$EXPECTED_XPAND_VERSION" ]]; then
        log_pass "xPand ${ver} found (certified)"
    else
        log_degraded "xPand ${ver:-unreadable} is not the certified ${EXPECTED_XPAND_VERSION}: runs are possible, not certifiable"
    fi
    return 0
}

check_psalter_installed() {
    log_info "Checking PSALTer installation..."

    local user_dir
    user_dir=$(wolframscript -code '$UserBaseDirectory' 2>/dev/null | tr -d '\n\r')
    local psalter_dir="${user_dir}/Applications/xAct/PSALTer"

    if [[ ! -d "$psalter_dir" ]]; then
        log_soft_fail "PSALTer not installed at ${psalter_dir}"
        log_info "  Run: ./scripts/install-psalter.sh"
        return 0
    fi

    local f
    for f in PSALTer.m Sources/ParticleSpectrum.m Sources/DefField.m; do
        if [[ ! -f "${psalter_dir}/${f}" ]]; then
            log_fail "  PSALTer install is incomplete: missing ${f}"
            return 1
        fi
    done
    log_pass "PSALTer found at ${psalter_dir}"
    PSALTER_PRESENT="true"

    # The Stage-1 exporter reads PSALTer's private symbols, so which revision is
    # installed is a correctness input rather than bookkeeping.
    local marker="${psalter_dir}/INSTALLED_COMMIT"
    if [[ ! -f "$marker" ]]; then
        log_soft_fail "  INSTALLED_COMMIT missing -- installed revision is unknown"
        log_info "  Re-run: ./scripts/install-psalter.sh --force"
        return 0
    fi

    local sha
    sha=$(head -1 "$marker" | tr -d '\n\r')
    if [[ ! "$sha" =~ ^[0-9a-f]{40}$ ]]; then
        log_fail "  INSTALLED_COMMIT does not start with a commit hash"
        return 1
    fi
    log_pass "  Installed revision: ${sha}"
    if [[ "$sha" != "$PSALTER_COMMIT" ]]; then
        log_warn "  Installed revision differs from the pin (${PSALTER_COMMIT})"
    fi
    return 0
}

# Check 8: PSALTer loads
check_psalter_loads() {
    log_info "Checking PSALTer package loading..."

    if [[ "$PSALTER_PRESENT" != "true" ]]; then
        log_info "  Skipped (PSALTer not installed)"
        return 0
    fi

    # Existence is tested with DownValues, not ValueQ: PSALTer defines its
    # functions through the StackSetDelayed wrapper, so they carry definition
    # rules rather than a value.
    local test_code='
    Needs["xAct`PSALTer`"];
    Print["PSALTER_VERSION=", xAct`PSALTer`Private`$Version[[1]]];
    Print["PSALTER_SYMBOLS_OK=",
      TrueQ[Length[DownValues[xAct`PSALTer`DefField]] > 0 &&
            Length[DownValues[xAct`PSALTer`ParticleSpectrum]] > 0]];
    Print["PSALTER_LOADED"];
    '

    local tmp result rc
    tmp=$(mktemp -d)
    set +e
    result=$(cd "$tmp" && QT_QPA_PLATFORM=offscreen \
        timeout "$PSALTER_LOAD_TIMEOUT" wolframscript -code "$test_code" 2>&1)
    rc=$?
    set -e
    rm -rf "$tmp"

    if [[ $rc -eq 124 ]]; then
        log_fail "PSALTer load timed out after ${PSALTER_LOAD_TIMEOUT}s"
        log_info "  The first load builds projection-operator tables and is slow;"
        log_info "  raise the budget with PSALTER_LOAD_TIMEOUT=<seconds>"
        return 1
    fi

    # Judged from the output, not the exit status: Wolfram Engine intermittently
    # segfaults during kernel shutdown, after producing correct results.
    if echo "$result" | grep -q "PSALTER_LOADED" &&
       echo "$result" | grep -q "PSALTER_SYMBOLS_OK=True"; then
        log_pass "PSALTer loads and exports DefField and ParticleSpectrum"
        local ver
        ver=$(echo "$result" | sed -n 's/^PSALTER_VERSION=//p' | head -1)
        [[ -n "$ver" ]] && log_pass "  Package version: ${ver}"
        return 0
    fi

    log_fail "PSALTer failed to load"
    echo "$result" | tail -10
    return 1
}

# Check 9: PSALTer DefField smoke test (headless PDF export)
#
# This is the check that catches the headless hang. DefField calls SummariseField
# (Sources/DefField.m:91), which exports a FieldKinematics<Field>.pdf through the
# Wolfram front end (Sources/DefField/SummariseField.m:85) with no guard and no
# time limit. The front end needs a Qt platform plugin whose libraries are all
# present; when none can be initialized, Qt aborts and the export call BLOCKS
# INDEFINITELY rather than failing. So a hang here is the diagnostic, and must be
# reported differently from a failure.
check_psalter_deffield_smoke() {
    log_info "Running PSALTer DefField smoke test (headless PDF export)..."

    if [[ "$PSALTER_PRESENT" != "true" ]]; then
        log_info "  Skipped (PSALTer not installed)"
        return 0
    fi

    if [[ "${QT_QPA_PLATFORM:-}" != "offscreen" ]]; then
        log_info "  Note: QT_QPA_PLATFORM is not 'offscreen' in this shell;"
        log_info "  this check sets it itself, but your own PSALTer runs must too"
    fi

    local smoke_script="${SCRIPT_DIR}/psalter_smoke.wl"
    if [[ ! -f "$smoke_script" ]]; then
        log_warn "PSALTer smoke test script not found: ${smoke_script}"
        return 0
    fi

    local tmp result rc
    tmp=$(mktemp -d)
    set +e
    result=$(cd "$tmp" && QT_QPA_PLATFORM=offscreen \
        timeout "$PSALTER_SMOKE_TIMEOUT" wolframscript -file "$smoke_script" 2>&1)
    rc=$?
    set -e
    rm -rf "$tmp"

    if [[ $rc -eq 124 ]]; then
        log_fail "PSALTer DefField smoke test HUNG (the front-end export blocked)"
        log_info "  Fix: export QT_QPA_PLATFORM=offscreen"
        log_info "  Cause: DefField exports a PDF unconditionally"
        log_info "  (Sources/DefField.m:91 -> Sources/DefField/SummariseField.m:85)"
        return 1
    fi

    if echo "$result" | grep -q "PSALTER SMOKE TEST PASSED" &&
       echo "$result" | grep -q "PSALTER_DEFFIELD_PDF=1"; then
        log_pass "PSALTer smoke test passed (DefField exported its PDF headlessly)"
        return 0
    fi

    log_fail "PSALTer smoke test did not complete successfully"
    log_info "  Last lines of output:"
    echo "$result" | tail -10
    return 1
}

# Check 10: PSALTer's Wolfram Function Repository dependencies
#
# PSALTer calls ResourceFunction["PolynomialDegree"] and
# ResourceFunction["LinearlyIndependent"] -- five call sites across
# ValidateLagrangian, UnresolvedPoleRow and the source-constraint null-space code
# (the latter on SUBKERNELS). Neither is declared anywhere in its README or install
# instructions. When either does not resolve, PSALTer does not stop: the result is
# non-Boolean, no null vector is ever appended, every source constraint is lost and
# every pseudo-determinant comes out zero -- the Tier-1 MISMATCH of GH #543 -- while
# the run completes and writes its .mx. That is why this check tests BEHAVIOR
# (values on the exact call shapes), on the master and on a fresh subkernel, and
# not merely whether a ResourceObject can be looked up.
#
# The genuine functions are supplied by scripts/psalter/register_resources.wl
# (GH #543); they are part of the certified configuration, so their absence is a
# hard failure under --require-psalter (the gate would otherwise produce a
# known-wrong answer) and a warning otherwise. This reverses the DEGRADED routing
# adopted for GH #549, when the resources were merely unavailable rather than a
# certified leg; the fix is one command, named in the message.
check_psalter_resources() {
    log_info "Checking PSALTer's Wolfram Function Repository dependencies (behavior, master + subkernel)..."

    if [[ "$PSALTER_PRESENT" != "true" ]]; then
        log_info "  Skipped (PSALTer not installed)"
        return 0
    fi

    local test_code='
    ok[expr_] := Quiet@Check[expr, $Failed];
    shapes := {ResourceFunction["LinearlyIndependent"][{{1, 0}, {0, 1}}],
               ResourceFunction["LinearlyIndependent"][{{1, 0}, {2, 0}}],
               ResourceFunction["PolynomialDegree"][rx^2 ry, {rx, ry}]};
    Print["RESOURCE_MASTER=", ToString[ok[shapes], InputForm]];
    Print["RESOURCE_UUIDS=", ToString[ok[{ResourceObject["LinearlyIndependent"]["UUID"], ResourceObject["PolynomialDegree"]["UUID"]}], InputForm]];
    Quiet@LaunchKernels[1];
    Print["RESOURCE_SUBKERNEL=", ToString[ok[First[ParallelEvaluate[{ResourceFunction["LinearlyIndependent"][{{1, 0}, {0, 1}}], ResourceFunction["LinearlyIndependent"][{{1, 0}, {2, 0}}], ResourceFunction["PolynomialDegree"][rx^2 ry, {rx, ry}]}]]], InputForm]];
    Print["RESOURCE_SUBKERNEL_UUIDS=", ToString[ok[First[ParallelEvaluate[{ResourceObject["LinearlyIndependent"]["UUID"], ResourceObject["PolynomialDegree"]["UUID"]}]]], InputForm]];
    CloseKernels[];
    Print["RESOURCE_CHECK_DONE"];
    '
    local tmp result
    tmp=$(mktemp -d)
    set +e
    result=$(cd "$tmp" && timeout 300 wolframscript -code "$test_code" 2>&1)
    set -e
    rm -rf "$tmp"

    if ! echo "$result" | grep -q "RESOURCE_CHECK_DONE"; then
        log_warn "Could not check the resource functions"
        return 0
    fi

    local expected='{True, False, 3}' bad=0 side
    for side in MASTER SUBKERNEL; do
        if echo "$result" | grep -qF "RESOURCE_${side}=${expected}"; then
            log_pass "  LinearlyIndependent / PolynomialDegree behave on the ${side,,}"
        else
            log_soft_fail "  LinearlyIndependent / PolynomialDegree do NOT behave on the ${side,,}: $(echo "$result" | grep -o "RESOURCE_${side}=.*" | head -1)"
            bad=1
        fi
    done
    # Provenance, not only behavior: a function fetched from the repository behaves
    # identically to the registered one and would pass the block above. The registered
    # objects carry FIXED identities (register_resources.wl), so identity is checkable
    # on both the master and a subkernel. Anything else means the local registry is
    # absent, corrupted, or silently unregistered -- the resolver drops a name whose
    # object fails ResourceObjectQ without a message (FindResource.m:76-83).
    local want_uuids uuids_master uuids_sub
    want_uuids="{\"${EXPECTED_RESOURCE_UUIDS%% *}\", \"${EXPECTED_RESOURCE_UUIDS##* }\"}"
    uuids_master=$(echo "$result" | grep -o 'RESOURCE_UUIDS=.*' | head -1 | cut -d= -f2-)
    uuids_sub=$(echo "$result" | grep -o 'RESOURCE_SUBKERNEL_UUIDS=.*' | head -1 | cut -d= -f2-)
    if [[ "$uuids_master" == "$want_uuids" ]]; then
        log_pass "  Resolved to the certified identities on the master"
    else
        log_soft_fail "  Master resolves to ${uuids_master:-<none>}, not the certified ${want_uuids}"
        bad=1
    fi
    if [[ "$uuids_sub" == "$want_uuids" ]]; then
        log_pass "  Resolved to the certified identities on the subkernel"
    else
        log_soft_fail "  Subkernel resolves to ${uuids_sub:-<none>}, not the certified ${want_uuids}"
        bad=1
    fi

    if [[ $bad -eq 1 ]]; then
        log_info "  PSALTer continues WITHOUT these and silently produces a wrong spectrum"
        log_info "  (GH #543). Register the genuine code with:"
        log_info "    wolframscript -file scripts/psalter/register_resources.wl"
        return 1
    fi
    return 0
}

# Check 11: Run smoke test
check_smoke_test() {
    log_info "Running xAct/xCoba smoke test..."
    
    local smoke_script="${SCRIPT_DIR}/xact_smoke.wl"
    
    if [[ ! -f "$smoke_script" ]]; then
        log_warn "Smoke test script not found: ${smoke_script}"
        return 0
    fi
    
    local result
    if result=$(timeout 120 wolframscript -file "$smoke_script" 2>&1); then
        if echo "$result" | grep -q "SMOKE TEST PASSED"; then
            log_pass "Smoke test passed"
            return 0
        else
            log_fail "Smoke test did not complete successfully"
            log_info "  Last lines of output:"
            echo "$result" | tail -10
            return 1
        fi
    else
        log_fail "Smoke test timed out or encountered an error"
        return 1
    fi
}

# Main verification function
main() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            --require-psalter) REQUIRE_PSALTER="true"; shift ;;
            --help|-h)
                awk 'NR==1 && /^#!/ {next} /^#/ {sub(/^# ?/, ""); print; next} {exit}' \
                    "${BASH_SOURCE[0]}"
                exit 0
                ;;
            *)
                log_fail "Unknown option: $1"
                exit 1
                ;;
        esac
    done

    echo ""
    echo "========================================"
    echo "Wolfram Engine, xAct/xCoba & PSALTer Verification"
    echo "========================================"
    echo ""
    
    # Run all checks.
    #
    # `|| true` on every check is deliberate -- one failure must not hide the
    # other nine -- but check 1 is the exception: if the mounted engine is not
    # there, every later check would start `wolframscript` anyway, and with the
    # pinned kernel missing that is a CLOUD evaluation. Slow, and it starts a
    # kernel the caller did not ask for. I-ONB hit exactly that while probing a
    # redirect with EXPECTED_ENGINE_DIR pointed at an empty path (#559): check 1
    # refused, and check 2 then ran a kernel against the real engine anyway.
    #
    # So this one check short-circuits: no engine, no kernel, report and stop.
    if ! check_wolfram_binary; then
        echo ""
        log_fail "Stopping: every remaining check needs the mounted engine, and"
        log_fail "  without it wolframscript evaluates in the CLOUD instead."
        log_fail "  Setup path: .devcontainer/docs/WOLFRAM_GUIDE.md"
        echo ""
        exit 1
    fi
    echo ""
    
    check_wolfram_activation || true
    echo ""
    
    check_applications_dir || true
    echo ""
    
    check_xact_installed || true
    echo ""
    
    check_xperm_binary || true
    echo ""
    
    check_xact_loads || true
    echo ""
    
    check_xpand_installed || true
    echo ""
    
    check_psalter_installed || true
    echo ""
    
    check_psalter_loads || true
    echo ""
    
    check_psalter_deffield_smoke || true
    echo ""
    
    check_psalter_resources || true
    echo ""
    
    check_smoke_test || true
    echo ""
    
    # Summary
    echo "========================================"
    echo "Verification Summary"
    echo "========================================"
    
    if [[ $ERRORS -eq 0 ]]; then
        if [[ $DEGRADED -eq 0 ]]; then
            log_pass "All checks passed!"
        else
            log_pass "All install checks passed"
            log_warn "${DEGRADED} capability degradation(s): the install works but is not certifiable"
        fi
        if [[ $WARNINGS -gt 0 ]]; then
            log_warn "${WARNINGS} warning(s) noted"
        fi
        # Exit 2 == usable but degraded. Callers that only need to exercise the
        # install (the Tier-1 gate) accept it; callers that need a certifiable
        # install treat any non-zero as refusal.
        [[ $DEGRADED -eq 0 ]] && exit 0 || exit 2
    else
        log_fail "${ERRORS} check(s) failed"
        if [[ $WARNINGS -gt 0 ]]; then
            log_warn "${WARNINGS} warning(s) noted"
        fi
        echo ""
        log_info "To fix issues, try:"
        log_info "  1. sudo ./scripts/install-wolfram-engine.sh"
        log_info "  2. ./scripts/activate-wolfram.sh"
        log_info "  3. ./scripts/install-xact-xcoba.sh"
        log_info "  4. ./scripts/install-psalter.sh"
        log_info "  5. export QT_QPA_PLATFORM=offscreen   # PSALTer exports PDFs from DefField"
        exit 1
    fi
}

# Run main function
main "$@"
