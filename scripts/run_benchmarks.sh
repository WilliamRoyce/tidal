#!/usr/bin/env bash
# Run all benchmark scripts listed in scripts/publication_manifest.yaml.
#
# Default: writes to benchmark_results/<host>-<date>/  (gitignored).
# --canonical: writes to benchmark_results/canonical/  (committed).
#
# The --canonical flag requires a clean git state so accidental local
# tweaks don't slip into the headline numbers cited in App C.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MANIFEST="${ROOT}/scripts/publication_manifest.yaml"

CANONICAL=0
for arg in "$@"; do
    case "$arg" in
        --canonical) CANONICAL=1 ;;
        -h|--help)
            sed -n '2,11p' "$0"
            exit 0
            ;;
        *)
            echo "unknown argument: $arg" >&2
            exit 2
            ;;
    esac
done

if [[ "$CANONICAL" -eq 1 ]]; then
    if [[ -n "$(git -C "$ROOT" status --porcelain)" ]]; then
        echo "refusing to write --canonical with dirty git state" >&2
        echo "commit or stash local changes first" >&2
        exit 1
    fi
    OUT_DIR="${ROOT}/benchmark_results/canonical"
else
    HOST="$(hostname -s 2>/dev/null || echo unknown)"
    DATE="$(date +%Y-%m-%d)"
    OUT_DIR="${ROOT}/benchmark_results/${HOST}-${DATE}"
fi

mkdir -p "$OUT_DIR"
echo "writing benchmark output to: $OUT_DIR"

# Extract benchmark_script paths from the manifest. The grep+sed pair is
# a minimal YAML reader for this fixed-shape file; no PyYAML dep at the
# orchestrator layer.
BENCHMARKS=$(
    grep -E '^\s+benchmark_script:' "$MANIFEST" \
    | awk -F': ' '{print $2}' \
    | sort -u
)

if [[ -z "$BENCHMARKS" ]]; then
    echo "no benchmark_script entries found in manifest" >&2
    exit 1
fi

FAILED=()
for script in $BENCHMARKS; do
    name="$(basename "$script" .py)"
    out="${OUT_DIR}/${name}.json"
    echo "==> $name"
    if uv run python "${ROOT}/${script}" --out "$out"; then
        echo "    wrote $out"
    else
        echo "    FAILED" >&2
        FAILED+=("$name")
    fi
done

if [[ "${#FAILED[@]}" -gt 0 ]]; then
    echo
    echo "benchmarks failed: ${FAILED[*]}" >&2
    exit 1
fi

echo
echo "all benchmarks complete: $OUT_DIR"
