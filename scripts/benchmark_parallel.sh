#!/bin/bash
# Benchmark parallel vs serial Wolfram derivation pipeline.
# Usage: bash scripts/benchmark_parallel.sh <theory.toml>
#
# Generates the WLS script, runs it in both parallel and serial modes,
# and compares timing + output correctness.

set -euo pipefail

TOML="${1:?Usage: $0 <theory.toml>}"
THEORY=$(basename "$(dirname "$TOML")")

echo "=== Benchmark: $THEORY ==="
echo "TOML: $TOML"
echo

# Generate the WLS script
WLS_BASE="/tmp/bench_${THEORY}.wls"
OUT_PARALLEL="/tmp/bench_parallel_${THEORY}.json"
OUT_SERIAL="/tmp/bench_serial_${THEORY}.json"

uv run tidal derive "$TOML" --dry-run > "$WLS_BASE" 2>/dev/null

# Create parallel version (as-is, just override output path)
WLS_PARALLEL="/tmp/bench_parallel_${THEORY}.wls"
sed "s|outputPath = \"[^\"]*\"|outputPath = \"$OUT_PARALLEL\"|" "$WLS_BASE" > "$WLS_PARALLEL"

# Create serial version: delete the entire parallel init block and
# keep $useParallel = False (header default).  This gives a true serial
# baseline without the 16s subkernel launch overhead.
WLS_SERIAL="/tmp/bench_serial_${THEORY}.wls"
sed "s|outputPath = \"[^\"]*\"|outputPath = \"$OUT_SERIAL\"|" "$WLS_BASE" | \
    sed '/=== Parallel Subkernel Initialization ===/,/After parallel init/d' | \
    sed '/Close parallel subkernels/,/Subkernels closed/d' \
    > "$WLS_SERIAL"

echo "--- PARALLEL run ---"
START_P=$(date +%s)
wolframscript -file "$WLS_PARALLEL" 2>&1 | \
    grep -E "TIMING|Parallel|Decomposing|term |component|Hamiltonian terms" | \
    sed 's/^/  /'
END_P=$(date +%s)
WALL_P=$((END_P - START_P))
echo "  Wall time: ${WALL_P}s"
echo

echo "--- SERIAL run ---"
START_S=$(date +%s)
wolframscript -file "$WLS_SERIAL" 2>&1 | \
    grep -E "TIMING|Serial|Decomposing|term |component|Hamiltonian terms" | \
    sed 's/^/  /'
END_S=$(date +%s)
WALL_S=$((END_S - START_S))
echo "  Wall time: ${WALL_S}s"
echo

echo "--- COMPARISON ---"
echo "  Parallel wall time: ${WALL_P}s"
echo "  Serial wall time:   ${WALL_S}s"
if [ "$WALL_P" -gt 0 ]; then
    SPEEDUP=$(python3 -c "print(f'{$WALL_S/$WALL_P:.2f}x')")
    echo "  Speedup: $SPEEDUP"
fi
echo

echo "--- OUTPUT CORRECTNESS ---"
if [ -f "$OUT_PARALLEL" ] && [ -f "$OUT_SERIAL" ]; then
    python3 -c "
import json, sys
with open('$OUT_PARALLEL') as f: p = json.load(f)
with open('$OUT_SERIAL') as f: s = json.load(f)

# Compare equations
peq = len(p['equations'])
seq_ = len(s['equations'])
print(f'  Equations: parallel={peq}, serial={seq_}', '✓' if peq == seq_ else '✗ MISMATCH')

# Compare hamiltonian terms
ph = sorted(p.get('canonical', {}).get('hamiltonian_terms', []),
            key=lambda t: (t['factor_a']['field'], t['factor_b']['field'], t['factor_a']['operator']))
sh = sorted(s.get('canonical', {}).get('hamiltonian_terms', []),
            key=lambda t: (t['factor_a']['field'], t['factor_b']['field'], t['factor_a']['operator']))
print(f'  H-terms:   parallel={len(ph)}, serial={len(sh)}', '✓' if len(ph) == len(sh) else '✗ MISMATCH')

errors = 0
for i, (pt, st) in enumerate(zip(ph, sh)):
    cdiff = abs(pt['coefficient'] - st['coefficient'])
    if cdiff > 1e-10:
        print(f'  ✗ H-term {i} coeff: {pt[\"coefficient\"]} vs {st[\"coefficient\"]} (diff={cdiff})')
        errors += 1
    if pt['factor_a'] != st['factor_a'] or pt['factor_b'] != st['factor_b']:
        print(f'  ✗ H-term {i} factors mismatch')
        errors += 1

# Compare equation coefficients
for i, (pe, se) in enumerate(zip(p['equations'], s['equations'])):
    if pe.get('field') != se.get('field'):
        print(f'  ✗ Equation {i} field mismatch: {pe.get(\"field\")} vs {se.get(\"field\")}')
        errors += 1
    pterms = sorted(pe.get('terms', []), key=lambda t: str(t))
    sterms = sorted(se.get('terms', []), key=lambda t: str(t))
    if len(pterms) != len(sterms):
        print(f'  ✗ Equation {i} term count: {len(pterms)} vs {len(sterms)}')
        errors += 1
    for j, (pt2, st2) in enumerate(zip(pterms, sterms)):
        if abs(pt2.get('coefficient', 0) - st2.get('coefficient', 0)) > 1e-10:
            print(f'  ✗ Eq {i} term {j} coeff mismatch')
            errors += 1

if errors == 0:
    print('  ✓ ALL OUTPUTS IDENTICAL')
else:
    print(f'  ✗ {errors} MISMATCHES FOUND')
    sys.exit(1)
" 2>&1
else
    echo "  ✗ Missing output files"
    ls -la "$OUT_PARALLEL" "$OUT_SERIAL" 2>&1 || true
    exit 1
fi

# Cleanup
rm -f "$WLS_BASE" "$WLS_PARALLEL" "$WLS_SERIAL" "$OUT_PARALLEL" "$OUT_SERIAL"
echo
echo "=== Benchmark complete ==="
