#!/usr/bin/env bash
# Phase E Atlas — multi-theory parallel launcher for one INTR slot.
#
# Designed to run INSIDE an interactive sapphire INTR allocation
# (booked via `bash scripts/hpc_shuttle.sh interactive`).  Launches all
# face chains for each requested theory in parallel via `srun --exclusive`,
# waits for global completion, writes BATCH_DONE sentinel.
#
# Usage (on the compute node):
#   bash run_atlas_slot.sh <theory_tag> [<theory_tag> ...]
#
# Theory tags: t4, t1, t5_0, t5_1, t5_2, t5_3
#
# Per-theory metadata (BSM coupling names, JSON, ndim, ranks/face) is
# encoded in the case statement below.  RANKS_PER_CHAIN is chosen so that
# total cores across all theories fits the 112-core sapphire node.

set -euo pipefail

TIDAL_ROOT="/rds/user/wr286/hpc-work/tidal"
CONFIG_ENV="${TIDAL_ROOT}/scripts/hpc_submit_drafts/v3_atlas/_atlas_config.env"

# shellcheck disable=SC1090
source "${CONFIG_ENV}"

cd "${TIDAL_ROOT}"

# Module + venv setup (nohup ssh runs a non-interactive non-login shell so
# nothing is auto-sourced). Mirrors what hpc_shuttle.sh attach does plus the
# MPI toolchain needed for pypolychord's PolyChord MPI ranks.
. /etc/profile.d/modules.sh 2>/dev/null || true
module purge                              2>/dev/null || true
module load rhel8/default-sar             2>/dev/null || true
module load python/3.11.0-icl             2>/dev/null || true
module load intel/2019.3.199              2>/dev/null || true
module load intel-oneapi-mpi              2>/dev/null || true

# shellcheck disable=SC1091
source "${TIDAL_ROOT}/.venv/bin/activate"

# Full path to tidal — srun's exec ranks don't inherit PATH reliably across
# the SLURM PMI handoff, so always invoke via absolute path.
TIDAL_BIN="${TIDAL_ROOT}/.venv/bin/tidal"

# Sanity check before launching any srun: confirm tidal is executable.
if [[ ! -x "${TIDAL_BIN}" ]]; then
  echo "ERROR: ${TIDAL_BIN} is not executable" >&2
  ls -la "${TIDAL_BIN}" >&2 || true
  exit 1
fi
echo "tidal binary: ${TIDAL_BIN}"
"${TIDAL_BIN}" --version 2>&1 | head -3 || true

: "${RESULTS_DIR:?RESULTS_DIR must be set by INTR slot setup (hpc_shuttle.sh attach exports it)}"

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <theory_tag> [<theory_tag> ...]"
  echo "  theory tags: t4 t1 t5_0 t5_1 t5_2 t5_3"
  exit 1
fi

# ─── Theory metadata: looked up per <theory_tag> ────────────────────────────
# Returns (via globals): JSON_PATH, NDIM, FACES, RANKS, BSM_NAMES, PARAMS_PIN
# Pin params keep the rest of the v3 JSON's BSM coupling space at zero.

lookup_theory() {
  local tag="$1"
  case "$tag" in
    t4)
      JSON_PATH="examples/data/torsion_gertsenshtein_nonminimal.json"
      NDIM=4
      FACES=8           # 2 * NDIM
      RANKS=6           # Slot S1 budget: 8 faces * 6 ranks = 48 cores
      BSM_NAMES="alpha1,alpha2,alpha3,delta1"
      PARAMS_PIN=""
      ;;
    t1)
      JSON_PATH="examples/data/dark_photon_plasma.json"
      NDIM=4
      FACES=8
      RANKS=6           # Slot S1: 48 cores
      BSM_NAMES="alpha3,xi,deltam,mA2"
      PARAMS_PIN=""
      ;;
    t5_0)
      JSON_PATH="examples/data/torsion_gertsenshtein_general_nonminimal.json"
      NDIM=5
      FACES=10
      RANKS=11          # Solo fresh: 10 * 11 = 110 cores (one chain per face, single-theory slot)
      BSM_NAMES="beta1,beta2,beta3,xi,delta1"
      PARAMS_PIN="--param chi=0.0 --param zeta1=0.0 --param zeta2=0.0 --param zeta3=0.0"
      ;;
    t5_1)
      JSON_PATH="examples/data/torsion_gertsenshtein_general_nonminimal.json"
      NDIM=6
      FACES=12
      RANKS=9           # Solo fresh: 12 * 9 = 108 cores
      BSM_NAMES="beta1,beta2,beta3,xi,delta1,chi"
      PARAMS_PIN="--param zeta1=0.0 --param zeta2=0.0 --param zeta3=0.0"
      ;;
    t5_2)
      JSON_PATH="examples/data/torsion_gertsenshtein_general_nonminimal.json"
      NDIM=8
      FACES=16
      RANKS=6           # bonus / phase β
      BSM_NAMES="beta1,beta2,beta3,xi,delta1,zeta1,zeta2,zeta3"
      PARAMS_PIN="--param chi=0.0"
      ;;
    t5_3)
      JSON_PATH="examples/data/torsion_gertsenshtein_general_nonminimal.json"
      NDIM=9
      FACES=18
      RANKS=6           # bonus / phase β
      BSM_NAMES="beta1,beta2,beta3,xi,delta1,chi,zeta1,zeta2,zeta3"
      PARAMS_PIN=""
      ;;
    *)
      echo "ERROR: unknown theory tag: $tag" >&2
      exit 2
      ;;
  esac
}

# ─── Build cubed-sphere sub-tile string ─────────────────────────────────────
# Per V3_PHASE_C_REFERENCE.md: cell directory layout is `<face_label>_tile<s_1>_<s_2>_..._<s_{N-1}>`.
# For M=1 single tile per face: sub = "1_1_..._1" with (N-1) ones.

build_sub_string() {
  local ndim="$1"
  local sub=""
  for ((i = 1; i < ndim; i++)); do
    if [[ -z "$sub" ]]; then sub="1"; else sub="${sub}_1"; fi
  done
  printf '%s' "$sub"
}

# ─── Common physics flags (uniform-field v3 geometry) ────────────────────────
COMMON_PHYS=(
  --baseline-formula 'sin(kappa*B0*t_end/2)**2'
  --soft-floor-noise 1.0
  --param kappa=1.0 --param "B0=${B0}"
  --grid-shape "${GRID}" "--bounds=0:${L}" --periodic
  --ic plane-wave --ic-component h_5 --ic-wavevector "${K_CARRIER}" --ic-amplitude "${H0}"
  --t-end "${T_END}" --snapshots "${SNAPSHOTS}"
  --measure conversion,peak_conversion --source h_5 --target a_1
  --method nested --sampler polychord --precision-criterion 0.01
)

# ─── Launch one face's chain in background ──────────────────────────────────
# Args: theory_tag face_idx
launch_face() {
  local tag="$1" face="$2"
  local sub
  sub="$(build_sub_string "${NDIM}")"
  # tidal sample --joint-prior writes results to <output>/<face_label>_tile<sub>/
  # automatically (psalter convention). So point --output at the survey dir
  # and let the cubed-sphere prior code create the face_label subdirectory.
  local outdir="${RESULTS_DIR}/atlas_${tag}"
  local logfile="${RESULTS_DIR}/atlas_${tag}/face_${face}.log"
  mkdir -p "${outdir}"

  local nlive=$((25 * NDIM))
  local num_repeats=$((2 * NDIM))

  # Optional resume: if TIDAL_READ_RESUME env var is set to 1, append
  # --read-resume so PolyChord picks up from the existing .resume
  # checkpoint in <outdir>/<face_label>_tile<sub>/_chains/tidal.resume.
  local resume_flag=()
  if [[ "${TIDAL_READ_RESUME:-0}" == "1" ]]; then
    resume_flag=(--read-resume)
  fi

  # Optional FRESH mode: if TIDAL_FRESH=1, delete the existing _chains/ dir
  # for this face so PolyChord starts from scratch (prevents
  # --read-resume from silently kicking in and avoids stale-state hangs).
  # Use case: prior resume attempts left chains stuck in re-verification.
  if [[ "${TIDAL_FRESH:-0}" == "1" ]]; then
    # face_label_tile<sub> dir; tidal sample writes into a child dir whose
    # name depends on the cubed_sphere face label — easier to wipe the whole
    # outdir/<face_label>* glob safely (only matches this theory's tiles).
    local face_label
    case "${face}" in
      1) face_label="01p" ;; 2) face_label="01m" ;;
      3) face_label="02p" ;; 4) face_label="02m" ;;
      5) face_label="03p" ;; 6) face_label="03m" ;;
      7) face_label="04p" ;; 8) face_label="04m" ;;
      9) face_label="05p" ;; 10) face_label="05m" ;;
      11) face_label="06p" ;; 12) face_label="06m" ;;
      *) face_label="" ;;
    esac
    if [[ -n "${face_label}" ]]; then
      rm -rf "${outdir}/${face_label}_tile${sub}"
    fi
  fi

  # shellcheck disable=SC2086
  srun --exclusive --ntasks="${RANKS}" --cpus-per-task=1 --mem-per-cpu=1500 --mpi=pmi2 \
    "${TIDAL_BIN}" sample "${JSON_PATH}" \
    --joint-prior "names=${BSM_NAMES};type=cubed_sphere;M=${M_TILES};face=${face};sub=${sub};r_lo=${R_LO};r_hi=${R_HI}" \
    --likelihood "P_max:maximize" \
    "${COMMON_PHYS[@]}" \
    ${PARAMS_PIN} \
    --nlive "${nlive}" --num-repeats "${num_repeats}" \
    --output "${outdir}" \
    "${resume_flag[@]}" \
    > "${logfile}" 2>&1 &
}

# ─── Main: launch all faces of all requested theories in parallel ────────────

echo "=== [$(date +%H:%M:%S)] ATLAS SLOT START — theories: $* ==="

PIDS=()
total_cores=0
for tag in "$@"; do
  lookup_theory "${tag}"
  echo "  ${tag}: ndim=${NDIM} faces=${FACES} ranks/face=${RANKS} (=${RANKS}*${FACES}=$((RANKS * FACES)) cores)"
  total_cores=$((total_cores + RANKS * FACES))
  for ((f = 1; f <= FACES; f++)); do
    launch_face "${tag}" "${f}"
    PIDS+=("$!")
  done
done

echo "=== [$(date +%H:%M:%S)] Launched ${#PIDS[@]} chains, total cores=${total_cores} ==="

# Wait for all chains
failures=0
for pid in "${PIDS[@]}"; do
  if ! wait "${pid}"; then
    failures=$((failures + 1))
  fi
done

echo
echo "=== [$(date +%H:%M:%S)] ATLAS SLOT DONE — ${#PIDS[@]} chains, ${failures} failed ==="
# BATCH_DONE sentinel for the external file-existence wait
touch "${RESULTS_DIR}/BATCH_DONE"
date
