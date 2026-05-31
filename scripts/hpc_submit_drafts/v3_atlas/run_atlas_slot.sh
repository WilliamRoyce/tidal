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

  # Optional per-theory RANKS override via TIDAL_RANKS_<TAG_UPPER>.
  # Use case: a small straggler-finish slot can bump ranks/face when only
  # a few chains run (e.g. TIDAL_RANKS_T5_1=28 with TIDAL_FACES_T5_1="11 12"
  # gives 2*28=56 cores per theory, fitting 4 chains in 112 sapphire cores).
  local tag_upper="${tag^^}"
  local ranks_env_var="TIDAL_RANKS_${tag_upper}"
  local ranks_override="${!ranks_env_var:-}"
  if [[ -n "${ranks_override}" ]]; then
    RANKS="${ranks_override}"
  fi
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

# ─── Face list override per theory ──────────────────────────────────────────
# Default: launch faces 1..FACES (whole atlas).  Set TIDAL_FACES_<TAG_UPPER>
# to a space-separated list to run only a subset — e.g.
#   TIDAL_FACES_T5_0="1 2 3 4 6"
# launches only the incomplete face indices and skips the converged ones.
# Use case: an earlier slot left some chains with inference.json (done) and
# some still incomplete; reschedule only the incomplete ones in this slot.

get_face_list() {
  local tag="$1"
  # bash ${var^^} → uppercase; e.g. t5_0 → T5_0
  local tag_upper="${tag^^}"
  local env_var="TIDAL_FACES_${tag_upper}"
  local list="${!env_var:-}"
  if [[ -n "${list}" ]]; then
    echo "${list}"
  else
    seq 1 "${FACES}"
  fi
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

  # Compute face_label (psalter convention: 01p=face 1, 01m=face 2, ...)
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
  local face_dir="${outdir}/${face_label}_tile${sub}"

  # Optional FRESH mode: if TIDAL_FRESH=1, delete this face's chain dir so
  # PolyChord starts from scratch (defends against stale .resume hangs).
  if [[ "${TIDAL_FRESH:-0}" == "1" && -n "${face_label}" ]]; then
    rm -rf "${face_dir}"
  fi

  # Auto-detect resume: if TIDAL_READ_RESUME=1 AND a .resume file ALREADY
  # exists for this face, append --read-resume. Skip --read-resume if no
  # checkpoint exists yet (PolyChord errors on missing resume; this lets us
  # mix RESUMED and FRESH theories in the same sequential slot).
  local resume_flag=()
  if [[ "${TIDAL_READ_RESUME:-0}" == "1" && -f "${face_dir}/_chains/tidal.resume" ]]; then
    resume_flag=(--read-resume)
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

# ─── Main: launch theories ───────────────────────────────────────────────────
# Default: all theories launched in parallel (max throughput when they fit).
# TIDAL_SEQUENTIAL=1: run theories one at a time — each theory's faces still
#   run in parallel within the theory, but the next theory only starts after
#   the previous theory's chains all return. Use this to (a) get full node
#   cores per theory, (b) maximize one theory's converge-rate, (c) avoid
#   the simultaneous-resume contention seen in slot 29878165.

echo "=== [$(date +%H:%M:%S)] ATLAS SLOT START — theories: $* (sequential=${TIDAL_SEQUENTIAL:-0}) ==="

failures=0

if [[ "${TIDAL_SEQUENTIAL:-0}" == "1" ]]; then
  # Sequential mode: launch one theory, wait for it, then next.
  for tag in "$@"; do
    lookup_theory "${tag}"
    echo "=== [$(date +%H:%M:%S)] starting theory ${tag}: ndim=${NDIM} faces=${FACES} ranks/face=${RANKS} ==="
    THEORY_PIDS=()
    for f in $(get_face_list "${tag}"); do
      launch_face "${tag}" "${f}"
      THEORY_PIDS+=("$!")
    done
    echo "=== [$(date +%H:%M:%S)] ${tag} launched ${#THEORY_PIDS[@]} chains; waiting ==="
    theory_fail=0
    for pid in "${THEORY_PIDS[@]}"; do
      if ! wait "${pid}"; then
        theory_fail=$((theory_fail + 1))
      fi
    done
    failures=$((failures + theory_fail))
    echo "=== [$(date +%H:%M:%S)] ${tag} done — ${theory_fail} chain failures ==="
    # Per-theory sentinel for the mid-slot monitor + downstream pipelines
    touch "${RESULTS_DIR}/${tag}_DONE"
  done
else
  # Parallel mode: launch all theories' chains at once, wait for everything.
  PIDS=()
  total_cores=0
  for tag in "$@"; do
    lookup_theory "${tag}"
    face_list_str="$(get_face_list "${tag}" | tr '\n' ' ')"
    n_faces=$(echo "${face_list_str}" | wc -w)
    echo "  ${tag}: ndim=${NDIM} faces=${n_faces} (${face_list_str% }) ranks/face=${RANKS} (=${RANKS}*${n_faces}=$((RANKS * n_faces)) cores)"
    total_cores=$((total_cores + RANKS * n_faces))
    for f in ${face_list_str}; do
      launch_face "${tag}" "${f}"
      PIDS+=("$!")
    done
  done
  echo "=== [$(date +%H:%M:%S)] Launched ${#PIDS[@]} chains, total cores=${total_cores} ==="
  for pid in "${PIDS[@]}"; do
    if ! wait "${pid}"; then
      failures=$((failures + 1))
    fi
  done
fi

echo
echo "=== [$(date +%H:%M:%S)] ATLAS SLOT DONE — ${failures} chain failures ==="
# BATCH_DONE sentinel for the external file-existence wait
touch "${RESULTS_DIR}/BATCH_DONE"
date
