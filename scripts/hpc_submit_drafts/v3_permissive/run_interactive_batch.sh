#!/usr/bin/env bash
# Interactive-session batch runner for v3 D2.x theories.
#
# Run this INSIDE an interactive sapphire allocation after attaching:
#   bash scripts/hpc_shuttle.sh interactive --ntasks 76
#   bash scripts/hpc_shuttle.sh attach <jobid>
#   # then on the compute node:
#   bash /rds/user/wr286/hpc-work/tidal/scripts/hpc_submit_drafts/v3_permissive/run_interactive_batch.sh
#
# Runs 5 remaining D2.x theories SEQUENTIALLY using all 76 ranks per theory.
# Validated timing (D2.0 sup, nlive=125, 76 ranks): 9 min.
# Expected wall times: D2.1(6p)~10min, D2.2(8p)~12min, D2.3(9p)~15min.
# Total: ~64 min — fits within 1hr INTR with the interactive session ending
# at walltime; any theory still running will be SIGTERM'd. D2.3 sup is last
# as it is least likely to complete; submit it separately if time runs out.
#
# Theories covered (in order):
#   1. D2.1 Barker sup   (6p, nlive=150, num_repeats=12)
#   2. D2.2 Shapiro amp  (8p, nlive=200, num_repeats=16)
#   3. D2.2 Shapiro sup  (8p, nlive=200, num_repeats=16)
#   4. D2.3 Full amp     (9p, nlive=225, num_repeats=18)
#   5. D2.3 Full sup     (9p, nlive=225, num_repeats=18)  ← may not reach in 1hr
#
# Output: ${RESULTS_DIR}/<theory_name>/ for each theory.
# After the session: bash scripts/hpc_shuttle.sh pull <jobid> --src <RESULTS_DIR>

set -euo pipefail

TIDAL_ROOT="/rds/user/wr286/hpc-work/tidal"

# ── Environment setup ────────────────────────────────────────────────────────
echo "=== [$(date +%H:%M:%S)] Setting up PolyChord environment ==="
LOCAL_SP="/tmp/site_packages_${SLURM_JOB_ID}"
mkdir -p "${LOCAL_SP}"
echo "Extracting venv tarball (~253 MB)..."
tar xf /home/wr286/venv_site.tar -C "${LOCAL_SP}"

source "${TIDAL_ROOT}/.venv/bin/activate"
export PYTHONPATH="${LOCAL_SP}/site-packages:${PYTHONPATH:-}"
export TIDAL_MODAL_BACKEND=jax

. /etc/profile.d/modules.sh
module load rhel8/default-icl   2>/dev/null || true
module load intel/2019.3.199    2>/dev/null || true
module load intel-oneapi-mpi    2>/dev/null || true

python -c "from pypolychord import run_polychord; import anesthetic; import mpi4py; print('venv OK')"

cd "${TIDAL_ROOT}"
export RESULTS_DIR="${TIDAL_ROOT}/hpc_results/${SLURM_JOB_ID}"
mkdir -p "${RESULTS_DIR}"
echo "RESULTS_DIR=${RESULTS_DIR}"
echo ""

MPI="mpirun -n ${SLURM_NTASKS}"

# Shared physics flags (all D2.x theories)
DATA="examples/data/torsion_gertsenshtein_general_nonminimal.json"
PHYS="--param kappa=1.0 --param B0=0.01
  --baseline-formula sin(kappa*B0*t_end/2)**2
  --soft-floor-noise 1.0
  --grid-shape 64 --bounds=0:50 --periodic
  --ic plane-wave --ic-component h_5 --ic-wavevector 2.0 --ic-amplitude 1e-2
  --t-end 10 --snapshots 2
  --measure conversion,peak_conversion --source h_5 --target a_1
  --method nested --sampler polychord --precision-criterion 0.01"

# ── 1. D2.1 Barker sup  (6p: β₁₋₃, ξ, δ₁, χ) ──────────────────────────────
echo "=== [$(date +%H:%M:%S)] 1/5: D2.1 Barker sup (nlive=150, 6p) ==="
# shellcheck disable=SC2086
${MPI} tidal sample ${DATA} ${PHYS} \
    --param zeta1=0 --param zeta2=0 --param zeta3=0 \
    --prior "beta1=arctan_uniform:-89:89" \
    --prior "beta2=arctan_uniform:-89:89" \
    --prior "beta3=arctan_uniform:-89:89" \
    --prior "xi=log_uniform:1e-3:1e3" \
    --prior "delta1=arctan_uniform:-89:89" \
    --prior "chi=arctan_uniform:-89:89" \
    --likelihood "P_max:minimize" \
    --nlive 150 --num-repeats 12 \
    --output "${RESULTS_DIR}/d21_barker_sup_v3"
echo "=== [$(date +%H:%M:%S)] 1/5 done ==="
echo ""

# ── 2. D2.2 Shapiro amp  (8p: β₁₋₃, ξ, δ₁, ζ₁₋₃) ─────────────────────────
echo "=== [$(date +%H:%M:%S)] 2/5: D2.2 Shapiro amp (nlive=200, 8p) ==="
# shellcheck disable=SC2086
${MPI} tidal sample ${DATA} ${PHYS} \
    --param chi=0 \
    --prior "beta1=arctan_uniform:-89:89" \
    --prior "beta2=arctan_uniform:-89:89" \
    --prior "beta3=arctan_uniform:-89:89" \
    --prior "xi=log_uniform:1e-3:1e3" \
    --prior "delta1=arctan_uniform:-89:89" \
    --prior "zeta1=arctan_uniform:-89:89" \
    --prior "zeta2=arctan_uniform:-89:89" \
    --prior "zeta3=arctan_uniform:-89:89" \
    --likelihood "P_max:maximize" \
    --nlive 200 --num-repeats 16 \
    --output "${RESULTS_DIR}/d22_shapiro_amp_v3"
echo "=== [$(date +%H:%M:%S)] 2/5 done ==="
echo ""

# ── 3. D2.2 Shapiro sup  (8p: β₁₋₃, ξ, δ₁, ζ₁₋₃) ─────────────────────────
echo "=== [$(date +%H:%M:%S)] 3/5: D2.2 Shapiro sup (nlive=200, 8p) ==="
# shellcheck disable=SC2086
${MPI} tidal sample ${DATA} ${PHYS} \
    --param chi=0 \
    --prior "beta1=arctan_uniform:-89:89" \
    --prior "beta2=arctan_uniform:-89:89" \
    --prior "beta3=arctan_uniform:-89:89" \
    --prior "xi=log_uniform:1e-3:1e3" \
    --prior "delta1=arctan_uniform:-89:89" \
    --prior "zeta1=arctan_uniform:-89:89" \
    --prior "zeta2=arctan_uniform:-89:89" \
    --prior "zeta3=arctan_uniform:-89:89" \
    --likelihood "P_max:minimize" \
    --nlive 200 --num-repeats 16 \
    --output "${RESULTS_DIR}/d22_shapiro_sup_v3"
echo "=== [$(date +%H:%M:%S)] 3/5 done ==="
echo ""

# ── 4. D2.3 Full amp  (9p: all couplings) ───────────────────────────────────
echo "=== [$(date +%H:%M:%S)] 4/5: D2.3 Full amp (nlive=225, 9p) ==="
# shellcheck disable=SC2086
${MPI} tidal sample ${DATA} ${PHYS} \
    --prior "beta1=arctan_uniform:-89:89" \
    --prior "beta2=arctan_uniform:-89:89" \
    --prior "beta3=arctan_uniform:-89:89" \
    --prior "xi=log_uniform:1e-3:1e3" \
    --prior "delta1=arctan_uniform:-89:89" \
    --prior "chi=arctan_uniform:-89:89" \
    --prior "zeta1=arctan_uniform:-89:89" \
    --prior "zeta2=arctan_uniform:-89:89" \
    --prior "zeta3=arctan_uniform:-89:89" \
    --likelihood "P_max:maximize" \
    --nlive 225 --num-repeats 18 \
    --output "${RESULTS_DIR}/d23_full_amp_v3"
echo "=== [$(date +%H:%M:%S)] 4/5 done ==="
echo ""

# ── 5. D2.3 Full sup  (9p: all couplings) — may not reach within 1hr ────────
echo "=== [$(date +%H:%M:%S)] 5/5: D2.3 Full sup (nlive=225, 9p) ==="
# shellcheck disable=SC2086
${MPI} tidal sample ${DATA} ${PHYS} \
    --prior "beta1=arctan_uniform:-89:89" \
    --prior "beta2=arctan_uniform:-89:89" \
    --prior "beta3=arctan_uniform:-89:89" \
    --prior "xi=log_uniform:1e-3:1e3" \
    --prior "delta1=arctan_uniform:-89:89" \
    --prior "chi=arctan_uniform:-89:89" \
    --prior "zeta1=arctan_uniform:-89:89" \
    --prior "zeta2=arctan_uniform:-89:89" \
    --prior "zeta3=arctan_uniform:-89:89" \
    --likelihood "P_max:minimize" \
    --nlive 225 --num-repeats 18 \
    --output "${RESULTS_DIR}/d23_full_sup_v3"
echo "=== [$(date +%H:%M:%S)] 5/5 done ==="
echo ""

echo "=== ALL DONE. Pull with: bash scripts/hpc_shuttle.sh pull ${SLURM_JOB_ID} --src ${RESULTS_DIR} ==="
