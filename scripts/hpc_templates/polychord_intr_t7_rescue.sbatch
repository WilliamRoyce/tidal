#!/bin/bash
#SBATCH -J {{JOB_NAME}}
#SBATCH -A {{ACCOUNT}}
#SBATCH -p icelake
#SBATCH --qos=intr
#SBATCH --nodes={{NODES}}
#SBATCH --ntasks={{NTASKS}}
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1500
#SBATCH --time={{TIME}}
#SBATCH -o {{REMOTE_ROOT}}/slurm_logs/slurm-%j.out
#SBATCH -e {{REMOTE_ROOT}}/slurm_logs/slurm-%j.out
#SBATCH --exclude=cpu-r-10,cpu-r-15

# T7 Complete-Even PGT rescue — t_end=1 (overflow-safe).
#
# v1 (29588680) used t_end=10 → modal solver expm overflows float64 at γ>71
# (max_eigenvalue × t_end > 708). 74% floor hits were numerical overflow, not
# physics. t_end=1 raises the threshold to γ>709 with zero quality trade-off:
# modal expm is O(1) in t_end, and A=P/P_GR is t_end-independent at machine
# precision (both P and P_GR scale as (κB₀t/2)²).
#
# ndim=18, nlive=450 (25×ndim), num_repeats=36 (2×ndim), t_end=1.
# Est. wall time: ~30-40 min at 76 ranks (38r amp + 38r sup in parallel).
#
# Theory: torsion_gertsenshtein_complete_even.json
#   Params: beta1-3 (torsion mass), xi (torsion kinetics), delta1 (Ricci×F),
#           zeta1-3 (∂T×F), chi1-10 (full R̃×∂T basis).
#
# Submit via:
#   bash scripts/hpc_shuttle.sh submit \
#     --template scripts/hpc_templates/polychord_intr_t7_rescue.sbatch \
#     --ntasks 76 --nodes 1 --name t7-rescue \
#     --account HANDLEY-SL3-CPU --cmd true
#
# After: bash scripts/hpc_shuttle.sh pull <jobid>

set -euo pipefail

TIDAL_ROOT="{{REMOTE_ROOT}}"
LOCAL_SP="/tmp/site_packages_${SLURM_JOB_ID}"
RESULTS_DIR="${TIDAL_ROOT}/hpc_results/${SLURM_JOB_ID}"

echo "=== Job ${SLURM_JOB_ID} on $(hostname), ntasks=${SLURM_NTASKS} ==="
echo "=== T7 Complete-Even PGT rescue (t_end=1, ndim=18, nlive=450) ==="
date

echo "--- Extracting site-packages tarball ---"
mkdir -p "${LOCAL_SP}"
tar xf $HOME/venv_site.tar -C "${LOCAL_SP}"

source "${TIDAL_ROOT}/.venv/bin/activate"
export PYTHONPATH="${LOCAL_SP}/site-packages:${PYTHONPATH:-}"

module load rhel8/default-icl 2>/dev/null || true
module load intel/2019.3.199 2>/dev/null || true
module load intel-oneapi-mpi 2>/dev/null || true

echo "--- Checking venv freshness ---"
if ! python -c "from pypolychord import run_polychord; import anesthetic; import mpi4py" 2>&1; then
    cat >&2 <<EOF
error: pypolychord, anesthetic, or mpi4py missing from the extracted venv.
Refresh: bash scripts/hpc_install_polychord.sh && bash scripts/hpc_refresh_venv_tar.sh
EOF
    exit 3
fi

mkdir -p "${RESULTS_DIR}"
cd "${TIDAL_ROOT}"

DATA="examples/data/torsion_gertsenshtein_complete_even.json"

PHYS_ARGS="--param kappa=1.0 --param B0=0.01
  --baseline-formula sin(kappa*B0*t_end/2)**2
  --soft-floor-noise 1.0
  --soft-floor-logl -15
  --grid-shape 64 --bounds=0:50 --periodic
  --ic plane-wave --ic-component h_5 --ic-wavevector 2.0 --ic-amplitude 1e-2
  --t-end 1 --snapshots 2
  --measure conversion,peak_conversion --source h_5 --target a_1
  --method nested --sampler polychord --precision-criterion 0.01"

# 18 priors — full T5 + chi1-10 basis
PRIOR_ARGS="--prior beta1=arctan_uniform:-89:89
  --prior beta2=arctan_uniform:-89:89
  --prior beta3=arctan_uniform:-89:89
  --prior xi=log_uniform:1e-3:1e3
  --prior delta1=arctan_uniform:-89:89
  --prior zeta1=arctan_uniform:-89:89
  --prior zeta2=arctan_uniform:-89:89
  --prior zeta3=arctan_uniform:-89:89
  --prior chi1=arctan_uniform:-89:89
  --prior chi2=arctan_uniform:-89:89
  --prior chi3=arctan_uniform:-89:89
  --prior chi4=arctan_uniform:-89:89
  --prior chi5=arctan_uniform:-89:89
  --prior chi6=arctan_uniform:-89:89
  --prior chi7=arctan_uniform:-89:89
  --prior chi8=arctan_uniform:-89:89
  --prior chi9=arctan_uniform:-89:89
  --prior chi10=arctan_uniform:-89:89"

echo "--- Launching T7 amp (38r) + T7 sup (38r) in parallel ---"
echo "    ndim=18, nlive=450, num_repeats=36, t_end=1, soft_floor_logl=-15"
echo "    Est. convergence: ~1125 dead pts ≈ 35 min"

# shellcheck disable=SC2086
mpirun -n 38 tidal sample ${DATA} ${PHYS_ARGS} ${PRIOR_ARGS} \
    --likelihood "P_max:maximize" \
    --nlive 450 --num-repeats 36 \
    --output "${RESULTS_DIR}/t7_amp_v2" \
    > "${RESULTS_DIR}/t7_amp.log" 2>&1 &
PID_A=$!

# shellcheck disable=SC2086
mpirun -n 38 tidal sample ${DATA} ${PHYS_ARGS} ${PRIOR_ARGS} \
    --likelihood "P_max:minimize" \
    --nlive 450 --num-repeats 36 \
    --output "${RESULTS_DIR}/t7_sup_v2" \
    > "${RESULTS_DIR}/t7_sup.log" 2>&1 &
PID_B=$!

echo "    T7 amp pid=${PID_A}  T7 sup pid=${PID_B} — waiting..."
wait ${PID_A}; RC_A=$?
wait ${PID_B}; RC_B=$?

echo "=== [$(date +%H:%M:%S)] Done — T7 amp rc=${RC_A}  T7 sup rc=${RC_B} ==="
echo "=== Pull: bash scripts/hpc_shuttle.sh pull ${SLURM_JOB_ID} ==="
