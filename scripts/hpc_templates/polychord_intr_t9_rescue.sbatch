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

# T9 Complete-Even + full ξ kinetic rescue — t_end=1 (overflow-safe).
#
# v1 (29596287) used t_end=10 → modal solver expm overflows float64 at γ>71
# (max_eigenvalue × t_end > 708). Floor-dominated results were numerical overflow,
# not physics. t_end=1 raises threshold to γ>709; modal expm is O(1) in t_end,
# A=P/P_GR is t_end-independent at machine precision.
#
# ndim=32 (xi11 vanishes in plane-wave reduction, confirmed 2026-05-23).
# nlive=800 (25×ndim), num_repeats=64 (2×ndim), t_end=1.
#
# Theory: torsion_gertsenshtein_complete_even_full_xi.json
#   Params: beta1-3 (torsion mass), xi1-16 (full (∇T)² kinetic basis, xi11 excluded),
#           delta1 (Ricci×F), zeta1-3 (∂T×F), chi1-10 (full R̃×∂T basis).
#
# ξ₁-ξ₁₆ all use arctan_uniform (NOT log_uniform — individual kinetic invariants
# don't have Barker's positivity constraint; positive/negative values both physical).
#
# Submit via:
#   bash scripts/hpc_shuttle.sh submit \
#     --template scripts/hpc_templates/polychord_intr_t9_rescue.sbatch \
#     --ntasks 76 --nodes 1 --name t9-rescue \
#     --account HANDLEY-SL3-CPU --cmd true
#
# After: bash scripts/hpc_shuttle.sh pull <jobid>

set -euo pipefail

TIDAL_ROOT="{{REMOTE_ROOT}}"
LOCAL_SP="/tmp/site_packages_${SLURM_JOB_ID}"
RESULTS_DIR="${TIDAL_ROOT}/hpc_results/${SLURM_JOB_ID}"

NDIM=32   # confirmed 2026-05-23: xi11 vanishes in plane-wave reduction (32/33 survive)
NLIVE=800  # 25 * NDIM
NREP=64    # 2 * NDIM

echo "=== Job ${SLURM_JOB_ID} on $(hostname), ntasks=${SLURM_NTASKS} ==="
echo "=== T9 Complete-Even + full xi rescue (t_end=1, ndim=32, nlive=800, nrep=64) ==="
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

DATA="examples/data/torsion_gertsenshtein_complete_even_full_xi.json"

PHYS_ARGS="--param kappa=1.0 --param B0=0.01
  --baseline-formula sin(kappa*B0*t_end/2)**2
  --soft-floor-noise 1.0
  --soft-floor-logl -15
  --grid-shape 64 --bounds=0:50 --periodic
  --ic plane-wave --ic-component h_5 --ic-wavevector 2.0 --ic-amplitude 1e-2
  --t-end 1 --snapshots 2
  --measure conversion,peak_conversion --source h_5 --target a_1
  --method nested --sampler polychord --precision-criterion 0.01"

# 32 priors — xi11 vanishes in plane-wave reduction and is excluded (all arctan_uniform)
PRIOR_ARGS="--prior beta1=arctan_uniform:-89:89
  --prior beta2=arctan_uniform:-89:89
  --prior beta3=arctan_uniform:-89:89
  --prior xi1=arctan_uniform:-89:89
  --prior xi2=arctan_uniform:-89:89
  --prior xi3=arctan_uniform:-89:89
  --prior xi4=arctan_uniform:-89:89
  --prior xi5=arctan_uniform:-89:89
  --prior xi6=arctan_uniform:-89:89
  --prior xi7=arctan_uniform:-89:89
  --prior xi8=arctan_uniform:-89:89
  --prior xi9=arctan_uniform:-89:89
  --prior xi10=arctan_uniform:-89:89
  --prior xi12=arctan_uniform:-89:89
  --prior xi13=arctan_uniform:-89:89
  --prior xi14=arctan_uniform:-89:89
  --prior xi15=arctan_uniform:-89:89
  --prior xi16=arctan_uniform:-89:89
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

NTASKS_PER_CHAIN=$(( SLURM_NTASKS / 2 ))

echo "--- Launching T9 amp (${NTASKS_PER_CHAIN}r) + T9 sup (${NTASKS_PER_CHAIN}r) in parallel ---"
echo "    ndim=32, nlive=800, num_repeats=64, t_end=1, soft_floor_logl=-15"

# shellcheck disable=SC2086
mpirun -n "${NTASKS_PER_CHAIN}" tidal sample ${DATA} ${PHYS_ARGS} ${PRIOR_ARGS} \
    --likelihood "P_max:maximize" \
    --nlive "${NLIVE}" --num-repeats "${NREP}" \
    --output "${RESULTS_DIR}/t9_amp_v2" \
    > "${RESULTS_DIR}/t9_amp.log" 2>&1 &
PID_A=$!

# shellcheck disable=SC2086
mpirun -n "${NTASKS_PER_CHAIN}" tidal sample ${DATA} ${PHYS_ARGS} ${PRIOR_ARGS} \
    --likelihood "P_max:minimize" \
    --nlive "${NLIVE}" --num-repeats "${NREP}" \
    --output "${RESULTS_DIR}/t9_sup_v2" \
    > "${RESULTS_DIR}/t9_sup.log" 2>&1 &
PID_B=$!

echo "    T9 amp pid=${PID_A}  T9 sup pid=${PID_B} — waiting..."
wait ${PID_A}; RC_A=$?
wait ${PID_B}; RC_B=$?

echo "=== [$(date +%H:%M:%S)] Done — T9 amp rc=${RC_A}  T9 sup rc=${RC_B} ==="
echo "=== Pull: bash scripts/hpc_shuttle.sh pull ${SLURM_JOB_ID} ==="
