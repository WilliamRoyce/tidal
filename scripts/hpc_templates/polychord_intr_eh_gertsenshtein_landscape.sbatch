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

# EH+Gertsenshtein landscape — GR+Maxwell+F^4 Euler-Heisenberg corrections.
#
# Theory: gertsenshtein_eh.json (derived 2026-05-25).
# Lagrangian: (1/κ²)R - (1/4)F²  + (ρ/8)(F²)² + (σ/8)(FF̃)²
# Both h_ij (graviton) and A_μ (photon) fields; B₀ background coupling.
#
# ndim=2 (ρ, σ). At the QED point: σ/ρ = 7/16 (Dunne/Adler).
# nlive=50 (25×ndim), num_repeats=4 (2×ndim), t_end=1.
# Est. wall time: ~5-10 min (trivially cheap 2D space).
#
# Note: pure-EM euler_heisenberg.json cannot test Gertsenshtein conversion
# (no graviton fields). This theory is the combined EH+Gert setup.
#
# Submit via:
#   bash scripts/hpc_shuttle.sh submit \
#     --template scripts/hpc_templates/polychord_intr_eh_gertsenshtein_landscape.sbatch \
#     --ntasks 8 --nodes 1 --name eh-gert-landscape \
#     --account HANDLEY-SL3-CPU --cmd true
#
# After: bash scripts/hpc_shuttle.sh pull <jobid>

set -euo pipefail

TIDAL_ROOT="{{REMOTE_ROOT}}"
LOCAL_SP="/tmp/site_packages_${SLURM_JOB_ID}"
RESULTS_DIR="${TIDAL_ROOT}/hpc_results/${SLURM_JOB_ID}"

NDIM=2
NLIVE=50   # 25 * NDIM
NREP=4     # 2 * NDIM

echo "=== Job ${SLURM_JOB_ID} on $(hostname), ntasks=${SLURM_NTASKS} ==="
echo "=== EH+Gertsenshtein landscape (t_end=1, ndim=2, nlive=50, nrep=4) ==="
date

echo "--- Extracting site-packages tarball ---"
mkdir -p "${LOCAL_SP}"
tar xf /home/wr286/venv_site.tar -C "${LOCAL_SP}"

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

DATA="examples/data/gertsenshtein_eh.json"

PHYS_ARGS="--param kappa=1.0 --param B0=0.01
  --baseline-formula sin(kappa*B0*t_end/2)**2
  --soft-floor-noise 1.0
  --soft-floor-logl -15
  --grid-shape 64 --bounds=0:50 --periodic
  --ic plane-wave --ic-component h_5 --ic-wavevector 2.0 --ic-amplitude 1e-2
  --t-end 1 --snapshots 2
  --measure conversion,peak_conversion --source h_5 --target a_1
  --method nested --sampler polychord --precision-criterion 0.01"

# 2 priors — rho and sigma (EH F^4 coefficients); arctan_uniform (both signs physical)
# QED point: sigma/rho = 7/16 ~ 0.4375 (Dunne 2012; Adler 1971)
PRIOR_ARGS="--prior rho=arctan_uniform:-89:89
  --prior sigma=arctan_uniform:-89:89"

NTASKS_PER_CHAIN=$(( SLURM_NTASKS / 2 ))

echo "--- Launching EH+Gert amp (${NTASKS_PER_CHAIN}r) + sup (${NTASKS_PER_CHAIN}r) in parallel ---"
echo "    ndim=2, nlive=50, num_repeats=4, t_end=1, soft_floor_logl=-15"
echo "    QED point: sigma/rho = 7/16; landscape surveys full (rho,sigma) plane"

# shellcheck disable=SC2086
mpirun -n "${NTASKS_PER_CHAIN}" tidal sample ${DATA} ${PHYS_ARGS} ${PRIOR_ARGS} \
    --likelihood "P_max:maximize" \
    --nlive "${NLIVE}" --num-repeats "${NREP}" \
    --output "${RESULTS_DIR}/eh_gert_amp_v1" \
    > "${RESULTS_DIR}/eh_gert_amp.log" 2>&1 &
PID_A=$!

# shellcheck disable=SC2086
mpirun -n "${NTASKS_PER_CHAIN}" tidal sample ${DATA} ${PHYS_ARGS} ${PRIOR_ARGS} \
    --likelihood "P_max:minimize" \
    --nlive "${NLIVE}" --num-repeats "${NREP}" \
    --output "${RESULTS_DIR}/eh_gert_sup_v1" \
    > "${RESULTS_DIR}/eh_gert_sup.log" 2>&1 &
PID_B=$!

echo "    EH+Gert amp pid=${PID_A}  sup pid=${PID_B} — waiting..."
wait ${PID_A}; RC_A=$?
wait ${PID_B}; RC_B=$?

echo "=== [$(date +%H:%M:%S)] Done — amp rc=${RC_A}  sup rc=${RC_B} ==="
echo "=== Pull: bash scripts/hpc_shuttle.sh pull ${SLURM_JOB_ID} ==="
