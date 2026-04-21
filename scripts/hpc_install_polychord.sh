#!/usr/bin/env bash
# hpc_install_polychord.sh — one-time PolyChord install on CSD3.
#
# Run AFTER `ssh csd3` has a live ControlMaster (authenticated session).
# Executes the install remotely, no interactive prompts needed.
#
# What it does:
#   1. Clones PolyChordLite to ${REMOTE_ROOT}/PolyChordLite (NOT
#      /tmp or .polychord_build — the sbatch template expects this path).
#   2. Builds libchord.so (Fortran via ifort) and _pypolychord.*.so
#      (C++ extension via gcc; CC=gcc is required because the intel
#      module pollutes CC with icx, which distutils can't introspect).
#   3. Injects the built extension into the project venv's
#      site-packages so `import pypolychord` works without any
#      per-job PYTHONPATH juggling.  The layout honours the compiled
#      rpath `$ORIGIN/pypolychord/lib`.
#
# Usage:
#   bash scripts/hpc_install_polychord.sh
#
# After this succeeds, jobs that use scripts/hpc_templates/polychord_intr.sbatch
# still rebuild libchord.so per job (fast, ~10 s) but pick up the
# installed Python package from the venv.

set -euo pipefail

readonly HOST="${HPC_HOST:-csd3}"
readonly REMOTE_ROOT="${HPC_ROOT:-/rds/user/wr286/hpc-work/tidal}"

# Verify ControlMaster is live (don't trigger fresh MFA prompt)
if ! ssh -O check "$HOST" 2>/dev/null; then
    cat >&2 <<EOF
error: no live SSH ControlMaster for '$HOST'.

Open one interactively first (password + TOTP):

    ssh $HOST

Then re-run this command.
EOF
    exit 2
fi

echo "==> installing PolyChord on ${HOST}:${REMOTE_ROOT}"

# Heredoc sent to remote: single-quoted `REMOTE_EOF` terminator would
# prevent client-side $REMOTE_ROOT interpolation; leave it un-quoted and
# escape the remote's own variable expansions with backslashes.
ssh "$HOST" 'bash -s' <<REMOTE_EOF
set -euo pipefail

. /etc/profile.d/modules.sh
module purge
module load rhel8/default-icl
module load python/3.11.0-icl
module load intel/2019.3.199 2>/dev/null || true
# MPI toolchain for PolyChord's --ntasks>1 runs via mpirun.  Required
# for the #292 parallelism work; pypolychord auto-detects MPI.COMM_WORLD
# and falls back to serial cleanly when mpi4py isn't importable.
module load intel-oneapi-mpi 2>/dev/null || true

REMOTE_ROOT="${REMOTE_ROOT}"
POLYCHORD_DIR="\${REMOTE_ROOT}/PolyChordLite"
VENV_SP="\$(\${REMOTE_ROOT}/.venv/bin/python -c 'import sysconfig; print(sysconfig.get_paths()["purelib"])')"

# 1. Clone or update PolyChordLite
if [[ -d "\${POLYCHORD_DIR}/.git" ]]; then
    echo "updating existing clone at \${POLYCHORD_DIR}"
    cd "\${POLYCHORD_DIR}" && git pull --ff-only || true
else
    echo "cloning PolyChordLite to \${POLYCHORD_DIR}"
    rm -rf "\${POLYCHORD_DIR}"
    git clone --depth 1 https://github.com/PolyChord/PolyChordLite.git "\${POLYCHORD_DIR}"
fi

cd "\${POLYCHORD_DIR}"

# 2. Activate venv + ensure setuptools + mpi4py are present.  mpi4py is
#    what lets pypolychord auto-detect MPI.COMM_WORLD at import time;
#    without it the run falls back to serial single-process mode.
source "\${REMOTE_ROOT}/.venv/bin/activate"
pip install --quiet setuptools
# Prefer the already-loaded intel-oneapi-mpi toolchain when building
# mpi4py from source (pip wheel may be built against a different MPI
# and segfault at MPI_Init — always compile against the cluster MPI).
pip install --quiet --no-binary=mpi4py mpi4py

# 3. Build libchord.so via ifort, WITH MPI for #292 parallel runs.  The
#    MPI=1 path picks up the intel-oneapi-mpi mpif90 wrapper from the
#    module loaded above and links libchord.so against libmpi_*.so.
echo "--- make all MPI=1 ---"
# Wipe any previous MPI=0 build products to avoid mixing.
make veryclean 2>&1 | tail -3 || true
make all MPI=1 2>&1 | tail -3
echo "libchord.so: \$(ls -lh lib/libchord.so)"
# Confirm the library is actually MPI-linked.
if ldd lib/libchord.so 2>/dev/null | grep -qi 'mpi'; then
    echo "  MPI-linked: yes (good — ranks will talk to each other)"
else
    echo "  MPI-linked: NO — the mpirun launch will behave as serial.  Fix: reload intel-oneapi-mpi before rerunning."
fi

# 4. Build _pypolychord.*.so via distutils.  The 'intel' module export
#    of CC=icx breaks distutils introspection ('icx -v' exits non-zero),
#    so unset CC/CXX in a subshell and let sysconfig's default (gcc)
#    drive the build.
echo "--- setup.py build_ext --inplace ---"
( unset CC CXX; python setup.py build_ext --inplace ) 2>&1 | tail -3
echo "_pypolychord: \$(ls _pypolychord*.so)"

# 5. Inject into the venv's site-packages.  Layout matches what
#    \`import pypolychord\` expects (see template header).
echo "--- installing into \${VENV_SP} ---"
rm -rf "\${VENV_SP}/pypolychord" "\${VENV_SP}"/_pypolychord*.so
cp -r pypolychord "\${VENV_SP}/"
cp _pypolychord*.so "\${VENV_SP}/"
mkdir -p "\${VENV_SP}/pypolychord/lib"
cp lib/libchord.so "\${VENV_SP}/pypolychord/lib/"

# 6. Verify
python -c 'from pypolychord import run_polychord; print("pypolychord installed OK")'
python -c 'from anesthetic import read_chains; print("anesthetic installed OK")'
REMOTE_EOF

echo "==> PolyChord installed on HPC (${HOST}:${REMOTE_ROOT}/PolyChordLite)"
