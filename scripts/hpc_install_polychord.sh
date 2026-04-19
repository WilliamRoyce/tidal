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

# 2. Activate venv + ensure setuptools is present
source "\${REMOTE_ROOT}/.venv/bin/activate"
pip install --quiet setuptools

# 3. Build libchord.so via ifort (MPI=0 for serial — MPI build is a
#    separate future task, see issue #269).
echo "--- make all MPI=0 ---"
make all MPI=0 2>&1 | tail -3
echo "libchord.so: \$(ls -lh lib/libchord.so)"

# 4. Build _pypolychord.*.so via gcc (CC override required — see header).
echo "--- setup.py build_ext --inplace ---"
CC=gcc python setup.py build_ext --inplace 2>&1 | tail -3
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
