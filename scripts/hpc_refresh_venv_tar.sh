#!/usr/bin/env bash
# hpc_refresh_venv_tar.sh — regenerate /home/wr286/venv_site.tar on CSD3.
#
# The sbatch template at scripts/hpc_templates/polychord_intr.sbatch
# extracts this tarball to /tmp per job because some CSD3 compute nodes
# have broken Lustre OSTs that hang Python imports from /rds
# (see scripts/hpc_install_polychord.sh and issue #269 session log).
#
# When the venv on HPC changes — e.g. after running
# hpc_install_polychord.sh, or after a manual `.venv/bin/pip install`
# on the login node — the tarball becomes stale.  Jobs then run against
# the old deps (no pypolychord, old numpy, etc.) and fail in obscure
# ways at import time.
#
# Run this script whenever you mutate the HPC venv.  It is idempotent
# and cheap (~2 s + tarball write to /home).
#
# Usage:
#   bash scripts/hpc_refresh_venv_tar.sh
#
# Preconditions:
#   - `ssh csd3` ControlMaster is live (authenticated).
#   - The HPC venv exists at ${HPC_ROOT}/.venv (usually
#     /rds/user/wr286/hpc-work/tidal/.venv).
#
# The tarball lives on /home (NFS, reliable) rather than /rds (Lustre,
# sometimes broken).  Per-job extraction to /tmp sidesteps both.

set -euo pipefail

readonly HOST="${HPC_HOST:-csd3}"
readonly REMOTE_ROOT="${HPC_ROOT:-/rds/user/wr286/hpc-work/tidal}"
readonly TAR_PATH="${HPC_VENV_TAR:-/home/wr286/venv_site.tar}"

if ! ssh -O check "$HOST" 2>/dev/null; then
    cat >&2 <<EOF
error: no live SSH ControlMaster for '$HOST'.

Open one interactively first (password + TOTP):

    ssh $HOST

Then re-run this command.
EOF
    exit 2
fi

echo "==> regenerating ${TAR_PATH} from ${HOST}:${REMOTE_ROOT}/.venv"

ssh "$HOST" 'bash -s' <<REMOTE_EOF
set -euo pipefail

REMOTE_ROOT="${REMOTE_ROOT}"
TAR_PATH="${TAR_PATH}"

# Resolve purelib inside the project venv so we tar the right dir even
# if the Python version changes (e.g. 3.11 -> 3.12 on a future upgrade).
SP="\$("\${REMOTE_ROOT}/.venv/bin/python" -c 'import sysconfig; print(sysconfig.get_paths()["purelib"])')"
if [[ ! -d "\${SP}" ]]; then
    echo "error: resolved site-packages path does not exist: \${SP}" >&2
    exit 1
fi

# Tar the site-packages/ directory (NOT its contents) so extraction
# reproduces the original parent layout at \${LOCAL_SP}/site-packages/.
# The sbatch expects exactly that shape — keep them in sync.
PARENT_DIR="\$(dirname "\${SP}")"
LEAF="\$(basename "\${SP}")"
echo "  site-packages at: \${SP}"
echo "  writing:          \${TAR_PATH}"

# Write to a temp path first, then rename — avoids leaving a half-written
# tarball if the caller ctrl-Cs mid-archive (jobs reading it mid-write
# would otherwise crash the next time this runs).
TMP="\${TAR_PATH}.new.\$\$"
tar cf "\${TMP}" -C "\${PARENT_DIR}" "\${LEAF}"
mv "\${TMP}" "\${TAR_PATH}"

echo "  size:             \$(du -sh "\${TAR_PATH}" | cut -f1)"

# Sanity check: the tar must contain pypolychord if install.sh has run.
# Flag (don't hard-fail) so a fresh venv that legitimately lacks it can
# still regenerate the tar.  Grep any entry under pypolychord/ since the
# leaf name may be any of {__init__.py, settings.py, polychord.py}.
# Disable pipefail for this check only: 'grep -q' closes the pipe as
# soon as it matches, which gives 'tar' SIGPIPE and with pipefail would
# fail the whole pipeline even on a successful match.
set +o pipefail
if tar tf "\${TAR_PATH}" | grep -q 'pypolychord/settings.py'; then
    echo "  pypolychord:      present (good)"
else
    echo "  pypolychord:      MISSING — run scripts/hpc_install_polychord.sh if you need it"
fi
set -o pipefail
REMOTE_EOF

echo "==> tarball refreshed on ${HOST}"
