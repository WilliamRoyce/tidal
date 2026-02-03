# Wolfram Engine - Complete Guide

## ✅ Current Status (Post-Setup)

**Wolfram Engine is now fully configured and persistent!** Both WolframKernel and WolframScript work automatically after container rebuilds.

## 🚀 Quick Usage

### Command Line

```bash
# WolframScript (simple)
wolframscript -code "2+2"
wolframscript -code "Integrate[x^2, x]"

# WolframKernel (direct)
/home/vscode/.local/wolfram/engine/14.3/Executables/WolframKernel -noprompt -run "Print[2+2]; Exit[]"

# Run script files
wolframscript -file mycode.wls
```

### Python Integration (Quickest Path)

Install the Wolfram Python client in the uv environment:

```bash
uv pip install wolframclient
```

Then use a `WolframLanguageSession` and evaluate xAct directly:

```python
from wolframclient.evaluation import WolframLanguageSession
from wolframclient.language import wlexpr

session = WolframLanguageSession()
session.evaluate(wlexpr('Needs["xAct`xTensor`"]'))
print(session.evaluate(wlexpr("$Version")))
print(session.evaluate(wlexpr("Integrate[x^2,x]")))
session.terminate()
```

If session startup fails, ensure the kernel is discoverable on `PATH` or pass the explicit kernel path:

```python
kernel = "/home/vscode/.local/wolfram/engine/14.3/Executables/WolframKernel"
session = WolframLanguageSession(kernel)
```

### xAct Tensor Package

**✅ xAct 1.3.0 is installed and ready for differential geometry computations!**

```bash
# xAct tensor calculations
wolframscript -code '<<xAct`xTensor`; DefManifold[M,4]; Print["4D spacetime dimension: ", Dim[M]];'
wolframscript -code '<<xAct`xCoba`; DefManifold[M,4]; DefChart[coord, M, {0,1,2,3}];'

# Test installation and compatibility fix
.devcontainer/tests/test-xact.wls
```

**Available xAct packages:**

- `xTensor`: Core tensor algebra ✅ Fully functional
- `xCoba`: Component calculations and coordinate charts ✅ Fully functional
- `xPerm`: Permutation groups ✅ **MathLink ENABLED - Full performance!**
- `xPert`: Perturbation theory ✅ Fully functional
- `Spinors`, `xIdeal`, `TexAct`, etc. ✅ Available

**Installation location:** `$WOLFRAM_USERBASE/Applications/xAct/` (persists across rebuilds)

**🎉 GLIBC Issue SOLVED!** xPerm now compiled from source with full MathLink support:

- Advanced permutation algorithms: ✅ Available
- Strong generating sets: ✅ Fast computation
- Stabilizer chain algorithms: ✅ High performance
- All xPerm functions: ✅ Native speed

**Maintenance:** Use `.devcontainer/scripts/build-xperm.sh` to rebuild xPerm after Wolfram updates.

## 🔍 Diagnostics

```bash
# Quick health check
bash .devcontainer/scripts/check-wolfram.sh

# Activation management
bash .devcontainer/scripts/wolfram-activation-manager.sh

# Check license
cat ~/.local/wolfram/userbase/Licensing/mathpass

# Check machine ID
cat /etc/machine-id
```

## 🔧 Persistence Architecture

### What's Mounted from Host:

- **Engine**: `~/.local/wolfram/engine/14.3` (6.7GB binaries)
- **License**: `~/.local/wolfram/userbase/Licensing/mathpass` (machine-specific licenses)
- **Activation Cache**: `~/.cache/Wolfram` (WolframScript tokens)
- **Machine ID**: `/etc/machine-id` (for license validation)

### Auto-Configuration:

- License symlinks created in all expected locations
- WolframScript.conf configured with correct kernel path
- Binary symlinks placed in PATH
- Activation data automatically restored from backup

### Backup Strategy:

- **Primary**: Activation tokens in mounted `~/.cache/Wolfram`
- **Backup**: Copy stored in `~/.local/wolfram/userbase/.activation_backup`
- **Recovery**: Automatic restoration via postCreateCommand on rebuild

## 🛠️ Troubleshooting

### If WolframScript stops working:

```bash
# Check if activation is present
bash .devcontainer/scripts/wolfram-activation-manager.sh status

# Restore from backup if needed
bash .devcontainer/scripts/wolfram-activation-manager.sh restore

# Re-activate if necessary (requires Wolfram ID)
wolframscript -activate
```

### If WolframKernel stops working:

```bash
# Check license files
bash .devcontainer/scripts/check-wolfram.sh

# Verify machine ID matches license
grep -o "^[a-f0-9]*" ~/.local/wolfram/userbase/Licensing/mathpass
cat /etc/machine-id

# If mismatch, check if on correct host machine
```

### Complete reset (if needed):

```bash
# Rebuild container (applies all mounts and auto-config)
# Ctrl+Shift+P → "Dev Containers: Rebuild Container"

# After rebuild, verify everything works
bash .devcontainer/scripts/check-wolfram.sh
```

## 📚 Understanding the Setup

### Two License Systems:

1. **mathpass** (offline): Machine-specific license file for WolframKernel
   - Located: `~/.local/wolfram/userbase/Licensing/mathpass`
   - Requires machine ID match
   - Works completely offline

2. **Cloud activation** (online): Wolfram ID tokens for WolframScript
   - Located: `~/.cache/Wolfram/WolframScript/`
   - Requires one-time online activation
   - Tokens persist via mount

### Why It Persists:

- **Host mounts**: Engine, license, and cache directories mounted from host
- **Machine ID**: Container uses host's machine ID for license validation
- **Automatic restore**: postCreateCommand restores activation on rebuild
- **Dual backup**: Both mounted cache and backup in userbase directory

## 🎯 Maintenance

**Normal operation**: No maintenance required! Everything works automatically.

**If activation is lost**: Run `bash .devcontainer/scripts/wolfram-activation-manager.sh restore`

**For new licenses**: Replace `~/.local/wolfram/userbase/Licensing/mathpass` on host

**Performance**: ~30 second startup vs ~10 minutes for full engine installation

## 🧹 Repo Hygiene (Never Commit)

Keep persistent Wolfram directories in `~/.local/wolfram/...` and `~/.cache/Wolfram/...` (outside the repo) so they never appear in git status. This repository already ignores license artifacts and installer bundles (see .gitignore for patterns like `*.mathpass`, `.Wolfram*`, `.WolframEngine*`, and `Wolfram*.{sh,run,tgz}`), so nothing sensitive should be committed.

## 🔄 Initial Setup (Historical Reference)

_This section documents the original setup process. Current users don't need to follow these steps._

### Prerequisites:

- Wolfram Engine 14.3 installed on host at `~/.local/wolfram/`
- Valid mathpass license file with machine-specific activation
- Wolfram ID account for cloud activation (if needed)

### Setup Process:

1. Configure devcontainer.json with required mounts
2. Add activation restoration to postCreateCommand
3. Create backup of activation tokens in mounted userbase
4. Test both WolframKernel and WolframScript functionality
5. Verify persistence across container rebuilds

### Mount Configuration:

```jsonc
"mounts": [
  "source=${localEnv:HOME}/.local/wolfram/engine/14.3,target=/home/vscode/.local/wolfram/engine/14.3,type=bind",
  "source=${localEnv:HOME}/.local/wolfram/userbase,target=/home/vscode/.local/wolfram/userbase,type=bind",
  "source=/etc/machine-id,target=/etc/machine-id,type=bind,readonly",
  "source=${localEnv:HOME}/.cache/Wolfram,target=/home/vscode/.cache/Wolfram,type=bind"
]
```

### Automatic Configuration (postCreateCommand):

- Create `.WolframEngine` symlink to userbase
- Generate `WolframScript.conf` with correct kernel path
- Link license file to all expected locations
- Create wolframscript binary symlink
- Restore activation tokens from backup if available

This setup provides a robust, maintenance-free Wolfram Engine environment that survives all container rebuilds while avoiding the overhead of repeatedly downloading and installing the 6.7GB engine.
