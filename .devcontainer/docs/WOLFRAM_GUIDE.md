# Wolfram Engine — the one setup path

**This file is the single source for Wolfram setup.** `.devcontainer/README.md`,
`.devcontainer/QUICKREF.md`, the root `README.md` and `scripts/README.md` all point here and
must not contradict it. Three paths used to disagree and none of them reached PSALTer (GH #559).

Required only for `tidal derive` — deriving linearized field equations from a Lagrangian.
Everything downstream of a JSON specification runs without any of this.

---

## The certified configuration

**Wolfram 14.3.0 × xAct 1.3.0 bundle × PSALTer v2.0.2 `bb45adb0` × the two Function Repository
resources registered locally under fixed UUIDs.**

"Set up correctly" has exactly one definition:

```bash
bash scripts/verify-wolfram-setup.sh --require-psalter    # exit 0
```

That asserts *provenance*, not just function: the engine running from its mount, the four xAct
package `$Version` strings, PSALTer at its pinned commit, and the two resources resolving to
their certified identities on both the master kernel and a subkernel.

It is **independent of any Wolfram Cloud login** — asserted both logged in and logged out. See
[Activation is not a cloud login](#activation-is-not-a-cloud-login).

> The certified xAct is the **code**, not the tarball label: xTensor 1.3.0, xPerm 1.2.4,
> xCore 0.6.10, xCoba 0.8.6, measured from the installed `.m` headers. Older documents said
> "1.2.1"; that was an installer default copied between files, never a measurement.

---

## First-time setup — seven steps

Roughly 30 minutes, most of it downloading. You need a free Wolfram ID
(<https://account.wolfram.com/login/create>) and about 8 GB of disk.

### 1. Download the installer

Get `WolframEngine_14.3.0_LIN.sh` from <https://www.wolfram.com/engine/> (a Wolfram account is
required to download) and put it in `third_party/`:

```text
third_party/WolframEngine_14.3.0_LIN.sh     # 1,750,578,010 bytes (~1.6 GiB)
```

### 2. Install the engine onto the mount

```bash
bash scripts/install-wolfram-engine.sh
```

Installs to `~/.local/wolfram/engine/14.3`, which is bind-mounted from your host, so it survives
container rebuilds. It needs no `sudo` (only the system-library step escalates, and the dev
container has already installed those). It expands to about 6.7 GB.

**Not `/usr/local`.** An engine there is wiped on every rebuild, and `verify-wolfram-setup.sh`
rejects a kernel outside the mount — because `wolframscript` answers from the *cloud* when its
pinned kernel is missing, and `1+1` passes there too. The installer now refuses `/usr/local`.

### 3. Activate, once, with your own Wolfram ID

```bash
wolframscript -activate          # or: bash scripts/activate-wolfram.sh
```

Interactive, one time. The engine writes `mathpass` into the mounted userbase, so you do **not**
re-activate after a rebuild. This is a license, not a cloud login.

If the container was created before the engine existed, wire it up now:

```bash
bash .devcontainer/scripts/setup-wolfram-links.sh
```

### 4. Install xAct

```bash
bash scripts/install-xact-xcoba.sh          # xAct 1.3.0, ~15.6 MB, into the userbase
```

This also rebuilds xPerm's MathLink binary — xAct's pre-built one needs a newer GLIBC than this
image has — and that is **all you need**: measured on a fresh userbase 2026-09-12, it detected
the dependency problem, recompiled with the engine's own `mcc`, and the result connected.
`build-xperm.sh` is a **fallback**, not a second step; see [xPerm](#xperm).

### 5. Install PSALTer and register its resources

```bash
bash scripts/install-psalter.sh
```

Installs PSALTer at `bb45adb0` and registers the two undeclared Function Repository
dependencies it needs. **Registration is not optional**: without it PSALTer does not fail — it
completes and writes a *silently wrong* spectrum (empty source constraints, zero
pseudo-determinants).

### 6. Install xPand — the FRW derivation engine

```bash
bash scripts/install-xpand.sh
```

Adds **only** `Applications/xAct/xPand/` (0.4.4, tarball pinned by sha256) and stamps it with an
`INSTALLED_VERSION`. Re-running it on a machine that already has 0.4.4 is a no-op; a different
version is refused rather than overwritten.

xPand carries the 3+1 split, the scalar–vector–tensor decomposition and the CAMB-named gauges
that the FRW derivation is built on (**D-C = O1′**, 2026-09-17; the evidence is
`docs/cosmology/perturbation_tooling.md`). Step 7 therefore reports **degraded** without it: the
spectrum branch still runs, FRW derivations do not.

### 7. Verify

```bash
bash scripts/verify-wolfram-setup.sh --require-psalter
```

Exit **0** = certified. Exit **2** = works but degraded (not certifiable). Exit **1** = broken.

---

## Supported hosts

Linux, macOS, or Windows **via WSL** — open the folder from the WSL side, not from a
Windows-side window. The three Wolfram mounts are rooted at the host's `$HOME`, and
`initializeCommand` creates them on the host before the container starts; it refuses loudly if
`$HOME` is unset, because the dev container spec leaves unset variables *blank*, which would
resolve the mount sources to the host filesystem root.

---

## What persists, and why

| Path | Kind | Holds |
| --- | --- | --- |
| `~/.local/wolfram/engine/14.3` | host bind | the engine itself (~6.7 GB) |
| `~/.local/wolfram/userbase` | host bind | `mathpass`, xAct, PSALTer, paclets |
| `~/.cache/Wolfram` | host bind | wolframscript cloud tokens (optional) |
| `/etc/machine-id` | host bind, read-only | license validation |
| `~/.claude` | named volume | Claude Code data |
| `~/.Wolfram` | named volume | the **resource registry** (`Objects/`), and `Logs/` |

`~/.Wolfram` is a container overlay directory: without the `wolfram-objects` volume the registry
is lost on every rebuild and nothing re-creates it, which is the silently-wrong-spectrum
condition above. The volume shadows the whole directory, so `Logs/` persists too.

**The first rebuild after this change starts with an empty registry.** That is exactly why
`postCreateCommand` ends with `scripts/psalter/ensure_registered.sh`, wrapped so it can never
abort container creation. Re-run it by hand at any time; it is idempotent.

`postCreateCommand` no longer runs `chown -R` over your home: recursive descent through a bind
mount rewrites ownership on the **host**, and walked the 6.7 GB engine tree every rebuild. If
the mounts ever come up unwritable (a host UID that does not match the container's `vscode`),
recover from inside the container with:

```bash
sudo chown -R vscode:vscode ~/.local/wolfram
```

---

## Activation is not a cloud login

Two different things are easy to confuse:

1. **`mathpass` — the license.** Written by `wolframscript -activate` into
   `~/.local/wolfram/userbase/Licensing/mathpass`. Machine-ID-bound, works entirely offline,
   and survives rebuilds because the userbase is mounted. `setup-wolfram-links.sh` links it into
   every location the engine consults once it exists.
2. **Cloud tokens — optional, and unused here.** `~/.cache/Wolfram/WolframScript/` holds Wolfram
   Cloud authentication. **Nothing in the TIDAL pipeline needs them**: `WolframScript.conf` pins
   `WOLFRAMSCRIPT_KERNELPATH` to the mounted kernel, so evaluation is local.

`.activation_backup` in the userbase backs up the **cloud tokens only** — never the license.
Earlier revisions of this guide described it as restoring "activation" on rebuild; it does not,
and it never needed to: `mathpass` persists because the userbase is a mount.

---

## xPerm

xAct ships a pre-built `xperm.linux.64-bit` needing a newer GLIBC than this image provides
(Debian GLIBC 2.36), so it has to be rebuilt. **Step 4 already does this** — you only need what
follows if `verify-wolfram-setup.sh` reports the xPerm external executable *not* connecting.

**What step 4 does, measured on a fresh userbase (2026-09-12):** `install-xact-xcoba.sh` finds
the shipped binary has unmet dependencies, resolves the engine's own MathLink compiler
(`<engine>/Executables/mcc` — reachable from `dirname "$(command -v wolframscript)"`, though
**not** after `readlink -f`, which lands in a directory with no `mcc`), recompiles, and the
result connects with **no wrapper**. That settles a question this guide previously left open.

**The fallback**, if that ever fails:

```bash
bash .devcontainer/scripts/build-xperm.sh
```

Preconditions: the engine installed (it needs `mprep` from the MathLink DeveloperKit), `gcc`
(installed if missing), and a set `TERM` — it colors its output and exits on an unset one, so
run it from a real terminal. It takes a different route from step 4: `mprep` on xAct's own
`xperm.tm`, then `gcc`, then a small wrapper putting the MathLink shared libraries on
`LD_LIBRARY_PATH`. That wrapper is what the *reference* machine runs, and it is why the
certified install has one where a fresh install does not — the `mcc` path was broken until
2026-09-12 and had never run.

It **reproduces the build**, not a byte-identical file: the wrapper it writes resolves `$HOME`
at run time, which is more portable than an older one you may find already installed.

Either way, success shows up in `verify-wolfram-setup.sh` as the xPerm external executable
connecting.

---

## Verifying and diagnosing

```bash
bash scripts/verify-wolfram-setup.sh                     # 11 checks
bash scripts/verify-wolfram-setup.sh --require-psalter   # the certification gate
tidal doctor                                             # the same diagnosis from the CLI
```

`.devcontainer/scripts/validate-setup.sh` and `check-wolfram.sh` now redirect here. Two of their
old checks have **no equivalent and are deliberately gone**: `wolframclient` (an optional Python
extra that nothing in TIDAL imports — the pipeline shells out to `wolframscript`) and the VS Code
extension check (which self-skips outside VS Code, and is owned by
`.devcontainer/scripts/install-extensions-final.sh`).

```bash
# what the license and machine ID say
cat ~/.local/wolfram/userbase/Licensing/mathpass
cat /etc/machine-id

# activation token backup/restore
bash .devcontainer/scripts/wolfram-activation-manager.sh status
```

`MATHEMATICA_HOME` may appear in your environment, set by the VS Code Wolfram extension rather
than by anything in this repository. **Nothing in TIDAL reads it**; a stale value there is
harmless and is not worth chasing.

### If a PSALTer run hangs instead of failing

**Set `QT_QPA_PLATFORM=offscreen`.** PSALTer exports a PDF through the Wolfram front end
unconditionally and without a time limit, starting at the *first field declaration*
(`DefField`) — seconds into any run. The front end needs a Qt **platform plugin** whose
shared libraries are all present; when none can be initialized, Qt aborts the front-end
process and the export **blocks indefinitely rather than erroring**. Measured: 20.001 s
against a 20 s cap versus **3.8 s** with the variable set.

```bash
QT_QPA_PLATFORM=offscreen wolframscript -file your_script.wls
```

**Do not rely on Inkscape to fix this, even though it appears to.** Of the nine platform
plugins Wolfram ships, `offscreen` was the only one with all its dependencies satisfied here;
`xcb` is missing five libraries. Installing Inkscape pulls in **`libwayland-egl1`**, which
makes the `wayland-egl` plugin loadable and fixes the hang *by accident*. That is luck:
`--skip-inkscape`, a slimmer image or another distribution puts you straight back at the
hang. PSALTer's own `$InkscapePath` is genuinely unused — it is only the dependency that
matters.

A hang is worse than an error here: it holds the single-license Wolfram lane indefinitely and
looks exactly like a slow theory. See `docs/cosmology/stage1_measurements.md` §2.3–2.4.

### If wolframscript answers but the answers look wrong

Check it is not evaluating in the cloud. The image ships a standalone
`/usr/bin/wolframscript`, so `command -v wolframscript` succeeding proves nothing about the
engine — `1+1` passes in the cloud too. The engine test is the kernel file:

```bash
ls -l ~/.local/wolfram/engine/14.3/Executables/WolframKernel
cat ~/.config/Wolfram/WolframScript/WolframScript.conf   # KERNELPATH must point at that file
bash .devcontainer/scripts/setup-wolfram-links.sh        # rewrites the conf and the links
```

### Complete reset

Rebuild the container (`Ctrl+Shift+P` → "Dev Containers: Rebuild Container"). The engine,
license and packages all live on mounts, so nothing is reinstalled; afterwards run
`bash scripts/verify-wolfram-setup.sh --require-psalter`.

---

## Usage

```bash
wolframscript -code "Integrate[x^2, x]"
wolframscript -file mycode.wls

# xAct
wolframscript -code '<<xAct`xTensor`; DefManifold[M,4]; Print[Dim[M]];'
```

Installed under `$WOLFRAM_USERBASE/Applications/xAct/` (~184 MB): xTensor, xPerm, xCoba, xPert,
xTras, Spinors, TexAct, PSALTer and more.

**One `wolframscript` session at a time** — the free license permits a single kernel. Never run
`tidal derive` in parallel.

---

## Repo hygiene (never commit)

Wolfram state lives in `~/.local/wolfram/...`, `~/.cache/Wolfram/...` and `~/.Wolfram/...`, all
outside the repository, so it never appears in `git status`. `.gitignore` additionally covers
license artifacts and installer bundles (`*.mathpass`, `.Wolfram*`, `.WolframEngine*`,
`Wolfram*.{sh,run,tgz}`).

---

## The mount configuration, for reference

```jsonc
"mounts": [
  "source=${localEnv:HOME}/.local/wolfram/engine/14.3,target=/home/vscode/.local/wolfram/engine/14.3,type=bind",
  "source=${localEnv:HOME}/.local/wolfram/userbase,target=/home/vscode/.local/wolfram/userbase,type=bind",
  "source=/etc/machine-id,target=/etc/machine-id,type=bind,readonly",
  "source=${localEnv:HOME}/.cache/Wolfram,target=/home/vscode/.cache/Wolfram,type=bind",
  "source=claude-code-data,target=/home/vscode/.claude,type=volume",
  "source=wolfram-objects,target=/home/vscode/.Wolfram,type=volume"
]
```

Container creation then: creates the host directories (`initializeCommand`, on the host), takes
ownership of the two named volumes, wires the engine and license
(`.devcontainer/scripts/setup-wolfram-links.sh`), restores Claude memory, reindexes sessions,
installs the Wolfram LSP paclets, and re-registers the PSALTer resources. Every Wolfram step is
guarded: a container whose engine is not installed yet still finishes creating, and prints
steps 1–3 above.
