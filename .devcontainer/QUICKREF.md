# .devcontainer Quick Reference

## 📁 Directory Map

| Folder      | Purpose         | Key Files                                                     |
| ----------- | --------------- | ------------------------------------------------------------- |
| `/docs/`    | Documentation   | WOLFRAM_GUIDE.md, XACT_TESTS.md                               |
| `/scripts/` | Utility scripts | check-wolfram.sh, install-extensions-final.sh, build-xperm.sh |
| `/tests/`   | Test suite      | test-all-xact.sh, test-xtensor.wls, test-xcoba.wls, etc.      |

## 🚀 Common Commands

### Run All Tests

```bash
.devcontainer/tests/test-all-xact.sh
```

### Check System Health

```bash
bash .devcontainer/scripts/check-wolfram.sh
```

### Manage Licensing

```bash
bash .devcontainer/scripts/wolfram-activation-manager.sh status
bash .devcontainer/scripts/wolfram-activation-manager.sh backup
bash .devcontainer/scripts/wolfram-activation-manager.sh restore
```

### Install Extensions

```bash
bash .devcontainer/scripts/install-extensions-final.sh
```

### Run Individual Tests

```bash
wolframscript .devcontainer/tests/test-xtensor.wls
wolframscript .devcontainer/tests/test-xcoba.wls
wolframscript .devcontainer/tests/test-xperm.wls
wolframscript .devcontainer/tests/test-xpert.wls
wolframscript .devcontainer/tests/test-integration.wls
```

## 📖 Documentation

- **Setup Guide**: `.devcontainer/docs/WOLFRAM_GUIDE.md`
- **Test Guide**: `.devcontainer/docs/XACT_TESTS.md`
- **Container Config**: `.devcontainer/README.md`

## 🔧 Scripts Reference

| Script                          | Purpose                    | Usage                                                                                |
| ------------------------------- | -------------------------- | ------------------------------------------------------------------------------------ |
| `build-xperm.sh`                | Compile xPerm from source  | `bash .devcontainer/scripts/build-xperm.sh`                                          |
| `check-wolfram.sh`              | Health check system        | `bash .devcontainer/scripts/check-wolfram.sh`                                        |
| `fix-xperm.sh`                  | Fix xPerm issues           | `bash .devcontainer/scripts/fix-xperm.sh`                                            |
| `install-extensions-final.sh`   | Install VS Code extensions | `bash .devcontainer/scripts/install-extensions-final.sh`                             |
| `wolfram-activation-manager.sh` | Manage licensing           | `bash .devcontainer/scripts/wolfram-activation-manager.sh [status\|backup\|restore]` |

## ✅ What's Working

- ✅ Wolfram Engine 14.3 (fully activated & persistent)
- ✅ xAct 1.3.0 with all packages (xTensor, xCoba, xPerm, xPert)
- ✅ xPerm MathLink (compiled from source, GLIBC compatible)
- ✅ Complete test suite (100% passing)
- ✅ VS Code extensions (auto-installing)
- ✅ Licensing persistence (through rebuilds)

## 📊 Test Status

```
✅ xTensor Core Algebra
✅ xCoba Coordinates
✅ xPerm Permutations (MathLink active)
✅ xPert Perturbations
✅ Full Integration (Schwarzschild computation)

Success Rate: 100% (5/5 tests passing)
```

---

**For detailed information**, see:

- Setup/troubleshooting → `.devcontainer/docs/WOLFRAM_GUIDE.md`
- Testing → `.devcontainer/docs/XACT_TESTS.md`
