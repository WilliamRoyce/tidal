# xAct Test Scripts Quick Reference

## 📋 Available Test Scripts

### Individual Package Tests

#### 🔹 **xTensor** - Core Tensor Algebra

```bash
.devcontainer/tests/test-xtensor.wls
```

**Tests:** Manifolds, metrics, tensors, index manipulations, covariant derivatives
**Duration:** ~30 seconds  
**Use case:** Verify abstract tensor algebra framework

#### 🔹 **xCoba** - Coordinate Calculations

```bash
.devcontainer/tests/test-xcoba.wls
```

**Tests:** Charts, metric components, coordinate transformations, Christoffel symbols
**Duration:** ~45 seconds  
**Use case:** Verify component-based calculations

#### 🔹 **xPerm** - Permutation Groups

```bash
.devcontainer/tests/test-xperm.wls
```

**Tests:** Basic permutations, group generation, MathLink performance, tensor symmetries
**Duration:** ~30 seconds  
**Use case:** Verify symmetry optimization (checks if MathLink is working)

#### 🔹 **xPert** - Perturbation Theory

```bash
.devcontainer/tests/test-xpert.wls
```

**Tests:** Parameter definitions, metric perturbations, systematic expansions, gauge theory
**Duration:** ~30 seconds  
**Use case:** Verify perturbative methods framework

### Integration Test

#### 🔹 **Full Integration** - Real GR Application

```bash
.devcontainer/tests/test-integration.wls
```

**Tests:** Schwarzschild solution with all packages working together
**Duration:** ~60 seconds  
**Use case:** Verify complete computational environment

### Master Test Runner

#### 🔹 **Complete Test Suite**

```bash
.devcontainer/tests/test-all-xact.sh
```

**Runs:** All individual tests + integration test with timing and error handling
**Duration:** ~3-4 minutes total  
**Use case:** Full system verification

## 🚀 Quick Commands

```bash
# Test everything at once
.devcontainer/tests/test-all-xact.sh

# Test specific functionality
.devcontainer/tests/test-xtensor.wls    # Abstract tensors
.devcontainer/tests/test-xcoba.wls     # Coordinates
.devcontainer/tests/test-xperm.wls     # Permutations (MathLink status)
.devcontainer/tests/test-xpert.wls     # Perturbations

# Full realistic application
.devcontainer/tests/test-integration.wls

# Rebuild xPerm if needed
.devcontainer/build-xperm.sh
```

## 📊 Expected Results

| Test            | Expected Output                       | Key Success Indicators                         |
| --------------- | ------------------------------------- | ---------------------------------------------- |
| **xTensor**     | Manifold creation, tensor definitions | ✅ Dim[Spacetime] = 4, metric signature        |
| **xCoba**       | Component calculations                | ✅ Schwarzschild & Minkowski metrics set       |
| **xPerm**       | MathLink connection status            | 🚀 "Connection established" (high performance) |
| **xPert**       | Perturbation parameters               | ✅ Multi-parameter expansions                  |
| **Integration** | Schwarzschild analysis                | ✅ Complete GR environment                     |

## 🔧 Troubleshooting

If tests fail:

1. Check Wolfram activation: `.devcontainer/check-wolfram.sh`
2. Rebuild xPerm MathLink: `.devcontainer/build-xperm.sh`
3. Verify package installation: `ls $WOLFRAM_USERBASE/Applications/xAct/`

The test scripts are designed to be:

- **Self-contained** - Each script loads its own packages
- **Educational** - Shows key functionality with explanations
- **Diagnostic** - Clearly indicates what's working/failing
- **Realistic** - Uses actual GR computations, not toy examples
