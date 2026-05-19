# Supervisor Meeting — 22 May 2026

**Period**: 8 May (last meeting) to 22 May 2026

---

## Corrections from the last meeting

- **Modal solver runtime.** In the 8 May meeting I quoted the per-evaluation cost of the modal solver as on the order of seconds. That was wrong: it is on the order of **milliseconds**.

- **Infinite vs localised setup.** I suggested the localised (wavepacket) setup as a natural next step without fully accounting for the cost. In a localised background the convolution means spatial Fourier modes do not decouple, so the evolution matrix is dense in k-space — significantly larger than the block-diagonal system of the infinite (plane-wave / periodic) setup. Individual likelihood evaluations are correspondingly slower, which compounds badly over the thousands of samples a nested-sampling run requires. The infinite setup was deliberately chosen for inference speed, not as an oversight.

---

## JAX as a modal-solver backend — what we tried, what we found

Acting on the suggestion to JIT-compile the modal solver, I built a JAX backend (`tidal/solver/modal_jax.py`). The code is correct but the performance result was the opposite of what we expected.

### Performance comparison on the production workload

Canonical case: `dark_photon_plasma`, N=32, 17 Fourier modes, 30×30 dynamical block after Schur elimination.

| Backend               | Median per call | Notes                           |
| --------------------- | --------------- | ------------------------------- |
| scipy `solve_modal`   | 9.2 ms          | baseline                        |
| JAX `solve_modal_jax` | 28 ms           | **3× SLOWER than scipy on CPU** |

Profile of the JAX path on the same case (post-JIT, steady state):

- `_build_evolution_matrices` (numpy Schur, shared with scipy): 2.1 ms
- JAX fused `expm` + `scan` (the part we hoped would win): **9.7 ms**
- numpy ↔ JAX boundary, IFFT, reconstruction: ~16 ms

For comparison, the same expm work done as 17 sequential `scipy.linalg.expm` calls runs in **2.6 ms**.

### Why JAX did not win

On CPU, `jax.scipy.linalg.expm` and `scipy.linalg.expm` both route through the same BLAS library (OpenBLAS / MKL). XLA compiles to the same C/Fortran kernels for the matrix multiplications inside Padé scaling-and-squaring. The Python loop overhead for 17 modes is ~85 μs — about 1 % of total compute. **"JIT compilation" eliminates Python interpreter overhead, not BLAS compute time**, so `vmap(expm)` ends up calling the same 17 BLAS routines, just batched, with no net win.

The regimes where JAX would actually beat scipy:

- **GPU**: fundamentally different compute path — XLA on GPU routinely gives 10–100× on similar workloads. Not available on CSD3 sapphire.
- **Larger grids**: N ≥ 192 in 1D, i.e. ≥ 100 Fourier modes, where XLA's dispatch overhead amortises. Our current production grids are N=64–128.

This is consistent with the published JAX-on-CPU benchmarks. The implementation is sound — the bottleneck is the hardware/library configuration we are running on.

### Productive correction

The ms-not-s correction above grounds the broader strategic picture:

- A PolyChord chain at O(10⁴) likelihood calls × ~30 ms/call ≈ **5 min** for a constant-coefficient theory. For Phase E position-dep theories via CVODE, ~10⁴ × 100 ms ≈ **17 min/chain**.
- These are not the multi-hour wall clocks the "seconds-per-call" framing implied.
- Aggressive per-call optimisation of the modal solver is **not the marginal-return investment we thought it was**. The bottleneck for campaign throughput is somewhere else (Wolfram derivation pipeline, inference-pipeline setup, possibly the analysis/plotting tail).

### Discussion questions

1. **Was the JAX suggestion grounded in a GPU expectation?** If we had GPU access anywhere (cluster other than CSD3 sapphire, cloud), the picture flips — we'd expect a real 10–100× from XLA on GPU. We did not investigate that path because CSD3 sapphire is CPU-only.
2. **Is there a problem-size regime we should be designing for** that we haven't reached yet? JAX wins above ~100 Fourier modes; if Phase F+ takes us to N=256+, re-enabling JAX becomes worth it.
3. **Given the ms-per-call cost** (not seconds), what is the realistic performance target for Phase E inference? At ~17 min/chain for position-dep theories via CVODE, is that acceptable, or does it justify the harder modal optimisations described below?

---

## Extending the parameter space: sampling over the base-theory coefficients

The Einstein-Cartan gravity term $\frac{1}{\kappa^2}\tilde R$ and the Maxwell kinetic term $-\frac{1}{4}F_{\mu\nu}F^{\mu\nu}$ — have been held fixed at their standard values while we sampled over the torsion and nonminimal couplings.

The effective-field-theory perspective says this is not obviously justified. If we regard the Lagrangian as a truncated expansion, integrating out heavy fields will in general shift the coefficients of _all_ operators at a given order, including the dimension-4 base terms.

The Maxwell coefficient is similarly renormalised by any loop involving charged matter or by EFT operators like the Euler–Heisenberg $F^4$ term. An $O(1)$ shift in the Maxwell coefficient could in principle amplify or suppress the Gertsenshtein channel in a way that mimics torsion structure.

**Concrete proposal:** add $\kappa^{-2}$ and the Maxwell coefficient as two additional free parameters in future inference runs, and check whether the nulls in the propagating-torsion sector survive marginalisation over these.

**Question for supervisors:** is there a physical argument for holding these fixed that we are missing, or would sampling over them be a straightforward extension worth including in the paper?
