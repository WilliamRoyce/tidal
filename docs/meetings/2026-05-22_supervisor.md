# Supervisor Meeting — 22 May 2026

**Period**: 8 May (last meeting) to 22 May 2026

---

## Corrections from the last meeting

- **Modal solver runtime.** In the 8 May meeting I quoted the per-evaluation cost of the modal solver as on the order of seconds. That was wrong: it is on the order of **milliseconds**.

- **Infinite vs localised setup.** I suggested the localised (wavepacket) setup as a natural next step without fully accounting for the cost. In a localised background the convolution means spatial Fourier modes do not decouple, so the evolution matrix is dense in k-space — significantly larger than the block-diagonal system of the infinite (plane-wave / periodic) setup. Individual likelihood evaluations are correspondingly slower, which compounds badly over the thousands of samples a nested-sampling run requires. The infinite setup was deliberately chosen for inference speed, not as an oversight.
