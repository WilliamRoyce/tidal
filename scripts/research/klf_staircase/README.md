# KLF staircase prototype (GH #473 stage 1 — research artifact)

Status: **checkpoint fired 2026-08-26 — gate (c) unachievable at float64;
no production code shipped.** Preserved because the implementation is
*verified correct* and is the instrument for both recorded fallbacks
(symbolic gauge-structure validation, GH arc-F issue; mixed-precision
structure pass, GH #470). Full evidence record: GH #473 (comment of
2026-08-26) and `docs/tex/pencil_engine.tex`.

## What this is

A faithful numpy port of the GUPTRI-style Kronecker-like staircase from
MatrixPencils.jl (Andreas Varga; `_preduceBF!`, `_preduce2!`,
`_preduce3!`, `_preduce4!` SVD variants + the `klf_right!` driver), with
a TIDAL-specific composition (`compose`/`klf_deflate`): trailing
(infinite + left) states forced identically to zero; right-block rows
solved min-norm via a nilpotent fixed point; the free (gauge) part =
span of the right block's polynomial null family; evolution generator,
manifold projector and pinned basis in the shape the pencil engine
expects.

## Verification status (all measured)

- Kernel equivalence `Q' A0 Z = M`, `Q' B0 Z = N` exact (~1e-15) on
  random complex windows, including every rank-deficient branch
  (`klf_dbg2.py`).
- Structural postconditions (E upper-triangular, B/C/D compressed
  shapes) exact on random deficient windows.
- Driver recovers a known scrambled KCF exactly: blocks L1 x3 + L2 x2 +
  finite x3 + J2 x2 + J1 x2 + L1' x2 + L2' give
  `nur=[5,2,0] mur=[5,5,2] nf=6` as theory predicts.
- Gate (a) synthetics (invertible-B, L1 orphan pair, duplicate rows
  #465, index-2 chain, left chain) at machine precision
  (`klf_prototype.gate_a` with `klf_deflate` swapped in).

Three bugs were found and fixed against the reference during
verification: a Q-update conjugation error in the paired Givens of
`_preduce3/_preduce4` (invisible on real data), a missing early-return
in `_preduce4` when C already has full column rank, and
tolerance-semantics (absolute tau on row-equilibrated pencils, not
`rtol*opnorm`).

## Why gate (c) failed (the finding, in one paragraph)

On the localized implicit-dynamical class (E.cal dual-Gaussian) the
gauge-determining content B0(x)^2 spans ~11 decades with no gap. Any
tau deep enough to preserve the physical structure leaves the
weakly-determined far-field directions *live* with noise growth rates
(maxRe +0.5..+4.5, fluctuating with tau and N) and nonzero cross
coupling (1e-4..0.3) into the physical block — certified (contract
~1e-10, physical {h_5, a_1} block exact to ~1e-13) but unusable at
t_end=40. Any tau shallow enough to classify those directions as free
lets the right-structure staircase march through the whole E-block and
consume the physical sector (nf ~ 10 << 64). Graded two-stage recursion
and two-sided (Ward-style) balancing do not move the pincer: the weak
pairs carry O(1) E-content through their velocity identities, so they
can never present as right (no-E) chains at any float64 tolerance. The
deflation contract certifies consistency with the tau-truncated
equations; on a determination continuum, consistency does not imply
rate fidelity in the weak directions.

## Files

- `klf_port.py` — the faithful port (kernels, driver, composition).
- `klf_prototype.py` — earlier clean-room prototype + `gate_a`
  synthetics (kept for the gate harness; its own staircase was refuted
  by gate (c) and superseded by the port).
- `klf_gate_c.py N TAU` — captures the E.cal pencil via the modal
  builder and runs `klf_deflate` (single-stage gate (c) probe).
- `klf_gate_c2.py N TAU1 TAU2` — graded two-stage recursion probe.
- `klf_dbg2.py` — kernel equivalence/structure harness on random
  windows.

Run from the repo root, e.g.
`uv run python scripts/research/klf_staircase/klf_gate_c.py 16 1e-8`.
Requires the E.cal spec `examples/data/gertsenshtein_ungauged_e_dual_gaussian.json`.
