# `tidalcosmo/background/` — the CAMB seam

> **Preliminary — planning-stage draft (H4, 2026-08-31).** This directory and its stated
> responsibility come from a design study written *alongside* the detailed investigations
> (H3 solver, H6 polology) rather than after them. It records why the boundary was drawn here so
> a later reader can weigh it — **it is expected to be revised or replaced.** Changing it does not
> require re-litigating H4.

**Responsibility.** Everything the engine is allowed to know about the background, and the only
place CAMB is called.

- `protocol.py` — the interface the engine sees. This is **H3's swap point**: whether the backend
  is stock CAMB, a patched CAMB, or DISCO-EB changes what sits behind this protocol and nothing
  else (`repo_reshape.md` §2.9).
- `camb_seam.py` — the verified calls of §2.3: `get_background_time_evolution` for `a`, `H`,
  `x_e`, `opacity`, `visibility`; `get_time_evolution` for initial conditions and untouched
  sectors **only**. Reading a standard mode that our sector couples to, while the coupling acts,
  is a correctness bug — see §2.2.
- `TabulatedBackground.py` + `.yaml` — **Seam A**: the optional `(a, ρ, P)` hook, off by default.

**Workstream.** WS5, with the Seam-A hook under WS0/#498. **Filled at.** M1 (`camb_seam.py`),
M2 (`TabulatedBackground`).

**Blocked on.** `TabulatedBackground` **cannot be built before the CAMB fork exists** — stock CAMB
offers only the pole-prone `set_w_a_table`. The fork is re-applied from the design of
`slegner/CAMB@2fb908af` onto the upstream `2.0.3` **tag** (H1 §7, GH #498), and it is a
deliverable in its own right, not a footnote. CAMB's repository carries no standard SPDX license
identifier, so **reading its license is a precondition**.

**Open question.** The gauge assertion of §2.8 lives at this seam: the spec's declared gauge is
passed as `frame=`, and a mismatch must fail loudly rather than produce a plausible wrong number.

---

**Standing design goal — conform to external conventions natively.** This package adopts the
naming, formats, gauge conventions and interchange of the tools it interoperates with — **CAMB
and PSALTer** — from the start. Legacy TIDAL notation has no claim here: there is no backward
compatibility to preserve, and no conversion layer to build. See
`docs/cosmology/repo_reshape.md` §2.8.
