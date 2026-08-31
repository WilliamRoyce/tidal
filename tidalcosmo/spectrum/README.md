# `tidalcosmo/spectrum/` — numerical polology — a gate on the sampling path

> **Preliminary — planning-stage draft (H4, 2026-08-31).** This directory and its stated
> responsibility come from a design study written *alongside* the detailed investigations
> (H3 solver, H6 polology) rather than after them. It records why the boundary was drawn here so
> a later reader can weigh it — **it is expected to be revised or replaced.** Changing it does not
> require re-litigating H4.

**Responsibility.** Deciding whether a coupling point is a *healthy* theory — free of ghosts and
tachyons — by screening the particle spectrum in vacuo.

**This is not an isolated corner.** The program document calls WS6 "independent"; it is not, and
`repo_reshape.md` §6 contradicts that deliberately. The spectrum answers exactly the question the
sampler needs answered *before* spending a forward solve: H2's dependency graph draws it as
`WS6 -.gates.-> O2` and `-.gates.-> O3`, and PSALTer's own `_extract/likelihood.py` treats the
spectrum as a likelihood ingredient. So this feeds `validity/` and the Cobaya prior surface, and
its output lands inside the likelihood loop.

Independent only in its **derivation** (Minkowski-only is correct and sufficient — the spectrum
screens in vacuo, and TorC used the same split) and in its **build order** (any time after M0).

**Workstream.** WS6 (#495). **Filled at.** M∥, wired into `validity/` before the first rung runs.

**Basis.** Barker's numerical polology (`psalter.tar.gz`, local; arXiv:2606.30785). D6 grants
explicit permission to copy — provenance in docstrings, attribution settled at publication. Not
clean-room. Start from the #360 plan, noting its Minkowski-first *project-wide* scope decision is
superseded.

**Design consequences.** Fast enough to run per sample (a performance requirement, per D5); a
verdict with a reason, not a boolean; rejection through the shared flagged-rejection mechanism.

**Gate.** Reproduce the Lin–Hobson–Lasenby inequalities exactly, and agree with the supplementary
materials implementation on pole masses and residues.

---

**Standing design goal — conform to external conventions natively.** This package adopts the
naming, formats, gauge conventions and interchange of the tools it interoperates with — **CAMB
and PSALTer** — from the start. Legacy TIDAL notation has no claim here: there is no backward
compatibility to preserve, and no conversion layer to build. See
`docs/cosmology/repo_reshape.md` §2.8.
