# New-package design: the Cobaya extension, from the goal backwards (H4)

**Executed:** 2026-08-31. **Status:** complete — design study, plus the scaffold it specifies.
**Workstream:** WS1 (#490), under umbrella #488.
**Companion artifacts:** the `tidalcosmo/` package tree (directories and READMEs only, no code).

## What this document is for

The program (`docs/COSMOLOGY_PROGRAM.md`) is building a **Cobaya extension**: a `Theory` class
that chains off upstream CAMB, evolves a candidate Lagrangian's perturbations as spectators on a
ΛCDM background, and returns CMB observables. WS1 decides what that package looks like.

The governing constraint is D3 plus the WS1 brief: **the existing framework is legacy, not a
template.** It is acknowledged as naive and gradually patched-on. The new package is designed
whole, from the goal backwards; **new code never imports old code**; useful capabilities are
*fully ported* — redesigned to fit, with docstrings, why-comments and issue references traveling
along — and the original stays in place as a test oracle until it is deliberately retired.

**The symbolic physics is not the problem; the wiring is.** The derivation pipeline works. What
does not work is its packaging: the forward model lives inside `tidal/cli/` (25,698 of ~67k
lines; `_simulate.py` alone is 3,126), four modules outside the CLI reach *back into* it for
private names — each suppressed with `# pyright: ignore[reportPrivateUsage]` — and the
configuration object is an `argparse.Namespace`, because the CLI *is* the config. So this port
is a restructuring of working physics, and what must survive the move is as much the
*explanation* as the code.

### Who this is written for

A theoretical physicist — field theory, Poincaré gauge theory, the Gertsenshtein effect — who is
**new to computational cosmology**. Nothing here explains Lagrangians, gauge symmetry, or the
physics of the couplings. Everything here explains the CMB-pipeline machinery: what a Boltzmann
code actually computes, what a transfer function is, what line-of-sight integration means, what
a "channel" is, and which conventions we must adopt from other people's codes. §0 is a glossary;
terms are also defined at first use. This document is meant to be a way into the field as much
as a design record.

---

## 0. Glossary

- **T, E, B** — the three CMB maps. **T** is the temperature anisotropy. **E** and **B** are the
  two parity components of *linear* polarization (E is curl-free, B is divergence-free).
  **V** is *circular* polarization, essentially absent in ΛCDM — which is exactly why it is a
  distinctive channel for new physics.

- **`C_ℓ`** — the angular power spectrum: how much variance a map has at angular scale ≈ 180°/ℓ.
  This is the quantity every CMB likelihood consumes. Everything upstream exists to produce it.

- **Scalar / vector / tensor (SVT)** — at linear order, perturbations on an FRW background split
  by how they transform under spatial rotations, and the three sectors **evolve independently**.
  Scalars are the density perturbations and gravitational potentials; vectors decay; **tensors
  are the transverse-traceless spin-2 part — gravitational waves.** This is the same irreducible
  decomposition familiar from field theory, applied to the perturbed metric. It is why a
  Boltzmann code can offer `WantScalars` and `WantTensors` as independent switches at all.

- **Channel** — *an observable pathway*, not one of our fields. A channel bundles "which modes
  are involved, how they evolve, and which `C_ℓ` they land in."

- **The tensor channel**, the worked example throughout. It is called that because the
  perturbations being evolved are the **tensor** sector: primordial gravitational waves. It is
  *observable* because a gravitational wave passing through the plasma at last scattering
  distorts the local photon distribution into a quadrupole, and Thomson scattering off a
  quadrupole produces linear polarization. The distinguishing fact is that at linear order
  **scalar perturbations cannot produce B-modes** — only tensors can (along with lensing and
  foregrounds, which are separately modeled). A primordial B-mode is therefore the clean
  fingerprint of a tensor perturbation. That is what makes "GW → B-mode" a *channel*: a specific
  set of modes, a specific evolution, landing in a specific spectrum (`C_ℓ^{BB}`) that other
  sectors cannot fake.

- **Transfer function** — the linear map from a primordial initial condition at wavenumber `k`
  to the value of some perturbation later. Boltzmann codes factorize the calculation as

  ```text
  C_ℓ = ∫ (dk/k) P_primordial(k) |Δ_ℓ(k)|²
  ```

  so the primordial physics sits in one factor and *all* the evolution and projection sits in
  `Δ_ℓ(k)`. This factorization is what makes it possible to replace one piece of the calculation
  without touching the rest.

- **LOS (line-of-sight) integration** — how `Δ_ℓ(k)` is obtained. Rather than evolving the
  anisotropy multipole by multipole to the present day, the code evolves a small number of local
  quantities, assembles them into a **source function** `S(k, η)` — "how much signal is emitted
  or altered here, now, along the path" — and then projects once:

  ```text
  Δ_ℓ(k) = ∫ dη  S(k, η)  j_ℓ(k(η₀ − η))
  ```

  integrating along the photon's line of sight from recombination to us. **Anything that happens
  to a photon en route enters as a term in `S`.** This is the single most important structural
  fact about a modern Boltzmann code, and most of §2 follows from it.

- **Gauge** — perturbations are defined relative to a chosen slicing of spacetime into surfaces
  of constant "time", and that choice is arbitrary. Different choices (**synchronous**,
  **Newtonian**) give literally different numbers for "the density perturbation" while all
  observables agree. It matters here because our equations and CAMB's must be **in the same
  gauge** before they can be combined. A gauge mismatch does not raise an error; it produces a
  plausible wrong number. §2.8 makes gauge a first-class, checked property of a derived theory.

- **`η` (conformal time)** — the working time variable throughout, `ds² = a(η)²[−dη² + dx²]`. It
  is CAMB's internal variable (`tau` in the Fortran, `t` in `camb.symbolic`), it is what TorC
  used, and it is what the program adopted. That agreement is free and load-bearing: there is no
  time-variable conversion anywhere at the CAMB seam.

---

## 1. Decisions: namespace, home, and the rename

### 1.1 A fresh top-level package, inside this repository

The new package is a **new top-level import root**, `tidalcosmo`, living as a sibling of
`tidal/` in this repository.

**Why a separate root rather than `tidal.cosmo`.** The primary reason is *separability*. A shared
namespace lets the two packages entangle — and the entire strangler-fig strategy depends on being
able to lift the legacy tree out cleanly at the end. Under a separate root, the rule "new code
never imports old code" is mechanically enforceable by a single hygiene test; under a shared one
it is not enforceable at all, because importing anything named `tidal.*` executes the legacy
`tidal/__init__.py` by definition. That is corroborated by measurement — `import tidal` costs
524 ms and pulls in matplotlib via [`tidal/__init__.py:24`](../../tidal/__init__.py) →
`.plot_pgf`, which would land inside every Cobaya `Theory` instantiation — but the measurement is
the supporting argument, not the reason.

**Why the same repository.** The git history, the local `literature/` corpus, and the issue trail
(#488 umbrella, #489–#495 workstreams, #497–#499) are load-bearing working context. A separate
repository forfeits all three, and makes "legacy as test oracle" clumsy exactly when it is needed
most. A later split, once the package stands on its own, remains available and costs little.

### 1.2 `tidalcosmo` is scaffolding, not an identity

The end state is: build `tidalcosmo/` alongside `tidal/`, migrate capability by capability,
delete the legacy `tidal/` tree, then **rename `tidalcosmo/` back to `tidal/`**. The name
`tidal` is already this project's own (`pyproject.toml` `name = "tidal"`,
`[project.scripts] tidal = "tidal.cli:main"`), so the swap is internal bookkeeping, not a
contest for a PyPI name.

Three rules follow, and they bind every session from here on:

1. **No cosmology-specific naming inside the tree.** `tidalcosmo/solver/` and
   `tidalcosmo/spectator/` rename cleanly; something like `tidalcosmo/cosmo_backend/` would read
   absurdly once the prefix is gone. Every directory in §2.11 survives the rename unchanged.
2. **Assume the final names.** The final console script is `tidal`; the final Cobaya dotted path
   is `tidal.SpectatorTheory`. Write docs and examples knowing that.
3. **Do not publish or circulate `tidalcosmo` externally.** No PyPI release, and no YAML shared
   outside the project naming `tidalcosmo.…`, before the rename — any external user's
   `theory: {tidalcosmo.SpectatorTheory: …}` breaks the day we rename. If early external use is
   wanted, do the rename first and publish once.

### 1.3 What the rename actually is, and how to keep the diff readable

`git mv tidalcosmo tidal` moves **just that one directory**, with everything inside it, and
stages the rename in one step (it is `mv` + `git rm` + `git add`). Nothing outside the directory
is touched.

Git stores content, not renames. A rename is **detected** at display time by comparing contents,
which means detection is exact only when a commit contains the move *and nothing else*: the
similarity is then 100%, and `git show` (and GitHub) render

```text
rename tidalcosmo/solver/magnus.py => tidal/solver/magnus.py (100%)
```

with **zero diff lines** — not a rewrite. Mix a content edit into the same commit and similarity
drops; a heavily edited file degrades to delete-plus-add and its history is visually severed.

So the cutover is **three commits, in this order**:

1. delete the last of legacy `tidal/`;
2. `git mv tidalcosmo tidal` — **pure move, zero content edits**;
3. string-replace `tidalcosmo` → `tidal` across imports, `pyproject.toml`, docs and YAML — a
   small readable diff against the already-moved paths.

Afterwards, `git log --follow <file>` traces a single file across the move, `git blame -C -C`
carries blame through it, and `git diff -M` (the default for `diff` and `show`) shows renames.
Verify step 2 was clean with:

```bash
git show --stat --find-renames=100% <sha>   # every path a rename; none add/delete
```

One caveat worth setting up for: rename detection gives up above `diff.renameLimit` files
(default 1000). `tidal/` currently holds 221 files across 13 directories, so there is ample
headroom — but if the tree has grown by then, raise the limit for that commit rather than letting
git silently fall back to delete-plus-add.

Step 2 **fails outright while `tidal/` still exists**. That is the safety property that makes the
capability-gated cutover mechanical rather than a judgment call.

---

## 2. How we interface with CAMB

### 2.1 The claim, and the assumption that licenses it

The question this whole section answers: **how do we change the way some standard perturbations
evolve — photon polarization, the tensor mode — without forking CAMB?**

**We never remove anything from inside CAMB. We supersede at the product level.** CAMB computes
its ΛCDM result as it always does. For the one channel our new sector touches, we compute that
channel ourselves and return *our* version through the Cobaya `Theory`. The likelihood
downstream sees a standard set of products and cannot tell the difference.

**What makes this legitimate rather than a fudge is the spectator assumption itself.** The
program's conditions are that the new sector **(a)** does not affect the background expansion
(negligible in the Friedmann equation) and **(b)** does not gravitationally disturb the standard
perturbations (negligible in the Einstein constraints). Because we assume both, the rest of
CAMB's calculation remains valid while we replace one channel. Drop (b) and the whole
construction is incoherent — we would have to be *inside* the Boltzmann hierarchy, which is
precisely the "gravitational sourcing is not reachable in the strict limit" boundary the program
already draws.

Note the deliberate asymmetry, because it is easy to misread. We **do** modify the propagation of
standard quanta — that is the entire point — through the **explicit coupling terms in the
Lagrangian**, kept at linear order. Condition (b) forbids only *gravity-mediated* sourcing via
the new sector's stress-energy. Coupling-mediated is not gravity-mediated.

**The spectator route is what licenses the no-fork architecture.** Those two decisions are not
independent choices that happen to sit well together; the second follows from the first.

### 2.2 What CAMB's standard variables may and may not be used for

An earlier draft of this study said CAMB's `get_time_evolution` supplies "the standard modes our
block couples to." That was wrong, and the correction is the sharpest point in the design: **if
torsion modifies how photon polarization evolves, then CAMB's photon polarization is wrong for
that channel.** We cannot read it out and use it as though the coupling were not there.

| Use of CAMB's standard variables | Valid? |
| --- | --- |
| **Background and thermal history** — `a`, `H`, `x_e`, `opacity`, `visibility`; the coefficients of our system | **Yes.** These are condition-(a) quantities, unaffected by assumption |
| **Initial conditions**, deep in the past, before the coupling becomes relevant | **Yes** — and this is how the coupled block gets started |
| **Channels our sector does not touch** — e.g. scalar TT when we modify only tensors | **Yes**, untouched by construction |
| **A standard mode our sector couples to, at times when the coupling acts** | **No.** That mode must be *inside our coupled block*, evolved by us, and CAMB's version discarded for that channel |

So **the coupled block is not "our fields". It is our fields *plus* every standard mode they
couple to, evolved together.** That is what "the coupled block" means in the program's
division-of-labor table, and it is worth stating explicitly because the shorter phrasing invites
exactly the error above.

### 2.3 What stock CAMB gives us — read from the source

The following were read from `camb/results.py` and `camb/model.py` on current upstream, not
inferred from documentation.

| Call | Returns | Our use |
| --- | --- | --- |
| `get_background_time_evolution(eta, vars)` | background and thermal history on a conformal-time grid; `background_names` = `x_e`, `opacity`, `visibility`, `cs2b`, `T_b`, `dopacity`, … | the coefficients of our `M(η, k)` |
| `get_time_evolution(q, eta, vars, frame=…)` | standard perturbation variables per `k` over `η`. `evolve_names` includes `a`, `etak`, `H`, `growth`, `v_photon`, `pi_photon`; `vars` accepts **sympy expressions built from `camb.symbolic`**; `frame=` selects the gauge | initial conditions and untouched-sector inputs **only** — see §2.2 |
| `set_custom_scalar_sources(sources, source_names, source_ell_scales, frame, code_path)` | accepts **sympy expressions**, generates and compiles CAMB Fortran (`code_path` dumps the generated `.f90`) | mechanism 1 of §2.4 |
| `get_cmb_unlensed_scalar_array_dict()` | **all** unlensed auto- and cross-spectra *including any custom sources*, keyed `TxT`, `TxE`, `custom_name_1xT`, … | a custom source's `C_ℓ` **and** its cross-spectra with T and E, for free |
| `get_cmb_transfer_data(tp='scalar'\|'vector'\|'tensor')` | `ClTransferData` with `delta_p_l_k` indexed by (source, `L`, `q`) | see the correction below |

Two questions this raised, answered directly:

**Does CAMB always compute all the standard perturbations, wasting time?** The compute switches
are coarse — `WantCls`, `WantScalars`, `WantTensors`, `WantTransfer`, `DoLensing` — and they
toggle **whole SVT sectors**, not individual species. There is **no switch that skips the photon
hierarchy**: the scalar system is monolithic. So if we take over the photon-polarization channel,
CAMB still evolves its own photons and that work is spent.

Two mitigations. Most of the scalar computation we *need anyway* — thermal history, the untouched
channels, the initial conditions. And **`WantTensors = False` genuinely does skip CAMB's tensor
computation**: precisely the tensor *sector* of §0, meaning the primordial gravitational-wave
contribution to the spectra (dominantly `C_ℓ^{BB}`, plus small tensor pieces of TT/EE/TE at low
`ℓ`). It touches nothing scalar. So a **tensor-channel takeover can switch its duplicate off
cleanly; a photon-polarization takeover cannot.**

*How much is actually wasted is deliberately left open.* It is not answerable by reading source —
it needs a measurement of a real forward model with and without the duplicated work. Recorded as
an M1 benchmark with its own issue, for the session that builds the seam. It feeds D5
(optimization as a standing first-class concern) and H3's 10×-versus-100× budget.

**Is `get_cmb_transfer_data` field-agnostic — will it hand us `Δ_ℓ(k)` for torsion?** **No.** It
returns *CAMB's own* sources, indexed by a `source` index over what CAMB computed: temperature,
E polarization, lensing potential. It has no notion of torsion and never will. Our new fields'
transfer functions are something **we** compute. What the call is genuinely good for is a
comparison oracle — CAMB's tensor `Δ_ℓ(k)` in the decoupled limit is exactly what our own tensor
machinery must reproduce, which makes it the natural first test of mechanism 2.

### 2.4 Three mechanisms — and where our own solver fits

These are not alternatives *to* the solver. They are **three different amounts of the calculation
we take over**, and the solver appears in two of them.

| # | Mechanism | Needs our solver? | When it applies |
| --- | --- | --- | --- |
| **1** | **Add a custom source.** Write the new term as a sympy expression in `camb.symbolic` variables; `set_custom_scalar_sources` compiles it into CAMB; spectra come back from `get_cmb_unlensed_scalar_array_dict` | **No.** No new dynamical variable exists, so there is nothing to integrate | the effect is a term in `S(k, η)` expressible in quantities CAMB already evolves |
| **2** | **Replace at transfer level.** We solve the coupled block ourselves to obtain **our** `Δ_ℓ(k)` for that channel, then form `C_ℓ` from it; CAMB's version of that channel goes unused | **Yes** — this is the solver's main job | a standard mode's *evolution* changes, but the way it is projected onto the sky does not |
| **3** | **Own the projection too.** Build `S(k, η)` from our solution and perform our own `∫dη S j_ℓ` | **Yes**, plus our own LOS integrator | the source function itself changes shape |

**"Replace at transfer level", unpacked.** Recall the factorization
`C_ℓ = ∫(dk/k) P_primordial(k) |Δ_ℓ(k)|²`. Under mechanism 2 the primordial factor is unchanged
and the projection is unchanged; only `Δ_ℓ(k)` — which encodes all the accumulated evolution —
differs. So we recompute that one factor and rebuild `C_ℓ` from it.

**"Its projection does not change", unpacked.** The projection is the geometric step: integrating
the source against `j_ℓ(k(η₀ − η))` along the line of sight. Saying it does not change means the
same integral with the same kernel still applies; only the thing being integrated is different.
When *that* stops being true — a new source shape, a different spin weight — mechanism 2 is
insufficient and we are in mechanism 3.

**`source_ell_scales`, and why it matters to us.** The parameter scales a custom source's
contribution at multipole `ℓ` by `√((ℓ+n)!/(ℓ−n)!)`, and CAMB's own docstring says `n = 2` is for
"a new polarization-like source". That factor is the standard spin-2 normalization which turns a
scalar-like source into one that projects onto **polarization** rather than temperature. The
relevance is direct: our headline targets — cosmic birefringence (E→B) and V-modes — *are*
polarization observables. CAMB therefore already supports adding a polarization-like source with
no Fortran work from us, which makes mechanism 1 the cheapest available route for any effect that
fits it.

**Where the solver fits.** Mechanism 1 needs no solver because nothing new is being integrated —
it is a term added to what gets projected. The moment the physics requires a **new dynamical
variable** to be evolved alongside the standard ones, mechanism 1 is out and the solver (WS3) is
the engine. That is the test, and it is a property of the physics, not a preference.

### 2.5 What we evolve, and the source functions we still owe

**We evolve** the coupled block's variables in conformal time `η`, per wavenumber `k`, under the
equations of motion our symbolic pipeline derives, with coefficients evaluated from CAMB's
background. Concretely, per channel:

- **tensor channel** — the two gravitational-wave polarizations `h_+`, `h_×` together with the
  torsion components they mix with;
- **photon-polarization channel** — the photon polarization variables together with the torsion
  components the Lagrangian couples them to.

In both cases, note §2.2: the standard modes are *inside* the block, not read from CAMB.

**We do not evolve `S`.** The source function is *constructed from* the solution, then projected.

**And constructing it is derivation work, per theory — not a lookup.** `S` is obtained by taking
the photon Boltzmann equation, integrating by parts to move the `ℓ`-dependence onto the Bessel
kernel, and grouping what remains into terms multiplying the visibility function. That is how
ΛCDM's temperature source ends up as the familiar Sachs–Wolfe + Doppler + integrated-Sachs–Wolfe
combination. The manipulation is specific to the equations you started from. Ours are not those
equations. Therefore:

- **each channel needs its own `S` derived** from *our* coupled system — "given the solved
  `(h, torsion)` history, what is the quantity integrated along the line of sight to produce a
  B-mode?";
- it is **not in the legacy code** (flat space has no line of sight to integrate along) and
  **not free from CAMB** (its sources are ΛCDM's);
- and by the project's standing rule it should be **derived symbolically in the Wolfram stage**,
  not hand-written in Python — which makes it a **WS2 export requirement** as much as a WS4
  observable, sitting alongside H2's eikonal reduction as a second thing WS2 owes WS4.

This is a substantial, currently unowned deliverable. It gets its own issue, and it should not be
estimated as small.

### 2.6 Forking: what we refuse, what we do, and why that is not a contradiction

Two statements below could look contradictory, so the resolution comes first:

> **The package's default install uses stock CAMB from PyPI. No fork, no Fortran build.**
> Exactly one feature — *optional, off by default, and needed only for the O1 plumbing gate* —
> wants a small patched CAMB. A user who never enables it never encounters a fork. **Nothing on
> the main physics path is forked.**

They are answers to two different questions, which is why both hold.

**Question 1 — do we fork CAMB in order to change how perturbations evolve? No.** This is the
architectural decision. The standard route in this literature is to fork the whole code:
`sfu-cosmo/MGCAMB` and `EFTCAMB/EFTCAMB` are complete CAMB forks — EFTCAMB adds `camb/eftcamb.py`
and `camb/eftcamb_tree.py` *inside* the `camb` package, plus new Fortran modules — and
`hi_class_public` is a CLASS fork. Adding a species that way means new evolution equations inside
the Boltzmann hierarchy, new initial conditions, and new entries everywhere the code sums over
species. It touches the core, and it never converges back to upstream. Four reasons we refuse it:

1. **A fork is not a patch you own; it is a copy of CAMB you own.** Every upstream release — new
   recombination, new neutrino treatment, a lensing bug fix — must be hand-merged through your
   changes, forever.
2. **The drift is empirically fast, and we have the measurement.** `slegner/CAMB`, forked for a
   paper published in 2025, was **180 commits behind upstream** when H1 audited it, carrying
   unrelated post-paper work, a breaking change to a public signature, a leftover debug `print`,
   and a stray settings file. One small feature, one paper, already past easy rebasing.
3. **It contradicts the program's own latest-upstream policy**, which exists so our results stay
   comparable with everyone else's and our users install a normal CAMB. A fork makes our
   package's install *be* a Fortran build of a stale cosmology code.
4. **The ratio is absurd.** Our new sector is a handful of extra linear equations per `k`.
   Forking to host them means adopting ~10⁵ lines of Fortran for ~10² lines of new physics, and
   it puts our symbolic pipeline's output on the far side of a Fortran compile step.

§2.4's three mechanisms are how we avoid all of that: the maintenance surface stays proportional
to the new physics.

**Question 2 — does any feature need a patched CAMB at all? Yes, exactly one, and it is
unrelated to the above.** Seam A feeds CAMB a tabulated `(a, ρ, P)` background, which stock CAMB
cannot do — it offers only `set_w_a_table`, the pole-prone route the TorC paper had to abandon
(H1 §4, §8). That is not "adding a species"; it is **adding one optional output to an existing
routine**. H1 §7 costed it: **two self-contained Fortran modules (~309 lines) plus ~115 lines of
Python surface**, branched from the released `2.0.3` **tag** rather than `master`, adding one
optional argument with no breaking API change. It is small enough to re-apply onto each new
upstream release, and it is an obvious candidate to contribute upstream — at which point our fork
disappears entirely. It is deliberately the same shape as Legner's change, which additionally
gives us their published work as a cross-check (H1 §8's free oracle at `ϖ_r = 0.8`).

So the four objections are objections to the *scale and permanence* of a change-the-hierarchy
fork. None of them bite on a ~400-line optional patch, off a release tag, that we intend to
upstream. And critically: **it is not on the critical path.** O1 is a plumbing gate; the physics
rungs need none of it.

### 2.7 The three seams

| Seam | What it reaches | Fork? | Where it lives |
| --- | --- | --- | --- |
| **A — background** | a tabulated `(a, ρ, P)` driving the expansion | **Yes** — the small optional-output patch of §2.6 (H1 §7 R2, GH #498) | `background/`, optional, off by default |
| **B — sources** | terms in `S(k, η)` expressible in quantities CAMB already evolves — mechanism 1 | **No** | `observables/rotation.py` |
| **C — new evolving fields** | the coupled block of §2.2 — mechanisms 2 and 3 | **No**, but we own the evolution and possibly the projection | `perturbations/` + `solver/` + `observables/los.py` |

Which seam an effect needs is decided by the physics — whether it requires a new dynamical
variable — not by preference, and the boundary moves: an effect that is a constant rotation sits
in Seam B, while the same effect with a time-dependent coefficient across recombination moves
into the LOS integral proper.

**Seam C is the design target.** Cases where `a(η)` drops out — conformal invariance in the
CS/photon sector, for instance — are rungs that happen to be cheap, **never architectural
drivers**. Per D7, cancellations are auto-detected fast paths, never scoping decisions. The
general `a(η)`-dependent path is the default everywhere, and WS2's existing verification gate —
*conformal fast path ≡ general path to machine precision* — is exactly the right relationship
between them.

### 2.8 Conventions: build to CAMB *and* PSALTer natively

> **This is a headline design goal, and later sessions should pick it up from here without being
> told.** We are rebuilding the symbolic layer, so it **adopts the naming, formats, gauge
> conventions and interchange of the tools it must interoperate with — CAMB and PSALTer — from
> the start.** Legacy TIDAL notation has no claim on the new package: there is no backward
> compatibility to preserve, and no reason to encode our own conventions and then translate.

An earlier draft of this study proposed a "variable dictionary" mapping our names onto CAMB's.
That is a conversion layer, and a conversion layer is what you build when you cannot change
either side. We *can* change our side. So we do.

#### CAMB

`camb/symbolic.py` defines ΛCDM's scalar linear perturbation equations in sympy, in covariant
notation, carrying an explicit map to the Fortran variables via `camb_var` and `camb_sub`
metadata — for example `rho_g` ↔ `grhog_t` with substitution `grhog_t/kappa/a**2`. It exposes
`a`, `H` (`adotoa`), every species density and pressure, `opacity`, `visibility` and `exptau` as
named symbolic functions.

Two things follow:

- **The time variable already agrees.** `camb.symbolic` uses `t` for conformal time, which is
  `tau` in the Fortran, which is TorC's `τ`, which is the program's `η`. No conversion layer at
  that seam — a free win worth stating so nobody builds one.
- **There is an existing codegen path to evaluate rather than duplicate.** `camb_fortran()`
  converts a sympy expression to CAMB source, and `compile_source_function_code()` compiles it to
  a function pointer at runtime. That is the supported route by which *a symbolic expression
  becomes something CAMB evaluates*, and it is the natural target for our Seam-B output. WS2/WS3
  should decide whether to emit into it rather than inventing a parallel path.

#### Gauge, as a concrete requirement on the symbolic work

This is the sharpest interop risk in the design, because **a gauge mismatch produces a plausible
wrong number, not an error**. Four requirements:

1. **The gauge is an explicit, named choice**, drawn from the set CAMB names — covariant,
   Newtonian, synchronous — declared as an *input* to the derivation, not left implicit in
   whatever the Wolfram stage happened to produce.
2. **The projection is applied symbolically, in Wolfram**, never patched in Python (standing
   project rule). `camb.symbolic`'s own covariant→Newtonian/synchronous projections are the
   reference for what the target conventions are.
3. **The choice is recorded in the spec** as first-class metadata, so it travels with the
   equations rather than living in someone's memory.
4. **The CAMB seam asserts agreement and refuses on mismatch.** `get_time_evolution` takes an
   explicit `frame=`; the seam passes the spec's declared gauge, and fails loudly if a caller
   tries to combine a spec derived in one gauge with CAMB variables requested in another. This is
   the honest-flags principle applied at the one place it matters most.

#### PSALTer

The supervisor's numerical polology package (`psalter.tar.gz`, local; D6 grants explicit
permission to copy) is a JAX package whose public API is

```text
extract(wxf_path) -> Theory
evaluate(theory, c, hyperparameters=...)
sample(theory, tiles, parallelism=..., ns=..., hyperparameters=...)
```

with `Theory`, `Sector`, `EvalResult`, `NestedSamplingResult` as its data types and `couplings`
as the parameter vector. Three concrete alignment items:

1. **Interchange format.** PSALTer ingests **WXF** (Wolfram Exchange Format —
   `psalter/_extract/wxf.py`), emitted by the Wolfram PSALTer. We emit JSON. Our exporter should
   either emit WXF alongside JSON, or adopt a schema PSALTer's `extract` can consume, so that
   `spectrum/` **interoperates instead of translating**.
2. **Vocabulary.** `Theory`, `Sector`, `couplings`, `tiles`. Adopt these names wherever they mean
   the same thing. This is not cosmetic: the cubed-sphere lineage of our legacy
   `RadialAngularPrior` is literally `psalter/_tile/geometry.py`, which is why that module's
   docstring says its face and tile conventions match the supervisor's reference implementation.
3. **Config idiom.** PSALTer uses frozen dataclass settings bundles — `Parallelism`,
   `NSSettings`, `LikelihoodHyperparameters` — over a central `config` module read **at call
   time, not import time**. That is precisely the typed-config architecture this document
   proposes for `tidalcosmo/config/`, so adopting the supervisor's own idiom is free consistency.

Also note `psalter/_extract/likelihood.py`: PSALTer already has a likelihood, which independently
confirms §6 — the spectrum belongs on the sampling path, not off to one side.

#### A structural note

PSALTer is *verb*-shaped (`psalter/<verb>.py` plus a private `_<verb>/` subpackage); SOLikeT is
*component*-shaped. These are compatible with the two-tier split of §2.11: the Cobaya surface
**must** be component-shaped, because Cobaya requires `<ClassName>.yaml` beside the class, while
the engine may be verb-shaped.

### 2.9 Not decided here: the integration target

**H3 owns the choice** between patching CAMB's Fortran, using DISCO-EB, and our own coupled-block
solver chained to unmodified CAMB. It is an open WS0 decision, named in H3's prompt, and this
document does not preempt it.

What this design *does* is make the tree survive either answer: the engine reaches a background
only through `background/protocol.py`, and receives `M(η, k)` from `perturbations/`. Swapping
what sits behind that protocol changes one directory, not the package shape.

H3's scope has widened since its prompt was written. H2 §0.1 showed that **O3 needs a second
engine kind**: a 100 GHz CMB photon makes ~10²⁶ oscillations over a Hubble time, and no
integrator steps through that, so the literature removes the carrier analytically and integrates
only the slowly varying amplitude, with patch averaging along the line of sight. `solver/`
therefore accommodates **two front-ends over a shared core**; H3 decides the internals.

### 2.10 Third-party codes: role, license, and what we take

| Code | Role for us | License | What we may do |
| --- | --- | --- | --- |
| `ohahn/DISCO-EB` | a candidate backend (H3's decision) and the **structural template for our engine-layer module naming** — `background`, `perturbations`, `thermodynamics`, `ode_integrators` | **GPL-3.0** | **Layout and naming inspiration is taken deliberately, and cited** — module organization is not copyrightable expression. Copying its *code* would impose GPL on our package; that stays a conscious decision routed to H3, never something that happens by accident |
| `adammoss/nanoCMB` | the **validation reference** — full LOS projection in ~1400 readable lines at sub-percent agreement with CAMB over `2 ≤ ℓ ≤ 2500`, which is exactly the program's stated O0 tolerance — and a readable implementation of the projection we must own for mechanism 3 | **MIT** | permissive; we may adapt code with attribution. Preferred use remains read-and-reimplement, with the source credited in the docstring |
| `psalter.tar.gz` (Barker) | WS6's basis, and the conventions source of §2.8 | supervisor's own; permission explicit under D6 | copy freely, provenance in docstrings, attribution settled with the supervisor at publication |

CAMB's own repository carries no standard SPDX identifier (`NOASSERTION`), so **reading its
license is a concrete precondition on #498**, before we fork and redistribute for Seam A.

### 2.11 The package tree

Two organizing principles, one per tier.

**The surface is organized by what a user names in their Cobaya YAML.** This is not a free
choice: Cobaya requires a component's defaults to live in `<ClassName>.yaml` **beside the module
defining the class** — `soliket/cosmopower/CosmoPower.yaml` sits next to `cosmopower.py`;
`mflike/TTTEEE.yaml` next to `mflike.py`. The directory is *addressed by the config file*.
Grouping by user-facing component is also what allows one class to carry several YAML presets,
the way `mflike/{TT,TE,EE,TTTEEE}.yaml` all sit beside one `mflike.py` — which is how our
channels should be expressed rather than as separate classes.

**The engine is organized by stage of the calculation** — background, perturbations, solver,
observables. No user ever names these in a YAML; they are called in sequence by the components
above. The naming follows DISCO-EB (§2.10) so that a cosmologist reading our tree already knows
what is in each file.

Every directory name below survives the §1.3 rename unchanged.

```text
tidalcosmo/
  __init__.py        version, and the re-export of the Cobaya components — nothing else. The
                     re-export is what gives the short dotted path (`tidal.SpectatorTheory:`),
                     which is how every reference package does it. The 524 ms lesson is "no
                     *unrelated* heavy imports" (matplotlib via a plotting helper), not "no
                     imports": re-exporting a component Cobaya will import anyway is fine.
  py.typed

  # ---- surface: what a user names in a YAML (Cobaya fixes this layout) ----
  spectator/         SpectatorTheory.py + SpectatorTheory.yaml + per-channel presets.
                     One class, many YAML presets (the mflike pattern).
  background/        protocol.py -- the only thing the engine may ask of a background, and
                     H3's swap point (§2.9). camb_seam.py -- the §2.3 calls.
                     TabulatedBackground.py + .yaml -- Seam A: optional, off by default,
                     and BLOCKED ON the CAMB fork (#498), which is a deliverable in its own
                     right, not a footnote.
  likelihoods/       per-rung likelihoods (birefringence Gaussian prior first)
  presets/           runnable ladder configurations (SOLikeT's presets/ pattern)

  # ---- engine: stage of the calculation (DISCO-EB naming, cited) ----
  config/            frozen dataclass settings bundles over a central config read at call
                     time -- PSALTer's idiom (§2.8). The CLI is an adapter over this; the
                     config object is never an argparse.Namespace.
  spec/              the equation-spec interchange contract, emitted in CAMB/PSALTer
                     conventions (§2.8), carrying the declared gauge as metadata
  derive/            WS2: .wls generation, the wolframscript driver, the ported Wolfram
                     modules
  coefficients/      symbolic coefficient -> numeric callable of (eta, k)
  perturbations/     assembly of the coupled block M(eta, k): our fields PLUS the standard
                     modes they couple to (§2.2)
  solver/            WS3: two front-ends over a shared core (§2.9); internals are H3's
  observables/       transfer.py (mechanism 2) · los.py (mechanism 3; nanoCMB as reference;
                     home of the per-channel source functions of §2.5) · rotation.py
                     (mechanism 1, via camb.symbolic)
  validity/          honest flags: flags.py · spectator.py (Delta-N_eff, |h|, growth-impact
                     monitor) · background_eom.py (PGT residual on the CAMB background)
  spectrum/          WS6 numerical polology -- a gate on the sampling path, not a corner (§6)
  diagnostics/       post-processing over Cobaya/anesthetic output
  cli/               thin: parse -> config -> library call. No physics.
```

**The one structural inversion relative to legacy**, stated plainly because it is the point of
the whole exercise: **the forward model is library code taking a typed config, and the CLI and
the Cobaya component are two thin callers of the same entry point.** Today the CLI *is* the
config object, with the physics living behind `argparse` — which is why four modules outside
`tidal/cli/` have to import private names back out of it. This is also the generalization of what
`tidal/measurement/_run_stages.py` learned the hard way under GH #454: *a stage used by more than
one caller lives in one place.*

---

## 3. The Cobaya surface

Our component is a `Theory` class — Cobaya's pluggable prediction component. The interface is
small: `get_requirements()` declares what it needs from other components, `must_provide()`
declares conditional requirements for the products it is asked for, `calculate()` does the work
and stores results in a state dict, and `get_X()` methods hand products downstream.

```python
class SpectatorTheory(Theory):
    def get_requirements(self):
        return {"CAMBdata": None}
    ...
```

**Wiring template.** H1 §3.1 documented the provider→consumer pattern from `slegner/cobaya`: a
component declaring `{"product": None}` in `must_provide`, pulling it via `self.provider.get_*()`
in `calculate`/`set`, and asserting the downstream state. The physics direction there is opposite
to ours, but the wiring shape is the same and is worth copying.

**Failure behavior.** H1 §3.3 flagged a patch worth stealing: `planck_clik.py` returning `−inf`
when any `C_ℓ` array contains NaN, instead of letting the likelihood crash. Exotic parameter
points *will* reach the likelihood inside a nested-sampling prior box. Our component fails the
same way — **a flagged rejection, not a crash** — and that mechanism is shared with the
`spectrum/` gate of §6 rather than reinvented per site.

**A caveat inherited from H1 §3.2.** Stock Cobaya's PolyChord interface does not give correct
evidences under non-uniform priors without the Ormondroyd patch (it passes the posterior on a
linearly rescaled hypercube rather than the likelihood through the true inverse CDF). Anyone
quoting a `log Z` from a run with non-uniform priors needs to know this.

**Packaging conventions**, taken from the reference packages:

- **Plain dotted-path reference, no entry-point registration.** Users write
  `theory: {tidal.SpectatorTheory: {...}}`, and `python_path:` in the input YAML locates a
  package that is not installed. This confirms the program's existing decision.
- **Conventional YAML keys** the reference packages all carry: `stop_at_error`, `extra_args`,
  `renames`, `python_path`. `renames` is genuinely useful to us for mapping sampled coupling
  names onto spec parameter names.
- **YAML must be declared package data.** LAT_MFLike ships a `MANIFEST.in` for exactly this. Our
  `[tool.setuptools.package-data]` needs `**/*.yaml`, or the defaults simply do not install and
  the failure is confusing.

---

## 4. The fate of `tidal/inference/` under D9

D9 migrates inference into the Cobaya ecosystem, which supplies priors, samplers (PolyChord
included) and output handling. Most of `tidal/inference/` (6,976 lines) is therefore superseded.

**Retired** — Cobaya covers it, and a second implementation is a liability:

| Module | Superseded by |
| --- | --- |
| `_prior.py` (604) | Cobaya's prior specification |
| `_mc.py` (299), `_nested.py` (277) | Cobaya's samplers |
| `_results.py` (512) | Cobaya's output / anesthetic |
| `_likelihood.py` (793) | the `Theory` + likelihood components |
| `_constraints.py` (207) | Cobaya external priors |
| `_visualize.py` (1,755) | GetDist and anesthetic |
| `_prior_stability.py` (204) | per-call validity flags — see below |

`_prior_stability.py` is worth a note, because its retirement is not a like-for-like swap. It
existed because nested-sampling chains contain only accepted samples, so points the pre-flight
guard rejected were invisible on a corner plot; it re-sampled the prior independently to
visualize the unstable region. In the new design, **validity is recorded per call as a flag with
a reason** (§6), so the information is present natively and the workaround is unnecessary.

**Ported:**

- `_importance.py` (1,024) — the KL divergence and Bayesian-model-dimensionality diagnostics
  (`D_KL`, `d_G`, marginal `D_KL` per parameter) become post-processing over anesthetic samples
  in `diagnostics/`. These answer "which parameters does the data actually constrain?", which no
  sampler gives you for free. The Handley (2019) and Handley et al. (2015) references travel with
  the code.

**Ported, deferred, and possibly better taken from PSALTer:**

- `_sphere.py` (398) + `RadialAngularPrior` — the cubed-sphere chart for a coupling-space prior.
  Worth having when coupling-sphere sampling is wanted, as a Cobaya external prior. But per §2.8
  the right move may be to **use `psalter/_tile/geometry.py`'s conventions directly** rather than
  port ours: our module's docstring already says its face indexing, sub-tile labeling and
  orientation match the supervisor's reference implementation, so adopting that implementation is
  more consistent than maintaining a parallel one.

**Informed-by only:** `_atlas.py` (870) — a figure script, not package surface.

---

## 5. Capability inventory: the cutover ledger

**Cutover is capability-gated and incremental.** Legacy stays runnable until each capability *we
intend to keep* has been ported and verified against frozen legacy output — but deletion happens
**per capability, as early as that capability allows**, not in one big-bang at the end. It is all
git-tracked; nothing is lost by removing a slice once its replacement is live. Only the *rename*
(§1.3) needs the whole legacy tree gone.

"Not needed" is a legitimate and cheap outcome. The bar is *the functionality we actually want to
use*, not everything that exists.

### 5.1 The eleven subcommands

| Subcommand | Verdict | Gate against frozen legacy output |
| --- | --- | --- |
| `derive` | **port** (WS2) | **semantic** equivalence — see §5.2 |
| `inspect` | **port** — the accessor layer is how specs are read at all (#401); `tests/data/spec_semantics.txt` is a committed corpus report | reports match, modulo the intended renaming |
| `validate` | **port**, narrowed to spec and stability validation | same verdicts on the example corpus |
| `simulate` | **re-implement, informed-by** — the cosmology forward model is per-`k` in `η`, not a grid PDE; `_simulate.py` (3,126 lines) is thesis-era in shape | not equality: WS3's own gates replace it |
| `measure` | **drop** — see §5.3 | — |
| `list`, `doctor` | **re-implement, informed-by** — small, trivially rewritten against the new config | none needed |
| `plot` | **drop** — see §5.3 | none |
| `sweep`, `sample` | **drop** — superseded by Cobaya under D9 | none |
| `analyze` | **drop** — chain analysis and the t-independence audit serve thesis-era campaign artifacts | none |

### 5.2 The `derive` gate is not a byte diff

The obvious gate would be `tidal inspect OLD --diff NEW` (exit 1 = a real physics change), which
already exists for precisely this comparison. **It is the wrong gate here**, because §2.8 has us
*deliberately* changing variable names, formats and conventions to match CAMB and PSALTer. A byte
diff would report every intended change as a failure.

Byte-equality was the right gate only under the assumption that we were preserving legacy
notation, and we are not — legacy notation has no claim on the new package.

**The gate instead:** the frozen legacy spec remains the oracle **for the physics**, and
equivalence is established by *careful comparison* — same equations, same coefficients, same
signs, under the new naming — recorded as **a written mapping committed alongside the fixtures**.
Slower and more human than a diff, and it is the honest gate for a deliberate reformatting.

The `--diff` tooling remains valuable *within* the new package, for detecting unintended change
between successive re-derivations. It is the cross-package comparison that it cannot do.

### 5.3 `measure` and `plot`, dropped

**`measure` — dropped, with two concepts re-implemented elsewhere.**

- **Energy.** `_energy.py` (1,422 lines) reconstructs a flat-space Hamiltonian from the E-L
  equations. There is no counterpart on FRW — no conserved energy of that kind in an expanding
  universe. What the spectator flag actually needs is `ρ_new/ρ_γ`, a stress-energy density on
  FRW, which is a *different quantity* computed fresh in `validity/spectator.py`. Consequently
  the known `t = 0.0` hardcode bug becomes moot rather than fixed.
- **Conversion.** `P(t) = E_target(t)/E_source(0)` is not an observable we are pursuing. This
  design does not need to settle what O3's observable *is* — that is WS4's business — only that
  the legacy energy-ratio measurement is not it. The `P_max`-versus-`P_final` blocker dissolves
  along with the code.
- **Gauge certificate.** Also dropped. It was an artifact of the halted symbolic-gauge arc
  (#477), and once gauge is handled consistently and checked at the CAMB seam (§2.8) there is
  nothing left for a per-result flag to warn about.

**`plot` — dropped, not re-implemented.** The ecosystem already covers it: **GetDist** for
posterior and triangle plots (Cobaya's standard companion) and **anesthetic** for
nested-sampling output. A third implementation of a solved problem is not worth maintaining.
Anything genuinely bespoke belongs in a one-off figure script, not in the package.

### 5.4 The symbolic surface — a port *and* a substantial extension

This is the real work, and the piece most at risk of being under-planned. The `.wl` modules are
the asset; `tidal/cli/_derive.py` is where TOML→`.wls` generation is fused to `argparse`. Porting
means **separating the generator from the CLI** — the same inversion as the rest of WS1, applied
to the single largest file in the codebase.

**It is not a lift-and-shift.** Legacy never had to handle FRW or genuine time dependence, so
beyond the port WS2 must *add*:

- conformal time as a first-class coordinate;
- `a(η)` and `H(η)` as **unspecified background functions** that survive component decomposition
  and reach the export symbolically;
- gauge projection into a CAMB-named gauge, performed in Wolfram (§2.8);
- emission in CAMB/PSALTer conventions, including the declared gauge as spec metadata;
- **per-channel source functions** (§2.5);
- a **second export form** for H2's eikonal reduction (O3).

The known `ExportJSON.wl:1638` `t`-filtering bug — which drops time-dependent terms from the
Hamiltonian export — must be fixed in the ported copy, but it is one bug inside a much larger
extension, not the extent of the work.

| Piece | Size | Verdict |
| --- | --- | --- |
| `tidal/wolfram/*.wl` — 7 files (`ComponentDecompose` 2,238, `ExportJSON` 1,740, `CommonUtilities` 1,238, `PerturbativeReduction` 396, `Linearize` 158, `GaugeFix` 118, `EulerLagrange` 94) | 5,982 | **port + extend** per above; largely CLI-independent, which is why this is tractable |
| `cli/_derive.py` (6,718) + `_derive_validate.py` (868) + `_wls_helpers.py` (219) | 7,805 | **port** — its own milestone; the generator comes out from behind `argparse` |
| `symbolic/json_loader.py` | 2,476 | **re-implement, informed-by.** Keep the *concepts* — operator vocabulary, kinetic-matrix handling, `implicit_dynamical_sector`, `dependency_closure`/`restrict_spec_dict` — not the schema, which changes under §2.8 |
| `symbolic/spec_query.py` (890), `sign_algebra.py` (1,160) | — | **port** near-verbatim: self-contained, needed early, #401 lineage. The soundness-over-coverage property (a verdict says `unknown` rather than guessing, and names the tactic that decided it) travels with its docstrings |
| `symbolic/latex.py` | 1,520 | **port and update.** Rendering a Lagrangian with its parameters made explicit is exactly how the polology literature presents theories; that is a presentation capability worth keeping, retargeted at the new conventions |
| `_eval_utils.py`, `_separable.py`, expression sandboxing (#406) | — | **port**; the allowlist-by-name security property must travel with its docstring |

**Cost basis:** the **v0.33.9 measured table** in `docs/tex/derivation_performance.tex`. The
per-theory derivation-timing headers in the TOMLs are declared untrustworthy there — treat them
as ceilings, not estimates.

### 5.5 The solver: no port presumed, and no structure inherited

`solver/modal.py` (5,525 lines) is **not ported and is not assumed to be an oracle.** It solves
`expm(M·t)` for constant `M` on a flat periodic grid; the cosmology problem is `M(η)` per `k`, in
two distinct regimes (§2.9). **H3 owns what replaces it.** This document reserves the `solver/`
directory and its seam — consume `M(η, k)` from `perturbations/`, return a transfer function —
and leaves the internals to H3.

The supporting files (`state.py`, `fields.py`, `_kinetic.py`, `kinetic_matrix.py`) hold reusable
*ideas* — the state slot layout, the `M⁻¹` handling — but they accreted over time and **their
structure is explicitly open to redesign.** The new package is not obliged to keep their file
boundaries; H3's architecture decides the shape.

What legacy contributes here is one durable **lesson**, not code: GH #367 and #379 — every
dispatch path must consume the kinetic matrix `M` identically, or cross-path regressions appear
silently. Two such regressions occurred. That becomes a **test contract** for whatever H3
designs, whatever its internal structure.

### 5.6 Remaining library-level rows

| Capability | Verdict |
| --- | --- |
| `solver/coefficients.py` | **re-implement, informed-by** — the multi-level cache concept survives; grid position-dependence gives way to `(η, k)` |
| `measurement/_stability.py` | **re-implement, informed-by** as a validity flag, **never a gate** (GH #454: the probe is a diagnostic and never blocks a run) |
| `_io.py`, `_sweep_results.py`, `_critical_field.py`, `_sensitivity.py`, `_posthoc_audit.py`, `_phase_e_transit.py`, `_resonance.py`, `_energy.py`, `_conversion.py`, `_mixing.py`, `_dispersion.py`, `_effective_mass.py`, operator kernels | **drop** — grid-era and survey-era concerns |
| `measurement/_run_stages.py` | **informed-by** — its rule (*a stage used by more than one caller lives in one place*) is the design principle of §2.11, not code to move |
| `cli/_console.py` `error_with_hint` | **port** — the repo-wide convention that user-facing errors carry actionable hints |

### 5.7 Documentation preservation, made checkable

Every **port** row carries its docstrings, why-comments and GitHub issue references across. That
is stated as an enforceable rule, not an intention:

- each ported module records its origin in the module docstring — `Ported from tidal/<path>
  (GH #NNN)` — and a **port manifest** maps new path → legacy path;
- **a test walks the manifest** and asserts that every `#NNNN` issue reference appearing in a
  legacy source file also appears somewhere in its ported counterpart. Explanations get lost
  silently; a lost issue number is detectable;
- where a port *deliberately* changes behavior — the `ExportJSON` `t` filter, the schema change
  of §2.8 — the ported code says so at the site, naming what it fixes;
- third-party inspiration is cited in the module docstring: DISCO-EB for engine layout, nanoCMB
  for the LOS projection, PSALTer for conventions and for WS6.

---

## 6. `spectrum/` sits on the sampling path

The program document describes WS6 as "independent". **It is not, and this document contradicts
that line deliberately** — flagged in §11 for the orchestrator.

The particle spectrum decides whether a coupling point is a *healthy* theory — free of ghosts and
tachyons. That is exactly the question the sampler needs answered **before** it spends a forward
solve on that point. H2's own dependency graph draws it that way (`WS6 -.gates.-> O2` and
`WS6 -.gates.-> O3`), and PSALTer's `_extract/likelihood.py` confirms it independently: the
reference implementation already treats the spectrum as a likelihood ingredient.

So `spectrum/` is a **fast gate feeding `validity/` and the Cobaya prior surface**: a cheap
in-vacuo screen that rejects sick points before the expensive `M(η)` integration. It is
independent only in its *derivation* (Minkowski-only is correct and sufficient — the spectrum
screens in vacuo, and TorC used the same split) and in its *build order* (it can be written any
time after M0). Its output lands squarely inside the likelihood loop.

Three design consequences:

1. it must be **fast enough to run per sample** — this is a performance requirement, not a
   nice-to-have, and it feeds D5;
2. its verdict is **a flag with a reason**, not a bare boolean — "rejected: ghost in sector X" is
   actionable, "False" is not;
3. its rejection path is **the same flagged-rejection mechanism** as the NaN guard of §3, not a
   separate concept.

---

## 7. Migration order — keyed to capability, not to rungs

This ordering is **deliberately independent of the observable-rung order**, which the orchestrator
has recorded as an open decision (H2 calls it a scientific-priority call, not a technical one).
M0–M4 build the shared spine that every rung needs; which observable is attempted first changes
what is *demonstrated* at M5, not what is *built* before it.

**Deletion is not a milestone.** It happens per capability, in the right-hand column, as soon as
that slice's replacement is live.

| M | Contents | Gate | Retire on completion |
| --- | --- | --- | --- |
| **M0** | **the scaffold** (this document's companion) plus `pyproject` second console script, extras, YAML package-data, the hygiene test, the CI lane | new lane green; no `import tidal` under the new root | the `drop` rows: `sweep`, `sample`, `analyze`, `plot` and their tests, once §5 records the decision |
| **M0.5** | **freeze the oracle** — run legacy across the example corpus and commit the outputs. **Before any porting begins** | fixtures committed and reproducible | — |
| **M1** | **O0** — `config/`, `background/camb_seam.py`, `spectator/` pass-through, `validity/flags.py` skeleton | sub-percent agreement with CAMB/nanoCMB, `2 ≤ ℓ ≤ 2500`; a standard ΛCDM posterior. **Also: benchmark the duplicated-scalar-compute overhead of §2.3** | **`tidal/inference/`** — outdated for a Cobaya interface (D9), goes as soon as the new path stands up |
| **M2** | **O1** — **the CAMB fork first** (re-apply off the `2.0.3` tag, #498; read CAMB's license), then `TabulatedBackground` against it | H1 §R1 gate, plus the free `set_w_a_table` cross-check at `ϖ_r = 0.8` | — |
| **M3** | **WS2** — `derive/` + `spec/` **port and extension** (§5.4), conformal time, CAMB/PSALTer conventions, per-channel source functions | semantic equivalence to the frozen specs (§5.2); de Sitter analytics; FRW-derived EOM → Minkowski EOM as `a → const` | `tidal/wolfram/`, `cli/_derive*.py`, `tidal/symbolic/` |
| **M4** | **WS3** — `solver/` per H3's design, two front-ends | H3's stated per-rung tolerances | `tidal/solver/` |
| **M5** | **WS4** — `observables/`, `validity/spectator.py`; the first rung attempted | validity flags on every run artifact | `tidal/measurement/`, `cli/_simulate.py`, `cli/_measure.py` |
| **M6** | whatever `port` rows remain | every §5 `port` row green | the rest of `tidal/`, `tests/`, the `tidal` console script |
| **M7** | **the rename** — trigger: `tidal/` is gone | the three-commit sequence of §1.3 | — · a single PyPI publish becomes possible after this |
| **M∥** | **WS6** — `spectrum/`; buildable any time after M0, wired into `validity/` before the first rung runs | Lin–Hobson–Lasenby inequalities reproduced exactly | — |

**Why early deletion is safe:** because the oracle is *frozen data* (§8), not live legacy code.
That is precisely what buys the freedom to retire `tidal/inference/` at M1 rather than waiting
for M6.

---

## 8. CI and test strategy for two coexisting packages

- **Two suites.** `tests/` (legacy) stays untouched and green; `tests_cosmo/` grows. Both in
  `testpaths`; coverage over both roots.

- **A hygiene test**, extending the existing `tests/test_repo_hygiene.py` pattern and adding
  `tidalcosmo/` to its `CHECKED_PREFIXES`: no `^(from|import) tidal\b` anywhere under
  `tidalcosmo/` **or under `tests_cosmo/`**. The second half matters most — a test that imports
  or shells out to legacy is exactly how the oracle becomes undeletable infrastructure.

- **Oracles frozen as committed data at M0.5, before any porting.** A script under
  `scripts/oracles/` — **the single place allowed to touch legacy** — runs the legacy CLI across
  the example corpus and commits the results under `tests_cosmo/data/oracles/`: derived spec
  JSONs, `tidal inspect --detail summary` reports, measured scalars, reference `C_ℓ` and transfer
  arrays. New-package tests assert against **those files**, never against a live import or
  subprocess. Deleting legacy then costs nothing and cannot silently break the suite.

- **The comparison is semantic, and the mapping is committed.** Per §5.2, the new naming differs
  from the old by design, so the fixtures ship with a written mapping recording how new
  corresponds to old. That document is part of the gate, not commentary on it.

- **A clean lint and type baseline.** The new tree starts with **no ruff per-file-ignore
  blanket** and strict pyright. The legacy ignore surface — dozens of per-path exemption blocks
  in `pyproject.toml` — is not inherited.

---

## 9. Packaging and documentation layout

One distribution during the transition, with **two console scripts** (`tidal` for legacy,
`tidalcosmo` for the new package — never mixed), and optional extras `camb`, `cobaya`, and a
`cosmo` aggregate. `[tool.setuptools.packages.find]` includes both roots;
`[tool.setuptools.package-data]` must ship `**/*.yaml` or Cobaya's defaults do not install (§3).

Documentation: narrative and program records in `docs/cosmology/`, technical documentation in
`docs/tex/`, runnable configurations in `tidalcosmo/presets/`.

**No PyPI release before M7** (§1.2).

---

## 10. Issues

**Still to open, at the milestone that needs them:** one per ledger row that warrants its own
ticket, plus M0–M7 as milestones. The `drop` rows get issues too — a recorded decision that
`sweep`/`sample`/`analyze`/`measure`/`plot` are superseded is what stops them being re-ported by
reflex two sessions from now.

Beyond the ledger rows, five came out of this study and had no home. All five are now filed:

| # | Finding | Section |
| --- | --- | --- |
| **#513** | Adopt CAMB and PSALTer conventions in the rebuilt symbolic layer — the headline design goal, including gauge-as-metadata and the WXF interchange question. Also records why the `derive` gate cannot be a byte diff | §2.8, §5.2 |
| **#514** | Derive the per-channel line-of-sight source functions — a substantial, previously unowned deliverable that WS2 owes WS4; #506 and #508 take `S` as given | §2.5 |
| **#515** | Benchmark the duplicated CAMB compute when we take a channel over — measurement, not estimation; M1 | §2.3 |
| **#516** | DISCO-EB is GPL-3.0 — decide layout-inspiration versus code-reuse explicitly, before any copying | §2.10 |
| *comment on #498* | CAMB carries no standard SPDX identifier; reading its license is a precondition on forking and redistributing | §2.10 |

H4's own outcome — the namespace decision, the ledger summary, and the two divergences of §11 —
is recorded as a comment on **#490**, the WS1 tracking issue that delegated the namespace
question here.

---

## 11. Divergences from the program document

Two, both deliberate, both for the orchestrator to fold in:

1. **WS6 is not independent** (§6). It gates O2 and O3 on the sampling path, per H2's own
   dependency graph and PSALTer's own structure. "Independent" is true only of its derivation and
   its build order.
2. **The integration target remains H3's open decision** (§2.9) — CAMB-Fortran versus DISCO-EB
   versus our own coupled-block solver. This design is built to survive any of the three; it does
   not preempt the choice.

And one addition worth propagating: **§2.8's convention-adoption goal** should appear in the
program document, because it constrains WS2's output format and it invalidates the obvious
byte-diff gate for `derive` (§5.2).

---

## See also

- `docs/COSMOLOGY_PROGRAM.md` — the operational record: decisions, observable ladder, workstreams
- `docs/cosmology/primer.md` — how a CMB pipeline computes `C_ℓ`; read first if new to the program
- `docs/cosmology/torc_pipeline_audit.md` — H1: the CAMB patch decision (#498) and the Cobaya
  wiring template referenced throughout §2 and §3
- `docs/cosmology/observable_ladder.md` — H2: the two solver kinds, and the rung dependency graph
- `docs/cosmology/spectator_route.md` — the spectator approximation that licenses §2.1
