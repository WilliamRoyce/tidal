# R-C — Cosmological perturbations: the equations we must produce, the methods, and the tools

> **STATUS: DISPATCHED — 2026-09-15.** R-1 merged (#568); D-A is recorded, and the sections that
> depended on it are updated below. Written 2026-09-13 at Wave-1 approval.

| | |
|---|---|
| **Issue** | **#567** (this memo) · #488 (umbrella) · #500 (FRW derivation mode), #501 (background-EOM residual), #504 (eikonal reduction) — the M3 work this informs |
| **Milestone** | M3 (it feeds WS2's design, not a Wave-1 build) |
| **Wave** | 1 — research, second of two memos; **second lane occupant, after R-1 has merged and its gate has been re-run** |
| **Wolfram lane** | **Yes — up to two days.** One `wolframscript` at a time, machine-wide. `ensure_registered.sh` and `verify --require-psalter` exit 0 **before** you start and **after** you finish; both outputs attached. |
| **Depends on** | R-1 merged ✅ 2026-09-15 (D-A recorded). The lane is free: the orchestrator's re-run of R-1's evidence finished 2026-09-15 with no kernel left running |
| **Owned paths** | `docs/cosmology/perturbation_tooling.md` (the memo) · `scripts/research/perturbations/` (install script, probes, reproductions) · **one additive install** into the real userbase: `~/.local/wolfram/userbase/Applications/xAct/xPand/` and nothing else there |
| **NOT owned** | anything under `tidalcosmo/` · `scripts/install-*.sh` · `scripts/psalter/**` · any existing package under `Applications/xAct/` (xTensor, xPerm, xCore, xCoba, xPert, PSALTer — read-only) · the design docs (report contradictions; do not edit) |

### Decision dependencies (the orchestrator updates these before marking READY)

| open decision | sections affected | what changes |
|---|---|---|
| **D-A** (R-1) — **recorded 2026-09-15: one committed Wolfram package run by a fixed driver**, theory passed as WXF data; no generated code; a `wolframclient` session was measured and cannot derive (`DefField` never returns inside one) | Part 4, "Recommendation" | Judge each tool on **whether it exposes callable functions that a committed package can call with data** (xPand and xPert do, per R-1 §2.6) — not on whether its API can be emitted as text, and not on session use. The requirements checklist (part 1) and the tool survey (part 3) do not change. |
| **D-B** (R-1) — **recorded: Option A′** | Part 1 (what the FRW derivation's *input* looks like) | The theory file is shared by both branches; the FRW branch will need a **separate solver-only block** (background fields, field → perturbation map, gauge choice) with **its own fingerprint** (`interfaces_decision.md` §3.8). Part 1 lists what that block must carry; do not design it. |

## Why this exists

M3 (WS2) must derive, for a user-supplied theory, the linear perturbation equations of its
new sector on an FRW background in a CAMB-named gauge, plus everything the ladder downstream
needs. Today the repository has: xPert installed (xAct's general perturbation-theory package,
`~/.local/wolfram/userbase/Applications/xAct/xPert/`, 1.0.6), xTras, xCoba; **xPand is not
installed** (the xAct package for cosmological perturbations on FLRW; the H7 niche analysis,
`COSMOLOGY_PROGRAM.md:302-310`, assessed it as a near-miss *as a product*, never as a
*component*); and no survey of which methods the literature uses or which other packages
exist. Legacy `tidal/wolfram/ComponentDecompose.wl` does a related job on flat and simple
curved backgrounds.

The user's instruction (round 5–6 of planning): this is **not** "does xPand work". It is:
**what equations must the FRW derivation produce, how does the literature derive them, what
packages exist — their documented scope and their actual uses — and how can we interface with
or adapt them**, so that M3 builds on existing tools rather than rewriting what exists, and
only builds what nothing provides.

Two physics facts that constrain the survey:

- **Standard sectors never enter our Lagrangian.** The user's Lagrangian specifies only the new
  sector and its couplings to the known fields it touches (metric/GW, photon). Neutrinos,
  baryons, CDM, dark energy are CAMB's; they reach our equations only as background functions
  and supplied source terms — the `S_std` in `h'' + 2ℋh' + k²h = S_std + S_new`
  (`solver_design.md:330-334`). A tool that insists on owning the whole Boltzmann system is
  answering a different question.
- **Covariant derivatives do not collapse to partials in PGT.** The connection carries
  torsion and is itself dynamical; a tool built for metric-only theories may or may not be
  able to represent a rank-3 torsion field and a torsionful connection. That representability
  is a required probe for every tool.

## Ordered reading list

1. `docs/COSMOLOGY_PROGRAM.md` — the register rows *Admissible theories*, *Two derivations*,
   *Convention carriage*, *Dispatch cadence*; the H7 niche paragraph (`:302-310`).
2. `docs/cosmology/spectator_route.md` §3–§4 (the double expansion; what is reachable).
3. `docs/cosmology/repo_reshape.md` `:228-233` (what the engine may ask of CAMB),
   `:329-348` (source functions are constructed from the solution, derived symbolically),
   `:453-469` (gauge as an explicit named input), `:824-833` (what WS2 must add).
4. `docs/cosmology/solver_design.md` §2 (the O2 equation, `:330-352`), §7.
5. `docs/cosmology/observable_ladder.md` §1–§2 (what O2 needs out of the derivation).
6. `docs/cosmology/interfaces_decision.md` — R-1's memo (exists once this is READY).
7. `tidal/wolfram/ComponentDecompose.wl` and `docs/tex/background_fields.tex` (what legacy
   already does, to compare against — not to build on).
7a. **`docs/cosmology/conventions.md` (canonical since 2026-09-15)** — our FRW equations use
   PSALTer's `(+,−,−,−)`, `ε₀₁₂₃ = +1` and **CAMB's perturbation-variable definitions**
   (§3: `etak = kη_s = −kη/2`, `ḣ_s = 6ḣ`, φ is the Weyl potential). **Ma & Bertschinger is
   `(−,+,+,+)`**: every sign difference in your reproductions must be attributed to convention or
   physics (§6), never left as "matches up to sign".
7b. `docs/cosmology/interfaces_decision.md` §2.6 (what D-A implies for tools), §3.8 (the
   solver-only block), §5.1 (R-1's pinned read of **SymBoltz.jl** at `3d1f20a3` — handed to you:
   its equations-in approach belongs in parts 2–3).
8. `literature/astro-ph_9506072/` (Ma & Bertschinger — your scalar reproduction target).

## What you produce

`docs/cosmology/perturbation_tooling.md`, in four parts.

### Part 1 — Requirements: the equation set, as a checklist

From what CAMB consumes and what the ladder needs (reading items 3–5): the linear perturbation
system for a new sector on FRW in a CAMB-named gauge (covariant, Newtonian, synchronous);
its coupling to the standard sectors it touches (tensor, photon), with those sectors entering
only as supplied source terms; the per-channel source functions `S(k,η)` (`repo_reshape.md:
329-348`); the background-EOM residual (#501 — the order-1 coefficient, which *is* the tadpole);
the eikonal-reduced form (#504) as a second export; conformal time as the coordinate; `a(η)`,
`ℋ(η)` as unspecified background functions. **One line per requirement with a source.** This
checklist is what every tool in part 3 is measured against.

### Part 2 — Methods in the literature

How cosmological perturbation equations are actually derived from an action, in practice:
the xPand/xPert line (Pitrou–Roy–Umeh and its successors); the EFT-of-dark-energy /
hi_class / EFTCAMB route; SymBoltz.jl's equations-in approach; CLASS/CAMB's hand-derived
hierarchies; PyTransport/CppTransport's symbolic core; and **any work deriving perturbations
for torsion or non-Riemannian sectors on FRW** — the papers `spectator_route.md` and
`birefringence_notes.md` cite are the starting list, not the end. For each: what it takes as
input, what it outputs, what it assumes about the connection, and what it would leave us to do.

### Part 3 — Tools

xPand, xPert, xTras, xCoba (installed), other xAct-based packages, SymBoltz.jl, the
hi_class / EFTCAMB code generators, Cadabra, and anything part 2 surfaces. For **each**:
documented scope and functionality (its manual and paper); **where it has been used** (citing
papers, dependent codebases); license; last release; compatibility with the certified bundle
(xAct 1.3.0 / xPert 1.0.6 / Wolfram 14.3.0); what it does and does not do against part 1's
checklist; and **whether a rank-3 torsion field and a torsionful connection are
representable**.

**xPand/xPert hands-on — mandatory:**

1. Install xPand **into the real userbase, additively, under `Applications/xAct/xPand/` only**
   (the certified fingerprint checks the four existing packages' versions, so an added package
   cannot change it; #563 already records the tree as non-pristine). Attach
   `verify --require-psalter` **before and after**.
2. Reproduce the FLRW tensor equation `h'' + 2ℋh' + k²h = 0` from the Einstein–Hilbert action
   at second order.
3. Reproduce the conformal-Newtonian scalar Einstein equations against Ma & Bertschinger.
4. Declare a torsion tensor as a perturbation on FLRW with `T̄ = 0`, apply xPand's SVT
   decomposition, and report what happens — representable, representable with work, or not.
5. Compare with what legacy `ComponentDecompose.wl` already does for the same field content.

Scripts committed under `scripts/research/perturbations/`; symbolic equality, not visual.

### Part 4 — Recommendation

Build-on / borrow-pieces / build-own, **per requirement in part 1**, with the cost of each,
and the interface each choice implies **under D-A as recorded**: our FRW derivation will be a
function in a committed Wolfram package, run by a fixed driver on WXF data — so say, per tool,
which of its functions that package would call, with what data, and what it would return.

## Success criteria — verified from artifacts

1. Part 1 is a checklist with a source per line, covering every item named above.
2. Part 2 cites each method's paper **and** code.
3. Part 3 covers every tool named in part 2 and every xAct-family package, each with scope,
   uses, license, release, compatibility, checklist coverage and the torsion probe.
4. The two xPand reproductions match to symbolic equality; the torsion probe has a definite
   answer with its script committed.
5. `verify --require-psalter` exit 0 after the session, attached; only `Applications/xAct/xPand/`
   added.
6. Part 4 recommends per requirement with cost — not "further work is needed".
7. Draft PR from the first commit; `ruff`, `cspell` over changed files; lane time ≤ 2 days.

## Scope fence

Nothing under `tidalcosmo/`. No change to `scripts/install-*.sh`. No FRW derivation of *our*
theories — the reproductions are standard GR results, chosen because their answers are known.
No modification of any existing xAct package. No PSALTer runs beyond the verify script.

## Working rules

- Worktree off `feat/cosmology-program`: `git worktree add /tmp/tidal-rc -b cosmo/rc-perturbation-tooling feat/cosmology-program`. **Never merge.**
- **Never version-bump, tag, or edit the changelog.**
- **Draft PR into `feat/cosmology-program` at your first commit**; CI is the gate; `CI <run-id>: <conclusion>`.
- Stay inside your owned paths; the one userbase write is `Applications/xAct/xPand/`.
- **Assertions get verified before they land.** "Tool X can represent torsion" is a script, not a sentence.
- **An argument for why something cannot be checked is not evidence about how the check comes out.**
- **A caveat that is a scope limit hides better than a doubt** — it must reach the memo's tables.
- **Read each tool's own known-issues** before concluding it cannot do something.
- One `wolframscript` at a time, machine-wide; the lane guard blocks a second.
- No environment-specific absolute paths in anything committed.
- Conventional commits; American English; no attribution trailers.

## If you find the design wrong

Amend a design-document error at the instruction site and report it. An architectural
contradiction: stop and report.

## Report back

Branch · PR · `CI <run-id>: <conclusion>` · each criterion with its output (the two
reproduction digests, the torsion-probe verdict, both verify outputs) · the part-4 table ·
amendments made · discoveries to route · what you could **not** test here and why.
