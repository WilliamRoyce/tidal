# M1b — The Cobaya `Theory` class, the veto wiring, and the ΛCDM identity gate

> **STATUS: HELD (outline) — Wave 2. Depends on D-B, I-532 and I-S1A-core.** Written
> 2026-09-13 at Wave-1 approval so nothing decided in Wave 1 is forgotten; the orchestrator
> updates the marked sections at each decision. **Not dispatched in Wave 1.** It becomes a
> full prompt, then READY, at Wave-2 planning.

| | |
|---|---|
| **Issue** | **#494** (M1b) · #515 (the benchmark) · #490 (the config layer this replaces) · #454 (the flag lesson) · #488 |
| **Milestone** | M1b |
| **Wave** | 2 |
| **Wolfram lane** | **NO.** Never start a kernel — a Cobaya run must never be able to (R-1's never-derive mechanism is what this class enforces). The complete list of things that would start one is in `I-532.md`'s header. |
| **Depends on** | D-B (R-1: the YAML block and the artifact handoff); I-532 merged (protocol, seam, flags, the veto-placement memo); I-S1A-core merged (the provisional config layer to replace) |
| **Owned paths** | `tidalcosmo/spectator/theory.py` (the `SpectatorTheory(Theory)` class) and `SpectatorTheory.yaml` · `tidalcosmo/__init__.py` (the re-export) · `tidalcosmo/config/` (**replacing** I-S1A-core's provisional layer, #490) · `tests_cosmo/test_theory*.py`, `tests_cosmo/test_config*.py` · `pyproject.toml` — the `cobaya` core promotion only (`:66-71`) · `benchmarks/` for #515 |
| **NOT owned** | `tidalcosmo/background/` (I-532's seam — call it, do not change it) · `tidalcosmo/validity/flags.py` beyond adding sites · `tidalcosmo/{derive,spectrum}/` · `docs/cosmology/veto_placement.md` (a decision record; a departure from it is a report, not an edit) |

### Decision dependencies (the orchestrator updates these before the Wave-2 rewrite)

| open decision | sections affected | what changes |
|---|---|---|
| **D-B** (R-1) | the YAML block's contents; `SpectatorTheory.yaml` defaults; how the artifact is located and its staleness detected | Path in YAML → `initialize()` opens it and checks the content hash against the theory block; embedded → the block carries the artifact and the hash check is internal; discovered → a documented search path with the same refusal. *(To be filled: R-1's recommendation and its worked YAML.)* |
| **I-532** | the seam calls; the flag type; **the veto wiring** | The memo's recommendation (external prior before any theory, or `calculate` returning `False`) is implemented as written; the policy object carries the group's severities. *(To be filled: the memo's verdict and the flag import path.)* |
| **I-S1A-core** | deliverable 3 — the config layer | Whatever shape the provisional layer took, this prompt **replaces** it with the D-B-format real one; every `EXPIRES-WITH: #490` guard is retired here. *(To be filled: the provisional layer's public surface.)* |

## Outline — what this prompt will ask for

1. **`SpectatorTheory(Theory)`** with `SpectatorTheory.yaml` defaults, `get_requirements`
   declaring only what the audited seam provides, `calculate` returning CAMB's own arrays in
   pass-through mode, and the artifact handoff per D-B — refusing at `initialize()` on a
   missing artifact or a content-hash mismatch, **never deriving**.
2. **The `tidalcosmo/__init__.py` re-export** — the cobaya-core trigger — with `cobaya`
   promoted from extra to core and the trigger comment at `pyproject.toml:66-71` updated.
3. **The veto wiring** exactly as `veto_placement.md` recommends, using I-532's flag type
   and policy object; one evaluator, never two code paths (#454).
4. **The ΛCDM-posterior identity gate** (`repo_reshape.md` §2.5): with the new sector off, a
   short Cobaya run's `logpost` at a fixed set of points equals a CAMB-only model's to
   machine precision — identity, not agreement.
5. **The #515 benchmark**: per-sample overhead of the pass-through Theory over bare CAMB,
   recorded with the machine fingerprint.
6. **Replacement of the provisional config layer** (#490): the D-B-format loader becomes the
   one owner of theory input; I-S1A-core's `EXPIRES-WITH: #490` guards are retired; the
   validator and convention mapping move with it unchanged.

## Success criteria — to be finalized at the Wave-2 rewrite

- A Cobaya YAML with our block runs `get_model` and evaluates `logpost` with no artifact
  work beyond a read; the refusal cases from R-1's prototype are tests here.
- Identity gate: `logpost` equal to a CAMB-only model at every test point, machine precision.
- The veto test from `veto_placement.md` reproduced in the real class: the ordering fact
  holds.
- `git grep 'EXPIRES-WITH: #490'` returns nothing after merge.
- #515 numbers recorded; `import tidalcosmo` ≤ 600 ms even with `cobaya` core.
- `CI <run-id>: success`; pyright; boundary test.

## Scope fence — provisional

No physics beyond pass-through. No solver, no stepper (M4). No likelihood of our own. No
change to the seam. No change to the veto decision — implement it.

## Working rules

The standard set from `docs/COSMOLOGY_PROGRAM.md` §"Prompt template" — worktree off
`feat/cosmology-program` (`git worktree add /tmp/tidal-m1b -b cosmo/m1b-cobaya-theory
feat/cosmology-program`); draft PR at the first commit; never merge, bump, tag or edit the
changelog; assertions verified before they land; a verification never writes into what it
verifies; guards carry `EXPIRES-WITH:`; no kernel, ever; no environment-specific absolute
paths; conventional commits; American English; no attribution trailers.

## Report back

Branch · PR · `CI <run-id>: <conclusion>` · each criterion with its output · the identity
gate's point set and values · the #515 table · what of the provisional layer survived and
what was retired · amendments made · discoveries to route.
