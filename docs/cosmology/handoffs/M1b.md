# M1b — The Cobaya `Theory` class, the veto wiring, and the ΛCDM identity gate

> **STATUS: HELD (outline) — Wave 2. Depends on I-532 and I-S1A-core** (D-B and a permanent
> theory loader recorded 2026-09-15 and written in below). Written
> 2026-09-13 at Wave-1 approval so nothing decided in Wave 1 is forgotten; the orchestrator
> updates the marked sections at each decision. **Not dispatched in Wave 1.** It becomes a
> full prompt, then READY, at Wave-2 planning.

| | |
|---|---|
| **Issue** | **#494** (M1b) · #515 (the benchmark) · #454 (the flag lesson) · #573 (the per-sample check) · #566 (R-1) · #488 |
| **Milestone** | M1b |
| **Wave** | 2 |
| **Wolfram lane** | **NO.** Never start a kernel — a Cobaya run must never be able to (R-1's never-derive mechanism is what this class enforces). The complete list of things that would start one is in `I-532.md`'s header. |
| **Depends on** | D-B ✅ (R-1, 2026-09-15); I-532 merged (protocol, seam, flags, the veto-placement memo); I-S1A-core merged (the permanent loader and the spectra store this class builds on) |
| **Owned paths** | `tidalcosmo/spectator/theory.py` (the `SpectatorTheory(Theory)` class) and `SpectatorTheory.yaml` · `tidalcosmo/__init__.py` (the re-export) · `tests_cosmo/test_theory*.py` · `pyproject.toml` — the `cobaya` and `wolframclient` core promotions only (`:66-71`) · `benchmarks/` for #515 |
| **NOT owned** | `tidalcosmo/background/` (I-532's seam — call it, do not change it) · `tidalcosmo/validity/flags.py` beyond adding sites · `tidalcosmo/{config,derive,spectrum}/` (I-S1A-core's permanent loader and store — build on them; a needed change is reported) · `docs/cosmology/veto_placement.md` (a decision record; a departure from it is a report, not an edit) |

### Decision dependencies (the orchestrator updates these before the Wave-2 rewrite)

| open decision | sections affected | what changes |
|---|---|---|
| **D-B** (R-1) — **recorded 2026-09-15: Option A′** | the YAML block; `SpectatorTheory.yaml`; locating and gating the spectrum | The block is `model:` (the theory, usually `!defaults ../theories/<Name>`) plus `spectra_path:`; priors only in top-level `params:`. **Couplings registered in `get_modified_defaults`** (the LAT_MFLike `Foreground` pattern, `mflike/foreground.py:215-225`); `initialize()` loads the theory with I-S1A-core's loader, finds the entry with its store and **raises `ComponentNotInstalledError`** (`camb.py:290-296` idiom) when missing, a hinted `LoggedError` when corrupt, stale or in other conventions; **`get_version()` returns `spectrum-<fingerprint>`**, so a resume against a different spectrum is refused by Cobaya itself; no `is_installed`/`install` hooks (`install.py:232-243, 374` never pass the block). Two lessons from R-1's run: `initialize_with_params` must **not** raise on *missing* couplings (it would preempt Cobaya's clearer requirement error — check only assigned non-couplings), and the unit likelihood `one` absorbs misspelled parameters, so typo tests need a real likelihood (§3.4). Worked YAML: `scripts/research/interfaces/runs/vector_one.yaml`; the prototype class: `proto/spectator_gate.py`. |
| **I-532** | the seam calls; the flag type; **the veto wiring** | The memo's recommendation (external prior before any theory, or `calculate` returning `False`) is implemented as written; the policy object carries the group's severities. *(To be filled: the memo's verdict and the flag import path.)* |
| **I-S1A-core** — **the loader is permanent (user, 2026-09-15)** | outline item 1 | This prompt **builds on** I-S1A-core's loader and store — it no longer replaces a provisional layer. *(To be filled: their import paths from I-S1A-core's report.)* |

## Outline — what this prompt will ask for

1. **`SpectatorTheory(Theory)`** with `SpectatorTheory.yaml` defaults, `get_requirements`
   declaring only what the audited seam provides, `calculate` returning CAMB's own arrays in
   pass-through mode, and the artifact handoff per D-B — refusing at `initialize()` on a
   missing artifact or a content-hash mismatch, **never deriving**.
2. **The `tidalcosmo/__init__.py` re-export** — the cobaya-core trigger — with `cobaya`
   promoted from extra to core and the trigger comment at `pyproject.toml:66-71` updated;
   **`wolframclient` promoted to core in the same change**, pinned `<2`, with a test that the
   sampling path never loads `wolframclient.evaluation`, `zmq` or `pkg_resources` (register row
   *`wolframclient` is core at M1b*, 2026-09-15).
3. **The veto wiring** exactly as `veto_placement.md` recommends, using I-532's flag type
   and policy object; one evaluator, never two code paths (#454).
4. **The ΛCDM-posterior identity gate** (`repo_reshape.md` §2.5): with the new sector off, a
   short Cobaya run's `logpost` at a fixed set of points equals a CAMB-only model's to
   machine precision — identity, not agreement.
5. **The #515 benchmark**: per-sample overhead of the pass-through Theory over bare CAMB,
   recorded with the machine fingerprint.
6. *(Removed 2026-09-15: the theory loader is permanent, so there is no provisional layer to
   replace. The Cobaya-facing modules keep R-1's static no-kernel test green.)*

## Success criteria — to be finalized at the Wave-2 rewrite

- A Cobaya YAML with our block runs `get_model` and evaluates `logpost` with no artifact
  work beyond a read; the refusal cases from R-1's prototype are tests here.
- Identity gate: `logpost` equal to a CAMB-only model at every test point, machine precision.
- The veto test from `veto_placement.md` reproduced in the real class: the ordering fact
  holds.
- The static no-kernel test covers `tidalcosmo/spectator/theory.py`; the twelve R-1 gate cases are tests of the real class.
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
