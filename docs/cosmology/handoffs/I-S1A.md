# I-S1A-core — The permanent theory loader, the spectra store, the driver launcher, and the Stage-2 reader

> **STATUS: HELD — depends on I-532 having merged (its flag type).** Its R-1 decisions — D-A,
> D-B, the conventions rule and a **permanent** theory loader — were recorded on 2026-09-15 and
> are written into the text below. Written 2026-09-13 at Wave-1 approval; rewritten 2026-09-15.
> The orchestrator fills in the flag type's import path and marks this READY after I-532 merges.
> **Do not dispatch while HELD.**

| | |
|---|---|
| **Issue** | **#527** (I-S1A) · #495 (WS6 umbrella) · #522 (linearity is ours), #523 (the harvest surface) · #566 (R-1, the decisions this builds) · #488 |
| **Milestone** | M-parallel |
| **Wave** | 1 — build, second; **after I-532 has merged** |
| **Wolfram lane** | **NO.** Never start a kernel — every test runs against committed fixtures and subprocess mocks. The complete list of things that would start one: `wolframscript`, `tidal derive`, every `examples/*/run.sh`, `scripts/psalter/{run_tier1_gate.sh,ensure_registered.sh,register_resources.wl,repro_543.wl,vector_smoke.wls,probe_*.wls}`, `scripts/verify-wolfram-setup.sh`, `scripts/run_wolfram_tests.sh`, `scripts/full_test.sh`, `scripts/research/interfaces/proto/{derive_proto,session_probe,session_probe_steps}.py`, the `/derive` skill (#555). |
| **Depends on** | R-1 ✅ (D-A, D-B, conventions, loader permanence — 2026-09-15); I-532 merged and its flag type's import path published |
| **Owned paths** | `tidalcosmo/config/` · `tidalcosmo/derive/` — **Python only, never `tidalcosmo/derive/wolfram/`** · `tidalcosmo/spectrum/` · `tests_cosmo/test_{config,derive,spectrum}*.py` · `tests_cosmo/fixtures/psalter/` (read-only) · `tests_cosmo/fixtures/theories/` (new) |
| **NOT owned** | `pyproject.toml` (I-532 owns it this wave; the orchestrator settles `wolframclient`'s placement before you start) · `tidalcosmo/{background,spectator,validity}/` · `tidalcosmo/cli/` (the `tidalcosmo derive` command is I-S1B's, which can run it end to end) · any `.wl` or `.wls` · the design docs |

### Decisions this prompt builds (recorded 2026-09-15)

| decision | what it means here |
|---|---|
| **D-A** — one committed Wolfram package run by a **fixed** `driver.wls`, theory passed as **WXF data** | Deliverable 3 is a *launcher* for that fixed driver. Python never writes Wolfram code. Symbol-level checks on operator text (names, whitelist, a coupling inside an operator) happen **in the package, on a held parse** (I-S1B) — not here. |
| **D-B** — Option A′: a theory YAML (`fields`, `lagrangian` coupling → operator, `derivation`) included into Cobaya run files with `!defaults`; a content-addressed spectra store; refusal at construction | Deliverable 1 loads it; deliverable 2 is the store the Cobaya Theory (M1b) will read. |
| **Conventions** — `docs/cosmology/conventions.md` is canonical: `(+,−,−,−)`, `ε₀₁₂₃ = +1` for both branches | The theory file is written in PSALTer's convention; nothing is converted. Every manifest records `conventions`; every read asserts them. |
| **Theory loader is permanent and schema-versioned** | Build it to last. **No `EXPIRES-WITH: #490`.** M1b builds on it. |
| **I-532's flag type** | Deliverable 4's verdict uses it. *(Import path to be filled at READY.)* |

## Why this exists

The two-stage spectrum architecture (`spectrum_design.md` §1; `stage1_engineering_plan.md`)
derives once with Wolfram/PSALTer and evaluates the result per sample with **no Wolfram at
sampling time**. R-1 (#566) settled the path a theory takes, and prototyped all of it under
`scripts/research/interfaces/` — measured, not argued: the prototype package reproduced the
reference wave operator `SameQ`; eleven malformed and injected inputs were refused before PSALTer
ran; twelve Cobaya gate cases refused and passed as designed. **Your job is to turn the Python half
of that prototype into production code**, fixing the one defect R-1 found in it, and to add the
Stage-2 reader the prototype did not need.

Four pieces, none of which starts a kernel:

1. **The theory loader.** Everything a user writes enters here. It must be exactly the format
   R-1 recommended and the user adopted — close-but-different from Cobaya's own conventions would
   be worse than none.
2. **The spectra store.** A derived spectrum costs minutes to hours and a Wolfram license; it is
   *data*, not cache. The store is what makes "`cobaya-run` can never start a kernel" true: the
   Theory only ever reads it.
3. **The driver launcher.** Carries three measured requirements (#543,
   `stage1_engineering_plan.md:534-563`): one kernel machine-wide; `QT_QPA_PLATFORM=offscreen` on
   every launch; and **never judge a run by its exit status** — 14.3 segfaults on shutdown after
   correct output, an uncaught `Throw` exits 0 (#561), and PSALTer has a bare `Quit[]` (#572).
4. **The WXF reader and Stage-2 contract dataclasses**, tested against the committed answer keys;
   `A[0]` represented in both forms as a named unknown (I-S1B answers it).

## Ordered reading list

1. `docs/COSMOLOGY_PROGRAM.md` — register rows **D-A**, **D-B**, **Conventions**, **Theory loader
   is permanent** (all 2026-09-15), *I-S1A split*, *`A[0]` owner*, *Dispatch cadence*.
2. **`docs/cosmology/interfaces_decision.md`** — §2.4 (what the package refuses, so you do not
   duplicate it), §2.6 and §3.7 (consequences for this prompt), §3.1 (the schema as run), §3.4 (the
   gate), §3.6 (store layout and lookup order), **§3.8 (schema stability — the gaps you must name,
   not fill, and the fingerprint defect you must fix)**.
3. **`scripts/research/interfaces/proto/theory_store.py`, `derive_proto.py`, `spectator_gate.py`,
   `run_gate_cases.py`** and `README.md` — the prototype you are productionizing. It is research
   code: port its design, not its shortcuts (hard-coded `14.3` engine path, `print` sentinels).
4. **`docs/cosmology/conventions.md`** §1, §4, §5.
5. `docs/cosmology/stage1_engineering_plan.md` §4.2 (`:481-500`), §4.4 with its 2026-09-15
   amendment, §5 (`:565-640`), the gates at `:753-764`.
6. `docs/cosmology/spectrum_design.md` §6.1, `:346-366`, `:610-613`.
7. `tests_cosmo/fixtures/psalter/PROVENANCE.md` (which artifact carries the published values).
8. I-532's report (the flag type's import path).
9. Pinned PSALTer sources: `Sources/DefField.m:11-24, 32, 58-88` (the 14 symmetry classes; no
   symmetry means `GenSet[]`), `ParticleSpectrum.m:74-81` (the private globals — enumerate, never
   hard-code a count).
10. Cobaya 3.6.2: `cobaya/yaml.py` (`yaml_load_file`, `!defaults`, duplicate keys rejected at
    `:104-113`).

## Deliverables

1. **`tidalcosmo/config/` — the permanent Option A′ theory loader.**
   - **Schema** (`schema` constant, versioned): `fields: {<name>: {indices, symmetry?, print_as?,
     print_source_as?}}`, `lagrangian: {<coupling>: <operator string>}`, `derivation:
     {max_laurent_depth ∈ {1,2,3}}`. `symmetry` is optional and, when present, is one of
     `{symmetric: [...]}` / `{antisymmetric: [...]}` over the field's own indices, validated so
     that (rank, type, slots) is one of **PSALTer's 14 registered classes** (`DefField.m:11-24`) —
     e.g. `{symmetric: [-a,-b]}` (`FieldKinematics.m:6`), torsion `{antisymmetric: [-b,-c]}`
     (`:18`); absent means PSALTer's `GenSet[]`.
   - **Loading** with Cobaya's own `yaml_load_file`, so `!defaults` expands exactly as `cobaya-run`
     expands it; from a theory file directly, or from a run file's single `theory: <component>:
     model:` block (the inline form uses the same schema and fingerprint).
   - **Validation is structural only**: required keys and types, unknown keys refused, the
     `symmetry` enum, `lagrangian` a non-empty mapping of string → string (so the legacy
     monolithic Lagrangian string — embed the string from
     `examples/torsion_gertsenshtein/theory.toml`, not the path — is refused with a hint that each
     coupling maps to its own operator). **Operator text is not parsed in Python**: names, the
     symbol whitelist and couplings-inside-operators are the package's checks on a held parse
     (`interfaces_decision.md` §2.4). One definition per check — a copied check is where two
     stages silently diverge.
   - **Hints** on every refusal, in the `error_with_hint` style (what is wrong, the legal values,
     an example).
   - **Fingerprint**: sha256 over a canonical form that **drops defaulted and absent keys** —
     fixing the prototype defect at `theory_store.py:75`, which merged `derivation` defaults into
     the hashed form so any new default would re-derive every stored spectrum — plus `schema`.
     Tests: (a) an explicit `max_laurent_depth: 1` and an omitted one fingerprint identically;
     (b) adding an absent `symmetry` leaves a Vector fingerprint unchanged; (c) any change to an
     operator changes it; (d) priors, likelihoods and sampler settings in a run file never enter
     it; (e) the fingerprint of `runs/vector_one.yaml` equals that of `theories/VectorTheory.yaml`
     (R-1 measured `57a91b35…c81d359` with the prototype's recipe; yours differs by design —
     record the new value and say why).
   - **Named gaps — represented nowhere, refused as unknown keys, listed in the module
     docstring with owners**: expansion geometry and derived fields (`F = dA`) → I-S1B;
     gauge-symmetry declarations → I-S1B (provisional); the solver-only background / linearization
     / gauge block → M3, with its own fingerprint.

2. **`tidalcosmo/spectrum/store.py` — the content-addressed spectra store (kernel-free, read and
   write side).** Layout `<store>/v<schema>/<fp[:2]>/<fp>/{spectrum.wxf, manifest.json}`; manifest
   = fingerprint, canonical theory, sorted couplings, `conventions`, spectrum sha256, generator
   versions, timings, host; lookup order `spectra_path` option → `TIDALCOSMO_SPECTRA_PATH`
   (`os.pathsep` list; first writable entry is the write target) → `$XDG_DATA_HOME/tidalcosmo/spectra`
   (default `~/.local/share`, **data not cache**, `interfaces_decision.md` §3.6; use the XDG
   variables directly — no `platformdirs`); atomic write (temp directory + rename on one
   filesystem); **verify on read**: fingerprint, sha256, `conventions` equal to `conventions.md`'s,
   schema, coupling roster — each mismatch a distinct, hinted error. Tests reproduce R-1's gate
   refusals at the store level (missing, stale, truncated, wrong conventions, second store on the
   path, moved store). **A static test** that walks the AST of `tidalcosmo/config/` and
   `tidalcosmo/spectrum/` and finds no `subprocess`, `multiprocessing`, `pexpect`,
   `wolframclient.evaluation` or `os.system`-style launcher — **shown failing** on
   `tidalcosmo/derive/` (the prototype's check, `run_gate_cases.py`, is the model).

3. **`tidalcosmo/derive/launcher.py` — runs the fixed driver.** Serializes the loaded theory plus
   fingerprint to `theory.wxf` with `wolframclient`'s **kernel-free** serializer; launches
   `wolframscript -file <driver.wls> <theory.wxf> <workdir> [--validate-only]` as one subprocess
   (the argv, sentinel lines `STAGE1 PASSED` / `STAGE1_ERROR=<tag>: <message>` and exit codes
   0/1/2 are the prototype's contract, `scripts/research/interfaces/wolfram/driver.wls`; I-S1B
   writes the production driver to it); the engine test is a **file**
   (`$HOME/.local/wolfram/engine/<series>/Executables/WolframKernel`, series discovered, never a
   hard-coded `14.3`; `command -v wolframscript` is satisfied by the image's cloud-only client,
   #565); `QT_QPA_PLATFORM=offscreen` set here; a **`pgrep` lane self-guard** before launch (the
   lane hook cannot see a Python-launched kernel); refusal unless `verify-wolfram-setup.sh
   --require-psalter` exited 0 (cached per process); timeout **off by default** (derivations can
   take hours); log capture; **success only when the sentinel is present, `spectrum.wxf` exists and
   its embedded fingerprint matches** — then the store's atomic write. Unit-tested with subprocess
   mocks covering: file-based engine guard, offscreen in the environment, lane busy, exit 0 with
   no artifact → failure, exit 143 with a good artifact and sentinel → success, a wrong embedded
   fingerprint → failure and nothing stored, verify refusal. **No kernel is ever started by this
   prompt.**

4. **`tidalcosmo/spectrum/` — the WXF reader and the Stage-2 contract dataclasses.** Reader with
   defensive key normalization (bare symbol = absent, `PROVENANCE.md:19-23`); dataclasses per
   `spectrum_design.md` §6.1 / `stage1_engineering_plan.md:605-620`: schema version, theory,
   PSALTer commit, `conventions` (asserted — a mismatch fails), ordered couplings, per-sector
   explicit `J^P` labels, the coupling-linear coefficient tensor `A[(n+1) × 3 × dim × dim]`, gauge
   rank, the massive/massless partition, source constraints, gauge generators, massless
   polarizations, massive spectrum, unitarity status as the **four-outcome** enum (closed form /
   provably impossible / absent / PSALTer halted), `z_degree`, `k4_guard`. `A[0]` representable in
   both forms with a named unknown **in a docstring** and a test exercising both. A test that reads
   `ParticleSpectrum.m:74-81` from the installed tree **when present** and compares the
   private-global set to the reader's, skipping cleanly otherwise (`EXPIRES-WITH:` the PSALTer
   pin). The verdict uses I-532's flag type *(import path filled at READY)*.

## Success criteria — verified from artifacts

1. `A23` block dims `(2,4,2)` with `2·1 + 4·3 + 2·5 = 24`; `Vector` blocks equal the published
   expressions from the committed fixture (state **which artifact** carried them,
   `PROVENANCE.md:87-93`); fixture digests asserted.
2. `scripts/research/interfaces/theories/VectorTheory.yaml` and `runs/vector_one.yaml` load, with
   equal fingerprints. Of the eleven files under `theories/bad/`, `duplicate_coupling.yaml` is
   refused at load (by Cobaya's own loader); **the other ten load** — underscore name, numeric-string
   key, undeclared symbol, coupling inside an operator, two expressions, zero operator, name
   collision, both injections, inexact number are symbol-level — and a test says they are the
   package's to refuse, so nobody later copies those checks into Python. Add your own structural
   fixtures (missing `fields`, a string `lagrangian`, an unknown key, an illegal `symmetry`).
3. The five fingerprint tests of deliverable 1; ≥ 1 test per loader and store message.
4. A convention mismatch in a manifest or a WXF fixture is a failing test.
5. The static no-kernel test passes on `config/` and `spectrum/` and is shown failing on `derive/`.
6. The launcher tests of deliverable 3.
7. `git grep 'EXPIRES-WITH: #490'` returns nothing in your paths.
8. `import tidalcosmo` ≤ 600 ms; pyright clean; `tests_cosmo/test_package_boundary.py` green.
9. Draft PR from the first commit; `CI <run-id>: success`; `ruff`, `ruff format --check`,
   `cspell`.

## Scope fence

No `.wl` or `.wls` of any kind; no production `driver.wls` (I-S1B). No `tidalcosmo derive` CLI
wiring (I-S1B). No Cobaya `Theory` class, no `get_modified_defaults` (M1b). No Stage-2
evaluator. No kernel. No `pyproject.toml` edits. **Do not fill the named schema gaps** —
represent none of them; an unknown key is refused. Do not parse operator strings in Python.

## Working rules

- Worktree off `feat/cosmology-program`: `git worktree add /tmp/tidal-is1a -b cosmo/is1a-stage1-python feat/cosmology-program`. **Never merge.**
- **Never version-bump, tag, or edit the changelog.**
- **Draft PR into `feat/cosmology-program` at your first commit**; CI is the gate; `CI <run-id>: <conclusion>`.
- Stay inside your owned paths.
- **Assertions get verified before they land; hypotheses do not have to be.**
- **A verification command must not write into what it verifies.**
- **Any guard you add carries its expiry** (`EXPIRES-WITH:`) — the private-globals guard expires with the PSALTer pin. The loader is permanent and carries none.
- **Follow the ecosystem's conventions exactly** — Cobaya's loader, Cobaya's error idioms, XDG's data directory. A convention that is close but subtly different is worse than none.
- **Read the tool's own known-issues before hypothesizing** (PSALTer's README known bug #1 is the one Wave 0 spent a week attributing to the engine).
- No kernel, ever (the list in the header).
- No environment-specific absolute paths in anything committed; run `tests/test_repo_hygiene.py` **after `git add`** — it scans tracked files.
- Conventional commits; American English; no attribution trailers.

## If you find the design wrong

Amend a design-document error at the instruction site and report it. An architectural
contradiction: stop and report.

## Report back

Branch · PR · `CI <run-id>: <conclusion>` · each criterion with its output · the new Vector
fingerprint and why it differs from the prototype's · where the `A[0]` unknown and the named
schema gaps are written (docstring paths, not only this report) · amendments made · discoveries to
route · suggested next step.
