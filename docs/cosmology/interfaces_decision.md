# Interfaces decision memo — user input, derivation, derived spectrum, Cobaya (R-1, #566)

> **Status: RECOMMENDATIONS, 2026-09-15 — not decisions.** D-A and D-B are recorded in
> `docs/COSMOLOGY_PROGRAM.md` by the orchestrator after the options go to the user.
> Evidence ledger and the planning revision log: `r1_planning_record.md`. Conventions
> reference: `conventions.md`. Prototypes and their reproduction commands:
> `scripts/research/interfaces/`.

<!-- cspell:words etak hdot mochi SUSY cymetric WSTP sympify covmat venvs MDMSM gravitymodel UFOs hiclass FlexibleSUSY yamls pexpect reweights platformdirs uₐuᵃ lmax ipynb DEVEL -->

## 0. Summary

**D-A — how Python drives Wolfram. Recommendation: one committed Wolfram package, run by a
fixed driver as one `wolframscript` subprocess per derivation; Python passes the theory as WXF
*data* and never writes Wolfram code.** Every mature Wolfram physics tool reads a model with
one generic package (FeynRules, SARAH, PSALTer, xPand, HiGGS, Hamilcar; none templates code
per model from another language), and established pipelines launch one process per derivation
(FlexibleSUSY, SARAH, PSALTer's own materials). Run here: the package reproduced the reference wave operator exactly (`SameQ`, §2.3), all
eleven malformed or malicious inputs were refused before PSALTer ran (§2.4), and the error handling #561 deferred now lives once in a fixed file.
The live alternative, the same package through a `wolframclient` session, was measured and
**did not complete**: `DefField` never returns inside a session (§2.5).

**D-B — what the user writes and how the derived spectrum reaches Cobaya. Recommendation:
two YAML files, one format (Option A′).** A theory file holds the physics — `fields`, and
`lagrangian` as a mapping from each coupling to its operator in PSALTer notation. Standard
Cobaya run files include it with Cobaya's own `!defaults`, and carry priors in top-level
`params:` exactly as for any theory parameter. The user runs `tidalcosmo derive run.yaml`
(needs Wolfram) and then `cobaya-run run.yaml` (does not), mirroring Cobaya's own
`cobaya-install` → `cobaya-run`. The derived spectrum lives in a content-addressed *data*
store, found by a fingerprint of the theory that the user never handles; the Cobaya Theory
refuses at construction with Cobaya's own `ComponentNotInstalledError` when it is missing,
stale, corrupt or in the wrong conventions (§3.4, twelve cases run). A run cannot start a
kernel: its import closure contains no launcher (checked statically, and shown able to fail).

**Convention carriage.** PSALTer and CAMB agree on (+,−,−,−); CAMB's Fortran is
signature-agnostic and defines no ε, so PSALTer's ε₀₁₂₃ = +1 is adopted; CAMB's
perturbation-variable definitions govern the seam. `spectrum_design.md:301, 312-313, 354`
are wrong about CAMB (§4). The rule, the pinned CAMB variable table and the legacy-comparison
rule are in `conventions.md`.

**M3's output location.** A spectra store inside the repository (same layout and manifests),
never `examples/data/`, which holds the frozen legacy oracle (§3.6).

## 1. What was verified here, and what was not

### 1.1 Run ledger (one Wolfram kernel at a time, 2026-09-15)

| step | what | kernel wall | result |
| --- | --- | --- | --- |
| lane | `ensure_registered.sh --no-verify` then `verify-wolfram-setup.sh --require-psalter` | 1 min 45 s (both) | both **exit 0**; verify: "All checks passed!", 1 warning ("Could not check xPerm binary dependencies", informational) |
| reference | `scripts/psalter/vector_smoke.wls`, unchanged, throwaway directory | 48 s | `PREFLIGHT PASSED`, `ParticleSpectrum` 26.1 s, 2 sectors |
| package, validate only | `derive_proto.py theories/VectorTheory.yaml --validate-only` | 9.1 s | `STAGE1 PASSED` |
| package, full | `derive_proto.py runs/vector_one.yaml` (theory via `!defaults`) | 55.6 s | stored; `ParticleSpectrum` 27.2 s; wave-operator digest `36115df6…d3dd76` |
| failure fixtures | eleven theories through the full driver | ≈ 10 s each | all refused, exit 1, nothing stored (§2.4) |
| injection control | `ToExpression` of the injected string, with and without `HoldComplete` | 2 s | without `HoldComplete` it prints `INJECTED`; with it, nothing runs (§2.4) |
| session, run 1 | `session_probe.py` (whole derivation, then timeout and `Quit[]`) | 12 min 39 s | derivation hung for its 600 s watchdog; later steps confounded (§2.5) |
| session, run 2 | `session_probe_steps.py bare` / `steps` | 25 s + 2 min 16 s | timeout ignored; `Quit[]` raises and loses state; `DefField` hangs (§2.5) |
| identity | `compare_waveoperators.wls` (+ a perturbed-input control) | 14 s + 25 s | **identical (`SameQ`)**; control proved-different (§2.3) |
| lane | `verify-wolfram-setup.sh --require-psalter` (after) | 2 min 8 s | **exit 0**, "All checks passed!", same single informational warning |

Total kernel time: **23 min 52 s** of the 30-minute budget, lane checks included (the 600 s
session watchdog is most of it). Logs are under gitignored
`third_party/interfaces_runs/` (reproduction commands in `scripts/research/interfaces/README.md`).

### 1.2 Not tested here, and why

| item | why not | what a test would add |
| --- | --- | --- |
| a torsion theory through the package (geometric objects expanded before PSALTer) | R-1 is Vector-only by scope; the expansion step is I-S1B's port | that operator whitelists and the post-Riemannian expansion compose |
| a real likelihood with CAMB | out of scope (no CAMB in the gate prototype) | ordering of a spectrum check relative to CAMB (I-532) |
| hi_class, EFTCAMB, CLASS, SymBoltz.jl, DISCO-EB built and run | pinned-source reads per the dispatch; their input formats are quoted from source (§5) | confirmation that documented formats are the real ones, and error behavior |
| per-sample stability hook inside Cobaya | I-532's decision; only evidence gathered (§3.5) | whether an external prior can locate the derived spectrum |
| CAMB notes PDF equations | read as extracted text only | machine-checked equations |
| why `DefField` does not return inside a `wolframclient` session | located, not diagnosed: kernel budget, and D-A does not depend on the cause | whether any session configuration could work (not needed by the recommendation) |

## 2. D-A — how Python drives Wolfram

### 2.1 Options

| | **A. committed package + fixed driver, subprocess** (recommended) | **B. same package, `wolframclient` session** | **C. Python-templated script per theory** (legacy) |
| --- | --- | --- | --- |
| Python writes | a WXF data file (kernel-free `wolframclient.serializers.export`) | calls over ZeroMQ | Wolfram source text |
| Wolfram code | `Stage1Proto.wl` (≈190 lines) + `driver.wls` (≈45 lines), fixed | same package | regenerated per theory; legacy: 5,764 of 6,756 lines of `tidal/cli/_derive.py` emit text |
| testability | package functions unit-testable (`tests/wolfram/test_gauge_fix.wls` pattern); fixtures in YAML | same, plus session plumbing | golden-text diffs only; legacy generated scripts are temp files, never kept (`_derive.py:6594-6613`) |
| where `Catch` + `Exit[1]` lives | once, in `driver.wls` | failures come back as objects, but a bare `Quit[]` raises and the next call silently runs in a fresh kernel without PSALTer (§2.5) | emitted into every script; #561 deferred it for that reason |
| `derivation_hash` | the fingerprint of the theory-defining keys + schema version; no script to hash | same | sha256 of generated text (changes with any template edit) |
| user-text safety | parse held, whitelist, then evaluate (§2.4) | same | legacy rewrote physics with regex (`_derive.py:193-265`) |
| readability | a physicist reads one package | same | reads a generator |
| precedent | FlexibleSUSY, SARAH, PSALTer materials, GUM (one link per run) | cymetric only | none among Wolfram physics tools |
| tested here? | **run**: identity (§2.3), 11 failures, injection control | **run: did not complete** — `DefField` never returns in a session (§2.5) | reference script only (it is itself a standalone script) |

### 2.2 Why A — evidence (full citations: `r1_planning_record.md` App. B2–B3)

- **Every mature Wolfram physics tool is one package reading a model**: FeynRules `LoadModel`
  (arXiv 1310.1921 lines 2678-2688), SARAH `Start["M"]` (1309.7223 lines 686-709), PSALTer's own
  models selected by `Get` (`SupplementalMaterials-2607@b49e9f1d:.../Models.m:3-18`).
- **Established invocation is one process per derivation, judged by exit status plus
  artifacts**: FlexibleSUSY runs `Get[...]; Quit[]` through `math` from make with `Quit[1]` on
  failure (`templates/module.mk.in:312-314`, `meta/Utils.m:344`).
- **PSALTer is hostile to a shared long-lived kernel**: bare `Quit[]` at
  `ConjectureInverse.m:37-39`; uncaught throws (`ParticleSpectrum.m:23`); it closes and relaunches
  the subkernel pool (`PreComputeComponents.m:32-34`).
- **`wolframclient` 1.4.0 (June 2024)**: Python-side timeout not passed through
  (`localsession.py:259-260`; issue #38 open) — **measured here: ignored** (§2.5); Linux
  startup/socket failures open (#17, #26, #33); `terminate()` left kernels running here.
  Its WXF (de)serializers are kernel-free (checked: no `wolframclient.evaluation`, `zmq` or
  `pkg_resources` imported) and round-trip `Fraction(-1,2)` to `Rational[-1, 2]`.
- **The polology JAX code reads PSALTer output from WXF without a kernel**
  (`SupplementalMaterials-2607:JAX/src/psalter/_extracting/wxf.py:25, 162`) — and replaced a
  regex-templated Julia exporter to do so (`JuliaExport.m:52-56`, `wxf.py:4-6`).
- **No design document argues for text emission**: `spectrum_design.md:392-394` sources it to
  legacy's generator; `architecture.tex` does not discuss it.

### 2.3 Identity: the package reproduces the reference wave operator

`wolfram/compare_waveoperators.wls`, reference first, in one kernel:

| label | source | wave-operator digest (sha256 of `BinarySerialize`) | stable twice | sectors | verdict |
| --- | --- | --- | --- | --- | --- |
| reference | `vector_smoke.wls` → `ParticleSpectrographVectorTheory.mx` | `36115df6579ece76cd70c72b4ce8336964b7aaa912d94060a914f158e4d3dd76` | True | 2 | — |
| package (`.mx`) | driver + package → `ParticleSpectrographTheory57a91b3528c2.mx` | `36115df6579ece76cd70c72b4ce8336964b7aaa912d94060a914f158e4d3dd76` | True | 2 | **identical (SameQ)** |
| package (`spectrum.wxf`) | the exported derived spectrum | `36115df6579ece76cd70c72b4ce8336964b7aaa912d94060a914f158e4d3dd76` | True | 2 | **identical (SameQ)** |

**The verdict is `SameQ`, not the digest, and it can fail**: the same script on a copy of
`spectrum.wxf` with one element perturbed by `Theta3/2` returns `proved-different` and exit 1.
The first element is `(Def^2*(-Theta1 + Theta2) - Theta3)/2` — the mass term is present, so the
caveat at `scripts/psalter/vector_smoke.wls:25-35` (association "missing the mass term") no
longer describes this install (#543 resolved). The session route is not in this table: it never
produced a wave operator (§2.5).

### 2.4 The package refuses bad input before PSALTer runs

All eleven fixtures (`scripts/research/interfaces/theories/bad/`) went through the full
driver: exit 1, **no store entry**, and a message naming the problem.

| fixture | refused by | message (abridged) |
| --- | --- | --- |
| `m_phi` coupling | name check | `BadName: coupling name "m_phi" must be letters and digits … (no underscores: m_phi parses as a pattern)` |
| undeclared `BField` | whitelist | `ForbiddenSymbol: … Global`BField` |
| key `'2'` | name check | `BadName: coupling name "2" …` |
| `Theta2` inside `Theta1`'s operator | whitelist | `CouplingInOperator: … contains the coupling(s) Theta2; each coupling must multiply its operator from outside (S = Sum c_i O_i)` |
| two lines, second starts `+` | one-expression check | `NotOneExpression` |
| operator `0` | zero check | `ZeroOperator` |
| coupling named `G` | context check | `NameCollision: coupling name G is already a symbol in xAct`PSALTer`` |
| `Print["INJECTED"]; …` | whitelist | `ForbiddenSymbol: System`CompoundExpression, System`Print` — `INJECTED` never printed |
| `Quit[]` | whitelist | `ForbiddenSymbol: System`Quit` — kernel kept running to report it |
| `0.5*…` | atom check | `InexactOrStringAtom … (use exact rationals like 1/2)` |
| duplicate `Theta1` key | Cobaya's YAML loader | `InputSyntaxError: Duplicate key Theta1` |

**The injection result is real, not an absence.** Control (same kernel build):
`ToExpression["Print[\"INJECTED\"]; 1"]` prints `INJECTED`;
`ToExpression[…, InputForm, HoldComplete]` returns `HoldComplete[Print["INJECTED"]; 1]` without
printing. **A second control settled an open question**: `SyntaxQ["-(1/2)*x\n+ y\n"]` is `True`,
but the held parse is `HoldComplete[(-(1/2))*x, Plus[y]]` — two expressions — so `SyntaxQ` alone
does not guarantee one expression; the structural check is required.

### 2.5 The session route, measured

**Run 1** (`proto/session_probe.py`): the same package, the same theory, a
`WolframLanguageSession` with `QT_QPA_PLATFORM=offscreen` confirmed in every kernel's
environment and the working directory set before PSALTer loads.

| step | observed |
| --- | --- |
| session start | ok, 5.7 s |
| `SetDirectory` → `Needs[PSALTer]` → `Get[package]` → `RunStage1` | **no result within the 600 s watchdog**; low CPU; no `spectrum.wxf`, no PSALTer PDF (a subprocess run of the identical package finished in 55.6 s) |
| `evaluate_wrap("Pause[8]; 1", timeout=2)` | no return within 60 s — **confounded**: the kernel was still busy with the step above |
| `evaluate("Quit[]")` | no return within 60 s — **confounded** for the same reason |

At the end of run 1, `session.terminate()` returned in 0.1 s but **left all five kernel processes
running**; the probe had to kill them — an orphaned kernel holds the single-license lane until killed.

**Run 2** (`proto/session_probe_steps.py`), each part in a fresh session, kernels killed between:

| part | step | observed |
| --- | --- | --- |
| bare (no PSALTer) | `evaluate_wrap("Pause[8]; 1", timeout=2)` | returned normally after **8.0 s** — the Python-side timeout is **ignored** (confirms wolframclient #38) |
| bare | `evaluate("Quit[]")` | **raised** `WolframKernelException: Kernel is not running anymore.` after 0.2 s |
| bare | next `evaluate("1+1")` | returned `2` after 6.4 s — a new kernel was started; all prior state is gone |
| bare | `terminate()` | returned at once; one kernel process remained and was killed |
| steps | `SetDirectory`, `Needs["xAct`PSALTer`"]`, `Get[package]` | ok (0.0 s, 4.7 s, 0.0 s) |
| steps | `RunStage1[..., "ValidateOnly" -> True]` | ok, couplings `Theta1, Theta2, Theta3` |
| steps | `DefField[ProbeField[-a], PrintAs -> "B", PrintSourceAs -> "k"]` alone | **no return within 120 s**; no PDF written |

**What this establishes.** Inside a `wolframclient` session on this install, PSALTer's `DefField`
does not return — while in a `wolframscript` subprocess it completes both in the full derivation
(55.6 s, §2.3) and in `verify-wolfram-setup.sh`'s `DefField` smoke test — so the session route
cannot derive at all here; the cause was not isolated further (both
runs had `QT_QPA_PLATFORM=offscreen` in every kernel's environment and no front-end process
appeared). **A planning hypothesis was wrong and is corrected here:** a bare `Quit[]` does not
hang a session — it raises, and the next call silently runs in a fresh kernel without PSALTer. The
timeout is ignored, and `terminate()` can leave kernels behind.

### 2.6 Recommendation for D-A, and consequences

**Recommend A.** Consequences:
- **I-S1B (exporter):** the exporter is a function in the committed package, called by the fixed
  driver; the spectrum WXF carries `schema`, `fingerprint`, `conventions`, sorted `couplings`,
  generator versions and the wave operator (prototype `Stage1Proto.wl:172-177`). #561's
  `Catch` + `Exit[1]` lives in the driver. Symbols must be created in `Global` (PSALTer derives
  contexts from `ToString`), the theory name is fingerprint-derived, and the artifact — not the
  exit status — is the verdict.
- **I-S1A-core (driver):** a Python launcher that self-guards the lane with `pgrep` (the hook
  cannot see Python-launched kernels), checks the pinned engine file (never `command -v`, #559),
  sets `QT_QPA_PLATFORM=offscreen`, writes atomically into the store.
- **R-C (D-C):** under this shape an FRW tool is judged on whether it exposes callable functions
  a committed package can call with data (xPand, xPert do), not on whether its API can be emitted
  as text.
- **Design-doc sites to amend (not edited here):** `spectrum_design.md:392-394`;
  `stage1_engineering_plan.md:502-538, 763, 803-805`.

## 3. D-B — what the user writes, and how the derived spectrum reaches Cobaya

### 3.1 Recommended input (Option A′), as run here

`scripts/research/interfaces/theories/VectorTheory.yaml` — the physics:

```yaml
# Conventions: signature (+,-,-,-), epsilon_0123 = +1 (PSALTer's).
fields:
  VectorField: {indices: [-a], print_as: A, print_source_as: j}
lagrangian:                   # S = ∫ Σᵢ θᵢ Oᵢ — key = coupling, value = its operator
  Theta1: '-(1/2)*CD[-a]@VectorField[-b]*CD[a]@VectorField[b]'
  Theta2: '(1/2)*CD[-a]@VectorField[a]*CD[-b]@VectorField[b]'
  Theta3: '-(1/2)*VectorField[-a]*VectorField[a]'
derivation: {max_laurent_depth: 1}
```

`scripts/research/interfaces/runs/vector_one.yaml` — the inference, standard Cobaya:

```yaml
theory:
  spectator_gate.SpectatorGate:          # production name: tidalcosmo.SpectatorTheory
    python_path: ../proto
    model: !defaults ../theories/VectorTheory     # Cobaya's own include
params:
  Theta1: {prior: {min: 0, max: 2}, ref: 1.0, proposal: 0.05, latex: \theta_1}
  Theta2: {prior: {min: 0, max: 2}, latex: \theta_2}
  Theta3: 0.5
likelihood:
  one:
sampler:
  evaluate:
```

Measured: Cobaya's file loader expands the include, and the fingerprint computed from the run
file equals the one computed from the theory file (`57a91b35…c81d359`). **The inline form**
(the same keys directly under `model:`) uses the same schema and fingerprint.

- **Why two files:** a large thesis Lagrangian lives in one place; many run files share it;
  `theories/` diffs are physics, `runs/` diffs are analysis. Precedent: SOLikeT's own run files
  (`examples/smooth/yamls/run_mflike.yaml:15-22` @`f6e00a9`: `theory: !defaults [theory_camb, …]`),
  Planck/MFLike `params: !defaults [...]` (`planck_2018_highl_plik/TT.yaml:14`;
  `mflike/TTTEEE.yaml:17`), CosmoSIS's separate values/priors files, GAMBIT's `!import`.
- **Why a mapping keyed by coupling, not a list:** Cobaya's include *replaces* lists on merge
  (tested); keyed terms merge term by term, and Cobaya's loader rejects duplicate keys
  (`yaml.py:104-113`), so "one operator per coupling" is enforced by YAML itself (fixture above).
- **Couplings are declared exactly as Cobaya users declare a theory's parameters**: priors,
  fixed values and labels only in top-level `params:` (`docs/likelihoods.rst:19`;
  `theories/camb/camb.py:33-43`); the Theory registers the couplings from the terms in
  `get_modified_defaults`, the pattern LAT_MFLike's `Foreground` uses for option-dependent
  parameters (`mflike/foreground.py:215-225`).
- **Operators are PSALTer notation**, the same text a PSALTer user types; conventions and input
  requirements in `conventions.md` §2, §4.

### 3.2 Options and trade-offs

| | **A′ theory YAML + run YAML** (recommended; inline A same schema) | **B run YAML names a Wolfram model file** (GUM/FeynRules pattern) | **D named operators from a fixed vocabulary** |
| --- | --- | --- | --- |
| formats a Cobaya user learns | none new | Wolfram + xAct declaration syntax | none new |
| couplings visible to Cobaya without a kernel | yes (errors at `cobaya-run`) | no (only after `derive` writes a manifest) | yes |
| linearity in couplings | by structure | separate Wolfram check | by structure |
| staleness safety | fingerprint of theory keys → missing entry → refusal | fingerprint of file bytes | as A′ |
| provenance | manifest + Cobaya `.updated.yaml` records `version: spectrum-<fp>` | same | same |
| testability | fixtures are YAML; package unit-testable | model files execute code unless parsed held | trivial |
| what M3 writes and where | theory YAMLs + a committed spectra store (§3.6) | model files + store | vocabulary entries + store |
| implication for I-S1B's exporter | reads a WXF theory association | must parse a `.wl` held, then extract couplings | a catalogue to maintain |
| precedent | SOLikeT, Planck/MFLike, CosmoSIS | GUM — because FeynRules/SARAH model libraries exist; PSALTer has no model-file format (`Models.m:3-18`) | hi_class (`include/background.h:20-27`), EFTCAMB (`09_EFTCAMB_main.f90:535-547`) |
| tested here? | **run** (§3.1, §3.4) | not run | not run |

A single monolithic Lagrangian string was dropped: term-level input is strictly better and was
already decided (`stage1_engineering_plan.md:483-496`).

### 3.3 What the user types

```bash
tidalcosmo derive runs/vector_one.yaml     # once per theory; needs Wolfram; minutes to hours
cobaya-run      runs/vector_one.yaml     # as often as needed; no Wolfram
```

Precedent: Cobaya's own `cobaya-install MyFile.yaml` then `cobaya-run` on the same file
(`docs/installation_cosmo.rst:28`). The derivation cannot be wired into `cobaya-install`: its
hooks never receive the component's YAML block (`install.py:232-243, 374`). Prototype flags
(`derive_proto.py`): `--test` (check only, never starts a kernel; mirrors `cobaya-install --test`),
`--force`, `--dry-run`, `--timeout` (default: no cap — legacy's 600 s is wrong for hours-long
derivations), `--validate-only`. `derive` also accepts a theory file directly. **Run** (no kernel
started by any of these): `--test` with the spectrum stored → `DERIVE_UP_TO_DATE`, exit 0; `--test`
against an empty store → `DERIVE_TEST=missing … searched <stores>`, exit 1; `--dry-run` → prints the
internal `wolframscript -file driver.wls …` command and target entry, exit 0; a second `derive`
without `--force` → `DERIVE_UP_TO_DATE`, exit 0.

### 3.4 The gate: sampling refuses without the matching derived spectrum (Cobaya 3.6.2, run)

`proto/spectator_gate.py` refuses in `initialize()`, which Cobaya calls during construction
(`component.py:400-401`), before any likelihood or sampler exists. Twelve cases, each from real
YAML files in a fresh interpreter (`proto/run_gate_cases.py`, transcript
`third_party/interfaces_runs/gate_cases.txt`):

| # | case | expected | observed (abridged) |
| --- | --- | --- | --- |
| 1 | no derived spectrum on the search path | refuse | `ComponentNotInstalledError`: "No derived spectrum for this theory (fingerprint 57a91b3528c2). Generate it with tidalcosmo derive &lt;your input&gt;.yaml on a machine with Wolfram (minutes to hours), or copy it from one. Searched: &lt;empty_store&gt;, &lt;home&gt;/.local/share/tidalcosmo/spectra" |
| 2 | operator edited after deriving | refuse | same message, fingerprint `b1f0656898a4` — a stale spectrum can never be loaded |
| 3 | `spectrum.wxf` truncated | refuse | `… is unusable: spectrum.wxf sha256 does not match its manifest (partial or edited file)` |
| 4 | manifest in legacy conventions | refuse | `… conventions {'signature': [-1, 1, 1, 1], 'epsilon0123': -1} != required {'signature': [1, -1, -1, -1], 'epsilon0123': 1}` |
| 5 | coupling `Theta2` missing from `params:` | refuse | Cobaya's own: `[model] *ERROR* Requirement Theta2 of spectator_gate.SpectatorGate is not provided by any component, nor sampled directly` |
| 6 | typo `Theta4` (real, non-absorbing likelihood) | refuse | Cobaya's own: `[model] *ERROR* Could not find anything to use input parameter(s) {'Theta4'}.` |
| 6b | same typo with the unit likelihood `one` | passes | **caveat:** `one` absorbs unclaimed parameters (`AbsorbUnusedParamsLikelihood`), so typo protection needs a real likelihood |
| 7 | everything matches | pass | `logpost = -1.386…`, component version `spectrum-57a91b3528c2…` |
| 8 | different priors, a different coupling fixed | pass | no re-derivation (fingerprint unchanged) |
| 9 | spectrum only in the second store of the search path | pass | found |
| 10 | resume a chain with a different theory | refuse | Cobaya's own: `Old and new run information not compatible! Resuming not possible!`; `.updated.yaml` recorded `version: spectrum-57a91b35…` |
| 11 | resume with the same spectrum in a moved store | pass | store location is not part of the physics |

**Static no-kernel check:** the AST of the Cobaya-facing modules (`spectator_gate.py`,
`theory_store.py`) imports no `subprocess`, `multiprocessing`, `pexpect` or
`wolframclient.evaluation` and calls no `os.system`-style launcher — PASS; the same check
applied to `derive_proto.py` and `session_probe.py` FAILS (controls), so it is a real gate. The
kernel-free WXF reader stays allowed.

**Two design corrections found by running it:** (i) the "compare assigned parameters with the
couplings" check (Planck precedent) must not raise on *missing* couplings — it runs before
Cobaya's own requirement check and replaced Cobaya's standard message with ours; it now raises
only on assigned non-couplings. (ii) Prototype runs that use `one` do not catch parameter typos
(case 6b).

### 3.5 The per-sample stability check — evidence for I-532 (not decided here)

The gate guarantees the derived spectrum is present, intact and matching before any sample.
Where the per-sample health check hooks in is I-532's (register row "Validity mechanism follows
group practice"). Evidence: TorC applied PSALTer's ghost/tachyon conditions **offline** and fixed
the couplings, putting only numerical-validity exclusions into Cobaya external priors
(`paper_Qtorsion.tex:180-188, 419, 632`; audit `:562`); polology treats health as a prior on the
couplings, offline, then reweights (`2606.30785` TeX `:341, 624, 820-840`). Cobaya's external
prior returning `-inf` is the only hook guaranteed to skip every theory and likelihood
(`model.py:650-667`) — but it receives only parameter values (`prior.py:553-555`), so it must
locate the derived spectrum itself. **There is no group precedent for a per-sample spectrum check
inside Cobaya.** Filed as #573.

### 3.6 Storage, and where M3 writes

- **A *data* store, not Cobaya's `packages_path` and not a cache.** `packages_path` is documented
  for installable "external packages … the original code itself, a cosmological dataset"
  (`docs/installation_cosmo.rst:4`); XDG defines the cache directory as "non-essential (cached)
  data", and platformdirs says cache "can be safely deleted without losing functionality"
  (`docs/explanation.rst:130-137`) — a spectrum that needs hours and a Wolfram license to
  regenerate is neither.
- **Layout (run):** `<store>/v<schema>/<fp[:2]>/<fp>/{spectrum.wxf, manifest.json}`; manifest
  holds fingerprint, canonical theory, couplings, conventions, spectrum sha256, generator
  versions, timings, host; write via temp directory + rename; verify sha256 on read.
- **Lookup order (run):** `spectra_path:` option → `TIDALCOSMO_SPECTRA_PATH` (colon list; first
  writable) → `$XDG_DATA_HOME/tidalcosmo/spectra`. Precedents: Cobaya's own override chain,
  Julia's depot stack, Snakemake's between-run cache, DVC's content-addressed store. On CSD3 the
  env var points at RDS, not `/home`.
- **Dependency:** `platformdirs` is not installed or declared; the prototype uses the XDG
  variables directly (routed to this wave's `pyproject.toml` owner).
- **M3:** re-derived specs go into a store committed inside the repository, with manifests, and
  **never** into `examples/data/` — the frozen legacy oracle M3 compares against
  (`repo_reshape.md:1011-1016`).

### 3.7 Recommendation for D-B, and consequences

**Recommend A′ with the gate and store above.** Consequences:
- **I-S1A-core (input model):** the theory schema is `fields` / `lagrangian` (coupling →
  operator) / `derivation`; loaded with Cobaya's file loader; fingerprint = sha256 of the canonical
  theory keys + schema; Python validates structure only, Wolfram validates symbols.
- **M1b (the Cobaya Theory):** `get_modified_defaults` coupling registration;
  `ComponentNotInstalledError` refusal; fingerprint as `get_version`; no `is_installed`/`install`.
- **I-532:** the per-sample hook question of §3.5, with its evidence.

## 4. Convention carriage

**Rule** (`conventions.md` §1): each quantity takes the convention of the tool that defines it.
Signature: PSALTer's (+,−,−,−), and CAMB agrees — `camb/symbolic.py:460`,
`docs/source/variables_guide.rst:223`, the CAMB notes ("uₐuᵃ = 1 signature"); CAMB's Fortran
equations are signature-agnostic (`equations.f90:2427, 2430-2431, 2442`). ε: PSALTer's +1
(`PSALTer.m:101`, `DefGeometry.m:69-72`); CAMB defines none. Perturbation variables: CAMB's,
pinned in `conventions.md` §3 — the sign hazards are `etak = kη_s = −kη/2`, `ḣ_s = 6ḣ`, φ as the
Weyl potential, and Hu et al.'s opposite Φ_N.

**Which terms change sign between conventions**: a scalar term's sign factor has three sources —
the number of inverse metrics (each flips under g → −g), how each field is defined relative to
the metric, and ε orientation for parity-odd terms. **Where a slip flips a ghost verdict**: the
parity factor multiplying every residue is defined by the signature (`spectrum_design.md:320-327`).
**The checks:** every manifest records `conventions` and the Theory refuses others (case 4); a
reference theory with a published verdict runs in CI; a `camb.symbolic` variable-mapping test at
the CAMB seam (I-532 / M1a). **M3's legacy comparison** attributes every difference to convention
or physics (`conventions.md` §6).

## 5. Comparables and pipelines

### 5.1 Cosmology ecosystem (group 1)

| package | status here | where the user writes the model (verbatim) | how precomputed data reach the sampler | staleness check |
| --- | --- | --- | --- | --- |
| Cobaya 3.6.2 | installed, run | `theory: camb: extra_args: …` / `params:` (`camb.py:33-43`) | `cobaya-install` → `packages_path`; `path:` | `version.dat` vs release (`InstallableLikelihood.py:104-116`); none on content |
| CAMB 2.0.4 | installed, read | Python/`.ini`; custom sources as sympy (`model.py:1111`) | `camb.symbolic` compiles Fortran at call time | in-memory cache keyed on source text (`symbolic.py:906-907`) |
| SOLikeT CosmoPower | **installed and run** (`soliket[emulator]==0.4.2`, throwaway venv, cobaya 3.6.2) | `theory: soliket.CosmoPower: network_path: … network_settings: tt: {filename: cmb_TT_NN}` (`cosmopower.py:49-62`) | path in YAML, loaded in `initialize()`; networks vendored in the package | **run:** vendored network → model builds; missing file → raw `FileNotFoundError: … does_not_exist_NN.pkl` (no Cobaya message, no hint); missing trained parameter → Cobaya's own `Requirement ns … not provided`; no identity or range metadata in the pickle (`cosmopower_NN.py:322-327`), so nothing ties a network to its theory |
| hi_class | source read `4f0aad6` | `gravity_model = propto_omega` / `parameters_smg = 1., 0., 0., 0., 1.` (`hi_class.ini:125-126`) | compiled catalogue | n/a |
| EFTCAMB | source read `16d9c4e` | `EFTflag = 2`, `AltParEFTmodel = 2` (`params_EFT.ini:14, 26`) | compiled catalogue | n/a |
| mochi_class | source read `2b0b16e` | `gravity_model = stable_params`, `smg_file_name = …` (`designer_fr.ini:883-884`) | table file | "Could not open file" only |
| CLASS | source read `64bbab7` | `output = tCl,pCl,lCl,mPk` (`explanatory.ini:72`) | n/a | n/a |
| DISCO-EB | source read `3aa88d2` | `param['w_DE_0'] = -0.99` (`nb_minimal_example.ipynb`) | JIT in process | n/a |
| SymBoltz.jl | source read `3d1f20a3` | `M = ΛCDM(lmax = 16)`; `System(eqs, τ; …)` (paper TeX `:632-645, 697-716`) | in-process compile | n/a |
| TorC | audit + TeX | YAML booleans `external_rhopa: False` (`slegner/cobaya@414a2e89:camb.yaml:26`) | per-sample provider Theory | none; no run config archived (audit §5.4) |

### 5.2 Symbolic → numeric pipelines (group 3)

| pipeline | generates code or passes data | what broke / lesson |
| --- | --- | --- |
| FeynRules → UFO → MadGraph | model file (Wolfram) in; UFO = generated Python module | text formats could not carry general structures (UFO 1108.2040 `:227-235`); kernel restart on model change (1310.1921 `:2688`); non-compliant vertices silently dropped (`:3779`) |
| SARAH → SPheno | model file in; Fortran generated | tadpole solving restricted by Mathematica `Solve` (1309.7223 `:995`) |
| FlexibleSUSY | spec file in; C++ generated by `math` subprocess | stamp-file staleness; `Quit[1]` on assertion |
| GUM (GAMBIT) | YAML names FeynRules/SARAH files; WSTP link per run | dry run still runs Mathematica (`gum.py:697-701`) |
| sympy codegen | `codegen((name, expr), "C89", …)` returns source *strings* (files only with `to_files=True`); `lambdify` builds Python source and `exec`s it (`sympy/utilities/codegen.py:1992`; `lambdify.py:182`) | `spectrum_design.md` §6.1 rules out sympy/`lambdify` on the per-sample path |
| `camb.symbolic` | sympy in; Fortran body generated, template hard-coded | unknown functions pass through (`allow_unknown_functions=True`, `symbolic.py:838-840`) |
| SymBoltz.jl | symbolic equations in; Julia code generated | slower C_ℓ than CAMB/CLASS (TeX `:954`) |
| numerical polology | Wolfram once → WXF data → JAX | symbolic "scales poorly … due to expression swell" (2606.30785 `:371`); identity = filename, no provenance (`wxf.py:165-172`) |
| legacy TIDAL | Python f-strings → `.wls`; `InputForm` strings back | bracket imbalance, underscore patterns, regex renames, silently dropped terms (`.claude/rules/wolfram.md`) |
| CppTransport, Cadabra | not read — no local sources; Cadabra known only via citations | — |

## 6. Findings routed to the orchestrator (not edited here)

- `spectrum_design.md:301, 354` mislabel CAMB as (−,+,+,+); `:312-313`'s "sign of det g" reason
  is wrong; `:392-394` models the generator on legacy.
- `stage1_engineering_plan.md:502-538, 763, 803-805` assume text emission.
- M3's written mapping needs a convention drift class (`tests_cosmo/data/oracles/README.md:80-110`).
- No β / polarization-handedness convention for O4.
- I-532: the per-sample hook evidence of §3.5.
- `platformdirs` dependency; legacy `derive` flags `--save-script`/`--dry-run`/600 s timeout.
- `R-1.md:38-40` cites `tidalcosmo/derive/README.md:12` for a claim that file does not make;
  `vector_smoke.wls:25-35` carries a caveat that expired with #543; `scientific_review.md` #28
  ("CAMB not installed") is stale.
- Citation slip in the planning record (non-authoritative archive): it gives
  `variables_guide.rst:133` for the synchronous line element; `:133` is the `etak` definition and
  the line element is `:223` (`conventions.md` has the correct pins).
- **The repo's pre-commit `check-yaml` hook rejects Cobaya's `!defaults` tag** ("could not determine a
  constructor for the tag '!defaults'"; run via `uvx --from pre-commit-hooks check-yaml`), so every
  Option A′ run file — and the deliberately invalid duplicate-key fixture — fails it. Pre-commit is
  not run in CI and no git hook is installed; if it is enforced, `.pre-commit-config.yaml` needs
  `--unsafe` or an exclude for Cobaya input files.
- `docs/README.md` needs rows for `conventions.md`, `interfaces_decision.md`,
  `r1_planning_record.md`.

GitHub issues filed: **#569** (`spectrum_design.md` CAMB convention errors), **#570** (M3 mapping
convention drift class), **#571** (O4 β / handedness convention), **#572** (PSALTer's bare
`Quit[]` and uncaught throws), **#573** (per-sample health check: no group precedent; external
prior cannot see the derived spectrum).

## 7. Evidence index

| artifact | where |
| --- | --- |
| planning record (evidence ledger, revision log) | `docs/cosmology/r1_planning_record.md` |
| conventions reference | `docs/cosmology/conventions.md` |
| prototypes and fixtures | `scripts/research/interfaces/` |
| logs, stores, transcripts | gitignored `third_party/interfaces_runs/` |
| derived spectrum (Vector) | fingerprint `57a91b3528c234e571dcf59a1e0b47d04b0bb063bcb96a95a2d365bc2c81d359`; `spectrum.wxf` sha256 `1777b64d48d950c900fe7e5e8373560e30413f11b1495cf581849ffdd8550f49`; wave-operator digest `36115df6579ece76cd70c72b4ce8336964b7aaa912d94060a914f158e4d3dd76` |
