> **ARCHIVE — verbatim record of the 2026-09-13…15 R-1 (#566) planning session plan file.**
>
> Preserved here because it lived only in `~/.claude/plans/` (not version-controlled;
> its backup directory is gitignored) and carries the full reasoning trail: every research
> pass, a pinned evidence ledger, and eighteen rounds of user correction that reshaped the
> original R-1 plan (Appendix A).
>
> **Not authoritative.** The living documents supersede it wherever they differ:
> `docs/cosmology/interfaces_decision.md` (the R-1 decision memo — options, evidence,
> recommendations), `docs/cosmology/conventions.md` (signature, ε, CAMB variable
> definitions, PSALTer input requirements, legacy-comparison rules) and
> `docs/COSMOLOGY_PROGRAM.md` (decisions register and wave board — D-A and D-B are recorded
> there by the orchestrator, not here). Everything below the rule is the plan file exactly as
> approved, with one redaction: §7 step 1 named the forbidden path patterns literally, which
> `tests/test_repo_hygiene.py` rejects, so that phrase is reworded as in `handoffs/R-1.md:186`.
> "You" in the text addresses the user who dispatched R-1.

<!-- cspell:words hdot Scal dgrho dgpi Challinor sympify mochi upvalues SUSY cymetric yamls hexic platformdirs cdef reweights venvs ayprime clxcdot birefring bispectrum shdw WSTP MDMSM hiclassy covmat Theo uₐuᵃ dxⁱdxʲ ipynb DEVEL -->

---

# R-1 (#566) — interfaces decision memo: plan

## TL;DR

- **Problem.** Four Wave-1 prompts are held on three open questions: how Python drives Wolfram (D-A), what the user writes and how the derived artifact reaches Cobaya (D-B), and how sign conventions are carried.
- **Research done (§2).** All of it read-only: source code, papers and issue trackers. It answers most of the memo.
- **Recommendations the memo will make (§3). The memo recommends; the user decides.**
  1. **Two YAML files, one format (Option A′).**
     - **A theory file** (`theories/VectorTheory.yaml`) holds the physics: fields, plus `lagrangian:` as a mapping from each coupling to its operator string in PSALTer's notation.
     - **Run files** hold standard Cobaya `params:`, `likelihood:` and `sampler:`, and include the theory with Cobaya's own `!defaults`.
     - This is the established pattern: SOLikeT run files, Planck/MFLike, CosmoSIS, GAMBIT. The include was verified on Cobaya 3.6.2.
     - A single-file inline form uses the same schema.
     - The live alternative, a Wolfram-code theory file (GUM/FeynRules style), is weighed and not recommended (§3.1).
  2. **The user runs two commands on the same file**, mirroring Cobaya's `cobaya-install` then `cobaya-run`: `tidalcosmo derive run.yaml` (needs Wolfram), then `cobaya-run run.yaml` (does not).
     - `derive` has `--test`, `--force`, `--dry-run` and `--timeout`, and a Python API mirror.
     - Underneath, **one committed Wolfram package** does all the work. Python passes data only, as WXF, and launches Wolfram as a background process. It never writes Wolfram code, and the user never types `wolframscript`.
  3. **Sampling refuses at construction** unless PSALTer's derived output for *this exact* theory block exists, found by a behind-the-scenes fingerprint the user never handles. This follows Cobaya's own `ComponentNotInstalledError` idiom.
     - **Storage** follows the standards for expensive derived data: a content-keyed store under the user data directory, movable by an env var (with shared read-only stores for HPC and collaborators), a manifest beside each file, and a checksum verified on read.
     - **Couplings** are declared exactly as Cobaya users already declare theory parameters: priors in top-level `params:`, the list advertised through the MFLike-style `get_modified_defaults` hook.
     - The per-sample stability check then reads that output.
  4. **Conventions:** each quantity uses the convention of the tool that defines it. For the metric signature, PSALTer and CAMB agree: both use (+,−,−,−). CAMB defines no ε, so we adopt PSALTer's. CAMB's perturbation-variable definitions are adopted verbatim.
- **The whole planning session is recorded, not just its conclusions.** This plan, including Appendix A (every correction and how the design changed), Appendix B (a pinned evidence ledger from every research pass) and Appendix C (what could not be verified), is archived verbatim as `docs/cosmology/r1_planning_record.md` in the **first commit**, following the `planning_session_record.md` precedent. Memory is updated with the project status and your standing preferences.
- **New durable output: `docs/cosmology/conventions.md`** (owned path added at your instruction). It is the single reference for signature, ε, CAMB's pinned variable definitions, PSALTer's input requirements, and the rules the M3 legacy comparison uses to separate convention differences from physics differences.
- **Still needs running (§4).** Prototypes that confirm the recommended design end to end, including its failure modes and an injection test, plus the two comparison runs R-1 requires.
- **Next step.** Approve, then start at §5.

---

## 1. What R-1 decides, in plain terms

A user does two things:

1. **Once per theory, a long symbolic run** (minutes to hours, needs Wolfram). PSALTer derives the wave-operator matrices, the *Stage-1 artifact*.
2. **Many times, sampling with Cobaya.** Per sample, a fast numerical check uses that artifact to decide whether the coupling point is healthy (no ghosts or tachyons). Then perturbations evolve on a CAMB background.

R-1 settles five interfaces around those steps:

- how Python starts the derivation (D-A);
- what the user writes;
- how sampling finds the artifact and **refuses** when it is missing or belongs to a different theory (D-B);
- conventions;
- where M3 writes re-derived specs.

Legacy design carries no weight here. It is cited only as a record of what went wrong.

---

## 2. Research findings (settled, with sources)

### 2.1 Conventions: who defines what, and where the tools compete

**Metric signature. No competition: PSALTer and CAMB agree.**

| Source | Signature | Evidence |
| --- | --- | --- |
| PSALTer (source and papers) | (+,−,−,−), ε₀₁₂₃ = +1 | `PSALTer.m:101` "West Coast signature (1,-1,-1,-1)"; `Sources/ReloadPackage/DefGeometry.m:69-72`; arXiv 2406.09500 `Manuscript.tex:102` |
| CAMB, declared | (+,−,−,−) | `camb/symbolic.py:460`; CAMB `docs/source/variables_guide.rst:133`; CAMB notes PDF ("Using the uₐuᵃ = 1 signature"); its formalism papers (astro-ph/9804301 line 308; astro-ph/9911481 line 157) |
| CAMB, Fortran at tag 2.0.4 (`a6de8cc5`) | signature-agnostic | `equations.f90:2427,2431,2442` are Ma–Bertschinger's synchronous-gauge equations term for term. The variables are unchanged under g → −g |
| Cobaya | none | a sampler; it defines parameter names, not physics conventions |

The docs match the code: CAMB's Fortran is written in variables that are the same under either signature, so no CAMB number depends on the choice.

**ε orientation. No competition.** CAMB defines none: a grep of every `.f90` and `.py` at the pinned tag finds no Levi-Civita or parity convention, and C_l^TB = C_l^EB = 0 is assumed. So we adopt PSALTer's ε₀₁₂₃ = +1. It matters because the programme's Lagrangians contain ε-odd couplings, for example the Nieh–Yan-type term `d14` in `examples/torsion_gertsenshtein/theory_parity_odd.toml:3-6`.

**Where CAMB *does* define things, we use CAMB's definition. "Wins" means over the alternatives that genuinely compete in the sources this programme uses:**
- **Legacy TIDAL**, which is (−,+,+,+) and uses its own variables.
- **Ma & Bertschinger and CLASS**, which use the synchronous-gauge h and η and the conformal-Newtonian ψ and φ. R-C validates against Ma–Bertschinger.
- **Hu et al.**, who use the opposite sign for Φ_N (CAMB `docs/ScalEqs.ipynb:278-285`).
- **SymBoltz**, which follows Ma–Bertschinger.
- **CAMB itself**, which exposes *two* families of names: covariant (`eta`, `hdot`, `sigma`) and Ma–Bertschinger-style (`eta_s`, `hdot_s`). They differ by stated factors (`symbolic.py:481`).

**How the conflict is resolved, without a run-time conversion.**
- Our FRW derivation and all CAMB-seam code are written in `camb.symbolic`'s **named variables**, which are CAMB's native "CDM" frame (`symbolic.py:403, 486`).
- Any formula taken from Ma–Bertschinger, CLASS or Hu et al. is rewritten **once, where it is imported**, using the relations CAMB itself states (e.g. `eta_s`↔`eta`, `hdot_s`↔`hdot`), with the source cited.
- `conventions.md` lists every competing definition beside CAMB's, so none is chosen by accident.

On the spectrum side there is **no competition**: all three Barker/PSALTer papers and the source agree on (+,−,−,−) and ε₀₁₂₃ = +1.

These are the real sign hazards at the CAMB seam. Every one is pinned:

| Quantity | CAMB definition | Source |
| --- | --- | --- |
| `etak` | k·η_s, with η_s = −η/2 (covariant η = −2·etak/k) | `symbolic.py:183-188`; `variables_guide.rst:223` |
| ḣ | h_s′ = 2kz; covariant hdot = h_s′/6 = kz/3 | `variables_guide.rst:252`; `symbolic.py:481` |
| φ | Weyl potential (Φ_N + Ψ_N)/2 | `symbolic.py:181, 439-440` |
| Φ_N, Ψ_N | signs as in Ma–Bertschinger; Hu et al. flip Φ_N | `symbolic.py:240-245`; CAMB `docs/ScalEqs.ipynb:278-285` |
| `grho`, `dgrho`, `dgq`, `dgpi` | a²κρ, a²κδρ, a²κq, a²κπ, with κ = 8πG | `equations.f90:1119-1125`; `constants.f90:49` |
| frame | "CDM" = synchronous gauge with covariant variable names | `symbolic.py:403, 486` |
| E₂ | related to the polarization tensor by a minus sign | CAMB notes PDF |
| time | conformal τ (= our η) | `repo_reshape.md:443-445` |

**Nothing in the repo documents these today.** Only fragments exist (`repo_reshape.md:248, 436-439`), and the definitions live only in the installed `camb/symbolic.py`. **R-1 therefore writes them into `docs/cosmology/conventions.md`** (§6), pinned to CAMB 2.0.4 / `a6de8cc5`, so the symbolic FRW derivation (R-C, I-S1B, M3) and every CAMB-seam test use one checkable set of definitions.

**Why other convention sources come up at all.** Nothing is imported from them. They matter only where the programme compares against them:

- **Ma & Bertschinger.**
  - R-C must reproduce the Newtonian-gauge scalar Einstein equations against it (`handoffs/R-C.md:69, 114`). Its signature is (−,+,+,+), so any mixed or raised-index component compared there flips sign.
  - It is also the (wrong) source of the CAMB row in `spectrum_design.md:301`.
- **Legacy TIDAL** ((−,+,+,+), ε₀₁₂… = −1, `tidal/wolfram/CommonUtilities.wl:30-34`).
  - M3 must show "same equations, same coefficients, same signs" against the frozen legacy oracle (`repo_reshape.md:787-790`).
  - The oracle includes `torsion_gertsenshtein_parity_odd`, whose ε-odd terms flip sign under the ε change.
  - **The comparison is of physics, not bytes.** The new package deliberately changes conventions, so its equations will not match legacy text. What must match is the physics, and every difference must be explained as either a convention change or a real change. That requires the written mapping to record a sign factor for each term. The factor has three sources:
    - **metric contractions:** each inverse metric flips under g → −g, so T_{λμν}T^{λμν} picks up (−1)³ (`docs/tex/pgt_stability_priors.tex:136-139`);
    - **field definitions:** how each field is defined relative to the metric, since a perturbation defined as δg_ab itself flips;
    - **ε orientation:** legacy ε₀₁₂… = −1 against PSALTer's +1.

    The parity-odd oracle specs are exactly where the ε factor appears. **R-1 states the rule and its sources, but does not claim the per-term table.** Deriving that table and checking it against the oracle is M3's job.
  - Today no document states this (`tests_cosmo/data/oracles/README.md:80-110` lists four drift classes, none of them conventions). `conventions.md` states the rule, with the flip factors, as the reference M3's mapping follows.
- **SymBoltz** and **Challinor–Lasenby.** No planned use of either as a convention source (only a survey item for R-C and a formalism citation). They drop out.
- **A genuine competing standard found, but not R-1's to settle.** Parity-odd CMB observables are on the ladder (O4 birefringence, EB/TB, `COSMOLOGY_PROGRAM.md:257`). No document fixes the sign convention for β or for polarization handedness.

### 2.2 How physics tools take a theory

- **Wolfram physics packages.** In FeynRules (`.fr`), SARAH (`Models/<M>/*.m`), PSALTer's supplementary models, xPand, Hamilcar and HiGGS, the user writes the model as Wolfram-language code that one generic package reads (arXiv 1310.1921 lines 3427–3582; arXiv 1309.7223 lines 1374–1437; `SupplementalMaterials-2607@b49e9f1d:.../Models.m:3-18`).
  - **Why HEP pipelines point to model files.** GAMBIT's GUM reads a YAML file that names a FeynRules or SARAH model file. It loads the model file "that the user has requested" and points users to the FeynRules model database (arXiv 2107.00030, lines 2795 and 3447). Those formats and model libraries predate GUM, and GUM exists to reuse them for many outputs (CalcHEP, micrOMEGAs, Pythia, SPheno). The model file is the *upstream tool's* format, not something GUM invented.
  - **PSALTer has no equivalent.** Its theory is a Wolfram expression typed into a live session, and its own supplementary materials choose a model by commenting out `Get` lines (`Models.m:3-18`). There is no model-file format and no model library. A separate model file for us would be a new format, not an adopted standard.
- **The repo has already decided on term-level input.** `stage1_engineering_plan.md:483-496` says: "Term-level, not one expression string … each term naming its coupling and its operator … Coupling per operator; reject, never auto-assign". The reason is that PSALTer requires the action to depend linearly on the couplings ("requires that all terms in the action be linearly parameterised by symbolic coupling coefficients", arXiv 2406.09500 `Manuscript.tex:130`) but never checks it. A non-linear coupling "does not raise — it yields a wrong … wave operator, silently" (`spectrum_design.md` §6.1 amendment). With one coupling per term, that requirement holds by construction.
- **Cobaya can already split one input across several YAML files.** Its `!defaults` directive loads other YAML files into the current one (`cobaya/yaml.py:71-101`). It works only when loading from a file, not from a Python dict. So a long or reusable theory definition can live in its own YAML file with no new format.
- **Cobaya users already write expressions as strings in YAML:**
  - `value: 'lambda logA: 1e-10*np.exp(logA)'` (`docs/src_examples/cosmo_basic/basic_camb.yaml:31`);
  - external priors and likelihoods as `lambda` strings (`prior.py:131-139`).
  - Cobaya turns these into functions with a plain `eval` and no validation (`tools.py:359-373`).
  - **What "unsafe loader" means, and why it matters little here.**
    - Cobaya reads YAML with PyYAML's full loader (`yaml.py:46` subclasses `yaml.Loader`). That loader obeys special tags such as `!!python/name:os.system`, so a malicious YAML file can make Python run code just by being opened. Tested here: `!!python/name:builtins.len` loads as the real `len` function, while the safe loader refuses it.
    - Together with the `eval` of `lambda` strings, this means Cobaya treats an input YAML like a script: **you only run YAML you trust**.
    - So our design adds no new risk by reading the same file. The Wolfram-side check on operator strings (§2.3) is mainly about **correctness**: catching `m_phi`, typos and non-linear couplings before PSALTer silently mis-derives. It also ensures nothing unintended executes inside the kernel.
- **Cosmology codes that take a Lagrangian.**
  - Only Hi-COLA does, putting SymPy strings in an `.ini` that are passed to `sympify`, which uses `eval`, with no validation (`Frontend/read_parameters.py:69`).
  - hi_class and EFTCAMB take a **catalogue name plus numbers**, confirmed from source (`include/background.h:20-27`; `09_EFTCAMB_main.f90:535-547`).
  - mochi_class takes catalogue names plus tables.
  - SymBoltz takes symbolic equations in Julia.

### 2.3 A Lagrangian string in YAML can be handled safely

The safe pattern: the **Wolfram side** parses the string *without evaluating it*, checks it, and only then evaluates. Python never edits the string.

- `ToExpression[s, InputForm, HoldComplete]` parses and wraps the result so nothing runs. `HoldComplete` also blocks upvalues, which xAct tensors have (Wolfram docs for `ToExpression` and `HoldComplete`).
- Then walk the held tree. Allow only the heads `Plus`, `Times`, `Power`, `Rational` and `Integer`, plus `CD`, the declared fields, the declared couplings and PSALTer's index symbols a…z. Refuse everything else.
- **An underscore name is caught automatically.** `m_phi` parses as `Pattern[m, Blank[phi]]`, and neither head is on the whitelist (Wolfram tutorial *OperatorInputForms*).
- **The kernel-free alternative** is `CodeParser`'s `CodeParse`, which returns a syntax tree with source positions (useful for error messages).
- **Established practice.** Mathematica StackExchange's canonical answer on safe `ToExpression` is exactly "parse with `HoldComplete`, then whitelist heads" (question 172553).

Four YAML traps, all tested with PyYAML 6.0.3 and Cobaya's loader:

- In a **plain** (unquoted) YAML scalar, ` #…` is silently stripped as a comment. Cobaya also expands `${HOME}`.
- A value starting with `@` or `[` fails to parse.
- Double quotes break `\[Theta]`.
- A `|` block whose lines begin with `+` or `-` may parse in Wolfram as *several* expressions.

**Fix:** the schema requires every operator to be a single-quoted string or a folded block scalar (`>-`), and the package refuses anything that is not exactly one expression. Per-term operators are short, so single quotes are the natural form.

### 2.4 How Python should invoke the package

**What we adopt, internally: the package runs Wolfram as a separate program.** This is plumbing inside the package, and **the user never types it**.

When the user runs `tidalcosmo derive runs/vector_planck.yaml` (the user-facing command is designed in §3.2), the package itself starts one background process, `wolframscript -file driver.wls theory.wxf <output-dir>`. That is the same way Python code calls any external program. The package waits for the process, then checks what it produced. The driver is a short committed file that loads the committed Wolfram package. One derivation is one process.

The alternative we are **not** adopting is a `wolframclient` session: a Wolfram kernel kept alive inside the Python process and sent commands one at a time.

**Established physics tools run one subprocess per derivation:**
- FlexibleSUSY runs `Get[...]` from make via `math`, with `Quit[1]` on failure (`templates/module.mk.in:312-314`, `meta/Utils.m:344`);
- SARAH is run the same way;
- PSALTer's supplementary materials use `math -run`.

A long-lived `wolframclient` session appears in only one physics code we found (cymetric).

**Why we do not drive the derivation through a `wolframclient` session:**
- **Its timeout is broken (#38, open).** The docstring says `timeout` raises, but `localsession.py:259-260` never passes it on.
- **Linux startup and socket failures are unresolved** (#17, #26, #33).
- **The last release was June 2024.** The `pkg_resources` fix is unreleased (#41, #47).
- **PSALTer would kill a session.** It calls a bare `Quit[]` (`ConjectureInverse.m:37-39`, exit 0, silent), throws uncaught (`ParticleSpectrum.m:23`), and relaunches its own subkernels (`PreComputeComponents.m:32-34`).
- **We would gain nothing.** The derivation is a one-off batch job, and a session only pays off for many repeated calls.

**Where `wolframclient` *is* used: WXF serialization only, never a kernel.**
- The WXF **reader** (`binary_deserialize`, as the polology code uses it) and the **writer** (`serializers.export`) are both kernel-free. Checked here: importing them loads no `wolframclient.evaluation`, `zmq` or `pkg_resources`.
- Exact rationals survive: `Fraction(-1,2)` round-trips to `Rational[-1, 2]`.
- So Python hands the theory block to Wolfram as WXF *data*, and Wolfram hands the artifact back as WXF.
- **No JSON anywhere.**

**Error signalling in the subprocess:**
- An uncaught `Throw` exits 0 (#561, probed in this repo).
- The committed driver wraps everything in `Enclose`/`Catch` and calls `Exit[1]` (Wolfram docs for `Exit` and `Enclose`). It is one fixed file, so #561's reason for deferring this ("it would change every generated script") no longer applies.
- Python still verifies the artifact positively, because of PSALTer's bare `Quit[]` and Wolfram 14.3's segfault on shutdown.

### 2.5 A generic package can take any theory: constraints from PSALTer's source

- **Field and coupling symbols must be created in `Global`, not in a package's private context.**
  - PSALTer turns the field's *name* into a context string (`DefField.m:37`), and `ValidateLagrangian.m:37` looks that context up.
  - A privately created symbol has a long qualified name, so the lookup fails. Worse, `DefField.m:52` silently creates the source tensor in a garbage context.
  - The package therefore creates symbols explicitly in `Global`.
- **The theory name must not collide with a PSALTer export** (`UpdateTheoryAssociation.m:8`).
- **PSALTer fixes its output directory when it loads** (`PSALTer.m:49-56`), so load it from the intended directory.
- **Indices are safe.** a…z live in PSALTer's *public* context (`DefGeometry.m:5-8`).
- **The theory name labels outputs only.** It never changes the wave operator (`ParticleSpectrum.m:41-66`).

### 2.6 How cosmology packages deliver artifacts and gate sampling

| Package | Artifact reaches sampler | Checks the artifact matches the inputs in use now? |
| --- | --- | --- |
| Cobaya/CAMB | `cobaya-install` into `packages_path`; `path:` | a version string only (`version.dat`) |
| SOLikeT CosmoPower | `network_path:` in YAML, loaded in `initialize()` | **no**; a missing file gives a raw `open()` error |
| mochi_class | `.dat` table path | only "could not open file" |
| TorC | per-sample provider Theory | **no** |
| Polology JAX | hard-coded WXF path | **no**; the theory's identity is its filename |

**Cobaya's gating idiom.**
- `cobaya-run` never installs anything; `run.py` has no install path.
- Components check during construction and raise `ComponentNotInstalledError` with "run `cobaya-install …`" (`InstallableLikelihood.py:36-57`).
- Wiring our derivation *into* `cobaya-install` fails, because `install()` never receives the YAML block (`install.py:232-243, 374`) and `--force` would re-run hours of work.
- SOLikeT breaks the idiom by installing at run time (`soliket/lensing/lensing.py:111-113`). We must never do that.

**What the missing staleness check means.**
- Cobaya confirms that a file *exists* and, for downloads, that its version is *new enough*.
- Nobody in the table confirms the file was produced from the theory being sampled *now*.
- Without that check, a user who edits their Lagrangian and forgets to re-derive gets stability verdicts for the **old** theory, with no error.
- We must build this check. There is no precedent to copy.

### 2.7 Where stability checks hook in, in the group's own work

- **TorC (arXiv 2507.09228).** PSALTer was used *offline* to justify coupling conditions: λ ≥ 0, μ < 0, (ν+2μ)(ν−μ) > 0 (`paper_Qtorsion.tex:180-188, 632`). The couplings were then fixed, and only cosmological parameters were sampled. Numerical-validity exclusions went in as **Cobaya external priors** (`:419`; audit `:562`).
- **Numerical polology (arXiv 2606.30785).**
  - The healthy region *is* a prior on the couplings, π(θ) (TeX `:341, 358`).
  - Nested sampling targets unitarity itself (`:820-840`), or health is applied as a hard post-hoc filter (`JAX/src/psalter/sample.py:169-195`).
  - Data enter afterwards as importance weights (`:624`; `Julia/src/run/reweighter.jl:4-6`).
  - There is no joint Cobaya run.
- **Cobaya's mechanics.**
  - An external prior returning `-inf` **skips every theory and likelihood by construction** (`model.py:650-667`).
  - `Theory.calculate` returning `False`, or raising, gives zero likelihood. It only saves computation if the check runs before CAMB, which depends on declaration order and is untested.
- **Already decided in the repo.**
  - Screening gate plus a designed-in soft prior (`spectrum_design.md:93, 837-842`).
  - Veto placement "inherited, not invented" (`COSMOLOGY_PROGRAM.md:215`), and **owned by I-532**.

**So is the suggestion an external prior returning `-inf` for unhealthy points?** Yes, that is where the evidence points:
- It is the only Cobaya hook guaranteed to run before, and skip, every theory and likelihood.
- It matches how TorC and polology treat health: as a prior.
- It carries both severities from one evaluator: `-inf` for a hard veto, a finite log-weight for the soft prior already designed in.

**What this evidence does *not* give us: a precedent for a per-sample spectrum check inside Cobaya.**
- **TorC's ghost/tachyon conditions never entered Cobaya.** PSALTer justified them offline, the couplings were fixed (λ = 0), and only cosmological parameters were sampled. What did go into Cobaya as external priors were *numerical-validity* exclusions (Ωm < 1 and the ϖr singularity), plus a NaN guard in their `planck_clik` fork.
- **Polology never runs inside Cobaya.** Its health check is an offline prior, applied afterwards by reweighting.
- So a coupling-level health check evaluated per sample inside a Cobaya run is **new ground** for this group. The design for it has to be justified from Cobaya's mechanics, not copied.

**This is recorded for the session that owns it, so it is not lost:**
- as a routed finding in the memo (§8);
- as a **GitHub issue** (§6) carrying this evidence (TorC `paper_Qtorsion.tex:180-188, 419, 632`; audit `:431-441, 562`; polology TeX `:341, 624, 820-840`; Cobaya `prior.py:553-555`, `model.py:650-667`);
- with the open questions spelled out below.

**The wrinkle, which I-532 must solve.**
- A Cobaya external prior is a function built from a YAML string, and it receives **only parameter values** (`prior.py:553-555`). It does not see the theory block or the loaded Stage-1 artifact.
- So the prior must locate the artifact itself, and the same fingerprint gate (§3.3) must apply to it.
- The formal decision stays with I-532, per the register.

### 2.8 What the polology paper says, in plain terms

- Symbolic algebra handles any quadratic theory, but the algebra blows up as theories grow (arXiv 2606.30785, line 371).
- So the paper uses computer algebra for **one bounded step**: extracting the wave-operator matrices, which are linear in the couplings (line 1217).
- Everything expensive is then done numerically at each sampled point (line 1223).

That is our Stage 1 / Stage 2 split. The lesson from their code is that handing off by file works, but it needs an identity check, and theirs has none.

### 2.9 Supporting facts

- **Legacy templating cost.** 85% of `tidal/cli/_derive.py` (5,764 of 6,756 lines) emits Wolfram text. No generated script is kept in the repo, the bracket validator only warns (`:6418-6430`), and the physics is rewritten by regex (`:193-265`).
- **CAMB extension.** A new perturbation sector needs a Fortran fork of the DarkEnergy slot (scalar-only, back-reacts; `DarkEnergyInterface.f90:9-25`), or integration outside CAMB from `get_time_evolution`. Routed to I-532 / M1a.
- **Licenses** of imitable designs: SOLikeT, CosmoPower-Cobaya, mochi_class and SymBoltz are MIT. cosmopower core, DISCO-EB and EFTCAMB core are GPL-3.

---

## 3. Recommendations the memo will make

### 3.1 What the user writes

**Recommended (Option A′): the theory definition in its own YAML file, included by each run file.**

**File 1: `theories/VectorTheory.yaml`, the physics.** Written once, fingerprinted, and re-derived only when it changes.

```yaml
# Conventions: signature (+,-,-,-), epsilon_0123 = +1 (PSALTer's)
fields:
  VectorField: {indices: [-a], print_as: A, print_source_as: j}
lagrangian:                   # S = ∫ Σᵢ θᵢ Oᵢ  — key = coupling, value = its operator
  Theta1: '-(1/2)*CD[-a]@VectorField[-b]*CD[a]@VectorField[b]'
  Theta2: '(1/2)*CD[-a]@VectorField[a]*CD[-b]@VectorField[b]'
  Theta3: '-(1/2)*VectorField[-a]*VectorField[a]'
derivation: {max_laurent_depth: 1}
```

**File 2: `runs/vector_planck.yaml`, the inference.** One of these per analysis. All of it is standard Cobaya.

```yaml
theory:
  tidalcosmo.SpectatorTheory:
    model: !defaults ../theories/VectorTheory   # Cobaya's own include; path relative to this file
params:                                          # exactly as for CAMB's ombh2
  Theta1: {prior: {min: 0, max: 2}, ref: 1.0, proposal: 0.05, latex: \theta_1}   # sampled
  Theta2: {prior: {min: 0, max: 2}, latex: \theta_2}                            # sampled
  Theta3: 0.5                                                                   # fixed
likelihood: {...}
sampler: {...}
```

**What the split gives you.**
- The Lagrangian lives in one place, however many couplings and terms it grows to.
- Many run files (different likelihoods, priors or samplers) share it without copies drifting apart.
- A diff of `theories/` is a diff of the physics; a diff of `runs/` is a diff of the analysis.
- Both files are YAML, so there is no second syntax to learn.

**Why this is the established pattern, not a local invention.** Every framework we checked separates the model definition from the run configuration and links the two by reference:
- **SOLikeT's own run files** are assembled from definition files: `theory: !defaults [theory_camb, theory_BandpowerForeground]` (`examples/smooth/yamls/run_mflike.yaml:15-22` @`f6e00a9`).
- **Planck and MFLike** share parameter definitions with `params: !defaults [...]` (`planck_2018_highl_plik/TT.yaml:14`, about 40 files; `mflike/TTTEEE.yaml:17`).
- **CosmoSIS** documents "three configuration files": a parameter file for sampler and pipeline, plus separate values and priors files ([parameter_files](https://cosmosis.readthedocs.io/en/latest/usage/parameter_files.html)).
- **GAMBIT** pulls model blocks into scan files with `!import` (`yaml_files/DMsimp.yaml:7-10`).
- **MadGraph** separates the UFO model directory from run cards.

**Verified behavior of Cobaya's include** (installed 3.6.2, in-memory tests with real Planck files as targets):
- **Placement.** `model: !defaults ../theories/VectorTheory` works on a sub-key, with sibling options such as `spectra_path:` kept in the run file.
- **Paths.** Resolved relative to the including file. `.yaml` is found automatically. `../` and absolute paths work, and nested includes resolve correctly.
- **Provenance.** Cobaya writes the *expanded* content into `[prefix].input.yaml` and `.updated.yaml` (`output.py:642-653`), so chain records stand alone. The file path is lost, but the fingerprint recorded as the component `version` restores it (§3.3).
- **MPI.** Only rank 0 reads includes (`input.py:143-145`).
- **Resume.** It refuses if the included Lagrangian changed (`output.py:574-585`).

**Pitfalls found, and how the design handles each:**
- **Lists are replaced, not merged**, when files are combined. This is why `lagrangian:` is a **mapping keyed by coupling name** rather than the earlier list of `{coupling, operator}` entries. Layered files (a base theory plus extra terms) then merge term by term. Cobaya's loader also rejects duplicate keys (`yaml.py:104-113`), so "one operator per coupling" is enforced by YAML itself.
- **No per-run override beside the include.** A duplicate key is rejected and `<<` merge keys are unsupported. So everything that defines the physics (fields, terms, derivation options) lives in the theory file, and the run file carries only non-physics options.
- **Includes work only when Cobaya loads from a file**, not from a Python dict or YAML string. So the Python API offers `load_theory(path)`, and `derive` loads from files.
- **Bad files crash with raw Python errors**, not Cobaya input errors (issue #463, and a top-level list). So our loader validates the theory file itself and reports a readable error.
- **`${YAML_ROOT}` is wrong after an include.** It is documented, and we do not use it.
- **No package-relative include paths** (issue #56, open since 2019). This matters only if example theories ship inside the package, and a helper covers it.

**The inline form (Option A)** puts the same `fields:`, `lagrangian:` and `derivation:` keys directly under `model:` in a single run file, for quick one-off runs. It uses the same schema and the same fingerprint (tested in §4), so moving between the two forms needs no conversion.

- **Two readers, one theory definition.** `tidalcosmo derive` and Cobaya read the same `model:` content. Each operator uses PSALTer's published notation and conventions: (+,−,−,−), ε₀₁₂₃ = +1.
- **Couplings are declared exactly the way Cobaya users already declare a theory's parameters.** This was verified against Cobaya's docs, source and in-memory experiments, not assumed.
  - **The standard.** Priors, fixed values, LaTeX labels and `renames` go in the **top-level `params:` block**, and only there.
    - `docs/likelihoods.rst:19`: "Likelihood parameters are specified within the ``params`` block".
    - CAMB's own docstring puts "any param that CAMB understands, fixed, sampled or derived" under top-level `params:` (`theories/camb/camb.py:33-43`).
    - **In plain terms: the theory code says *which* parameters exist; the user's `params:` block says *what to do* with them.**
      - The theory code tells Cobaya "I need values called Theta1, Theta2, Theta3 at every sample". Cobaya uses that list to check the user supplied them, catch misspellings, and pass the right numbers in.
      - The *choices* about those parameters (prior range, fixed or sampled, starting point, label) belong to the user's analysis, so they live in `params:` and never inside the theory code or theory file.
      - **Why it matters here:** it is why the theory file holds coupling *names* only, and why changing a prior never changes the physics, the fingerprint, or the need to re-derive. CAMB works the same way: CAMB knows it needs `ombh2`, and your YAML says how to sample it.
  - **The mechanism, copied from the closest precedent.** In LAT_MFLike's `Foreground`, the parameter set depends on a user option (`mflike/foreground.py:215-225`). It builds its parameter list in the `get_modified_defaults` classmethod from that option, and users still write priors in top-level `params:`. Our Theory does the same:
    - it reads the couplings from the keys of `lagrangian:`;
    - it registers each coupling as a required parameter (`defaults["params"][coupling] = None`);
    - it also returns them from `get_can_support_params`, as `DEVEL.rst:150` asks for option-dependent parameters;
    - it stays non-agnostic.
  - **Errors are Cobaya's own, verified by experiment.**
    - A coupling missing from `params:` gives "Requirement Theta1 of … is not provided by any component, nor sampled directly".
    - A misspelled or extra parameter gives "Could not find anything to use input parameter(s) {'Theta2'}".
    - Users see nothing new. Fixing `Theta1: 0.3` holds it constant, and a `lambda` reparameterization works as usual.
  - **One extra check, with precedent.** In `initialize_with_params`, compare Cobaya's assigned `input_params` with the couplings from the terms, and raise on mismatch. The Planck likelihood does exactly this (`base_classes/planck_clik.py:123-135`).
  - **Pitfalls found by experiment, all avoided:**
    - relying on `get_can_support_params` alone: an omitted coupling **passed silently with a `nan` likelihood**;
    - a user `params:` key inside the theory block: it silently dropped a coupling;
    - priors or LaTeX beside each term: no precedent in Cobaya;
    - making the Theory agnostic: it swallows typos, and conflicts with CLASS.
  - **Changing a prior or fixing a coupling never forces a re-derivation.** Only the terms define the physics.
  - A related question stays open in the register: whether constants that are never sampled belong in the coupling vector, owned by I-S1B. The memo records how this layout interacts with it.
- **Linearity in the couplings holds by construction.** Every operator hangs off exactly one coupling key. Three things are refused with a hint: a key that is not a valid coupling name (e.g. a bare number), an operator that itself contains a coupling symbol, and a coefficient that is a function of a coupling. This carries out the decision already recorded at `stage1_engineering_plan.md:483-496`.
- **`derivation: {max_laurent_depth: 1}` means the following.** It passes PSALTer's `MaxLaurentDepth` option. That option sets how deep PSALTer searches for massless poles in the propagator: `1` finds ordinary quadratic 1/k² poles, the potentially healthy massless modes, while `2` or `3` also searches for quartic and hexic poles, which are always sick but expensive to find (`PSALTer.m:82`; arXiv 2406.09500 `Manuscript.tex:132`). It changes the artifact, so it is part of the fingerprint. It defaults to 1 and is omitted unless a theory needs a deeper search.
- **Torsion theories: does writing `CD` (or a torsionful derivative) create new coupling symbols PSALTer would demand? No, provided the package, not the user, does the expansion.**
  - PSALTer requires the Lagrangian to be "quadratic in the perturbed fields … and linear in the couplings" (`PSALTer.m:80`). That is *linear parameterization*, not one coupling per monomial.
  - Legacy's torsion theory wrote geometric invariants, `(1/kappa^2) R̃ + alpha1 I1 + … + b5 R̃²` (`examples/torsion_gertsenshtein/theory.toml:3-6`). The pipeline expanded R̃ with the post-Riemannian decomposition, which splits the Riemann–Cartan curvature into Riemannian curvature plus torsion terms. That produces many monomials, **all multiplied by the same coupling**, so the result is still linear. That is why it worked in the HPC campaigns.
  - For the new package, a term's operator may therefore be a **composite expression in geometric objects**: the torsionful curvature, the torsion T, the field strength F. The package expands it and linearizes about the vacuum before PSALTer. This is the ported decomposition step at `stage1_engineering_plan.md:510-516`, which I-S1B owns.
  - Two consequences for the input design. The operator whitelist must include those geometric objects. And a coefficient that is a *function* of a coupling, like legacy's `1/kappa^2`, must be refused; Einstein–Hilbert gets its own generic coupling (`stage1_engineering_plan.md:497-500`).
  - R-1's prototype uses the Vector theory only, so the expansion path is **not tested here**. The memo says so in its option table.
- **Safety** comes from §2.3: each operator is parsed held, whitelisted, and must be exactly one expression.
**Option B, live alternative: the separate theory file is Wolfram code, not YAML** (the GUM/FeynRules pattern).

The run file looks almost the same as A′. `params:`, `likelihood:` and `sampler:` stay in it exactly as before; the earlier short snippet only omitted them. **The difference is what the theory file contains and who can read it.**

```yaml
# runs/vector_planck.yaml — identical to A′ except the model line
theory:
  tidalcosmo.SpectatorTheory:
    model: ../theories/VectorTheory.wl        # a path to Wolfram code, not an include
params: {Theta1: {prior: {min: 0, max: 2}}, Theta2: {...}, Theta3: 0.5}
likelihood: {...}
```

```wolfram
(* theories/VectorTheory.wl — Wolfram-language model file *)
DefConstantSymbol[Theta1]; DefConstantSymbol[Theta2]; DefConstantSymbol[Theta3];
DefField[VectorField[-a], PrintAs -> "A", PrintSourceAs -> "j"];
Lagrangian := -(1/2)*Theta1*CD[-a]@VectorField[-b]*CD[a]@VectorField[b] + ...;
```

| | A′ (recommended) | B |
| --- | --- | --- |
| theory file language | YAML | Wolfram Language |
| how the run file uses it | Cobaya's `!defaults` pastes the content in, so Cobaya sees fields and couplings | a path string; the content is opaque to Cobaya |
| couplings known to Cobaya and Python without a kernel | yes, so errors surface at `cobaya-run` time | no, only after `derive` has written a manifest |
| linearity in couplings | enforced by structure (one operator per coupling key) | needs a separate Wolfram-side check |
| syntax a Cobaya user must learn | none new | Wolfram and xAct declaration syntax |
| precedent | SOLikeT run files, Planck/MFLike, CosmoSIS values files | GUM, but only because FeynRules/SARAH model files and libraries already existed; PSALTer has no such format (§2.2) |

- **B's one real strength:** a Wolfram-fluent author writes exactly what they would type in a PSALTer session.
- **Why it is not recommended:** it adds a second language to the user's workflow, which is the opposite of "as smooth as possible to pick up", and it hides the couplings from Cobaya.

**Also presented, not recommended: Option D, named operators from a fixed vocabulary** (no expression strings). It is fully kernel-free to validate, but it is a closed catalogue: the hi_class/EFTCAMB design (§2.2), where every new operator needs a package release.

*(A single monolithic Lagrangian string was considered and dropped: term-level input is strictly better and was already decided at `stage1_engineering_plan.md:483`.)*

### 3.2 D-A — what the user runs, and how the derivation runs underneath

**What the user types: two commands on the same file, mirroring Cobaya's own install-then-run.**

```bash
tidalcosmo derive runs/vector_planck.yaml     # once per theory; needs Wolfram; minutes–hours
cobaya-run      runs/vector_planck.yaml     # as often as you like; no Wolfram needed
```

- **The precedent is Cobaya's own two-step.** Its docs have `cobaya-install MyFile.yaml` followed by `cobaya-run` on the same file (`docs/installation_cosmo.rst:28`; `docs/cosmo_basic_runs.rst:70`). Users already know "prepare with one command, run with another, same YAML".
- **The name `derive`** already belongs to this repo: legacy has `tidal derive`, and the port plan lists `derive` (`repo_reshape.md` §5.1). It becomes `tidal derive` at the rename (`repo_reshape.md:143`). It is one subcommand of the package's single console script, matching the existing `tidal`/`tidalcosmo` entry points (`pyproject.toml:79-83`).
- **Why not piggyback on `cobaya-install`?** Its hooks never receive the component's YAML block (`install.py:178-196, 334-338, 372-374`), so it cannot know which Lagrangian to derive. For the same reason the Theory defines no `is_installed`/`install`; `cobaya-install` then just reports it as built-in (`install.py:320-324`).
- **`derive` also accepts a theory file directly** (`tidalcosmo derive theories/VectorTheory.yaml`), useful before any run file exists. Both inputs give the same fingerprint.
- **The Python API mirrors Cobaya's `install()`/`run()`:** `from tidalcosmo.derive import derive`, then `derive(yaml_or_dict, test=False, force=False, timeout=...)`, returning or raising. The CLI and the Cobaya Theory are two thin callers of this one entry point (`repo_reshape.md:611-612`). The Theory calls only the check half.

**Flags, each with one meaning, taken from the closest precedents:**

| flag | meaning | precedent |
| --- | --- | --- |
| `--test` | check only, **never starts a kernel**: reports which theories are derived and up to date. Exit 0 if all are, 1 if any need deriving | `cobaya-install --test` ("Just check whether components are installed", `install.py:758-762`); `make -q` |
| `--force` | re-derive even if an up-to-date derived spectrum exists | Cobaya's `-f/--force` ("Force re-installation of apparently installed packages") |
| `--timeout SECONDS` | wall-clock cap. The legacy default of 600 s is wrong for hours-long derivations; the memo proposes no cap by default, with progress logged | legacy `tidal derive --timeout` |
| `--dry-run` | print what *would* happen (resolved theory, fingerprint, store location, the internal command) without starting Wolfram | `dvc repro --dry`, `make -n` |

Legacy's `--save-script` has no equivalent, because nothing is generated any more. When the derived spectrum is already current, `derive` says so and stops, as legacy's cache-hit message does.

**What `cobaya-run` says when the derived spectrum is missing.** It raises Cobaya's own `ComponentNotInstalledError`, worded on Cobaya's own messages. Compare CAMB's "Could not find CAMB … To install it, run `cobaya-install camb`" (`camb.py:291-296`) and "The data for this likelihood is not up to date … run `cobaya-install … --upgrade`" (`InstallableLikelihood.py:54-55`):

- *missing:* "No derived spectrum for this theory (fingerprint `3fa9…`). Generate it with `tidalcosmo derive <your input>.yaml` on a machine with Wolfram (minutes to hours), or copy it from a machine that has it."
- *store searched:* the message lists the store directories it looked in.

Cobaya gives a component its options, not the input filename, so the message uses a `<your input>` placeholder, as Cobaya's own messages name the component rather than the file. On the `tidalcosmo` CLI side, the same failures use the repo's `error_with_hint` (CLAUDE.md).

**How the derivation runs underneath (internal flow):**
1. `tidalcosmo derive` reads the YAML with Cobaya's own file loader, so `!defaults` includes expand exactly as `cobaya-run` will expand them.
2. It extracts the theory block and writes it as **WXF data** with `wolframclient`'s kernel-free serializer.
3. It launches `wolframscript -file` on a **committed fixed driver**, as one subprocess.
4. The driver `Import`s the WXF and calls the **committed package**. The package validates the block (§2.3, §2.5), declares couplings and fields in `Global`, evaluates the Lagrangian, runs `ParticleSpectrum`, and exports the Stage-1 WXF with an identity and conventions block.
5. `Enclose` plus `Exit[1]` handle failure, and Python verifies the artifact positively.

Alternatives the memo presents:
- **The same package through a `wolframclient` session.** Not recommended, for the reasons in §2.4, which the prototype measures.
- **Python-templated scripts.** The status quo being replaced; costs in §2.9.

**The fingerprint replaces `derivation_hash`.**
- It is a hash of the *theory-defining* keys: fields, couplings, Lagrangian and derivation options, canonically serialized.
- It also covers the package's schema version.
- Priors, likelihoods and sampler settings are excluded, so changing them never forces a re-derivation.
- The artifact *records*, but does not gate on, the PSALTer commit and the Wolfram and xAct versions.

### 3.3 D-B — how the PSALTer output reaches Cobaya, and the gate

**Terminology.** "The artifact" means **PSALTer's symbolic output for one theory**, exported by our package as a WXF file: the wave-operator coefficient tensors plus spectrum metadata. It is what the per-sample numerical stability check reads. It is produced once by `tidalcosmo derive` and read at the start of every sampling run. Below it is called the **derived spectrum**.

**The user never handles the fingerprint.** It is computed behind the scenes from the theory block, used to find the derived spectrum, and recorded for provenance. It appears only in error messages and in Cobaya's run records.

**Where the derived spectrum lives: standards and verdict.**

- **Not Cobaya's `packages_path`.** Cobaya documents it as holding "external packages … the original code itself, a cosmological dataset" (`docs/installation_cosmo.rst:4`; `docs/input.rst:31`), i.e. things `cobaya-install` can fetch. A user-specific, Wolfram-generated file does not match that meaning. We borrow only two things from it: its override order (input YAML > environment variable > default) and its cluster advice to use scratch, not home (`installation_cosmo.rst:55`).
- **It is *data*, not *cache*.**
  - The XDG spec defines the cache directory as "non-essential (cached) data".
  - platformdirs says "If it can be regenerated, use cache … Cached data can be safely deleted without losing functionality" (`docs/explanation.rst:24-28, 130-137`).
  - A derived spectrum takes hours and a single-seat Wolfram license to regenerate, and **cannot** be regenerated at all on an HPC node or a collaborator's machine without Wolfram. Deleting it loses functionality.
  - HF, astropy and pooch put large files in the cache directory only because they can always re-download them. We have no remote to re-download from.
- **Precedents for expensive, content-keyed, portable derived data:**
  - Snakemake's between-run output cache, keyed on hashed inputs and located by an env var (`docs/executing/caching.rst`);
  - DVC's content-addressed store, `files/md5/ab/cdef…`;
  - sympy2c, PyCosmo's cache of symbolically generated code, keyed by an md5 of the symbolic module;
  - Julia's depot stack: the first entry writable, later entries read-only (`JULIA_DEPOT_PATH`);
  - astropy's `export_download_cache`/`import_download_cache`, for moving files to offline machines.

**Recommended design:**
- **Default store:** `platformdirs.user_data_dir("tidalcosmo")/spectra`, i.e. `~/.local/share/tidalcosmo/spectra`.
- **Lookup order** (Cobaya's override order, plus a read-only search path in the Julia/XDG style):
  1. an optional `spectra_path:` key in the theory block;
  2. an environment variable holding a colon-separated list of stores, the first writable and the rest read-only (e.g. a shared group store on RDS);
  3. the default.
- **Layout:** `<store>/v<schema>/<fp[:2]>/<fp>/`, holding `spectrum.wxf` and `manifest.json`.
- **Manifest contents:**
  - the fingerprint and the canonicalized theory block it was computed from;
  - derivation options and schema version;
  - the generator stack (Wolfram 14.3.0, xAct fingerprint, PSALTer `bb45adb0`, package version);
  - a sha256 of the WXF, the wall time, the host and a timestamp.
- **Writes are atomic:** build the entry in a temporary directory inside the store, then rename it into place. This is numba's approach, and it is safe for many MPI ranks reading at once on cluster filesystems.
- **Reads verify** the sha256 before use, as pooch does. A partial copy fails loudly instead of loading silently.
- **Sharing and HPC.** Each `<fp>/` directory is self-contained, so copying it with rsync works with the project's existing `hpc_shuttle.sh push` flow. Also provide `tidalcosmo spectra export`/`import` bundle commands, following astropy. On CSD3 the env var points at `~/rds/hpc-work/...`, because `/home` "should not be used for I/O from running jobs" (CSD3 user guide).

**Provenance in Cobaya's own records.** The Theory reports the fingerprint as its component `version`. Cobaya then writes it into `[prefix].updated.yaml`, and its resume check refuses to continue a chain against a different derived spectrum (`output.py:576-631`). Whether the store path can be kept out of Cobaya's resume comparison is unverified, and is tested in §4.

**The gate.** A refusal at construction, before any sample. The Theory's `initialize()` does four things, and none can start a kernel:

1. **Recompute the fingerprint** from the theory block in use, which is a cheap hash.
2. **No derived spectrum at that fingerprint** in any store: raise `ComponentNotInstalledError` (Cobaya's idiom, §2.6) with a message naming the fingerprint and the command to run on a Wolfram machine, `tidalcosmo derive run.yaml`. Editing an operator changes the fingerprint, so a stale derived spectrum **can never be loaded**.
3. **The sha256 check fails, or the manifest's conventions ≠ (+,−,−,−), ε₀₁₂₃ = +1, or its schema is unknown:** refuse.
4. **The couplings in the manifest differ from the couplings in the theory block:** refuse. Parameter-level mismatches are already caught by Cobaya's own errors (§3.1).

A **static test** also asserts that the Cobaya-facing modules import no `subprocess`, no `wolframclient.evaluation` and no `os.system`, so a Cobaya run *cannot* start a kernel.

**This gate is what makes the numerical stability check possible.** Once sampling starts, the matching derived spectrum is guaranteed present and intact, and the per-sample check reads it. Where that check hooks in is **I-532's decision**. §2.7's evidence points to a Cobaya external prior.

Alternatives the memo presents:
- **Cobaya's `packages_path`.** Rejected for the documented-meaning reason above.
- **The XDG cache directory.** Rejected: the file cannot be regenerated without Wolfram.
- **A project-local store only, DVC-style.** Live, and supported through `spectra_path:`, but it is not the default, because it would duplicate spectra across projects.
- **User pastes a hash into YAML.** Rejected: the fingerprint is behind-the-scenes plumbing.
- **Derivation as Cobaya's `install()`.** Rejected: `install()` never sees the YAML block.

**M3's output location.** Re-derived specs go into a spectra store *inside the repository* (a committed store directory, set with `spectra_path:` in the example YAMLs), using the same layout and manifests. They never go into `examples/data/`, which holds the frozen legacy oracle M3 compares against (`repo_reshape.md:1011-1016`). Committed stores also let CI and collaborators run examples without Wolfram, as `tests_cosmo/fixtures/psalter/` already does for two WXF files. The licensing position on committing derived outputs follows that fixture's `PROVENANCE.md`.

### 3.4 Convention carriage

**The rule.**
- Each quantity takes the convention of the tool that defines it.
  - Signature and ε: PSALTer. CAMB agrees on signature and defines no ε.
  - Perturbation variables, gauge, time and units: CAMB's definitions (table in §2.1).
  - Parameter names: Cobaya's, via `renames`.
- Every artifact records `conventions`, and every consumer asserts it.
- Nothing is converted at run time, so #513 is kept.

**The memo documents CAMB's variable definitions** with pinned citations, as the normative definitions for the symbolic FRW derivation (R-C, I-S1B, M3).

**Checks that catch a slip.**
- **Spectrum side.** A reference theory with a published healthy/ghost verdict runs the whole pipeline in CI, and a flipped sign flips the verdict. The prototype theory is the massive vector.
- **CAMB side.** A test that our variables match `camb.symbolic`'s `camb_var`/`camb_sub` mapping. Recommended to I-532 and M1a.
- **M3 oracle.** The written mapping must record the (−,+,+,+) → (+,−,−,−) flip and the ε₀₁₂… = −1 → +1 flip, including for the parity-odd oracle specs. Routed.

---

## 4. What still needs running, and why reading cannot settle it

**Wolfram lane.** Strictly serial, about 10 minutes expected against a 30-minute cap.

1. **Reference.** Run `scripts/psalter/vector_smoke.wls` unchanged. It is itself a standalone script, so this is also R-1's "script shape". Its wave operator is the baseline.
2. **Recommended design, end to end.** A theory file (Option A′, `lagrangian:` keyed by coupling) is included from a run file, loaded with Cobaya's file loader, converted to WXF, then run through the driver and package to the Stage-1 WXF. Its wave operator must `SameQ` the reference. This also confirms the §2.5 context analysis empirically.
3. **Failure and safety cases.** Each must exit non-zero, leave no artifact, and print a clear message. All of them fail in validation, before `ParticleSpectrum`:
   - an underscore coupling (`m_phi`);
   - an undeclared symbol inside an operator;
   - a `lagrangian:` key that is not a valid coupling name (e.g. a bare number);
   - a duplicated coupling key, which Cobaya's loader must reject;
   - a coupling appearing inside an operator, so the action is not linear in the couplings;
   - an operator string holding two expressions;
   - a zero Lagrangian;
   - a theory name colliding with a PSALTer export;
   - **an injection string** such as `Print["INJECTED"]` or `Quit[]` inside an operator, which must be refused and **never execute**.
4. **Session run** (approved). The same package through a `WolframLanguageSession`, completing R-1's three-shape criterion.
   - Measures what `Quit[]` does to a session (simulating PSALTer's bare `Quit[]`).
   - Measures whether the Python-side timeout is ignored.
   - Self-guards the single-kernel lane with `pgrep`, because the hook cannot see `wolframclient`.

**No kernel.**

5. **Cobaya gate prototype.**
   - **Core, ≤40 lines** (R-1's criterion 3): `cobaya.likelihoods.one`, no CAMB. It shows refusal on a missing derived spectrum, refusal on a fingerprint mismatch, then a pass. It uses a store in a temporary directory with the manifest layout from §3.3.
   - **Extended cases** in a separate script:
     - edited operator → refuses, because the fingerprint changed;
     - corrupted `spectrum.wxf` → refuses on sha256;
     - wrong conventions in the manifest → refuses;
     - a coupling missing from `params:` → Cobaya's own "Requirement … not provided" error;
     - an extra parameter → Cobaya's own "Could not find anything to use" error;
     - a changed prior or a fixed coupling only → passes, with no re-derivation;
     - two read-only stores on the search path → found in the second;
     - fingerprint reported as the component `version` → a resume against a different spectrum is refused, and moving the store path alone does not block a resume (or the memo records that it does).
   - **Static no-kernel import test.**
6. **Option A′ include check, with real files** (the in-memory tests used Planck files as stand-ins). Using `theories/VectorTheory.yaml` and `runs/vector_one.yaml` with `model: !defaults ../theories/VectorTheory` plus a sibling option, confirm:
   - `cobaya-run --test`'s loader expands the include;
   - `tidalcosmo derive`'s loader sees identical values;
   - the fingerprint is identical to the inline form;
   - `.input.yaml` holds the expanded theory;
   - a layered pair of theory files merges term by term (the mapping form);
   - a malformed theory file gives our readable error, not a raw `AttributeError`.

   Any failure is recorded in the memo, and A′ falls back to an explicit `theory_file:` path key.
7. **CosmoPower check** (as dispatched). Install SOLikeT's CosmoPower in a throwaway venv and run it with a missing network file, to confirm the error behavior read from source.

---

## 5. Execution steps

1. **Worktree, planning record, and PR. Do this first, so no planning evidence can be lost.**
   - `git worktree add /tmp/tidal-r1 -b cosmo/r1-interfaces feat/cosmology-program`.
   - **Archive this plan file verbatim** as `docs/cosmology/r1_planning_record.md`, including Appendices A–C. Follow the `planning_session_record.md` precedent: an `ARCHIVE — verbatim record of the 2026-09-13…15 R-1 planning session plan file` header, a note that it is **not authoritative**, and pointers to the living documents that supersede it (`interfaces_decision.md`, `conventions.md`, `COSMOLOGY_PROGRAM.md`). Nothing is paraphrased; the header is the only addition.
   - Commit it together with the `scripts/research/interfaces/README.md` skeleton.
   - Run the hygiene test on the archive, since it quotes paths.
   - Open a **draft PR into `feat/cosmology-program`** immediately.
   - **Update memory in the same sitting:**
     - the project status record (`project_cosmology_program.md`: R-1 planned, the design direction and where the record lives);
     - **feedback memories** for the standing preferences in Appendix A: research fully during planning; no attachment to legacy; reason from comparable packages' conventions; never a subtly different convention; one familiar format; plain explanations; record the full reasoning trail;
     - then `bash .devcontainer/scripts/sync-claude-memory.sh backup`.
2. **Wolfram lane.**
   - First `bash scripts/psalter/ensure_registered.sh`, then `bash scripts/verify-wolfram-setup.sh --require-psalter`. Both outputs are captured and both must exit 0.
   - Then run items 1–4, strictly serial.
   - Every launch sets `QT_QPA_PLATFORM=offscreen` and uses a throwaway working directory. No run is judged by exit status alone.
   - Run `verify --require-psalter` again at the end.
3. **In parallel.** Items 5–7, and drafting the memo from §2–§3.
4. **Memo and conventions reference.**
   - Write `docs/cosmology/conventions.md` first. It holds the research the memo cites.
   - Write `docs/cosmology/interfaces_decision.md`, ≤~500 lines, tables over prose.
   - File the GitHub issues (§6 item 3).
   - Include a "tested here?" column on every option table.
   - Every scope cut is a row with its pinned revision and what a run would have added.
5. **Verify (§7), push, and report back.**

## 6. Deliverables

0. **`docs/cosmology/r1_planning_record.md`** (new; owned path added at your instruction). The verbatim archive of this plan file:
   - the full findings with sources (§2);
   - the recommendations and alternatives (§3);
   - the chronological revision log of every correction and question, and what research and decision followed (Appendix A);
   - the pinned evidence ledger from every research pass (Appendix B);
   - the list of what was not verified (Appendix C).

   The memo and `conventions.md` cite it for evidence. Future sessions read it instead of re-researching or re-proposing rejected designs. **The orchestrator adds its row to `docs/README.md`** (not an owned file), next to `planning_session_record.md`.
1. **`docs/cosmology/interfaces_decision.md`:**
   - findings with sources;
   - input Options A′ (recommended) and inline A, B and D, with worked YAML, the A′-vs-B comparison table, and the precedents for separating theory files from run files;
   - the user-facing `derive` command, its flags, the Python API, and the `cobaya-run` refusal wording, each with its precedent;
   - the D-A and D-B recommendations and alternatives;
   - the trade-off table: user learning cost, staleness safety, provenance, testability, what M3 writes and where, implications for I-S1B's exporter, tested here?;
   - the convention-carriage rule and its checks, summarized with a pointer to `conventions.md`;
   - prototype results and transcripts.
2. **`docs/cosmology/conventions.md`** (new; owned path added at your instruction). The single durable reference, written so R-C, I-S1B, I-532, M1a and M3 can cite it directly, with every entry pinned to a source line:
   - **ownership rule:** which tool defines which convention;
   - **signature and ε:** PSALTer's, with sources, and CAMB's agreement and absence of ε;
   - **CAMB perturbation variables:** the full definition table (§2.1), pinned to CAMB 2.0.4 / `a6de8cc5`, including the sign hazards;
   - **PSALTer input requirements:** quadratic in fields, linear in couplings, one coupling per term, no coupling functions like `1/kappa^2`, plus the symbol-context rules (§2.5);
   - **legacy and Ma–Bertschinger comparisons:** the three sources of per-term sign factors, and the rule that M3's written mapping attributes every difference to convention or physics;
   - **open items:** CAMB's authoritative notes read only through `pdftotext`; the per-term flip table left to M3; no β or polarization-handedness convention yet for O4.
3. **GitHub issues** for each discovered defect, following CLAUDE.md's proactive-issue rule and deduplicated first with `gh issue list -S`:
   - the `spectrum_design.md` convention errors;
   - the M3 mapping convention requirement;
   - the missing β and handedness convention;
   - PSALTer's bare `Quit[]`;
   - **the per-sample stability-check hook** (for I-532), with the TorC and polology evidence from §2.7, the external-prior wrinkle, and the open questions: how the prior locates the derived spectrum, how the fingerprint gate applies to it, and hard vs. soft severity. This is so the question is remembered when that work starts.

   Issue numbers go in the report.
4. **`scripts/research/interfaces/`:**
   - README with a status line and a contents table with a provenance column;
   - prototype package, driver and theory YAMLs;
   - failure and injection fixtures;
   - comparison script;
   - Cobaya gate prototype;
   - fetch routes (payloads stay in gitignored `third_party/`).
5. **Draft PR** into `feat/cosmology-program` with green CI.
6. **Report back** in protocol format: branch · PR · `CI <run-id>: <conclusion>` · each success criterion with its output · the recommendations · routed findings (§8) · what could not be tested.

**Nothing** is written under `tidalcosmo/`, `pyproject.toml`, `uv.lock`, `.venv`, `scripts/psalter/`, other design docs, or the Wolfram userbase. No merge, bump, tag or changelog edit.

## 7. Verification before each push

1. `git add` first, then `uv run pytest tests/test_repo_hygiene.py`. It scans tracked files for the container clone path, a user home directory and the Claude project slug, so transcripts are scrubbed at print time.
2. `uv run ruff check` on changed Python.
3. Local `cspell` as a pre-check. New terms go in an in-file `cspell:words` directive, since `cspell.json` is not owned.
4. `git diff --stat` shows owned paths only.
5. `gh run list` for the head commit gives `CI <run-id>: <conclusion>`.

## 8. Findings routed to the orchestrator (reported, not edited)

**`spectrum_design.md`:**
- `:301` labels CAMB (−,+,+,+) with a Ma–Bertschinger citation, and `:354` repeats it.
- `:312-313`: "sign of det g" is wrong (det g = −1 in both signatures; the real mechanism is (−1)ⁿ from raised indices, as in `docs/tex/pgt_stability_priors.tex:136-139`).
- The two-signature premise of §4.3/§4.5 collapses to one.
- `:392-394` models the generator on legacy.

**`stage1_engineering_plan.md`:** §4.3/§4.4 and gate rows `:763` and `:803-805` assume text emission.

**M3 oracle mapping:** it must attribute every difference to convention or physics, recording the signature, field-definition and ε sign factors, including for the parity-odd oracle specs (`tests_cosmo/data/oracles/README.md:80-110`). The rule is stated in `conventions.md`.

**CAMB variable definitions:** documented nowhere in the repo until now. They now live in `conventions.md`, for the orchestrator to make canonical.

**Dependency (for this wave's `pyproject.toml` owner):** `platformdirs` is **not installed** and not declared (checked: import fails; absent from `pyproject.toml` and `uv.lock`). The storage design needs it, or an equivalent ~10-line XDG lookup. R-1's prototypes use the XDG environment variables directly, so no install is needed here.

**I-S1B (torsion input):** operators must be able to reference geometric objects that the package expands before PSALTer. Coefficients that are functions of couplings must be refused.

**O4 owners:** no convention for the sign of β or for polarization handedness.

**I-532:**
- The stability-check hook evidence from §2.7. No group precedent exists for a per-sample spectrum check inside Cobaya: TorC fixed couplings offline, and polology reweights offline.
- The external-prior wrinkle: a prior sees only parameter values.
- The CAMB extension constraint.
- All of this is also filed as an issue.

**Legacy `tidal derive` flags:** `--save-script` and `--dry-run` (print the generated script) have no meaning once nothing is generated. The 600 s `--timeout` default is wrong for hours-long derivations.

**PSALTer behaviors, for I-S1B:** a bare `Quit[]` at `ConjectureInverse.m:37-39`, uncaught throws, `Message` silenced on subkernels, and published `.wxf` files that come from an unpublished build.

**Stale items:**
- `R-1.md:38-40` cites `tidalcosmo/derive/README.md:12` for a claim that file does not make.
- `vector_smoke.wls:25-35` appears to carry a caveat that expired with #543.
- `scientific_review.md` #28 ("CAMB not installed") is stale.

## 9. Risks

1. **The recommended design's wave operator differs from the reference.** A three-valued in-kernel comparison (proved equal / proved different / could not decide). If it differs, stop, report, and do not predict an outcome.
2. **The held-parse whitelist rejects a legitimate xAct construct**, for example an unanticipated head from `CD` or `@`. The prototype expands the whitelist from the actual parse of the reference Lagrangian and records it.
3. **YAML folding changes the expression.** `>-` joins lines with spaces. Test with the multi-line reference Lagrangian.
4. **A transcript trips the hygiene test.** Scrub at print time and test after `git add`.
5. **A comparable install damages the repo venv.** Throwaway venvs only.
6. **The memo exceeds 500 lines.** Use tables, and reference full logs by path.
7. **`ensure_registered.sh` fails.** It is the first step; its output is attached; exit 2 from the verifier means degraded, not failed.
8. **Planning evidence is lost before it is written down.** Everything below lives only in this plan file and the session transcript. → §5 step 1 archives this file, including Appendices A–C, as the **first commit**, before any prototype work.

---

## Appendix A — How the plan changed during planning, and why

Chronological. Each entry gives what was proposed, the user's correction or question, what research followed, and the resulting decision. This is the reasoning trail future sessions need, so they do not re-propose rejected designs.

| # | Proposed | User feedback | Research done | Outcome |
| --- | --- | --- | --- | --- |
| 1 | A bake-off: build three Wolfram "shapes" and pick whichever matches; a JSON input file; CAMB's signature left as a hypothesis to fetch later | "research can be done in planning mode"; "why depend on downloading things?"; explain "shapes" and "one input, three carriages" | Read CAMB's installed source and literature | CAMB's declared signature is (+,−,−,−), the same as PSALTer's (§2.1). Research moves *before* implementation |
| 2 | Keep the legacy generator as the baseline | "we must avoid being modeled on legacy… start afresh"; "committed package is what we need" | PSALTer source hazards; the templating cost numbers | A committed Wolfram package. Legacy is cited only as a failure record |
| 3 | ε orientation left open | "normalize against PSALTer; CAMB likely has no ε" | Grep of the whole CAMB tree at the pinned tag | Confirmed: CAMB has no ε. Adopt PSALTer's ε₀₁₂₃ = +1 |
| 4 | CAMB signature checked in docs only | "explore the code too, to understand formulas and conventions" | CAMB Fortran at 2.0.4, CAMB notes, the Challinor/Lewis papers | The Fortran is signature-agnostic (Ma–Bertschinger's equations in invariant variables). The docs match the code |
| 5 | "Test all three shapes" | "shouldn't conventions from similar packages already inform the decision?" | Wolfram ecosystem (FeynRules, SARAH, PSALTer materials, xPand, HiGGS, Hamilcar, GUM, FlexibleSUSY), `wolframclient` issues | Reason from precedent: run as a subprocess, no session. Prototypes confirm the choice; they don't decide it. The session run is kept, approved as measured evidence |
| 6 | Stability gate stated abstractly | "gate sampling on the PSALTer output"; "how do TorC and polology hook stability in?" | Cobaya rejection routes; TorC TeX and audit; polology code | Refuse at construction via `ComponentNotInstalledError`. The per-sample hook is I-532's, recorded as an issue |
| 7 | The YAML names a Wolfram model file (GUM pattern) | "avoid a new file format; easy to pick up" | Cobaya `lambda`-string precedent; safe held parsing; YAML traps | The Lagrangian is written in YAML |
| 8 | Dropped the model file on the user's steer | "I didn't ask to drop it — evaluate the standard" | GUM paper motivation; `stage1_engineering_plan.md:483-496`; Cobaya `!defaults` | A Wolfram model file is standard only where a model-file ecosystem exists. PSALTer has none, so it stays a live alternative (Option B). Term-level input |
| 9 | One monolithic Lagrangian string | "discount this option, A far better" | — | Dropped |
| 10 | "Couplings declared once, next to their term", with optional labels beside terms | "is this how Cobaya users define parameters? A subtly different convention would be very bad" | Cobaya docs, source, in-memory `get_model` experiments; LAT_MFLike; Planck | Priors and labels go only in top-level `params:`; the MFLike `get_modified_defaults` mechanism. Found: the naive route gives a silent `nan` |
| 11 | Store derived spectra in `packages_path` | "look for standards on where reusable generated data is stored" | Cobaya docs, XDG, platformdirs, pooch, HF, astropy, numba, sympy2c, DVC, Snakemake, Julia depot, CSD3 | Not `packages_path` (documented as installable packages). Data, not cache. A content-keyed user data store with an env-var search path |
| 12 | Inline theory block, with A′ as an aside | "separate files seem best, may declutter big thesis Lagrangians" | Framework precedents (SOLikeT, CosmoSIS, GAMBIT, MadGraph, Hydra); `!defaults` tests | **A′ recommended.** `lagrangian:` becomes a mapping keyed by coupling, because lists are replaced on merge |
| 13 | "`derive` launches `wolframscript -file …`" | "does the user type this? follow package conventions" | Cobaya CLI family; GUM, FlexibleSUSY, CONNECT, CosmoPower, DVC, make | The user types `tidalcosmo derive run.yaml`, then `cobaya-run run.yaml`, mirroring Cobaya's install-then-run. Flags from precedents |
| 14 | "CAMB's definition wins" | "wins over what? competing standards?" | Competing definitions enumerated | Use `camb.symbolic` variables; convert anything imported once, with a citation |
| 15 | Conventions noted in the memo only | "document all this; where should conventions live?" | — | New `docs/cosmology/conventions.md` (user chose) |
| 16 | M3 mapping gap | "of course — compare physics, not bytes; explain every difference" | Oracle README drift classes | The mapping attributes each difference to convention or physics, with three sign-factor sources |
| 17 | Torsion concern: does `CD` force new couplings? | Asked | PSALTer usage text; legacy torsion TOML | Linear *parameterization*, not one coupling per monomial. The package expands geometric objects. Coupling functions are refused |
| 18 | This plan's evidence lived only in the session | "record all investigations, evidence, decisions and feedback" | Precedent: `docs/cosmology/planning_session_record.md` | Archive this plan verbatim as `docs/cosmology/r1_planning_record.md` (Appendices A–C included) as the first commit |

**Standing preferences this session showed.** These are saved to memory in §5.
- Do all read-only research in planning, before proposing.
- Feel no attachment to legacy design.
- Reason from established conventions of comparable packages, not bake-offs.
- Never introduce a convention subtly different from what users already expect.
- Keep the user's workflow to one familiar format.
- Explain terms plainly.
- Record the full reasoning trail for future sessions.

---

## Appendix B — Evidence ledger

Pinned sources for claims not fully cited in the body. Local paths are repo-relative. `$USERBASE` means Wolfram's `$UserBaseDirectory`. The PSALTer install is `$USERBASE/Applications/xAct/PSALTer` at `bb45adb0`.

### B1. CAMB (cmbant/CAMB tag 2.0.4 = `a6de8cc59124c8bbe3924f1c546f96cef73dcabb`; the installed Python files are byte-identical)

**Coded Einstein equations** (`fortran/equations.f90`):
- `:2365` `adotoa = sqrt(grho/3)`
- `:2427` `z = (0.5_dl*dgrho/k + etak)/adotoa`
- `:2430` `sigma = (z + 1.5_dl*dgq/k2)`
- `:2431` `ayprime(ix_etak) = 0.5_dl*dgq`
- `:2442` `clxcdot = -k*z`
- `:2830` φ, `:2876` φ̇, `:2879` σ̇, `:3279,3284` tensors
- These translate term for term to Ma–Bertschinger via `etak = kη_s`, `h_s' = 2kz`, `dgrho = κa²δρ`.

**Definitions and conventions:**
- `:38-39` "Equations and relation to synchronous gauge variables documented in the notes: https://cosmologist.info/notes/CAMB.pdf"
- `:1119-1125` `grho = a^2 kappa rho`, `dgq`, `dgpi`; `constants.f90:49` `kappa = 8·π·G`
- CAMB notes (dated 2026-07-17): "Using the uₐuᵃ = 1 signature"; `etak = kη_s`, with η_s = −η/2; "parity symmetric ensemble, so ClTB = ClEB = 0". Read via pdftotext, so the maths is garbled and only unambiguous text was quoted.
- `docs/source/variables_guide.rst:133` ds² = a²(τ)[dτ² − (δᵢⱼ+hᵢⱼ)dxⁱdxʲ]; `:223`, `:252`
- `docs/ScalEqs.ipynb:278-285` defines the metric sign choice and gives Hu et al.'s opposite Φ_N sign.

**Formalism papers:**
- Challinor & Lasenby, astro-ph/9804301 TeX: `:308` "(+---) metric signature"; `:511` η₀₁₂₃ = −√(−g)
- astro-ph/9911481 `:157` (+---)
- astro-ph/0406096 `:1105-1108` "the signature where uₐuᵃ = 1"
- astro-ph/0203507 `:24761`

**No ε convention in code.** Grepping all `.f90`/`.py` at the tag for `epsilon|levi|parity|birefring|EB|TB|chiral|antisym` finds only machine epsilon, `do_parity_odd` (lensing bispectrum slices), pseudo-Cℓ mask matrices and Bessel parity relations.

**Extension points:**
- `DarkEnergyInterface.f90:9-25` (`num_perturb_equations`, `PerturbationEvolve`); `equations.f90:621-622` slot reservation; `DarkEnergyQuintessence.f90:256-257, 273-274` example ODE.
- Python cannot add physics: `baseconfig.py:858-879`; `dark_energy.py:198` "need to define a new derived class in Fortran".
- Limits: one DarkEnergy object replacing Λ; back-reacts via `dgrho`/`dgq` (`equations.f90:2421-2422`); scalar-only (`derivst` has no hook).
- `set_custom_scalar_sources` is read-only over existing variables (`symbolic.py:917-936`).
- No-fork route: `get_background_time_evolution` + `get_time_evolution(q, eta, vars)` (`results.py:557`). The docstring warns potentials "may be numerically unstable far outside the horizon".

**TorC's fork (slegner/CAMB vs upstream master):**
- ahead 18, behind 98, 24 files;
- adds `is_no_mod_w`/`is_no_mod_P` and `P_de`;
- `BackgroundDensityAndPressure(this, grhov, a, grhov_t, w, P)`;
- an early-abort print "tauend larger than all switch error encountered";
- changes background pressure only.

**SymBoltz.jl v1.7.0 `3d1f20a3`:**
- (−,+,+,+) in the conformal Newtonian gauge: `docs/src/conventions.md:7`, `src/models/metric.jl:4,11-12`.
- Equations at `src/models/gravity.jl:38-39` match Ma–Bertschinger.
- Paper typesetting disagrees with the code at TeX `:1075`, `:1094`, and in `extended_models.md` (`3*g.Φ` vs `3*D(g.Φ)`).
- Scalar modes only; conformal Newtonian only (TeX `:1073`).

### B2. PSALTer source (installed `bb45adb0`, v2.0.2)

- **`ValidateLagrangian.m` (39 lines).** Refuses:
  - a zero Lagrangian (`:20`);
  - an `UnknownCoupling` that is not `ConstantSymbolQ` (`:32`);
  - `epsilonG` is exempted (`:35`);
  - an `UnknownField` via `Names["xAct`PSALTer`"<>ToString@#<>"`*"]` (`:37`);
  - `NonQuadraticFields` via `ResourceFunction["PolynomialDegree"]` (`:38`).
- **Never thrown anywhere:** `NonLinearCouplings` (`:6`) and `ParityOdd` (`:10`).
- **Dead code:** `ValidateNeglect` is called at `ParticleSpectrum.m:30` and defined nowhere. `Method` is declared (`:22,116`) and never read (#521).
- **`DefField.m`:**
  - `:37` `FieldContext = "xAct`PSALTer`"<>ToString@InputField<>"`"`;
  - `:50-53` `DefTensor` in the caller's context, and `ToExpression["Source"<>ToString@InputField]`;
  - `Begin[FieldContext]` `:57` … `End[]` `:92`;
  - wrapped in `Catch` (`:32-33`), so `DefField` errors are contained;
  - `SummariseField.m:26,85` read the ambient `Context[]`.
- **`PreComputeComponents.m:32-34`** calls `CloseKernels[]` then `LaunchKernels[$ProcessorCount]`: PSALTer owns the kernel pool.
- **`ParticleSpectrum.m`:**
  - Options (`:22,116`): `TheoryName->False, MaxLaurentDepth->1, Neglect->{}, MasslessSpectrum->True, AspectRatio->Landscape, ShowPropagator->True, Method->"Easy"`.
  - No `Catch` (`:23`).
  - Eight `$Local*` globals (`:74-81`).
  - Constructors use the literal `"xAct`PSALTer`Private`ClassName"` (`:41-66`).
  - The PDF uses `UsingFrontEnd@Export` (`:82`).
  - `UpdateTheoryAssociation` is called twice through `MapThread` (`:87-101`).
- **`UpdateTheoryAssociation.m:8`** `Symbol@Name` resolves at runtime; `:14` `DumpSave[…"ParticleSpectrograph"<>Name<>".mx"]`.
- **`ValidateTheoryName.m:6-8`** throws `WrongTheoryName` on a non-string, so the default `False` throws uncaught.
- **`PSALTer.m`:**
  - `:15-20` redefines `Print`/`Message` to `Null` when `$KernelID≠0`;
  - `:11-12` `Off@Solve::fulldim`, `Off@General::shdw`, never restored;
  - `:26` the `BeginPackage` dependencies;
  - `:49-56` freezes `$WorkingDirectory`;
  - `:80` `ParticleSpectrum::usage` ("quadratic in the perturbed fields … linear in the couplings");
  - `:82` `MaxLaurentDepth::usage`;
  - `:101` the West Coast signature.
- **`DefGeometry.m`:** `:5` `DefManifold[M4,4,IndexRange[{a,z}]]` read in the public context (via `PSALTer.m:160-161`, `ReloadPackage.m:15,19,24`); `:69-72` metric and ε component values.
- **`ConjectureInverse.m:37-39`** has a bare `Quit[]` when the null-space lengths differ.
- **README "Known bugs" 1** (`:305`): gauge symmetries are sporadically not identified, because of runtime RNG; re-running usually fixes it.
- **Supplementary materials** (`wevbarker/SupplementalMaterials-2607@b49e9f1d`):
  - `Models.m:3-18` selects a model by uncommenting a `Get`;
  - `VectorTheory.m` is the whole-model example;
  - `JuliaExport.m:52-56` has regex bugs (`"I" -> "im"`, and a `Sqrt` regex that breaks on nesting).
  - **The published `.wxf` files have more keys than 2.0.2 writes**, with unevaluated private symbols, so they came from an unpublished build. The writer code was not found.
- **Polology JAX code** (same repo):
  - `JAX/src/psalter/sample.py` `sample(wxf_path, …)`;
  - `_extracting/wxf.py:25-26,162-172`: `binary_deserialize`, theory identity from the filename stem, contexts stripped, couplings sorted by name, no provenance;
  - `sampling.py:106-110`: uniform prior on the sphere;
  - `sample.py:169-195`: healthy filter;
  - `_measuring/likelihood.py:22-23`: classifier;
  - `Julia/src/run/reweighter.jl:4-6`: post-hoc reweighting.

### B3. Wolfram ecosystem and `wolframclient`

**FeynRules** (arXiv 1310.1921):
- model file syntax at `:3427-3582`; `LoadModel` at `:2678-2688`;
- "the model must be reloaded … kernel must be quit" (`:2688`);
- restrictions are irreversible (`:1976`);
- non-compliant vertices silently discarded (`:3779`).

**UFO:**
- 1108.2040 `UFO_Paper.tex:227-235` explains why a generic format was needed.
- UFO 2.0, 2304.09883 `1_intro.tex:1,3`; `2_general.tex:65-70` optional version signature.

**SARAH** (1309.7223): model at `SARAH4.tex:1374-1437`; `Start`/`MakeSPheno` at `:686-709`; tadpoles solved on the Mathematica side (`:995`); first gauge-invariant contraction taken by default (`:478`).

**xPand** (1302.6174) `:1340-1359, 1407-1411`. **HiGGS** (2206.00658) `:677-680, 819, 963-973`. **Hamilcar** (2512.25007) listings.

**GUM** (`GambitBSM/gambit_2.0@a4742ac9`, arXiv 2107.00030):
- `:2795` loads the requested FeynRules/SARAH model over WSTP;
- `:3447` points to the FeynRules model database;
- `gum/Tutorial/MDMSM.gum`;
- `math_package.cpp:66-75`: one WSTP link per run.

**FlexibleSUSY** `@913bbbf1`: `templates/module.mk.in:312-318` (subprocess, stamp file, "Error: The code generation failed!"); `meta/Utils.m:344` `Quit[1]`.

**cymetric** uses `WolframLanguageSession` per call and notes a consumer bug in 1.1.6.

**Wolfram docs:**
- `Exit[n]` passes exit codes;
- `Enclose` catches `Confirm*`;
- `Check` ignores messages switched off by `Off`/`Quit`;
- `wolframscript` falls back to cloud evaluation if no local kernel is found;
- `Throw::nocatch` exits 0 (probed in #561).

**`wolframclient` 1.4.0 (PyPI 2024-06-04; commits to 2026-03-02):**
- #47: release planned with Mathematica 15;
- #41/#42 fixed on master only;
- `kernelcontroller.py:105-132, 391-413` (`STARTUP_TIMEOUT` 20 s, `TERMINATE_TIMEOUT` 3 s);
- `localsession.py:33-34` docstring vs `:259-260` implementation (the timeout is not passed).

**Issues:**
- #38 (open) timeout ignored;
- #29 no interrupt from Python, use `TimeConstrained`;
- #34 busy kernel cannot be aborted;
- #30 hang on kernel death (fixed 1.1.7);
- #17 (open) socket startup failures on Ubuntu;
- #26 (open) hangs from missing front-end libraries;
- #33 (open) socket failures under MPI;
- #31 pool launch failures;
- #28 kernel not found (permissions);
- #37 parser version error on 13.3;
- #39 no entitlement licenses.

**Kernel discovery** (`utils/environment.py:62-76`) searches only `/usr/local/Wolfram/{Desktop,Mathematica,WolframEngine}`, plus `WOLFRAM_INSTALLATION_DIRECTORY`.

**Footprint (checked this session):** importing `serializers.export` and `deserializers.binary_deserialize` loads no `wolframclient.evaluation`, `zmq` or `pkg_resources`. `Fraction(-1,2)` round-trips to `Rational[-1, 2]`.

### B4. Safe parsing and YAML

**Wolfram docs:**
- `ToExpression[input, form, h]` wraps the head before evaluating;
- `HoldComplete` blocks upvalues, while `Hold` does not;
- `SyntaxQ`/`SyntaxLength`;
- `Cases` evaluates the right-hand side for matches, so wrap them in `HoldComplete`;
- symbols are created on parse (`$NewSymbol`);
- multi-line input may split into several expressions;
- `x_y` is `Pattern[x, Blank[y]]`;
- `CodeParser` `CodeParse` exists from 12.2.

**Precedent:** Mathematica StackExchange q172553 (the whitelist answer by Theo Tiger, id 172631).

**Cobaya:** `tools.py:359-373` `eval` of `lambda` strings, where `scope = globals()` includes `os`; `replace_optimizations` regex at `:321-341`. `yaml.py:46` `ScientificLoader(yaml.Loader)`; tested: `!!python/name:builtins.len` loads the function.

**PyYAML 6.0.3 test results:**
- **Parsed as intended:** plain, quoted, `|`, `|-` and `>`, including leading `+`/`-` continuation lines.
- **Errors:** a plain value starting with `@`, `` ` ``, `%`, `[`, `'a'*b`, `- x`, `*x`, `!x` or `|x`; `{L: CD[-a]@V[-b]}` inside a flow mapping.
- **Silently wrong:**
  - `a*b #c` drops the comment;
  - `{a}` becomes a mapping and `&x` becomes `None`;
  - `${HOME}` is expanded (Cobaya, plain scalars only);
  - `1e2` becomes `100.0`;
  - `yes` becomes `True`.
- **`\[Theta]`:** fine except in double quotes.

**YAML 1.2.2:** `@` and `` ` `` are reserved; ": " and " #" are forbidden in plain scalars.

**Other tools:**
- Hi-COLA `Frontend/read_parameters.py:67-69, 216, 388` (sympify/eval, no validation);
- the SymPy `parse_expr` eval warning (`sympy_parser.py:919-921`);
- ALOHA regex, then PLY, then eval (`create_aloha.py:222, 263`);
- Cadabra has a LaTeX pre-processor.

### B5. Cobaya 3.6.2 mechanics (upstream `b76b6fed`; docs at v3.6.2 = `899f30a`)

**Rejection routes:**
- `prior.py:131-139` (string eval), `:155-165` (input parameters only), `:551-555` (only named arguments passed), `:712-719, 765-772`;
- `model.py:607-608, 650-667`: a `-inf` prior skips likelihoods;
- `theory.py:132, 282-297` and `model.py:447-452`: `calculate` False, or an exception, sets `-inf`;
- `camb.py:101-104, 747-765, 1019-1047`: CAMB errors mean zero likelihood;
- `mcmc.py:713-714`: prior rejections counted separately;
- component ordering `model.py:789-810, 1338` (dependency-free components in list order; untested whether a check runs before CAMB).

**Install vs run:**
- `ComponentNotInstalledError` defined at `component.py:824`, raised at `InstallableLikelihood.py:36-57`, `DataSetLikelihood.py:48-56`, `camb.py:290-296`;
- `run.py` has no install path; flags `run.py:198-259`;
- `install.py:117-482` flow, with `is_installed`/`install` receiving only path and flags (`:232-243, 374`); options `:751-819`;
- SOLikeT self-installs at runtime (`soliket/lensing/lensing.py:111-113`), the counterexample.

**Parameters:**
- docs: `likelihoods.rst:19, 131`, `input.rst:18`, `DEVEL.rst:142-162`, `prior.py:399-406`, `cosmo_basic_runs.rst:65-67`, `theories_and_dependencies.rst:124-133`;
- routing: `model.py:1155-1224, 1308-1318, 852-857, 913-918, 1001-1005`;
- hooks: `theory.py:208-225`, `component.py:241-242, 333-338`, `input.py:279-281, 361-365, 436-454, 506-528`;
- MFLike `foreground.py:215-225`; Planck `planck_clik.py:123-135`; plik `TTTEEE.yaml:14`;
- classy is agnostic (`boltzmannbase.py:41-42`);
- SOLikeT's `get_can_support_parameters` is misnamed and a no-op (`cosmopower.py:286, 372`).

**In-memory experiments this session:**

| Mechanism | Coupling missing | Extra parameter |
| --- | --- | --- |
| `get_can_support_params` only | silent `nan` | "Could not find anything to use …" |
| `get_requirements` or `get_modified_defaults` | "Requirement … not provided" | "Could not find anything to use …" |

Also: a block-level `params:` silently dropped a coupling; duplicate terms made the last default win; fixed values and overrides behaved normally.

**`!defaults`:**
- `yaml.py:71-101, 104-113, 116, 121, 150-151`;
- tests A–L: sub-key with siblings OK; list merge left to right; scalar-then-dict crash (issue #463); lists replaced; duplicate key rejected; `<<` merge unsupported; `${YAML_ROOT}` wrong after include; no string or dict input; missing-file message; absolute paths OK; top-level list crashes; class default replaced wholesale;
- output expansion `output.py:642-653`; resume `:574-587, 592-631`; MPI root-only `input.py:143-145`;
- grids: `gridconfig.py:192-235, 319`, `tests/test_cosmo_grid.yaml`;
- issue #56 (no package-relative include); zip-safe caveat `component.py:279, 299`.

### B6. Cosmology comparables (SHAs, last commits, licenses)

- **SOLikeT** `f6e00a9` (2026-09-01, MIT): `cosmopower.py:49-62, 120, 128-133, 147, 286-300`; `cosmopower_NN.py:218-222, 321-327`.
- **cosmopower_cobaya** `simonsobs@e5a06c5` (2024-08-29, MIT). **cosmopower** `7cac5e7` (2024-12-21, GPL-3).
- **hi_class** `emiliobellini@4f0aad6` (2026-05-22, no LICENSE): `include/background.h:20-27`, `gravity_models_smg.c:48-468, 575`, `background_smg.c:1238-1243`, `hiclassy.pyx:587-589`.
- **mochi_class** `2b0b16e` (2025-11-12, MIT): README `:13`; `designer_fr.ini:858-884`; `gravity_models_smg.c:116, 149-188, 208, 255-260, 895, 1017`.
- **EFTCAMB** `16d9c4e` (2026-04-02, GPL-3 core): `09_EFTCAMB_main.f90:455-547`; `008p0_Horndeski.f90:857-885`; `HorndeskiExample.yaml:14-30`; `params_EFT.ini:7-86`.
- **CLASS** `64bbab7` (no LICENSE): `explanatory.ini` flat keys; `input.c:140-143, 185-195` (no include).
- **DISCO-EB** `3aa88d2` (2025-10-15, GPL-3). **Hi-COLA** `97f38dd` (no license).
- **TorC:**
  - `paper_Qtorsion.tex:168, 180-188, 395-419, 632, 656`;
  - audit `docs/cosmology/torc_pipeline_audit.md:431-441, 562-571, 743`;
  - `slegner/cobaya@414a2e89:camb.yaml:25-26`, `camb.py:832-834, 920-922`;
  - `slegner/CAMB@2fb908af:DarkEnergyPressure.f90:34, 149, 152`.
- **Local literature already present:** TorC, PSALTer ×3, CAMB2, CLASS II, nanoCMB, ABCMB, SymBoltz, DISCO-DJ, class_rot, anesthetic.
- **Not present:** hi_class, EFTCAMB, cosmopower, SOLikeT, and the Cobaya paper.

### B7. Storage standards

- **Cobaya:** `installation_cosmo.rst:4, 55, 73-92, 120`; `input.rst:31, 103, 107`; `conventions.py:122-123`; `tools.py:1062-1071, 1100-1103, 1188-1207`; `autoselect_covmat.py:17, 54-63`; `planck_2018_CamSpec_python.py:194-208`; `output.rst:29-30`.
- **XDG** basedir spec 0.8. **platformdirs** `c5ef1edb` `docs/explanation.rst:24-40, 130-137`, `unix.py:42, 60`. **CACHEDIR.TAG** spec.
- **Precedent details:**
  - joblib `memory.py:886-890, 983-987, 1046`;
  - pooch `utils.py:77-103`, `core.py:310-313`;
  - astropy `data.rst:41, 219-231`, `ASTROPY_CACHE_DIR`;
  - HF `manage-cache.md:27, 75`;
  - numba `envvars.rst:466-472`, `caching.rst:77-87`;
  - JAX `compilation_cache.py:207-210`;
  - sympy2c 2.5.1 `utils.py:114-138`, `compiler.py:436-453`;
  - Julia `DEPOT_PATH`;
  - Wolfram `$UserBaseDirectory`/`LocalObject`/`$CacheBaseDirectory`;
  - MadGraph `import_ufo.py:69-100, 357-371`;
  - DVC `internal-files.md:75`;
  - Snakemake `docs/executing/caching.rst`.
- **HPC:** CSD3 io_management (50 GB home quota; no job I/O from home; RDS not backed up); `docs/hpc_workflow.md:43-45`; `scripts/hpc_shuttle.sh:53-76`.
- **`platformdirs` is not installed or declared** (checked this session).

### B8. CLI precedents

- **Cobaya:** `cobaya-3.6.2.dist-info/entry_points.txt` (no `cobaya-doctor`); `install.py:125, 320-324, 366-369, 722, 758-762`; `run.py:236-238`; `__main__.py` fallback to run; grid `--dryrun`.
- **GUM:** `gum.py:44-45, 697-701, 997-1009`; Tutorial README `:84, 96`.
- **FlexibleSUSY:** README `:88-90`. **CONNECT:** README `:88, 92`. **CosmoPower:** training notebooks only. **SymBoltz:** `getting_started.md:48, 58-59`.
- **DVC:** `status`, `repro --dry`. **Snakemake:** `-n`. **make:** `-n`/`-q` semantics.
- **Repo:** `pyproject.toml:79-83`; `tidal/cli/__init__.py:71-127` legacy `derive` flags; `tidalcosmo/cli/__init__.py`, `README.md`; `_console.py::error_with_hint`; `scripts/oracles/freeze_legacy_oracle.py:798-810` `--check`/`--staleness`; `repo_reshape.md:143, 611-612`.

### B9. Legacy record

- `tidal/cli/_derive.py`: 6,756 lines, 53 of 69 top-level functions emit text (5,764 lines, 85.3%), 2,621 string-literal lines.
- Largest emitters `:2205-2870, 4518-5045, 5763-6267, 5060-5549`.
- Escaping hazards `:305, 316, 329, 361, 5797-5799, 6235`; regex rename `:193-265`; `Hold` power workaround `:1024-1030`; bracket check only warns `:6418-6430`; temp script unlinked `:6594-6613`; `pkill -f WolframKernel` `:6469`; exit-0 probe `:6650-6659`; test `tests/test_cli.py:2296-2347`.
- **Zero** generated `.wls` are tracked; `examples/polar_kg/run.sh:3` references a missing `.wls`.
- Failure catalogue: `.claude/rules/wolfram.md` (`:30-106`), `docs/tex/troubleshooting.tex:261-263`; return path `tidal/wolfram/ExportJSON.wl:333, 487, 1638`; `tidal/symbolic/_eval_utils.py:40-83, 119-133, 447-465`.
- `tests/wolfram/test_harness.wl` and `test_gauge_fix.wls:72, 101` show package functions are unit-testable.

### B10. Repo planning facts that constrain R-1

**What each held prompt waits on:**
- R-C: D-A only (`R-C.md:3-6, 20-23`).
- I-532: D-B plus conventions (`I-532.md:3-6, 20-24`).
- I-S1A-core: D-A, D-B, carriage and I-532's flag type (`I-S1A.md:3-6, 20-25`).
- I-S1B: D-A, D-B, D-C, carriage and I-S1A-core (`I-S1B.md:20-26`).
- M1b: D-B, I-532 and I-S1A-core (`M1b.md:3-6, 20-24, 31`).

**Register rows:** `COSMOLOGY_PROGRAM.md:184, 189, 194, 203, 208, 213-217, 219`; wave board `:707-728`; Wave-1 orchestrator commitments `:938-954`; carried-forward `:1041`.

**Design-doc sites a D-A verdict touches:** `spectrum_design.md:392-394`; `stage1_engineering_plan.md:502-538, 606, 763, 771, 803-805`.

**Flag-vs-gate contradiction:** `tidalcosmo/validity/README.md:28-29` vs `stage1_engineering_plan.md:786-791` (I-532's to resolve).

**Parity-odd scope:**
- O4 is on the ladder (`COSMOLOGY_PROGRAM.md:257, 211`; `observable_ladder.md:336-339, 367, 568`);
- ε-odd operators in `examples/torsion_gertsenshtein/theory_parity_odd.toml:3-31`;
- these specs sit in the frozen oracle manifest.

**Oracle:** `repo_reshape.md:787-790, 1011-1033, 1073-1083`; drift classes `tests_cosmo/data/oracles/README.md:24-27, 80-110`.

**`vector_smoke.wls`:** 75 lines; harvest from `Global`VectorTheory` (`:66`); stale caveat `:25-35` (see `stage1_engineering_plan.md:658-665`, #543 expired 2026-09-11).

**Lane guard:** `.claude/hooks/wolfram-guard.sh:7-21, 53-66, 79-86` is fail-open for `wolframclient`, `.py` wrappers and compound commands. `scripts/psalter/run_tier1_gate.sh:76-83` self-guards.

**Verifier exit codes:** `scripts/verify-wolfram-setup.sh:20-24, 99-105` (0/1/2; `--require-psalter` promotes soft-fails). `ensure_registered.sh:23-66`.

---

## Appendix C — Not verified, or could not be read (carried into the memo's "tested here?" columns)

**CAMB:**
- The CAMB notes PDF maths could not be read.
- η₀₁₂₃ sign in Challinor 2000 and Lewis 2004 only grepped.
- slegner/cobaya not re-diffed.
- Whether EFTCAMB "User defined" options need Fortran edits.

**Wolfram behavior that needs a kernel:**
- Multi-line `ToExpression`.
- `$Pre`/`$PreRead` inside `ToExpression`.
- `Level[..., HoldComplete, Heads -> True]` non-evaluation.
- Which xAct/PSALTer symbols carry upvalues.
- Script exit status after `Abort[]`.
- The code that writes the polology WXF was not found.

**Outside tools:**
- sympy2c cache key internals and license.
- SARAH's actual `SM.m` (quoted from the paper instead).
- FeynRules website model files.

**Cobaya:**
- Whether a component option can be excluded from the resume comparison.
- Whether a dependency-free check Theory runs before CAMB.

**Physics and design:**
- Whether a per-sample spectrum check can locate the derived spectrum from inside an external prior.
- The torsion expansion path is not exercised by R-1's Vector prototype.
