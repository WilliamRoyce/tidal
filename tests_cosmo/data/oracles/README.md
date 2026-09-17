# The frozen legacy oracle

Committed output from legacy `tidal`, captured once so the new package can be checked
against the old physics **without ever running legacy from a test**. New-package tests
assert against these files. Nothing here is edited by hand.

Produced by [`scripts/oracles/freeze_legacy_oracle.py`](../../../scripts/oracles/) —
the one place in the repository allowed to touch legacy. Regenerate with:

```bash
uv run python -m scripts.oracles.freeze_legacy_oracle          # rewrite
uv run python -m scripts.oracles.freeze_legacy_oracle --check  # compare only
```

GH #525 · milestone M0.5 · design `docs/cosmology/repo_reshape.md` §7, §8.

## The gate these support is semantic, not a byte diff

The obvious gate would be `tidal inspect OLD --diff NEW`. **It is the wrong one here.**
#513 has the new package emit CAMB and PSALTer conventions *natively* rather than
translating, so names, formats and conventions change deliberately — a byte diff would
report every intended change as a failure (§5.2).

The frozen spec is the oracle **for the physics**: same equations, same coefficients, same
signs under the new naming, established by careful comparison and recorded as **a written
mapping committed alongside these fixtures**. That mapping is *part of the gate, not
commentary on it*. **It does not exist yet** — WS2 owes it at M3.

## What is frozen

49 theory/spec pairs × 4 artifact kinds, plus `manifest.json` — 197 files. (**49, not the original 46**: the four theories excluded at M0.5 for want of the Wolfram lane were derived on 2026-09-12 (#547); three succeeded in under a minute each and are now pairs, and `theory_radial` cannot be derived at all — the manifest's single exclusion carries its measured reason.)

| kind | path | what it is |
| --- | --- | --- |
| spec JSON | `specs/<id>.json` | copied byte-for-byte from `examples/data/`; never regenerated |
| sign table | `inspect/<id>.summary.txt` | `tidal inspect --detail summary` — the proven signs |
| index structure | `inspect/<id>.families.txt` | `tidal inspect --families` |
| verdicts | `validate/<id>.txt` | `tidal validate` **and** `tidal validate --stability` |

`<id>` is the spec stem, so `cmp specs/x.json examples/data/x.json` is a one-line check.

**Not frozen**, removed by the 2026-09-06 amendment to §8: *measured scalars*, because
`simulate` is re-implemented against WS3's own gates and `measure` is a `drop` row whose
energy and conversion quantities have no FRW counterpart; and *reference `C_ℓ`/transfer
arrays*, because legacy is flat-space and cannot produce a `C_ℓ` at all. Those are M1a's.

## Why the specs are *copied* rather than referenced

Not for preservation — that reason was checked and is false. All 48 spec JSONs are
git-tracked and cannot be lost, and `examples/` appears in no retire column: M6 retires
`tidal/`, `tests/` and the `tidal` console script, not the examples. Deleting legacy does
not touch them.

They are copied because **a copy is a pin and a reference is a moving target**:

1. **M3's own output location is undecided.** M3 re-derives this corpus under the new
   conventions (#513), and nothing says where the new specs are written. If they land in
   `examples/data/`, the old spec is overwritten and the reviewer of the §5.2 mapping
   loses one side of the comparison at exactly the moment it is needed.
2. **#397 makes drift concrete rather than hypothetical.** 18 of these specs are
   known-defective under an **open** issue. Against a referenced corpus, re-deriving them
   would leave the frozen verdicts below describing specs that no longer produce them — an
   oracle silently inconsistent with itself, which is the worst failure mode a gate has.

### What this oracle is the oracle *for*, and three drift classes it carries (2026-09-11)

The frozen spec is **legacy as committed** — a pin of what the derivation produced when each
file was last derived — not "legacy as of the version that froze it". The corpus is
**mixed-vintage** by construction (files last derived between v0.31 and v0.53), and
`--check` compares fresh `inspect`/`validate` runs **over the committed spec** against the
committed reports: it detects an edited spec or a changed reader, and it is **structurally
blind** to "the committed spec is older than what `derive` would produce today" (#554).
That blind spot has a license-free detector, `freeze_legacy_oracle.py --staleness`
(`derivation_hash` against the hash of the driver script `tidal derive` would generate here —
see `scripts/oracles/README.md`; as of 2026-09-11 it reports **0 current / 40 stale / 6
without a hash**, because every committed spec predates the generator's last text change on
2026-08-18 — vintage, not drift), and it is why re-deriving in place from
`examples/*/run.sh` now needs `FORCE_DERIVE=1`.

Three drift classes are known and are **representational, not physical**. The §5.2 mapping
must not treat their absence in the new package as a regression, and a port that reproduces
the *physics* while dropping them is correct:

1. **#397 — 18 specs carry the stale `a_0` sign** (the section below). A defect the port must
   not carry forward; their frozen verdict is a specification, not a target.
2. **The `coupling` block is present in 46 of the 49 frozen specs and is no longer emitted**
   (the three added on 2026-09-12 are fresh derivations, so they do not carry it — which is
   itself the clearest statement of the class).
   `efd18a5b` (#403/#404, 2026-08-17) deleted `ExportJSON.wl`'s `coupling` output as
   write-only dead data — nothing in `tidal/` reads it (`git grep '\["coupling"\]'` → zero) —
   and verified the loaded model byte-identical with and without it. A fresh derive omits
   it; the committed files keep it. `coupled_scalars` is the case #554 noticed.
**A fifth class arrived with R-C (2026-09-16), and it is a *defect in the physics*, not a
representation (#582).** Legacy's post-Riemannian rewrite
(`_derive.py::_wls_torsion_curvature_decomposition`) uses an expression that is not a
contortion and carries the opposite sign to xTensor's connection-difference tensor, so **every
`R̃` torsion term it produced differs from the correct identity — the quadratic mass terms
included.** The frozen specs derived from `examples/graviton_torsion/theory.toml`,
`examples/torsion_dark_photon/theory.toml` and `examples/torsion_gertsenshtein/theory.toml`
(and the FV variant where it uses `R̃`) therefore encode wrong torsion terms. **M3 must
attribute that difference to this defect rather than report a port regression**, and must port
the corrected identity (`stage1_engineering_plan.md` §4.3 item 2). **The user decided on
2026-09-16 that legacy is not corrected or re-derived and that no thesis-era number is
re-checked** — open by decision, not by oversight; evidence and the correct identity are in
#582 and `docs/cosmology/perturbation_tooling.md` §5, §7.1.

**A fourth class arrived with the #547 theories (2026-09-12), and it is a *defect in the
reader*, not in the spec.** `de_sitter_kg` fails plain `tidal validate` on the #394
volume-element check, which demands a **spatial** first-derivative term whenever `sqrt|g|` is
non-constant. A measure that depends on **time only** — every FRW background — generates a
*temporal* one instead, and this spec carries it (`first_derivative_t(phi_0): -dSH`, i.e.
`phi_tt + H phi_t - laplacian(phi) + a^2 m^2 phi = 0`, the analytic result). The check had only
ever seen spatial measures (`polar_kg`, `spherical_kg_1d`, both of which pass). Filed as
**#560**; the frozen `validate` fixture pins the *current* reader, so it will change when #560
is fixed — that is a reader change, which is exactly what `--check` exists to catch, and the
mapping must not read it as a physics regression. Related: `de_sitter_kg`'s
`canonical.volume_element` is `e^{2Ht}` where `sqrt|g| = e^{3Ht}` for its metric, and the
velocity-squared `hamiltonian_terms` carry `1/(2 e^{2Ht})`. The **equations are verified
correct**; the energy-side normalization is **not** verified, so treat this fixture as an
oracle for the equations and not for energy (noted on #547 and #560, and it matters because
`de_sitter` is WS2's designated oracle).

3. **`component_metadata` (`tensor_head`/`tensor_rank`/`tensor_indices` on `fields[]`) is
   absent from 11 of the 49** that predate `78374c1`: `conformal_kg_static`, `navier_cauchy_2d`,
   `graviton_torsion`, `massive_3form`, `scalar_vector_coupling`, and six
   `torsion_gertsenshtein*` variants. The port's field metadata is richer than these by
   design.

**A full re-freeze was considered and rejected (2026-09-11).** It would take `--stability`
coverage from 6/49 to 1/49 — five of the six passes depend on `metadata.parameters`, which
#232 stopped injecting — three theories cannot be re-derived at all on the current pipeline
(#321 ×2 including the main thesis theory, #402 ×1), and it costs 8–12 hours of the single
kernel to strip data nothing reads. The pin stays; the classes above are the contract.

Cost: 7.0 MB raw, ~0.3 MB compressed, which is what git stores.

**Considered and rejected:** recording each spec's git blob SHA instead of copying. It
gives the same pin, but it puts a `git cat-file` step between a reviewer and the evidence
months later, and it degrades under shallow clones and exported archives — the same reason
`tests/data/gertsenshtein_ungauged_pre397.json` is vendored rather than read from history.

## The standing rule

> **If `tidal/` or `examples/data/` changes, re-run `scripts/oracles/` in the same commit.**

An action rather than a prohibition: "never edit legacy" is unenforceable, and the point
is that the fixtures never describe a corpus that has moved underneath them.
`--check` is the detector.

## Provenance

Frozen from `feat/cosmology-program` at `ca5d524d`, with:

| | |
| --- | --- |
| legacy `tidal/` tree | `d6a0c5f9` |
| `examples/data/` tree | `1eb0b451` |
| `tidal --version` | 0.53.0 |
| Python | 3.11.2 |
| NumPy | 2.3.4 |
| platform | Linux x86\_64 |

Tree hashes rather than a commit SHA: they identify exactly what was frozen and do not
change when an unrelated file is committed. This file is hand-written and never
regenerated, which is why the provenance lives here — anything varying with `HEAD`, the
machine or the clock inside a *generated* file would break the byte-for-byte gate on the
next re-run.

## The scrub

Every captured stream passes through this, in order. Each layer is a no-op that the next
can still catch.

1. **A constructed environment**, never inherited: `NO_COLOR=1`, `TERM=dumb`,
   `COLUMNS=80`, `LC_ALL=C.UTF-8`, `TZ=UTC`, `PYTHONHASHSEED=0`, `MPLBACKEND=Agg`, and
   BLAS threads pinned to 1. A developer's `FORCE_COLOR` or `COLUMNS` therefore cannot
   reach committed data. `NO_COLOR` is explicit because `tidal/cli/_console.py` caches its
   color decision from a TTY check on `sys.stderr` in a module global.
2. **Relative paths from the repository root.** Neither `inspect` nor `validate` resolves
   its argument, and both print only the basename on success, so no absolute path reaches
   an artifact by this route.
3. **ANSI escape stripping** — a no-op under (1), asserted anyway.
4. **Path substitution**, `<repo>` and `<venv>` only. Deliberately *not* the home
   directory or `sys.base_prefix`: the latter is `/usr`, so substituting it would rewrite
   a stray stdlib path into a portable-looking token and the guard in (6) would never see
   it. A path from outside the repository means something leaked that this scrub did not
   anticipate, and that aborts the run.
5. **`PYTHONWARNINGS=ignore`.** A decision, not tidiness: Python's warning format embeds
   an absolute path *and* a legacy line number, so an unsuppressed warning would fail
   repo-hygiene and churn on any unrelated legacy edit. **The consequence is that a
   warning legacy emits is not in the oracle.** Exactly one spec trips this —
   `torsion_gertsenshtein_exact` under `--families`, which emits
   `UserWarning: [json_loader] time_order > 2 fields present without [perturbation];
   loaded for read-only use.` Recorded here because it is otherwise invisible.
6. **Trailing-whitespace normalization** on the text artifacts: CRLF→LF, per-line trailing-space removal, exactly one terminating newline. `inspect` ends on a bare `print()`, so its
   stdout finishes on a blank line; left alone, pre-commit's `end-of-file-fixer` would
   rewrite the committed file and the next run would rewrite it back. **So the text
   artifacts are not literally raw stdout.** The spec JSONs are untouched.
7. **A guard that aborts the whole run** on any absolute path, ANSI escape, timestamp,
   timing line, hostname, username, or non-UTF-8 byte. Non-UTF-8 matters specifically
   because `tests/test_repo_hygiene.py` silently *skips* a file it cannot decode, so such
   an artifact would escape the absolute-path check rather than fail it.

A spec legacy cannot load is recorded as a status row and the run continues; a scrub
violation aborts. The asymmetry is deliberate — the failure being prevented is committing
a machine-specific path, not missing a row.

## What the verdicts actually say

Measured across the 49 pairs (regenerated 2026-09-12):

| | |
| --- | --- |
| `tidal validate` (plain) | **29** pass, **20** fail |
| `tidal validate --stability` | **6** pass, 43 fail — of which **42** are "cannot evaluate symbolic coefficient" |
| `tidal inspect --detail summary` | 47 exit 0; `gertsenshtein_dipolar_centered` and `torsion_gertsenshtein_exact` exit 1 (both `sign_algebra: unsupported node Subscript`) |
| `tidal inspect --families` | all 49 exit 0 |

### The 18 #397 specs — "same verdicts" is the wrong gate for these

18 of the 19 plain-`validate` failures are the same verdict:
`Field 'a': components 'a_0' and 'a_1' have opposite effective self-laplacian`. That is
**#397**, which is **open**: 25 committed JSONs predate the #381 sign fix (`49fcd59`,
2026-05-25, which re-derived only five) and encode `∂²ₜa₀ = −∇²a₀`, a temporal-only
tachyon confirmed empirically — the committed spec diverges where the re-derived one
conserves energy exactly.

Three things are easy to conflate. **#381** fixed the *derivation*. **#401**'s semantic
accessors fixed how a spec is *read* — naive scans had produced six confidently-wrong
diagnoses — and could not change committed data. **#397** is the remaining *data* defect.
So `validate` rejecting these is correct behavior, not a misreading.

**Consequence for M3, and this is the important part.** Legacy is migrated, not repaired:
nobody re-derives these 18 in `examples/data/`, so their frozen verdict stays `fail`
permanently. That makes it **a specification of a defect the port must not carry
forward**, not a verdict to reproduce. A new package that *passes* these 18 has fixed the
bug; treating that as a regression would either block M3 on a gate it ought to fail or
paper over the difference silently.

The 19th failure is `torsion_gertsenshtein_exact`, which fails to load on `time_order > 2`
without a `[perturbation]` section — unrelated to #397. It is also the one spec where
`inspect --detail summary` exits 1 while `--families` succeeds, because `--detail` is
missing from the query-flag list that relaxes the v6 strict guard
(`tidal/cli/_inspect.py:497`).

### `--stability` is a weak oracle on its own, and that is recorded rather than fixed

42 of the 49 `--stability` captures are **not stability verdicts**. They record that the
spec has free symbolic parameters with no committed defaults, so the tachyon check never
runs. Only **5** specs carry `metadata.parameters`; those plus one needing none are the
6 that pass. As a port gate, this corpus exercises the stability path on
**6** specs and the error path on 42 (plus one load failure).

No `--param` values were invented to improve this. Choosing them would be choosing
physics, and would make the oracle a record of *our* parameter choices rather than of
legacy's behavior. The lesson is carried forward instead as a requirement on the **new**
package's spec contract: parameter defaults must travel with the spec. Legacy is not
patched for it.

## The corpus, and what is not in it

Enumerated as: every `examples/*/*.toml` whose `[output].path`, resolved relative to the
TOML's own directory, ends in `.json` — **50** TOMLs, of which **49** have a committed
spec (46 until 2026-09-12; see #547).

> The recipe originally given in §8 was `examples/*/theory*.toml`. It matches 48 files and
> silently drops `examples/curved_spacetime/conformal_static.toml`, giving 45. Amended at
> both instruction sites; see GH #525.

**The four exclusions**, named rather than silently missed. Each has no committed spec,
and deriving one needs a `wolframscript` kernel — one machine-wide, held by another
session in Wave 0:

| theory | expected spec |
| --- | --- |
| `examples/curved_spacetime/de_sitter.toml` | `examples/data/de_sitter_kg.json` |
| `examples/gertsenshtein/theory_dipolar.toml` | `examples/data/gertsenshtein_dipolar.json` |
| `examples/gertsenshtein/theory_dipolar_centered.toml` | `examples/data/gertsenshtein_dipolar_centered.json` |
| `examples/gertsenshtein/theory_radial.toml` | `examples/data/gertsenshtein_radial.json` |

**Two adjacent categories**, so the arithmetic audits cleanly:

- **Five `sweep_*.toml`** files under `coupled_scalars/` and `coupled_scattering/` are not
  theories; their `[output].path` is a directory, which is what the `.json` test excludes.
- **Two spec JSONs have no theory TOML** — `cylindrical_kg_1d.json` and
  `gw_plane_wave_1d.json`. They are live (used by `tests/test_solver_constraint_solve.py`
  and `tests/test_implicit_dynamical_sector.py`, and listed in
  `tests/data/spec_semantics.txt`), but with no TOML they are outside the pair predicate
  and are therefore not frozen. That means two spec shapes — cylindrical coordinates and a
  1-D plane wave — are absent from this oracle.

Six specs carry no `metadata.derivation_hash`, so their provenance is weaker than the
rest: `conformal_kg_static`, `navier_cauchy_2d`, `massive_3form`, `scalar_vector_coupling`,
`spherical_kg_1d`, `torsion_gertsenshtein_exact`.

## Retirement

**These fixtures have an expiry.** They exist to gate M3's port. Once M3's §5.2 mapping is
written and recorded, they retire with the rest of legacy rather than being maintained —
see §7's M7 row, which also retires `tests_cosmo/test_package_boundary.py` and the
shell-out check in `tests_cosmo/test_oracles.py`. A frozen oracle nobody remembers to
retire is the mirror image of a guardrail nobody remembers to delete.
