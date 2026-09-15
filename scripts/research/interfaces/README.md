# `interfaces/` — research artifacts behind the R-1 interfaces decision

**Status: prototypes run 2026-09-15; results in `docs/cosmology/interfaces_decision.md`.**

Supporting material for `docs/cosmology/interfaces_decision.md` (R-1, #566), which recommends
how Python drives the Wolfram derivation (D-A), what a user writes and how the derived
spectrum reaches Cobaya (D-B), and how sign conventions are carried
(`docs/cosmology/conventions.md`). The full planning reasoning trail — evidence ledger,
revision log of user corrections, what could not be verified — is archived verbatim in
`docs/cosmology/r1_planning_record.md`. Nothing here runs in the pipeline; it exists so the
sessions that build the interface (I-S1A-core, I-S1B, M1b) inherit measured evidence rather
than re-deriving it. **Research code, not production code.**

## Contents

| path | what it is | provenance |
| --- | --- | --- |
| `theories/VectorTheory.yaml` | the recommended theory-file schema (Option A′): `fields`, `lagrangian` as coupling → operator, `derivation` | ours; physics from `SupplementalMaterials-2607` `Models/VectorTheory.m` (Barker) |
| `theories/bad/*.yaml` | eleven failure and injection fixtures (underscore name, undeclared symbol, coupling inside an operator, two expressions, zero operator, name collision, `Print`/`Quit[]` injection, inexact number, numeric key, duplicate key) | ours |
| `runs/vector_one.yaml` | a standard Cobaya run file that includes the theory with `model: !defaults ../theories/VectorTheory` | ours |
| `wolfram/Stage1Proto.wl` | the committed Wolfram package prototype: parse operators held, whitelist, declare in `Global`, `ParticleSpectrum`, export `spectrum.wxf` | ours |
| `wolfram/driver.wls` | the FIXED driver (never generated), with `Catch` + `Exit[1]` | ours |
| `wolfram/compare_waveoperators.wls` | in-kernel `SameQ` of wave operators across routes, three-valued on failure | ours |
| `proto/theory_store.py` | theory loading with Cobaya's file loader, fingerprint, content-addressed data store | ours |
| `proto/derive_proto.py` | prototype of `tidalcosmo derive` (`--test`, `--force`, `--dry-run`, `--validate-only`, `--timeout`) | ours |
| `proto/spectator_gate.py` | the Cobaya Theory that refuses at construction without a matching derived spectrum | ours |
| `proto/run_gate_cases.py` | twelve refusal/pass cases for the gate, plus the static no-kernel import check | ours |
| `proto/session_probe.py` | the `wolframclient` session measurement (package through a session; timeout; `Quit[]`) | ours |

## How to reproduce

All Wolfram steps are strictly serial (one kernel, machine-wide). Outputs go to gitignored
`third_party/interfaces_runs/`.

```bash
bash scripts/psalter/ensure_registered.sh --no-verify
bash scripts/verify-wolfram-setup.sh --require-psalter
cd scripts/research/interfaces
S=../../../third_party/interfaces_runs
# reference: the unchanged standalone script, run from a throwaway directory
(mkdir -p $S/reference && cd $S/reference && QT_QPA_PLATFORM=offscreen wolframscript -file ../../../scripts/psalter/vector_smoke.wls)
python proto/derive_proto.py runs/vector_one.yaml --store $S/store --log-dir $S/logs/derive_good
for f in theories/bad/*.yaml; do python proto/derive_proto.py $f --store $S/failures/store; done
python proto/session_probe.py theories/VectorTheory.yaml $S/session
QT_QPA_PLATFORM=offscreen wolframscript -file wolfram/compare_waveoperators.wls \
  reference=$S/reference/ParticleSpectrographVectorTheory.mx \
  package=$S/store/vr1-proto-1/57/<fingerprint>/spectrum.wxf
python proto/run_gate_cases.py $S/store $S/gate_scratch    # no Wolfram
```

## Rules this directory follows

- **Route, not payload.** Nothing upstream is vendored; throwaway virtual environments for
  comparable packages live under gitignored `third_party/`; the repository `.venv` is never
  modified.
- **One Wolfram kernel at a time.** Every launch sets `QT_QPA_PLATFORM=offscreen`, runs from a
  throwaway working directory, and is judged by the artifacts and sentinel lines it produced —
  never by its exit status alone. Python-launched kernels self-guard with `pgrep`, because the
  lane hook cannot see them.
- **No machine-specific paths** in anything committed; transcripts are scrubbed at print time
  (`tests/test_repo_hygiene.py` scans tracked files).
