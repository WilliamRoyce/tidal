# `interfaces/` — research artifacts behind the R-1 interfaces decision

**Status: IN PROGRESS (2026-09-15).** Planning complete and archived; prototypes not yet run.

Supporting material for `docs/cosmology/interfaces_decision.md` (R-1, #566), which
recommends how Python drives the Wolfram derivation (D-A), what a user writes and how the
derived spectrum reaches Cobaya (D-B), and how sign conventions are carried. The full
planning reasoning trail — evidence ledger, the revision log of user corrections, and what
could not be verified — is archived verbatim in `docs/cosmology/r1_planning_record.md`.
Nothing here runs in the pipeline; it exists so the sessions that build the interface
(I-S1A-core, I-S1B, M1b) inherit measured evidence rather than re-deriving it.

## Contents

| path | what it is | provenance |
| --- | --- | --- |
| `README.md` | this file | ours |

*(Rows are added as each prototype lands: the committed Wolfram prototype package and
driver, the theory/run YAML fixtures, the failure and injection fixtures, the wave-operator
comparison script, the `wolframclient` session measurement, the Cobaya gate prototype, and
the `!defaults` include check.)*

## Rules this directory follows

- **Route, not payload.** Anything fetched from upstream is downloaded by a committed script
  into gitignored `third_party/` at a pinned revision — the rule in
  `scripts/research/psalter_stage1/README.md`. Throwaway virtual environments for
  comparable packages live there too; the repository `.venv` is never modified.
- **One Wolfram kernel at a time, machine-wide.** Every Wolfram launch sets
  `QT_QPA_PLATFORM=offscreen`, runs from a throwaway working directory, and is judged by the
  artifacts and sentinel lines it produced — never by its exit status alone.
- **No machine-specific paths** in anything committed; transcripts are scrubbed at print
  time (`tests/test_repo_hygiene.py` scans tracked files).
