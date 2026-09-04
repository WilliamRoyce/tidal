# `psalter_stage1/` — research artifacts behind the H8 Stage-1 study

Supporting material for `docs/cosmology/stage1_engineering_plan.md` (H8, WS6/#495).
Nothing here runs in the pipeline; it exists so the session that implements Stage 1
inherits the evidence rather than re-deriving it.

## Contents

| path | what it is | provenance |
| --- | --- | --- |
| `fetch_reference_sources.sh` | pinned downloader for the upstream Barker sources the study cites | ours |
| `wxf_decode/pub_wxf.py` | minimal WXF reader: PSALTer spectrograph export → per-sector `M0/M1/M2` as SymPy matrices, keyed on `xAct`PSALTer`Def` | **ours** (H6, 2026-08-30/09-03) |
| `wxf_decode/tar_wxf.py` | the same for the `psalter.tar.gz` WIP snapshot's variant layout | **ours** (H6) |

The two decode scripts were written during the H6 design study and lived only in `/tmp`,
where they would have been lost on the next container rebuild. They are the basis for
`spectrum_design.md` §6.1 — the observation that the released positional spin labeling is
wrong on the release's own `A23Theory` fixture came out of running them. They are preserved
here verbatim as research code; the production reader that Stage 2 consumes is a separate,
tested module (`tidalcosmo/spectrum/`, per the study §5), not these files.

## Why the upstream sources are fetched, not vendored

PSALTer and its supplemental materials are GPL-3.0-or-later; TIDAL is MIT. Reading and
adapting them is authorized (decision D6 — the author's permission is explicit, and every
borrowed function records its origin in a docstring). *Redistributing* them from this
repository is the trigger for the license question recorded as a release blocker on #495.
So this directory commits the route and not the payload:

```bash
bash scripts/research/psalter_stage1/fetch_reference_sources.sh
```

downloads into `third_party/psalter_reference/` (gitignored) at pinned revisions —
PSALTer `bb45adb0` (v2.0.2), plus the specific files from `SupplementalMaterials-2506b`
(the TorC/CTEG companion) and `-2607` (numerical polology) that the study analyses.
