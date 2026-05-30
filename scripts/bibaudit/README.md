# bibaudit — citation audit tooling

Read-only helpers for the manuscript citation audit. **None of these scripts edit
`references.bib` or any `.tex` file.** They fetch canonical records, compare, and
*propose* exact edits; all actual changes go through the Edit tool (exact string
replacement that fails loudly on mismatch) so a parser bug can never silently
corrupt the bibliography. See `manuscript/planning/citation_audit_tracker.md` for
live progress and `/home/vscode/.claude/plans/working-on-the-manuscript-*.md` for
the full plan.

## Modules

| file | role | network |
|------|------|---------|
| `bibio.py` | tolerant field extractor (handles `{}`/`""`, nested braces, accents); tracks `% === section ===` banners | no |
| `normalize.py` | field normalization + `classify()` → MATCH/MINOR/ENRICHED/SUBSTANTIVE | no |
| `journal_map.json` | journal full↔abbrev equivalence (comparison aid only) | no |
| `fetch_canonical.py` | INSPIRE → Crossref → Google Books, cached + rate-limited | yes |
| `compare.py` | Task 1 driver: classify all 231, write `cache/report.json` + proposed edits | yes |
| `extract_cites.py` | Task 2 worklist: every `\cite{}` grouped by key, with claim context | no |
| `audit_status.py` | progress dashboard: parses `%@audit`/`%@cite` markers → X/231, Y/192 | no |

## Usage

```bash
# from repo root, always via -m so the scripts.bibaudit package resolves
uv run python -m scripts.bibaudit.fetch_canonical        # prime the cache
uv run python -m scripts.bibaudit.compare                # classify + report
uv run python -m scripts.bibaudit.compare --summary-only # re-print last report
uv run python -m scripts.bibaudit.extract_cites          # Task 2 worklist
uv run python -m scripts.bibaudit.audit_status           # progress
```

`cache/` (gitignored) holds raw API responses + `report.json` + `cites.json`,
making every run idempotent and resumable. Pass `--refresh` to force re-fetch.

## Verdicts

- **MATCH** — fields agree after normalization.
- **MINOR** — only cosmetic differences (journal abbrev, page dash, brace/case,
  volume-numbering convention). Safe to normalize.
- **ENRICHED** — canonical has fields the local entry lacks (fill from cache).
- **SUBSTANTIVE** — author surname / title token-set / year differ → **STOP**,
  resolve with the user (the mis-attribution class).
- **NOTFOUND / NOID** — has an identifier but no record / has no identifier
  (books, pre-arXiv papers) → manual review.

## Marker schema (written to references.bib via Edit, above each entry)

```bibtex
%@audit v=1 date=YYYY-MM-DD src=<url> verdict=MATCH
%@cite v=1 status=SUPPORTS
%! purpose: <what this paper backs>
%! used: results.tex (§4.1), theory.tex (§2.3)
%! checked: literature/<id> [abstract+§3], YYYY-MM-DD
@article{key, ...}
```
