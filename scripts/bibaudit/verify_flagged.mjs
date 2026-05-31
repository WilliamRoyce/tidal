export const meta = {
  name: 'verify-flagged-deep',
  description: 'Deep-completion: re-verify every flagged citation (MISMATCH/WEAK/PARTIAL/UNVERIFIABLE) against the best available source (OA-PDF full text / INSPIRE+S2 abstract / Google Books / secondary), per-site, with adversarial re-check',
  phases: [
    { title: 'Load', detail: 'list flagged-source records' },
    { title: 'Verify', detail: 'per-key deep adjudication against gathered source' },
    { title: 'Recheck', detail: 'adversarial second opinion on MISMATCH/WEAK' },
  ],
}

const VERDICT = {
  type: 'object',
  required: ['key', 'status', 'purpose', 'depth_reached', 'source_checked', 'confidence', 'sites'],
  properties: {
    key: { type: 'string' },
    status: { type: 'string', enum: ['SUPPORTS', 'PARTIAL', 'WEAK', 'MISMATCH', 'UNVERIFIABLE'] },
    purpose: { type: 'string' },
    depth_reached: { type: 'string' },
    source_checked: { type: 'string' },
    confidence: { type: 'string', enum: ['high', 'medium', 'low'] },
    evidence_quote: { type: 'string' },
    sites: {
      type: 'array',
      items: {
        type: 'object',
        required: ['file', 'line', 'role', 'site_status'],
        properties: {
          file: { type: 'string' }, line: { type: 'integer' },
          role: { type: 'string', enum: ['EXISTENCE', 'ATTRIBUTION', 'SPECIFIC'] },
          site_status: { type: 'string', enum: ['SUPPORTS', 'PARTIAL', 'WEAK', 'MISMATCH', 'UNVERIFIABLE'] },
          evidence: { type: 'string' },
        },
      },
    },
  },
}

function prompt(fname) {
  return `Re-verify ONE flagged manuscript citation as thoroughly as the available source allows.
This is the deep-completion pass: the goal is to move UNVERIFIABLE -> a real verdict and to confirm
or refute each PARTIAL/WEAK/MISMATCH per-site.

STEP 1 — read scripts/bibaudit/cache/flagged_src/${fname}.json
  Fields: key, final_status (current), sites[] (each \\cite with {file,line,role,claim}), and the
  gathered source: fulltext_paths[] (local .tex or oa.txt full text — READ THESE FIRST if present),
  abstract, tldr, gbooks (book preview), oapdf_url.

STEP 2 — read the richest source available, in order:
  (a) fulltext_paths[] — Read/grep the full text for the specific formula/number/claim of each site.
  (b) else the abstract / tldr / gbooks text in the json.
  (c) if the json source is thin AND the claim is specific, you MAY WebFetch the oapdf_url or
      WebSearch for the paper's content / a review describing it (secondary source) — cite what you used.
  Do NOT just restate the abstract; for ATTRIBUTION/SPECIFIC sites, find the actual passage.

STEP 3 — adjudicate EACH site by role:
  EXISTENCE/CONTEXT: topic match in title/abstract -> SUPPORTS.
  ATTRIBUTION/SPECIFIC: the exact claim must appear in the source -> SUPPORTS (quote it); if the source
  is silent -> WEAK or UNVERIFIABLE; if it contradicts -> MISMATCH (quote it).

STEP 4 — return the verdict for key: status = worst across sites (but a clean topic match for pure
  EXISTENCE cites is SUPPORTS, not UNVERIFIABLE). Fill purpose, depth_reached (e.g. deep:oa-pdf,
  abstract, secondary:<source>), source_checked (what you actually read), confidence, evidence_quote
  (for non-SUPPORTS), and per-site verdicts with quotes. Be exact and skeptical — this is the audit's
  final safety net.`
}

function recheckPrompt(fname, v) {
  return `A pass flagged this citation as ${v.status}. ADVERSARIAL re-check: find the STRONGEST evidence
the paper DOES support each claim; if you cannot, confirm with a quote. Read
scripts/bibaudit/cache/flagged_src/${fname}.json, then the deepest source (fulltext_paths, else
WebFetch oapdf_url / WebSearch the paper). Return the same verdict schema with your honest re-assessment.`
}

phase('Load')
const listed = await agent(
  'Read directory scripts/bibaudit/cache/flagged_src/ (ls or Glob). Return {files:[...]} of the '
  + '.json basenames WITHOUT extension.',
  { label: 'load-flagged', phase: 'Load', schema: { type: 'object', required: ['files'], properties: { files: { type: 'array', items: { type: 'string' } } } } },
)
const FILES = (listed && listed.files) || []
log(`flagged keys to re-verify: ${FILES.length}`)

const results = await pipeline(
  FILES,
  (f) => agent(prompt(f), { label: `verify:${f}`, phase: 'Verify', schema: VERDICT }),
  (v, f) => {
    if (v && ['MISMATCH', 'WEAK'].includes(v.status)) {
      return agent(recheckPrompt(f, v), { label: `recheck:${f}`, phase: 'Recheck', schema: VERDICT })
        .then((r) => ({ ...v, recheck: r })).catch(() => v)
    }
    return v
  },
)
const clean = results.filter(Boolean)
const tally = {}
for (const r of clean) { const s = (r.recheck && r.recheck.status) || r.status; tally[s] = (tally[s] || 0) + 1 }
log(`re-verified ${clean.length}/${FILES.length}; final tally: ${JSON.stringify(tally)}`)
return { n: clean.length, tally, results: clean }
