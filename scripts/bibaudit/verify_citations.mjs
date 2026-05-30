export const meta = {
  name: 'verify-citations',
  description: 'Task 2: verify each cited paper supports its claim (layered read), with adversarial re-check of weak/mismatch verdicts',
  phases: [
    { title: 'Load', detail: 'list per-key worklist files' },
    { title: 'Verify', detail: 'one agent per key, layered read vs offline TeX / cached abstract' },
    { title: 'Recheck', detail: 'adversarial second opinion on MISMATCH/WEAK/PARTIAL' },
  ],
}

// Structured verdict the verify agents must return.
const VERDICT = {
  type: 'object',
  required: ['key', 'status', 'purpose', 'depth_reached', 'source_checked', 'confidence', 'sites'],
  properties: {
    key: { type: 'string' },
    status: { type: 'string', enum: ['SUPPORTS', 'PARTIAL', 'WEAK', 'MISMATCH', 'UNVERIFIABLE'] },
    purpose: { type: 'string', description: 'one line: what this paper is cited for' },
    depth_reached: { type: 'string', description: 'cheap | abstract | abstract+intro | deep:<section>' },
    source_checked: { type: 'string', description: 'file/path actually read, e.g. literature/2406.12826 [abstract,intro]' },
    confidence: { type: 'string', enum: ['high', 'medium', 'low'] },
    evidence_quote: { type: 'string', description: 'load-bearing quote for non-SUPPORTS (else empty)' },
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

function verifyPrompt(safe) {
  return `You are verifying ONE manuscript citation for a physics paper (torsion / Poincaré gauge
theory + the Gertsenshtein graviton–photon effect). Goal: does the cited paper actually
SUPPORT the claim it is attached to in the prose?

STEP 1 — read the worklist record:
  Read scripts/bibaudit/cache/worklist/${safe}.json
  It has: key, tier (A/B/C/D/E), eprint, doi, ident, in_notes, and sites[] — each site is a
  \\cite call with {file, line, claim (the enclosing sentence), role (EXISTENCE/ATTRIBUTION/SPECIFIC)}.

STEP 2 — locate the paper's content (by tier), cheapest layer first:
  tier A: ls literature/<eprint>/  then Read the main .tex (title + abstract + intro + conclusion).
          You MAY deep-read/grep the body for a specific formula/number ONLY if step 3 needs it.
  tier B or C: Read scripts/bibaudit/cache/abstracts/<ident>.json (title + abstract). This is the
          only offline source — do NOT attempt network. If a SPECIFIC/ATTRIBUTION claim cannot be
          adjudicated from the abstract, the right verdict is UNVERIFIABLE (not MISMATCH).
  tier D/E: no abstract — judge only from the bibliography fields + claim wording.
  If in_notes is true, also grep manuscript/planning/literature_content_notes.md for the eprint or
  key as a PRIOR (a hint, not gospel — still confirm against the source).

STEP 3 — judge support, per site, by role (LAYERED — do not over-read):
  EXISTENCE/CONTEXT site (e.g. "see also", "as detected by", grouped landmark cite): if the paper's
      title/abstract plainly matches the topic -> SUPPORTS at depth "cheap". Do NOT deep-read these.
  ATTRIBUTION / SPECIFIC site (a formula, number, result, or "first derived by"): confirm the exact
      claim appears in the source. If confirmed in abstract -> SUPPORTS. If the abstract is silent and
      tier A -> deep-read the body for it. If silent and only an abstract exists -> UNVERIFIABLE.
      If the source CONTRADICTS the claim -> MISMATCH (quote the contradicting passage).

STEP 4 — return the structured verdict for key from the worklist:
  status = worst across sites (MISMATCH > WEAK > PARTIAL > UNVERIFIABLE > SUPPORTS, but a clean
  topic-match for pure EXISTENCE cites is SUPPORTS/low-confidence, NOT UNVERIFIABLE).
  Fill purpose (one line), depth_reached, source_checked (the actual path you read), confidence,
  evidence_quote (for any non-SUPPORTS), and sites[] with a per-site verdict + short evidence.
  Be skeptical and exact: this is the mis-attribution safety net.`
}

function recheckPrompt(safe, v) {
  return `A first pass flagged a manuscript citation as possibly unsupported (status=${v.status}).
Your job is ADVERSARIAL: find the STRONGEST evidence the paper DOES support the claim; if you
cannot, confirm the problem with a quote. Do not anchor on the first verdict.

Read scripts/bibaudit/cache/worklist/${safe}.json for the key, tier, sites[] and claims.
Then read the deepest available source: tier A -> the FULL body of literature/<eprint>/*.tex
(grep for the specific term/formula/number); tier B/C -> scripts/bibaudit/cache/abstracts/<ident>.json.
No network. Return the same structured verdict; set status to your honest re-assessment
(SUPPORTS if you found real support; keep MISMATCH/WEAK with a quote if not; UNVERIFIABLE if no
deep source exists to settle a SPECIFIC claim).`
}

// 1. load the list of per-key worklist basenames (one agent reads the dir; avoids huge args)
phase('Load')
const listed = await agent(
  'Read the directory scripts/bibaudit/cache/worklist/ (use ls or Glob). Return {keys: [...]} where '
  + 'each element is a .json filename WITHOUT the .json extension. Return ALL of them (expect ~193).',
  { label: 'load-worklist', phase: 'Load', schema: { type: 'object', required: ['keys'], properties: { keys: { type: 'array', items: { type: 'string' } } } } },
)
const SAFE = (listed && listed.keys) || []
log(`worklist: ${SAFE.length} keys to verify`)

// 2+3. verify each key (parallel), then adversarially re-check weak/mismatch verdicts
const results = await pipeline(
  SAFE,
  (safe) => agent(verifyPrompt(safe), { label: `verify:${safe}`, phase: 'Verify', schema: VERDICT, agentType: 'Explore' }),
  (v, safe) => {
    if (v && ['MISMATCH', 'WEAK', 'PARTIAL'].includes(v.status)) {
      return agent(recheckPrompt(safe, v), { label: `recheck:${safe}`, phase: 'Recheck', schema: VERDICT, agentType: 'Explore' })
        .then((r) => ({ ...v, recheck: r }))
        .catch(() => v)
    }
    return v
  },
)

const clean = results.filter(Boolean)
const tally = {}
for (const r of clean) tally[r.status] = (tally[r.status] || 0) + 1
log(`verified ${clean.length}/${SAFE.length}; status tally: ${JSON.stringify(tally)}`)
return { n: clean.length, tally, results: clean }
