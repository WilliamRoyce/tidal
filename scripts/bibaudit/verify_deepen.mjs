export const meta = {
  name: 'verify-citations-deepen',
  description: 'Task 2 deepen-pass: deep full-text re-verification of citations that were UNVERIFIABLE/WEAK/PARTIAL from the abstract alone (now fetched into literature/)',
  phases: [
    { title: 'Load', detail: 'list deepen worklist files' },
    { title: 'Deep-verify', detail: 'full-text read of each fetched paper' },
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
  return `DEEP full-text re-verification of one manuscript citation. A first (abstract-only) pass
could not settle its specific claims; the paper's full TeX source has now been fetched locally.

STEP 1 — read scripts/bibaudit/cache/worklist_deepen/${fname}
  It has: key, eprint, litpath (a local directory with the paper's full .tex), and sites[]
  (each \\cite site with {file, line, claim, role}).

STEP 2 — DEEP read the full text at litpath:
  ls the litpath dir, Read the main .tex (and grep across *.tex in it) for the SPECIFIC
  formula / number / theorem / attribution each ATTRIBUTION or SPECIFIC site claims. You now
  have the body, not just the abstract — actually find (or fail to find) the claimed content
  and quote it. For EXISTENCE sites a topic match is enough.

STEP 3 — return the structured verdict for key:
  status SUPPORTS if the full text corroborates the claims; PARTIAL if some sites only;
  WEAK if topically related but the specific claim is not established; MISMATCH if the text
  contradicts/does-not-concern the claim (quote it); UNVERIFIABLE only if the body genuinely
  does not address it. Fill depth_reached="deep:<section/line>", source_checked (the file you
  read), confidence, evidence_quote, and per-site verdicts with quotes. Be exact and skeptical.`
}

phase('Load')
const listed = await agent(
  'Read the directory scripts/bibaudit/cache/worklist_deepen/ (ls or Glob). Return {files: [...]} '
  + 'with each .json filename WITHOUT the .json extension.',
  { label: 'load-deepen', phase: 'Load', schema: { type: 'object', required: ['files'], properties: { files: { type: 'array', items: { type: 'string' } } } } },
)
const FILES = (listed && listed.files) || []
log(`deepen worklist: ${FILES.length} keys for full-text re-verification`)

const results = await parallel(
  FILES.map((f) => () =>
    agent(prompt(f), { label: `deepen:${f}`, phase: 'Deep-verify', schema: VERDICT, agentType: 'Explore' })),
)
const clean = results.filter(Boolean)
const tally = {}
for (const r of clean) tally[r.status] = (tally[r.status] || 0) + 1
log(`deep-verified ${clean.length}/${FILES.length}; tally: ${JSON.stringify(tally)}`)
return { n: clean.length, tally, results: clean }
