# Reports for upstream package authors

<!-- cspell:words xPand xMAG xBrauer -->

What is here is written for other people — package authors — and nothing here is filed or sent
by a lane session. The user decides; the drafts only have to be right. The internal evidence
behind each report stays in `../xpand_upstream_issue.md`, `../xmag_upstream_issue.md` and
`../xbrauer_upstream_issue.md`.

| file | what it is | destination | state |
| --- | --- | --- | --- |
| `xmag_issue2_followup.md` | a follow-up comment correcting three errors in the filed issue | [xMAG #2](https://github.com/THelpin/xMAG/issues/2) | **posted 2026-09-18**, in the wording from before the contraction issue existed; a one-line comment linking xBrauer_Bundle #3 is still to post |
| `xbrauer_issue_contraction.md` | a new issue: the contraction and separation regressions | [xBrauer_Bundle #3](https://github.com/THelpin/xBrauer_Bundle/issues/3) | **filed 2026-09-18**, identical to this draft |
| `xpand_email.txt` | plain-text email body (Outlook-safe: no markdown) | `pitrou@iap.fr` — xPand has no issue tracker | ready |
| `xpand_item{1,2,3,4}_*.wls` | the email's attachments, one per item | attached to that email | ready |
| `claims.md` | every claim in every report, how it was checked, the result | for the user, not sent | — |
| `check_snippets.py` | runs every snippet as a reader would and compares its outputs | — | — |

Filed by the user on 2026-09-18: [xBrauer_Bundle #3](https://github.com/THelpin/xBrauer_Bundle/issues/3) (from `xbrauer_issue_contraction.md`, checked before filing), the
follow-up comment on xMAG #2, and, checked after the fact (see `claims.md`):
[xMAG #2](https://github.com/THelpin/xMAG/issues/2) (three errors, corrected by the follow-up) and
[xBrauer_Bundle #2](https://github.com/THelpin/xBrauer_Bundle/issues/2) (correct as filed).

## The rule for anything in this directory

A report is ready only when every row of its table in `claims.md` reads **verified** and
`python3 check_snippets.py` prints `all as expected`. The checker runs markdown snippets line by
line in one fresh kernel, as a reader working down the page would, and compares every
`(* expected *)` comment with what actually comes back; attachments must print only `ok` and
`reproduced`. It was confirmed to reject the snippet filed in xMAG #2, which is the check that
would have stopped that error before filing.

Formats: GitHub text is raw markdown, handed over in a code block so the source is copied, not
the rendering. Email is plain text with code in attached scripts, never in the body.
