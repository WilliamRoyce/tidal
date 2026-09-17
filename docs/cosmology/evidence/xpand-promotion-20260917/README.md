# xPand promotion — the gate, watched failing (2026-09-17)

D-C (O1′) makes xPand part of the certified configuration, so `verify-wolfram-setup.sh` gained
its own xPand check. A gate is not real until it has been watched failing, so it was:

| log | state of the userbase | expected | observed |
| --- | --- | --- | --- |
| `pass.log` | certified (xPand 0.4.4 present) | exit 0, PASS | `xPand 0.4.4 found (certified)`, `All checks passed!`, exit 0 |
| `missing.log` | `Applications/xAct/xPand/` moved aside | exit 2, DEGRADED | `xPand not installed … FRW derivations cannot run`, exit 2 |
| `wrong-version.log` | restored, `$Version` edited to `0.4.3` | exit 2, DEGRADED | `xPand 0.4.3 is not the certified 0.4.4: runs are possible, not certifiable`, exit 2 |
| `re-certified.log` | restored | exit 0, PASS | `All checks passed!`, exit 0 (one informational xPerm-ldd warning, as always) |

**The tree was restored byte-exactly**, verified against R-C's own manifest rather than by eye:
all 26 files under `Applications/xAct/xPand/` match `third_party/perturbations_runs/manifest/
after.new.sha256` (`xPand.m` = `aa374fdc…`).

The check reads the version from `xPand.m`'s header, so it needs no kernel and cannot be
confused by a loaded session — the same tactic as the four-package xAct fingerprint, which it
deliberately does **not** join: many documents quote that string as the xAct bundle.
