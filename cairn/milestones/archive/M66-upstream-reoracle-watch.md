# M66: Upstream re-oracle watch (weekly CRAN-current CI) — done 2026-07-17

**Goal.** Mechanize DESIGN §2's track-CRAN-current posture (D-031, IP8/IP9): a weekly CI job
reruns the full check against fresh CRAN psych/lavaan and files a GitHub issue on failure, so an
upstream break (like lavaan 0.7, unnoticed until the 0.1.1 release walk) surfaces on a clock.

**Outcome.** New `.github/workflows/re-oracle.yaml` — weekly cron (Mon 06:00 UTC) +
`workflow_dispatch` with a `dry_fail` boolean. Single ubuntu-latest / R-release job installing
from **current CRAN** (`use-public-rspm: false`, not lagging RSPM binaries). Records psych/lavaan
versions to the log/step-outputs, then `check-r-package`. On any failure a `gh`-CLI step opens —
or comments on, never duplicates — one `upstream-watch`-labeled issue with the versions + run
link. A workflow-level `concurrency: { group: re-oracle }` serializes runs so the issue
read-then-write can't race a duplicate. `.github/`-only; no package surface.

**Verification.** End-to-end on the branch via a temporary push trigger (reverted before review):
green run recorded psych 2.6.5 / lavaan 0.7.2; two forced-fail runs opened then commented on one
issue (#69, no duplicate), closed as evidence. 3-lens review + scorer: one below-threshold finding
(F1, concurrency race, 72/100) fixed at the owner's election at the merge gate (the group above).

**Key decisions.** Notification = `gh` CLI (preinstalled, no third-party action). No matrix /
R-devel (push-PR CI + rhub cover those; this hunts upstream drift). 60-day cron auto-disable
documented in the header, not worked around.

**PR:** [#70](https://github.com/jmgirard/ackwards/pull/70) · squash `3ac3c64`. **Owner
follow-up:** one confirming `workflow_dispatch` now the file is on `master`.
