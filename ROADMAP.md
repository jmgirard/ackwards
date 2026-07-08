# ROADMAP.md — pending work (candidates, deferrals, release tail)

Forward-looking counterpart to
[`MILESTONES.md`](https://jmgirard.github.io/ackwards/MILESTONES.md).
`MILESTONES.md` is the source of truth for **completed** milestones;
this file captures the **intent and source notes** for *not-yet-built*
work, so their context survives across planning sessions. Keep this file
scoped to *pending* work only — when an item ships, delete its section
and rely on its `MILESTONES.md` entry.

**There are currently no scheduled milestones.** The items below are
either owner-side release mechanics or unscheduled ideas that each need
a scoping discussion before a `/plan-milestone` run. The completed M41
review’s reference results (verified-clean list, findings ledger,
defaults audit, arbitrary-constant inventory) previously lived here;
they were folded into the M41 entry in `MILESTONES.md` (audit appendix)
when the M41→M44 arc closed on 2026-07-03.

------------------------------------------------------------------------

## Owner-only release tail (from M49, still pending)

Not a coding milestone — maintainer actions to finish the 0.1.0 CRAN
submission:

- win-builder / R-hub remote checks, then the interactive
  `devtools::submit_cran()` (emails the maintainer to confirm).
- Until CRAN **accepts** 0.1.0, the `install.packages("ackwards")` line
  in the README points at a package not yet on CRAN — update the release
  date and any “on CRAN” phrasing when it lands.
- If CRAN **bounces** 0.1.0, a patch branches from the `v0.1.0` tag (the
  `0.1.0.9000` dev bump on `master` does not block that).

------------------------------------------------------------------------

## Unscheduled ideas (each needs a scoping discussion before /plan-milestone)

- **[`comparability()`](https://jmgirard.github.io/ackwards/reference/comparability.md)
  and
  [`boot_edges()`](https://jmgirard.github.io/ackwards/reference/boot_edges.md)
  engine/basis extensions** (deferred from M46/M47, DESIGN
  §14.35/§14.36) — **feasible but demand-gated; keep off the schedule
  until asked.** Feasibility verdict (M49 Phase A):
  - **ESEM
    [`comparability()`](https://jmgirard.github.io/ackwards/reference/comparability.md)**
    is the plausible one: `2 * n_splits` (~20) lavaan hierarchies per
    call, minutes with the M26 cached-sample-stats + `future.apply`
    machinery; the matching / anchoring code is already engine-agnostic.
    Main work is per-half convergence handling (Invariant 7).
  - **ESEM
    [`boot_edges()`](https://jmgirard.github.io/ackwards/reference/boot_edges.md)**
    is painful: `n_boot` (~1000) × (k_max − 1) WLSMV fits — tens of
    minutes to hours even parallel; users would have to cut `n_boot`,
    undermining the percentile CIs, and resample non-convergence /
    Heywood cases inflate the dropped-replicate count.
  - **Polychoric basis** (both verbs) is the least realistic — not just
    slow: a bootstrap resample can drop an entire response category,
    changing the threshold structure mid-replicate (a drop-or-merge
    policy question — a genuine statistical wrinkle, not just
    performance). Mirror
    [`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)’s
    Pearson/Spearman-only scope; screen on Pearson, fit the final model
    polychoric.

  Net: only
  ESEM-[`comparability()`](https://jmgirard.github.io/ackwards/reference/comparability.md)
  is worth reconsidering if demand appears; the EFA-on-Pearson screen
  covers the common case for each verb.

------------------------------------------------------------------------

**How to maintain this file:** when a future milestone is scoped,
capture its intent and source notes here (raw review questions + banked
decisions — inputs to planning, not a finished spec); when it ships,
delete its section and rely on its `MILESTONES.md` entry. Completed
review results belong in the reviewing milestone’s `MILESTONES.md`
entry, not here. Keep this file scoped to *pending* work only.
