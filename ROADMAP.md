# ROADMAP.md — planned milestones (none pending)

Forward-looking counterpart to [`MILESTONES.md`](MILESTONES.md). `MILESTONES.md` is the source of
truth for **completed** milestones; this file captures the **intent and source notes** for the
*not-yet-built* milestones, so their context survives across planning sessions.

**There are currently no pending milestones.** The M31–M40 arc — the correctness sweep (M31–M32),
the documentation/UX epic (M33–M39), and the deferred code/viz asks (M40) — has shipped in full. The
last pending item, M40's three spin-off asks, was completed on **2026-07-01**: the ordinal
correlation-comparison viz (dodged bar chart) and the Forbes pruned-level italic axis labels shipped
as code/viz; the ordinal `categorical` convenience flag was **declined** as redundant (`cor =
"polychoric"` already auto-selects WLSMV). See M40's `MILESTONES.md` entry.

For deferred / out-of-scope work that is *not* scheduled, see `DESIGN.md` §14 ("Decisions resolved &
remaining") and the "Out of scope for now" list in `CLAUDE.md` — e.g. bootstrap CIs on skip-level
edges, EAP scoring for ordinal ESEM, oblique rotation, higher-order SEM. None of these is a planned
milestone; each would need its own scoping discussion before a `/plan-milestone` run.

## Provenance

The M31–M40 arc was decomposed from a page-by-page pkgdown-site review the owner did on
**2026-06-30**: ~90 doc/function notes across the README and eight vignettes, grouped into milestones
correctness-first. Origin session transcript: `9a5dc6bd-40ca-4f66-9f11-5bf0c0a4e19a.jsonl` (under
`~/.claude/projects/-Users-jmgirard-GitHub-ackwards/`). The roadmap was first committed with M31's
`## Current focus` update (commit `28513da` on the `m31-correctness-sweep` branch). Review notes and
banked decisions were retired into each milestone's `MILESTONES.md` entry as it shipped, per the
maintenance rule below.

**How to maintain this file:** when a future milestone is scoped, capture its intent and source notes
here (raw review questions + banked decisions — inputs to planning, not a finished spec); when it
ships, delete its section and rely on its `MILESTONES.md` entry. Keep this file scoped to *pending*
work only.
