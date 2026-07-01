# ROADMAP.md — planned milestones (M40)

Forward-looking counterpart to [`MILESTONES.md`](MILESTONES.md). `MILESTONES.md` is the source of
truth for **completed** milestones; this file captures the **intent and source notes** for the
*not-yet-built* milestones, so their context survives across planning sessions.

> **Renumbering note (2026-07-01):** M38 was inserted as a new *code* milestone
> (`missing = "fiml"` for PCA/EFA); the former "M38 — Narrative & remaining prose" was renumbered to
> **M39**. The M31–M39 documentation/UX epic is now **complete** — M39 (narrative & prose) shipped;
> its section was removed from this file per the maintenance rule (see its `MILESTONES.md` entry).
> **M40** was spun off *from* M39's planning: three code/visualization asks that surfaced in the
> M39 review were deliberately excluded from that doc-only milestone and parked here. **M40 is the
> only pending milestone left.**

M35 (autoplot & visualization) has shipped; see its `MILESTONES.md` entry. M40's pruned-level-label
styling (item 3 below) builds on M35's orientation-aware `autoplot()` and its `compress_levels`
label logic; note the M35 sign-propagation fix means primary-parent edges are now always
non-negative (this is why M39 dropped the README/intro red-arrow explanation — the example plot has
no red arrow).

**Read this before running `/plan-milestone 40`.** The one-line index in `CLAUDE.md` is a lossy
pointer; the driving rationale and the raw review notes that shaped M40 live here.

## Provenance

This epic was decomposed from a page-by-page pkgdown-site review the owner did on **2026-06-30**:
~90 doc/function notes across the README and eight vignettes, grouped into milestones
correctness-first. Origin session transcript: `9a5dc6bd-40ca-4f66-9f11-5bf0c0a4e19a.jsonl` (under
`~/.claude/projects/-Users-jmgirard-GitHub-ackwards/`). The roadmap was first committed with M31's
`## Current focus` update (commit `28513da` on the `m31-correctness-sweep` branch). The M31–M39 epic
has since shipped in full; only the M40 spin-off remains below. Review notes that were resolved along
the way are recorded in the relevant milestone's `MILESTONES.md` entry, not here.

**How to use / maintain this file:**
- The notes below are the owner's *raw review questions* plus the *banked decisions* from the
  decomposition discussion. They are inputs to planning, not a finished spec — the concrete
  file-level plan and acceptance criteria are produced at `/plan-milestone` time.
- When a milestone ships, delete its section here and rely on its `MILESTONES.md` entry. Keep this
  file scoped to *pending* work only. (The entire M31–M39 epic has shipped; its sections were
  removed accordingly. Only M40 — carved out of M39's planning — remains.)

---

## M40 — Deferred code/visualization asks (spun off from M39)

- **Type:** code + viz · **Depends on:** M22 (`cor` arg), M35 (autoplot rendering) · **Status:**
  pending, unplanned. Run `/plan-milestone 40` before starting.

M39 was scoped **doc-only** (owner decision, 2026-07-01). Three asks from the same pkgdown review
imply real code or non-trivial in-vignette visualization, so they were carved out here rather than
smuggled into a prose milestone. Each is independent; M40 can ship any subset.

**1. Ordinal `categorical` convenience flag (API — needs owner sign-off first).**
- Source note: "would it be better to have a binary flag of `categorical` to switch from pearson to
  polychoric or MLR to WLSMV?"
- The idea: a single `categorical = TRUE/FALSE` argument on `ackwards()` that flips the correlation
  basis (pearson → polychoric) *and* the ESEM estimator (ML/MLR → WLSMV) together, instead of the
  user setting `cor` (and, for ESEM, letting the estimator auto-switch) themselves.
- **Open design question (do not decide unilaterally at plan time — this is an API/defaults change,
  §9/§14 territory):** the package already has explicit `cor` and (auto-resolved) estimator
  handling. A `categorical` flag partially *duplicates* `cor = "polychoric"` and could create two
  ways to say the same thing (what wins if `categorical = TRUE` and `cor = "pearson"` are both
  passed?). Decide whether it's a genuine ergonomic win or redundant surface area before building.
  Likely wants an `AskUserQuestion` at plan time.

**2. Ordinal correlation-comparison visualization (in-vignette viz).**
- Source note: "i dont love showing the subsetting and rounding code here and the matrices are hard
  to compare visually. could we hide the code and maybe make a dodged bar chart showing the
  lower-triangle correlations mapping edge to fill? or would a gt showing the lower-triangle
  correlations in long format but with cols for pearson and polychoric (and maybe delta) be best?"
- Target: the two raw `round(x$r[1:5, 1:5], 2)` matrix chunks in `vignettes/ackwards-ordinal.Rmd`
  (the "Correlation matrices" subsection). M39 left them untouched.
- Two candidate designs floated by the owner: (a) a **dodged bar chart** of the lower-triangle
  item-pair correlations with fill = basis (pearson vs polychoric); (b) a **gt long-format** table,
  one row per item pair, columns pearson / polychoric / Δ. Pick one at plan time (the gt Δ-table
  mirrors the loadings/edges comparison tables already in that vignette, so it's the lower-risk,
  more consistent option; the bar chart is more visual). Viz-only — `ggplot2`/`gt` are already in
  Suggests; **no package-code or dependency change** required for this item.

**3. Forbes pruned-level axis-label styling (`autoplot()` code).**
- Source note: "should a pruned level (level 4 in this example) still have its level label on the
  left? maybe include but put its level label in parentheses or make it italic to denote its status?"
- Requires touching `autoplot.ackwards()` label rendering (M35 territory), so it is code, not prose —
  hence deferred out of doc-only M39. When `drop_pruned = TRUE` / `compress_levels`, mark a fully
  pruned level's axis label to denote its status (parenthesized or italic). Coordinate with the
  existing orientation-aware layout and the `compress_levels` label logic already in `autoplot.R`.

**Already resolved in M39 (do not reopen):** every prose/formatting ask from the intro, suggest_k,
ordinal, forbes, and README review — including the `sim16`-vs-`bfi25` contrast (now in the suggest_k
vignette), the Waller-citation nudge removal, the tetrachoric/WLSMV clarifications, the verbatim-span
headings, the structural-table gt highlighting, the visible `redundancy_r` chunk, and the
Lorenzo-Seva reference. See M39's `MILESTONES.md` entry.
