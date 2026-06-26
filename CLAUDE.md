# CLAUDE.md — `ackwards`

Operating manual for AI-assisted development of this package. Read `DESIGN.md` (repo root) first
and treat it as the **source of truth** for all design decisions; this file covers *how we work*,
not *what we're building*. When this file and `DESIGN.md` disagree, `DESIGN.md` wins for design and
this file wins for process — and flag the conflict.

## What this is

`ackwards` is an R package implementing Goldberg's (2006) bass-ackwards method and modern
descendants (PCA / EFA / ESEM engines) for hierarchical structural analysis. Extract solutions from
1..k factors, then characterize the hierarchy via between-level factor-score correlations. Full
rationale, contracts, object spec, and resolved defaults are in `DESIGN.md`.

## Current focus — Milestone 1

Build the core path only; do not start later engines until this is green.

**Scope:** package skeleton → PCA engine → algebra `compute_edges()` → result object →
`print` / `tidy` / `glance`.

**Acceptance criteria:**
- `ackwards(psych::bfi[1:25], k = 5, method = "pca")` returns a valid `ackwards` object.
- PCA cross-level correlations **match `psych::bassAckward()`** on the same data within tolerance
  (snapshot test).
- The algebra-vs-scores cross-check test passes for the PCA engine.
- `print`, `tidy(what=...)`, `glance` implemented with cli formatting.
- `devtools::check()` is clean (0 errors, 0 warnings; notes triaged).

If a step needs a design decision not covered in `DESIGN.md`, **stop and ask** rather than guessing.

## Invariants — do not violate without flagging

These encode hard-won reasoning from the design phase. Changing them is a design decision, not a
refactor.

1. **One edge path.** All between-level correlations go through `compute_edges()`. Use the algebra
   (`W'RW`, standardized) when scoring is linear; materialize scores only when nonlinear (EAP) or
   when the user asks. **Always** standardize by real score SDs `sqrt(diag(W'RW))` — never assume
   unit variance (Bartlett/oblique scores are not unit-variance).
2. **Keep the cross-check.** Retain the scores route even where algebra is the default, and keep the
   test asserting they agree within tolerance for linear engines.
3. **Light core, heavy opt-in.** The object always carries loadings/variance/fit/weights/edges/
   lineage/`R`/meta. `scores`, raw `fits`, raw `data` are NULL by default and recomputable.
4. **Sign alignment anchors to the primary parent**, not "all positive" (that's impossible).
5. **Lineage lives in edges, never in IDs.** `m{k}f{j}` are stable labels; parentage is in the edge
   structure.
6. **Loud defaults.** Announce consequential auto-choices via cli (e.g., the ordinal-detection
   warning). Advise loudly; never switch basis silently.
7. **Convergence is data, not an error.** A non-converging level warns + is skipped; the object
   still builds to the deepest converged level. Never let one bad level abort the run.

## Resolved defaults (see `DESIGN.md` §9, §14)

Orthogonal CF (`cfT`, ≈ varimax) rotation · `cor = "pearson"` with ordinal-detection warning ·
`tenBerge` scoring (on the active basis) · Forbes extension **off** · `k` required ·
sign `align = TRUE` · `scores`/`keep_fits` stored = `FALSE`. Don't change these silently.

## Dependencies (see `DESIGN.md` §12)

Keep `Imports` lean: `stats`, `methods`, `cli`, `rlang`. Everything else (`psych`/`GPArotation`,
`lavaan`, `clue`, `ggraph`+friends, `EGAnet`/`paran`, `future`) goes in `Suggests`, gated by
`rlang::check_installed()`. **Do not add to `Imports` without flagging it.** **No Rcpp** — profile
first; the heavy compute already lives in compiled deps (§3).

## Dev workflow

R >= 4.1 (native pipe `|>` and `\(x)` lambdas allowed). Standard devtools loop:

```r
devtools::load_all()      # load for interactive testing
devtools::document()      # regenerate roxygen docs + NAMESPACE after any roxygen change
devtools::test()          # run testthat suite
devtools::check()         # full R CMD check
styler::style_pkg()       # format
lintr::lint_package()     # lint
```

Scaffolding helpers: `usethis::use_r()`, `use_test()`, `use_package()`, `use_testthat(3)`,
`use_github_action("check-standard")`. Use testthat 3e, roxygen2 for all exported functions
(document the *why* of each default, runnable `@examples`, `@seealso` cross-links).

## Definition of done (every change)

- Tests written/updated and passing; new behavior has a test.
- `devtools::document()` run if roxygen changed; NAMESPACE committed.
- `devtools::check()` clean.
- Styled and linted.
- Public-facing change reflected in NEWS.md and (if user-visible) the relevant `@examples`/vignette.

## Git

- Default branch is **`master`**.
- **Do not touch** the `legacy` branch or the `v0-legacy` tag — they preserve the pre-AI code.
- Small, focused commits with imperative messages (e.g., `Add PCA engine and level contract`).
- Don't force-push `master`. Don't commit data, credentials, or large binaries.

## Ask-first / guardrails

- Ambiguity in `DESIGN.md` → ask; don't invent a design decision.
- Adding an `Imports` dependency, introducing Rcpp, or changing a resolved default → flag for
  approval first.
- Touching `git history`, `legacy`, tags, or anything destructive → confirm first.
- Prefer wrapping established engines (`psych`, `lavaan`, `GPArotation`) over reimplementing
  numerics.

## Out of scope for now

ESEM/EFA engines, ordinal/polychoric path, Forbes extension, and visualization
(`ba_layout`/`autoplot`) are **later milestones** (`DESIGN.md` §15). Don't build ahead of Milestone 1
unless asked.
