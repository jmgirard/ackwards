# CLAUDE.md — ackwards

Operating manual for AI-assisted development of this package. Read
`cairn/DESIGN.md` first and treat it as the **source of truth** for all
design decisions; this file covers *how we work*, not *what we’re
building*. When this file and `cairn/DESIGN.md` disagree,
`cairn/DESIGN.md` wins for design and this file wins for process — and
flag the conflict.

## What this is

`ackwards` is an R package implementing Goldberg’s (2006) bass-ackwards
method and modern descendants (PCA / EFA / ESEM engines) for
hierarchical structural analysis. Extract solutions from 1..k factors,
then characterize the hierarchy via between-level factor-score
correlations. Full rationale, contracts, object spec, and resolved
defaults are in `DESIGN.md`.

**Note:** Forbes (2023) footnote 3 cites this package
(`github.com/jmgirard/ackwards`) as the reference implementation of the
extended bass-ackwards approach. Fidelity to the paper’s algorithm is
the baseline contract for anything Forbes-related; additive enrichments
are acceptable but the default output must reproduce Forbes’s examples
exactly. **This contract is test-backed** (M44 + M53):
`tests/testthat/test-forbes-fidelity.R` reproduces the paper’s three
simulation studies *and* her 155-variable AMH applied example
(`k_max = 10`) against expected values computed with Forbes’s own
reference implementation (OSF `pcwm8`; provenance in each fixture). AMH
matched to 1.3e-14 across all 45 level-pairs. The AMH matrix ships as a
CC-BY fixture (M53 — Forbes licensed it CC-BY 4.0; declared in
`LICENSE.note`). M53 also reproduced her redundancy chase exactly (54/54
components): her `ChaseCorrPaths()` uses the **direct/skip-level**
correlation to each ancestor level, which `prune("redundant")` now
adopts as its default (`redundancy_criterion = "direct"`) — the pre-M53
adjacent-hop walk diverged on 7/54 AMH components because correlation is
non-transitive (they agree on the shallow simulations). See M53 in
`MILESTONES.md`.

## Invariants — do not violate without flagging

These encode hard-won reasoning from the design phase. Changing them is
a design decision, not a refactor.

1.  **One edge path.** All between-level correlations go through
    [`compute_edges()`](https://jmgirard.github.io/ackwards/reference/compute_edges.md).
    Use the algebra (`W'RW`, standardized) when scoring is linear;
    materialize scores only when nonlinear (EAP) or when the user asks.
    **Always** standardize by real score SDs `sqrt(diag(W'RW))` — never
    assume unit variance (Bartlett/oblique scores are not
    unit-variance).
2.  **Keep the cross-check.** Retain the scores route even where algebra
    is the default, and keep the test asserting they agree within
    tolerance for linear engines.
3.  **Light core, heavy opt-in.** The object always carries
    loadings/variance/fit/weights/edges/ lineage/`R`/meta. `scores`, raw
    `fits`, raw `data` are NULL by default and recomputable.
4.  **Sign alignment anchors to the primary parent**, not “all positive”
    (that’s impossible).
5.  **Lineage lives in edges, never in IDs.** `m{k}f{j}` are stable
    labels; parentage is in the edge structure.
6.  **Loud defaults.** Announce consequential auto-choices via cli
    (e.g., the ordinal-detection warning). Advise loudly; never switch
    basis silently.
7.  **Convergence is data, not an error.** A non-converging level
    warns + is skipped; the object still builds to the deepest converged
    level. Never let one bad level abort the run.
8.  **Oracle-backed numerics.** Every numeric result is verified against
    ≥2 independent oracle *types* (published/closed-form, an independent
    package, or seeded simulation); **no unsourced or unreproducible
    reference value ships**. Committed fixtures carry a structured
    `provenance` attr naming their `data-raw/` generator (guarded by
    `test-oracle-provenance.R`); the full catalogue of every oracle in
    the suite, classified by type, is `cairn/ORACLES.md`. Live
    independent-impl oracles (`psych`/`lavaan`, recomputed at test time)
    are the stronger form — do not freeze them into fixtures unless they
    become expensive/network-bound (M57). *(This is the interim home for
    the oracle principle; fold it in as a numbered DESIGN IP/GP when the
    design-interview pass runs.)*

## Resolved defaults

High-stakes defaults (varimax rotation, `cor = "pearson"`, `tenBerge`
scoring, WLSMV for ordinal ESEM, `k_max` required, `align_signs`,
`keep_scores`/`keep_fits = FALSE`, `redundancy_phi` auto-resolve,
`redundancy_criterion = "direct"`) and their full rationale live in
`cairn/DESIGN.md` §9 + §14. Announce auto-resolved choices loudly
(Invariant 6); never change these silently.

## Dependencies

Dependency rationale and the full Imports/Suggests split live in
`cairn/DESIGN.md` §12. Hard rules: `psych` is the only heavy **Imports**
dep (engine substrate for PCA/EFA + polychoric); everything else is
guarded **Suggests**; **no Rcpp**; **do not add to `Imports` without
flagging**. Current `Imports`: `cli`, `generics`, `psych`, `rlang`,
`stats`, `utils`.

## Dev workflow

R \>= 4.1 (native pipe `|>` and `\(x)` lambdas allowed). Standard
devtools loop:

``` r

devtools::load_all()      # load for interactive testing
devtools::document()      # regenerate roxygen docs + NAMESPACE after any roxygen change
devtools::test()          # run testthat suite
devtools::check()         # full R CMD check
styler::style_pkg()       # format
lintr::lint_package()     # lint
```

**Efficiency (don’t re-run the suite needlessly).** The suite is
parallel since M48 (`Config/testthat/parallel`, slowest files first):
~27s with `TESTTHAT_CPUS=8`, ~81s serial — **always prefix suite runs
with `TESTTHAT_CPUS=8`** (testthat defaults to 2 workers without it).
`check()` still runs the full suite *and* examples, and
[`covr::package_coverage()`](http://covr.r-lib.org/reference/package_coverage.md)
runs the suite *again* — so `test()` → `check()` → `coverage()` at one
gate executes the suite ~3×. Vignette rebuild is now cheap: seven of the
eight vignettes are **precomputed** (CRAN 2026-07-05 flagged the old
~317s vignette rebuild) — authored as `vignettes/*.Rmd.orig`, knitted
ahead of time into static `*.Rmd` (results + `vignettes/assets/` figures
baked in) by `vignettes/precompute.R`; only `ackwards-interpret.Rmd`
stays live (needs full bfi25 for the IPIP-label lesson). **After editing
any `*.Rmd.orig`, re-run `Rscript vignettes/precompute.R` and commit the
regenerated `*.Rmd` + `vignettes/assets/`** — the `.orig` and the script
are `.Rbuildignore`d, so an un-regenerated edit ships stale output.
Iterate with **targeted** `devtools::test(filter = "<x>")` /
[`testthat::test_file()`](https://testthat.r-lib.org/reference/test_file.html);
run failing tests **once** in a way that shows the details (capture
`res <-` or use a non-silent reporter — never silent-then-rerun); skip
the vignette rebuild during mid-work checks with
`check(vignettes = FALSE)`; and at the final gate run
`Rscript tools/dod-gate.R` — it executes the whole DoD sequence (check →
coverage → style → lint → pkgdown) serially in one process with sensible
`TESTTHAT_CPUS`, and prints/exits on any failure. Never run two
package-touching R processes concurrently. Two transcript-mined
anti-patterns to avoid (M48): a bare `devtools::load_all()` in its own
`Rscript` call does nothing persistent (each `Rscript` is a fresh
process; `test()`/`check()` load the package themselves), and
`cd <repo> && …` compound commands trigger avoidable permission prompts
— use absolute paths. In tests, reuse the `cached()` fit memo
(`tests/testthat/helper-data.R`) instead of refitting identical
[`ackwards()`](https://jmgirard.github.io/ackwards/reference/ackwards.md)
objects — but never for reproducibility/serial-vs-parallel oracles (a
cached second call asserts nothing) or fits wrapped in condition
expectations; and treat cached objects as read-only — rebinding a
returned value is copy-on-modify safe, but environment-bearing
components (e.g. lavaan fits under `keep_fits = TRUE`) are shared by
reference across cache hits.

Scaffolding helpers: `usethis::use_r()`, `use_test()`, `use_package()`,
`use_testthat(3)`, `use_github_action("check-standard")`. Use testthat
3e, roxygen2 for all exported functions (document the *why* of each
default, runnable `@examples`, `@seealso` cross-links).

## Definition of done (every change)

- Tests written/updated and passing; new behavior has a test.
- `devtools::document()` run if roxygen changed; NAMESPACE committed.
- **New/removed exported object → pkgdown reference updated.** Whenever
  `NAMESPACE` gains or loses an `export()`/`S3method()`, or a dataset is
  added under `data/`, update the `reference:` list in `_pkgdown.yml` to
  match, then verify with
  [`pkgdown::check_pkgdown()`](https://pkgdown.r-lib.org/reference/check_pkgdown.html).
  This is the exact check the pkgdown GitHub Action runs — an exported
  topic missing from the reference index fails that workflow (and only
  that workflow; local `R CMD check` won’t catch it). Run it at the DoD
  gate unconditionally (it’s sub-second and gated on
  `rlang::is_installed("pkgdown")`); it costs nothing and closes a
  recurring CI-break.
- `devtools::check()` clean (run **once** at the gate — it subsumes
  `devtools::test()`; see the efficiency note under *Dev workflow*).
  Coverage checked once, not per sub-step.
- Styled and linted.
- The whole gate above is one command: `Rscript tools/dod-gate.R` (M48)
  — check → coverage → style → lint → pkgdown, serial, one process,
  non-zero exit on any failure.
- Public-facing change reflected in NEWS.md and (if user-visible) the
  relevant `@examples`/vignette.
- **Milestone/hotfix tracking is cairn’s now** (see the *Project
  tracking (cairn)* section below): status lives in `cairn/ROADMAP.md`,
  task/work-log detail in the milestone file, decisions in
  `cairn/DECISIONS.md`, finished-milestone history in
  `cairn/milestones/archive/` + git. Pre-cairn milestone history
  (M1–M53) is entombed in `cairn/legacy/MILESTONES.md`. Do not
  re-introduce a milestone index or status slot in this file.

## Git

cairn’s git & merge-approval model (`skills/shared/tracking-rules.md`,
loaded when a cairn skill fires) governs branching, PRs, and the
review/merge-approval gate. **This repo’s default branch is `master`,
not `main`** — read every “main” in cairn’s rules and CLAUDE section as
`master` here. Repo-specific facts cairn does not know:

- **Branch protection & CI gating (owner decision, 2026-07-01) —
  standing override of cairn’s default gate:** `master` is deliberately
  **not** branch-protected. Non-release milestones have historically
  **merged on local-green** (`devtools::check()` 0/0/0 + tests +
  coverage + style/lint) with CI treated as an after-the-fact signal —
  *not* gated on the ~8–15 min CI matrix (pure latency for a solo
  pre-CRAN repo where local `check()` already ran). Don’t use
  `gh pr merge --auto` (silently no-ops without required checks on an
  unprotected branch) and don’t synchronously watch or
  background-poll CI. This **consciously overrides cairn’s default
  “never merge pending CI”** (tracking-rules) and holds until the
  package has real users/collaborators — at which point reconsider
  branch protection. Revisit deliberately; don’t let it lapse silently.
- **Release / CRAN-submission milestones are the exception to the
  above**: wait for the *full* green CI matrix (macOS/Windows/Ubuntu ×
  release/devel/oldrel) before merging, because CRAN runs exactly that
  matrix and will reject platform failures the local macOS `check()`
  can’t see.
- **Do not touch** the `legacy` branch or the `v0-legacy` tag — they
  preserve the pre-AI code.
- Small, focused commits with imperative messages. Don’t force-push
  `master`. Don’t commit data, credentials, or large binaries.

## Ask-first / guardrails

- Ambiguity in `DESIGN.md` → ask; don’t invent a design decision.
- Adding an `Imports` dependency, introducing Rcpp, or changing a
  resolved default → flag for approval first.
- Touching `git history`, `legacy`, tags, or anything destructive →
  confirm first.
- Prefer wrapping established engines (`psych`, `lavaan`, `GPArotation`)
  over reimplementing numerics.

## Out of scope

EAP scoring (declined M28), oblique rotation, and higher-order SEM /
Schmid-Leiman are out of scope — rationale in `cairn/DESIGN.md` §2 +
§14. (Bootstrap skip-level edge CIs shipped M47 as
[`boot_edges()`](https://jmgirard.github.io/ackwards/reference/boot_edges.md);
no longer deferred.)

## Project tracking (cairn)

This repo uses the cairn plugin. **Before acting on any request,
classify it and route** — the tracking rulebook only loads once a cairn
skill fires, so starting work in plain conversation silently bypasses
the work tiers and the git model. Classify first:

- **Trivial** (no runtime surface — typo, comment, tracking edit):
  commit directly to master.
- **User-visible bug**: invoke `/hotfix`.
- **New work, a design decision, or more than one sitting**: invoke
  `/milestone-plan` (then `/milestone-implement` → `/milestone-review`).
- **Status, “what’s next”, or unsure which tier**: invoke `/milestone`.
- **Never implement code on master** outside a milestone/hotfix branch;
  nothing reaches master without the user’s explicit approval at the
  review gate.

Whenever the request is anything but trivial, invoke the skill *first*
so the full rulebook (the plugin’s `skills/shared/tracking-rules.md`)
and its conduct load — do not reconstruct the rules here from memory.
All project state lives under `cairn/` (**Architecture → DESIGN · Status
→ ROADMAP · Tasks → milestone files · Decisions → DECISIONS · History →
archive + git**); never record status or TODOs in this file. Claude’s
persistent memory never holds project state; `cairn/` files win any
conflict.
