# ROADMAP.md — planned milestones (M37–M38)

Forward-looking counterpart to
[`MILESTONES.md`](https://jmgirard.github.io/ackwards/MILESTONES.md).
`MILESTONES.md` is the source of truth for **completed** milestones;
this file captures the **intent and source notes** for the
*not-yet-built* milestones in the M31–M38 documentation/UX epic, so
their context survives across planning sessions.

`sim16` (the M33 dataset — 1000×16 continuous, known 1→2→4 hierarchy)
has shipped; see its `MILESTONES.md` entry and
[`?sim16`](https://jmgirard.github.io/ackwards/reference/sim16.md).
M37/M38’s dependency on M33 is satisfied.

M34 (the
[`prune()`](https://jmgirard.github.io/ackwards/reference/prune.md)
verb) has shipped; see its `MILESTONES.md` entry. M38’s forbes-vignette
dependency on “final
[`prune()`](https://jmgirard.github.io/ackwards/reference/prune.md)
API + artifact naming” is satisfied — canonical rule name is
`"artifact"` (`"artefact"` accepted as an alias).

M35 (autoplot & visualization) has shipped; see its `MILESTONES.md`
entry. The M38 forbes-vignette display note about pruned-level labels
can build on M35’s orientation-aware
[`autoplot()`](https://jmgirard.github.io/ackwards/reference/autoplot.md);
the sign-propagation fix means primary-parent edges are now always
non-negative.

**Read this before running `/plan-milestone N` for any milestone
below.** The one-line index in `CLAUDE.md` (“Remaining milestones in the
epic”) is a lossy pointer; the driving rationale and the raw review
notes that shaped each milestone live here.

## Provenance

This epic was decomposed from a page-by-page pkgdown-site review the
owner did on **2026-06-30**: ~90 doc/function notes across the README
and eight vignettes, grouped into milestones correctness-first. Origin
session transcript: `9a5dc6bd-40ca-4f66-9f11-5bf0c0a4e19a.jsonl` (under
`~/.claude/projects/-Users-jmgirard-GitHub-ackwards/`). The roadmap was
first committed with M31’s `## Current focus` update (commit `28513da`
on the `m31-correctness-sweep` branch). M31 and M32 have since shipped;
several review notes were resolved there and are marked **(✔ resolved in
M31/M32)** below so future planning doesn’t re-litigate them.

**How to use / maintain this file:** - The notes below are the owner’s
*raw review questions* plus the *banked decisions* from the
decomposition discussion. They are inputs to planning, not a finished
spec — the concrete file-level plan and acceptance criteria are produced
at `/plan-milestone` time. - When a milestone ships, delete its section
here and rely on its `MILESTONES.md` entry. Keep this file scoped to
*pending* work only.

------------------------------------------------------------------------

## M37 — Engines vignette

- **Type:** doc-heavy · **Depends on:** M31, M32, M33 · **Status:**
  pending.

**Banked decisions (from the M33 post-review):** - **Frame `sim16` as
the *idealized* case, not the norm.** Its
[`suggest_k()`](https://jmgirard.github.io/ackwards/reference/suggest_k.md)
criteria all agree on `k = 4` because the planted signal is
clean/strong; real data doesn’t. Whenever a vignette shows `sim16`’s
consensus, contrast it explicitly with `bfi25`, where the same six
criteria span `k = 4`–`6` (PA-PC 5, PA-FA 6, MAP 5, VSS-1 4, VSS-2 5, CD
6). Teaching arc: `sim16` = “watch the method recover a known
structure”; `bfi25` = “reason under criterion disagreement.” (The
[`?sim16`](https://jmgirard.github.io/ackwards/reference/sim16.md) doc
already states this; carry it into the prose.) Decision: keep `sim16`
clean rather than muddying it for realism — `bfi25` already carries the
realism lesson. - **Deferred M33 nice-to-have:** add an Invariant-2
algebra-vs-scores cross-check on `sim16` (a clean continuous fixture) if
the engines vignette materializes scores — folded here from the M33
review.

**Source notes (engines vignette):** - “the ‘what it models’ row should
be the same for efa and esem, no? … aim for parallelism (either common
variance for both or common variance via psych/lavaan if the package
underlying it matters).” - “instead of ‘chi’ should this say chi-squared
in html symbols?” - “should we add rows to the table about which
correlation types it supports or what estimation methods it
uses/supports? i would want to remove the ‘WLSMV for ordinal’
parenthetical from the loadings SEs row, that is confusing.” - “is it
true that esem is better than efa for ordinal data?” - “should we add
that esem unlocks ML/MLR/FIML for gaussian data? the jump straight to
WLSMV makes it sound like esem is just for ordinal data, which i dont
think we want to imply.” - “but does bad fit at k = 3 mean that we
shouldnt trust the edges connecting k=2 to k=3?” - “im not really sure
how to interpret the autoplot(x_sem), especially with only k=2 and k=3
being included. since we know k~5 is where we should go to, why isnt
that included here?” - “this vignette is very long. should we shorten or
split it somehow?” - “is plan() not re-exported by future.apply?” /
“would it be more readable to just load library(future) and not use
colon notation here?” - “if we dont end up removing cutoffs, this
vignette should cite hu and bentler. do we need any other citations for
this vignette?” *(cutoffs arg was kept — M32 — so the Hu & Bentler
citation applies.)* - “if a user used psych::corFiml() and passed that
to ackwards with engine = pca or efa, would that be a way to smuggle
FIML into those engines?” *(open Q — relates to correlation-matrix
input, M22.)* - **(✔ resolved in M31/M32):** index column names
(`chi`/`dof`/`p_value` → `statistic`); ESEM p-values NA; `glance` BIC
NA; the `*_meets` columns; whether to keep the `cutoffs` arg; and
`cor="polychoric"` + non-WLSMV estimator behavior (now guarded). Don’t
reopen these.

------------------------------------------------------------------------

## M38 — Narrative & remaining prose (intro, suggest_k, ordinal, forbes, README)

- **Type:** doc · **Depends on:** M33, M34 (both shipped) · **Status:**
  pending.

**Banked decision (from the M33 post-review):** the suggest_k vignette
is the natural home for the `sim16` (idealized, criteria converge on 4)
vs. `bfi25` (realistic, criteria span 4–6) contrast — see the M37 note
above for the numbers and framing. Don’t present `sim16`’s clean
consensus as typical.

**README prose:** - “should the step 3 explanation explain why m2f2 -\>
m3f2 gets a negative sign arrow?” - “should the step 4 output avoid
indexing to stay simpler? can we present this more straightforwardly?” -
“i dont think we need to encourage people to cite waller (2007) here.
they probably wont even know if they are using that algebra unless they
dive deep into our docs.”

**Intro vignette prose:** - “the end of the note referencing efa vs esem
seems out of place here.” - “should we explicitly say `print(sk)`? the
sk is lost in the `#>` preamble.” / “step 3: should we use print(x)
instead of x again?” - “is suggest_k suggesting k or max_k? seems
confusing.” - “maybe this would be cleaner to read if cut = 0.5.” - “we
should include a note about why se and cis are NA here and that they
could be estimated had we used the esem engine.” *(the machine behavior
— SE/CI as NA + signal — was settled in M31; this is the prose note.)* -
“in step 6, is the standardizing warning relevant or distracting
here?” - “i dont think we need to show the keep_scores = TRUE part
here.” *(M36 shipped `augment(append = FALSE)` and already replaced the
intro’s base-R score indexing with it; the remaining ask — whether to
drop the `keep_scores = TRUE` demo chunk entirely — is still M38’s.)* -
“should we explain briefly why we use orthogonal rotations within each
level in bassackwards? and specifically why we use varimax (and how that
is equivalent to crawfer in case people see that reference in kim and
eaton and get confused … like i did).” *(see the rotation-rationale
memo.)* - **(✔ resolved in M31/M32):** “no warning this time” mismatch;
drifted hardcoded percentages (now dynamic); 5-vs-6 criteria count;
`variance_pct` 0–100 → 0–1. Prose should reflect the fixes.

**suggest_k vignette prose:** - “VSS limitation has an extra space in
‘cross- loadings’.” - “use print(sk) here.” / “what does the check vs
star convention mean?” - **(✔ resolved in M31/M32):** 5-vs-6 criteria;
the [`set.seed()`](https://rdrr.io/r/base/Random.html)/`fa.parallel`
reproducibility question (doc confirmed correct — not a bug to report);
`k_max` meaning across functions (kept in both, roxygen-disambiguated).

**Ordinal vignette prose:** - “i dont like the em dash and examples
here, makes it seem like likert items with 3 or 8 wouldnt apply. what
about binary items with 2?” - “i dont love showing the subsetting and
rounding code here and the matrices are hard to compare visually. could
we hide the code and maybe make a dodged bar chart showing the
lower-triangle correlations mapping edge to fill? or would a gt showing
the lower-triangle correlations in long format but with cols for pearson
and polychoric (and maybe delta) be best?” - “does WLSMV really use the
polychoric matrix as its base? would it be better to have a binary flag
of categorical to switch from pearson to polychoric or MLR to WLSMV?” -
“do we need to estimate tetrachoric separately or does polychoric
collapse into that perfectly?” - “is this note saying we can or cannot
trust the scores in downstream analysis when data are ordinal?” -
“references should be alphabetical in APA style.” - **(✔ resolved in
M31):** missing p-values in the table.

**Forbes vignette prose:** - “i dont like how it looks when a heading
starts with a verbatim span - can we reorganize to avoid that?” -
“should a pruned level (level 4 in this example) still have its level
label on the left? maybe include but put its level label in parentheses
or make it italic to denote its status?” *(display detail — coordinate
with M35; M34 did not touch rendering.)* - “in these gt tables, the true
values need to be colored or something - too hard to visually find.” -
“the code is hidden in the redundancy_r chunk, only the output is
displayed in a massive block. was that intentional? looks confusing.” -
“lorenzo-seva citation mentioned but not realized at the end.” -
*Depended on M33 (guaranteed Tucker finding) and M34 (final
[`prune()`](https://jmgirard.github.io/ackwards/reference/prune.md)
API + artifact naming) — both shipped; the vignette’s code chunks were
mechanically updated to the piped
[`prune()`](https://jmgirard.github.io/ackwards/reference/prune.md) API
in M34 itself (see its `MILESTONES.md` entry), but the deeper
prose/formatting notes above are still M38’s to resolve.*
