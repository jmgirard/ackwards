# ROADMAP.md — planned milestones (M34–M38)

Forward-looking counterpart to [`MILESTONES.md`](MILESTONES.md). `MILESTONES.md` is the source of
truth for **completed** milestones; this file captures the **intent and source notes** for the
*not-yet-built* milestones in the M31–M38 documentation/UX epic, so their context survives across
planning sessions.

`sim16` (the M33 dataset — 1000×16 continuous, known 1→2→4 hierarchy) has shipped; see its
`MILESTONES.md` entry and `?sim16`. M37/M38's dependency on M33 is satisfied.

**Read this before running `/plan-milestone N` for any milestone below.** The one-line index in
`CLAUDE.md` ("Remaining milestones in the epic") is a lossy pointer; the driving rationale and the
raw review notes that shaped each milestone live here.

## Provenance

This epic was decomposed from a page-by-page pkgdown-site review the owner did on **2026-06-30**:
~90 doc/function notes across the README and eight vignettes, grouped into milestones
correctness-first. Origin session transcript: `9a5dc6bd-40ca-4f66-9f11-5bf0c0a4e19a.jsonl` (under
`~/.claude/projects/-Users-jmgirard-GitHub-ackwards/`). The roadmap was first committed with M31's
`## Current focus` update (commit `28513da` on the `m31-correctness-sweep` branch). M31 and M32 have
since shipped; several review notes were resolved there and are marked **(✔ resolved in M31/M32)**
below so future planning doesn't re-litigate them.

**How to use / maintain this file:**
- The notes below are the owner's *raw review questions* plus the *banked decisions* from the
  decomposition discussion. They are inputs to planning, not a finished spec — the concrete
  file-level plan and acceptance criteria are produced at `/plan-milestone` time.
- When a milestone ships, delete its section here and rely on its `MILESTONES.md` entry. Keep this
  file scoped to *pending* work only.

---

## M34 — Pruning verb: extract `prune()` from `ackwards()`

- **Type:** code (breaking, pre-CRAN) · **Depends on:** none · **Blocks:** M38 (forbes vignette)
- **Status:** pending.

**Banked design decisions (from the roadmap discussion — these are settled, not open):**
- Pruning is already **flag-only, never removes** (`prune.R` writes `pruned`/`prune_reason`
  annotations into `x$prune$nodes`, preserving Invariant 5). That makes it architecturally a
  post-hoc annotation pass — the shape of a pipeable verb.
- **Extract a `prune()` verb:** `ackwards(...) |> prune(...)`. The five prune args (`prune`,
  `redundancy_r`, `redundancy_phi`, `min_items`, `orphan_r`) **leave `ackwards()` entirely** and
  live on the new verb. Lets users re-prune with new thresholds *without re-extracting* (re-running
  ESEM is the expensive path M26 optimized).
- **Clean move, no deprecation shim** — pre-CRAN, no users, so the old args don't need to forward.
- **Manual + mixed pruning:** `prune(x, "redundant", manual = c("m4f3","m4f4"))` runs the auto rules
  then unions in user-named nodes flagged `prune_reason = "manual"`. Stays flag-only; `.drop_pruned_nodes()`
  (M8) already handles display. Invariant 5 untouched.
- The **`prune = "artefact"` vs `"tucker"` naming** decision moves *here* (the whole surface is being
  rebuilt), rather than staying in the M32 naming pass. Consider permitting either as aliases.
- **DESIGN.md update required:** §5 and §14 currently frame pruning as *inside* the Forbes extension;
  moving it to a standalone verb is a genuine design change and must be reflected there.

**Source notes (forbes vignette + design fork):**
- "is 'artefact' a better option name for prune than 'tucker'? any value in permitting either as
  equal alternatives?"
- "still not clear what im supposed to do with the artefactual pruning, is this automatic and youre
  just explaining what happens under the hood? let's clarify"
- "where do artefactual factors come from? overextraction? misspecified data?"
- design fork: "wondering if it is better to keep the pruning stuff folded into ackwards() ... or
  have that be a separate follow-up step that you pipe the output of ackwards into. similarly, if
  the pruning decisions are ... user determined, we need some way for them to manually prune factors
  not just have it done automatically (and probably allow a mix of automatic and defined prunes)."

---

## M35 — autoplot & visualization

- **Type:** code + viz/README · **Depends on:** none · **Status:** pending.

**Banked decision:**
- **`ggsave` cannot be re-exported** — doing so would drag `ggplot2` from Suggests into Imports,
  violating the dependency guardrail. Instead, add a docs *section* that calls `ggplot2::ggsave()`.

**Source notes (README + visualization vignette + intro layout):**
- README: "should autoplot legend define what dashed vs solid line means?"
- viz: "why doesnt linetype get a legend and why is m2f2 -> m3f2 red?"
- viz: "i would prefer to use color not colour."
- viz: "mono switching the meaning of dashed linetype seems like a recipe for confusion. is that
  really a good idea? what if negative was one of those double-lines instead of dashed? then maybe
  it could pair with dashed?" *(design Q — linetype encoding of sign)*
- viz: "here you are using show_r before defining it, not good teaching. probably drop that example
  and replace with a text reference and pointer below."
- viz: "heading says annotate with |r| but we are annotating with r (no absolute value transformation)."
- viz: "many readers wont know the R convention of L after integers, so just say r_digits = 1."
- viz: "do we need this diagnostic plot section here? seems redundant with suggest_k vignette unless
  we customize it in new ways."
- viz: "we need a section at the end here about how to export the plots using ggsave, which we might
  want to re-export?" *(see banked decision above)*
- intro: "should we add an arg in autoplot to make it plot left-to-right instead of top-to-bottom?
  might be useful for slideshows or posters that are wider than tall." *(also fixes the intro text
  that says "k=1 to the left / k=5 to the right" while the layout is top-to-bottom.)*

---

## M36 — Interpretation functions (`augment`, `top_items`)

- **Type:** code + interpret vignette · **Depends on:** none · **Status:** pending.

**Code changes (the two functional asks):**
- **`augment()` scores-only option** — optionally return only the factor scores. Also removes the
  base-R indexing (`names(scored)[26:40]`, `[26:40]` for k5) the owner disliked in the intro.
- **`top_items()` enhancements:** (a) if a column carries a variable label attribute (from
  `labelled`/`haven`), display that label instead of the bare colname (e.g. "How much do you like
  parties?" not `bfi1`); (b) a **group-by-item** mode that lists factors per item (inverse of the
  default group-by-factor), to make cross-loadings legible.

**Interpret-vignette prose (rides with M36):**
- "let's adjust cut to 0.5 for output clarity and parsimony ... then no need to duplicate with
  another cut, could just mention in text that it is adjustable."
- "can probably drop the last part of this section about top_items returning an object with $data —
  that's detail better saved for reference docs."
- "probably remove the comment about why polychoric."
- "should we give advice on or examples of naming higher level factors? could reference metatraits
  in big five theory or spectra in hitop."
- "the where to go next section points to visualization but that is not next in the order of the
  menu ... maybe remove that or make it a note rather than a full section."

---

## M37 — Engines vignette

- **Type:** doc-heavy · **Depends on:** M31, M32, M33 · **Status:** pending.

**Source notes (engines vignette):**
- "the 'what it models' row should be the same for efa and esem, no? ... aim for parallelism (either
  common variance for both or common variance via psych/lavaan if the package underlying it matters)."
- "instead of 'chi' should this say chi-squared in html symbols?"
- "should we add rows to the table about which correlation types it supports or what estimation
  methods it uses/supports? i would want to remove the 'WLSMV for ordinal' parenthetical from the
  loadings SEs row, that is confusing."
- "is it true that esem is better than efa for ordinal data?"
- "should we add that esem unlocks ML/MLR/FIML for gaussian data? the jump straight to WLSMV makes
  it sound like esem is just for ordinal data, which i dont think we want to imply."
- "but does bad fit at k = 3 mean that we shouldnt trust the edges connecting k=2 to k=3?"
- "im not really sure how to interpret the autoplot(x_sem), especially with only k=2 and k=3 being
  included. since we know k~5 is where we should go to, why isnt that included here?"
- "this vignette is very long. should we shorten or split it somehow?"
- "is plan() not re-exported by future.apply?" / "would it be more readable to just load
  library(future) and not use colon notation here?"
- "if we dont end up removing cutoffs, this vignette should cite hu and bentler. do we need any other
  citations for this vignette?" *(cutoffs arg was kept — M32 — so the Hu & Bentler citation applies.)*
- "if a user used psych::corFiml() and passed that to ackwards with engine = pca or efa, would that
  be a way to smuggle FIML into those engines?" *(open Q — relates to correlation-matrix input, M22.)*
- **(✔ resolved in M31/M32):** index column names (`chi`/`dof`/`p_value` → `statistic`); ESEM
  p-values NA; `glance` BIC NA; the `*_meets` columns; whether to keep the `cutoffs` arg; and
  `cor="polychoric"` + non-WLSMV estimator behavior (now guarded). Don't reopen these.

---

## M38 — Narrative & remaining prose (intro, suggest_k, ordinal, forbes, README)

- **Type:** doc · **Depends on:** M33, M34 · **Status:** pending.

**README prose:**
- "should the step 3 explanation explain why m2f2 -> m3f2 gets a negative sign arrow?"
- "should the step 4 output avoid indexing to stay simpler? can we present this more straightforwardly?"
- "i dont think we need to encourage people to cite waller (2007) here. they probably wont even know
  if they are using that algebra unless they dive deep into our docs."

**Intro vignette prose:**
- "the end of the note referencing efa vs esem seems out of place here."
- "should we explicitly say `print(sk)`? the sk is lost in the `#>` preamble." / "step 3: should we
  use print(x) instead of x again?"
- "is suggest_k suggesting k or max_k? seems confusing."
- "maybe this would be cleaner to read if cut = 0.5."
- "we should include a note about why se and cis are NA here and that they could be estimated had we
  used the esem engine." *(the machine behavior — SE/CI as NA + signal — was settled in M31; this is
  the prose note.)*
- "in step 6, is the standardizing warning relevant or distracting here?"
- "i dont think we need to show the keep_scores = TRUE part here." *(dovetails with M36's augment
  scores-only, which removes the base-R indexing.)*
- "should we explain briefly why we use orthogonal rotations within each level in bassackwards? and
  specifically why we use varimax (and how that is equivalent to crawfer in case people see that
  reference in kim and eaton and get confused ... like i did)." *(see the rotation-rationale memo.)*
- **(✔ resolved in M31/M32):** "no warning this time" mismatch; drifted hardcoded percentages (now
  dynamic); 5-vs-6 criteria count; `variance_pct` 0–100 → 0–1. Prose should reflect the fixes.

**suggest_k vignette prose:**
- "VSS limitation has an extra space in 'cross- loadings'."
- "use print(sk) here." / "what does the check vs star convention mean?"
- **(✔ resolved in M31/M32):** 5-vs-6 criteria; the `set.seed()`/`fa.parallel` reproducibility
  question (doc confirmed correct — not a bug to report); `k_max` meaning across functions
  (kept in both, roxygen-disambiguated).

**Ordinal vignette prose:**
- "i dont like the em dash and examples here, makes it seem like likert items with 3 or 8 wouldnt
  apply. what about binary items with 2?"
- "i dont love showing the subsetting and rounding code here and the matrices are hard to compare
  visually. could we hide the code and maybe make a dodged bar chart showing the lower-triangle
  correlations mapping edge to fill? or would a gt showing the lower-triangle correlations in long
  format but with cols for pearson and polychoric (and maybe delta) be best?"
- "does WLSMV really use the polychoric matrix as its base? would it be better to have a binary flag
  of categorical to switch from pearson to polychoric or MLR to WLSMV?"
- "do we need to estimate tetrachoric separately or does polychoric collapse into that perfectly?"
- "is this note saying we can or cannot trust the scores in downstream analysis when data are ordinal?"
- "references should be alphabetical in APA style."
- **(✔ resolved in M31):** missing p-values in the table.

**Forbes vignette prose:**
- "i dont like how it looks when a heading starts with a verbatim span - can we reorganize to avoid that?"
- "should a pruned level (level 4 in this example) still have its level label on the left? maybe
  include but put its level label in parentheses or make it italic to denote its status?"
  *(display detail — coordinate with M34/M35.)*
- "in these gt tables, the true values need to be colored or something - too hard to visually find."
- "the code is hidden in the redundancy_r chunk, only the output is displayed in a massive block.
  was that intentional? looks confusing."
- "lorenzo-seva citation mentioned but not realized at the end."
- *Depends on M33 (guaranteed Tucker finding) and M34 (final `prune()` API + artefact naming).*
