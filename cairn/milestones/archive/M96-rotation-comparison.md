# M96: Varimax and oblimin side by side on bfi25 in the engines vignette

**Status:** done (2026-10-01, PR #105 https://github.com/jmgirard/ackwards/pull/105)

**Goal:** If users fit an oblique rotation in place of varimax, show them on one real dataset what changes and what stays the same.

**Outcome:** The rotation section of the engines vignette now fits EFA to `na.omit(bfi25)` twice, on the polychoric basis with `k_max = 5`. One fit uses the default varimax, and one uses oblimin with `seed = 1`. They replace the `sim16` oblique example. The section shows five things:
- The oblimin `factor_cor` table, largest .34 for m5f3 and m5f4.
- A one-to-one match by signed Tucker congruence (`psych::factor.congruence()`). The smallest is .94, and level 5 swaps m5f2 and m5f3.
- The two primary-parent trees joined after the match, which are equal.
- The above-cut secondary edges: none for varimax, and three for oblimin with `beta` .016, .0015, and .16.
- Both `autoplot()` figures.

A paragraph points at the oblimin fit's printed advisory, and "How to decide" has its own heading. `tests/testthat/test-vignette-rotation.R` guards the four findings and the two primary `beta` values over 1. A NEWS entry was added. No R, man, or manuscript file changed.

**Decisions:** none cross-cutting. The plan gate kept the comparison in the engines vignette, not a new article. It chose EFA on the polychoric basis over PCA.

**Review:** one pass, and all six criteria passed with fresh evidence. Three lenses gave 14 findings, and none showed a criterion failing. The prior-review lens found no prior-review evidence. Eight were fixed before merge. The .94 match is now a close match under the .95 `prune()` default. A second `beta` over 1 is named. The congruence keeps its sign. The `cut_show` wording names the default. The m4f3 → m5f4 `beta` is explained. The .0015 test pin has a tolerance. The advisory has a pointer, and "How to decide" has a heading. Five were rejected and one was noted. The M87 lesson gained the section-splice rebuild. Merged on local green under the repo's standing non-release rule.
