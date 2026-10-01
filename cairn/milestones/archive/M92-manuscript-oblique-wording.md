# M92: Manuscript wording for the oblique rotation option

**Status:** done (2026-09-30, PR #101 https://github.com/jmgirard/ackwards/pull/101)

**Goal:** The manuscript describes the shipped `rotation` option: varimax by default, oblique rotations as an option, and edge algebra that holds without orthogonality.

**Outcome:** Two passages of `manuscript/manuscript.qmd` changed. The methods passage says the `W'RW` closed form is exact for any fixed linear scoring weights, oblique included. It cites Waller's oblique form. It explains the varimax default: with the default scores an edge equals the ancestor's unique contribution, and varimax matches Goldberg (2006) and Forbes (2023). The Discussion scope paragraph says `r` is a total correlation under either rotation and defines `beta`. It states the oblique cost: a primary parent can be a factor that only correlates with the real parent, and the `cut_show` and `prune()` thresholds were set under varimax. The old "orthogonal only" and "correlated factors confound" claims are gone. The Discussion banner comment no longer says "author-owned stub". No package code changed.

**Decisions:** none. The wording follows RR01, D-034, and D-036.

**Review:** all five criteria passed on fresh evidence, before and after the review fixes. The search found no forbidden claim, the passages have no em dash, citations resolve, and the render exited 0 with PDF and docx. Package check was clean. Three-lens fan-out: 15 findings, none failing a criterion. Nine were fixed before merge. Four were gaps: the oblique cost, the varimax rationale, a definition of `beta`, and the parent level. Four were prose faults: the paragraph flow, an unsupported "therefore", "total correlation" reading as oblique-only, and an "X, not Y" sentence. The ninth was the stale banner comment. Four older claims went to one candidate row: "exact" under a polychoric R, the fallback sentence, the Waller credit, and the AI-use dates. Two were rejected as package-help detail. These were GPArotation for some rotations and `beta` as `NA` for a singular level. Merged on local green under the repo's standing non-release rule.
