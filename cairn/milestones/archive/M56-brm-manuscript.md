# M56: BRM manuscript (reproducible Quarto scaffold + worked example) — done 2026-07-12

**Goal.** Stand up a reproducible, in-repo Quarto/apaquarto manuscript for a *Behavior Research
Methods* tool paper on `ackwards`, rendering to PDF + Word, with worked examples computed live from
the package and the verifiable prose seeded (scholarly framing left to the author).

**Outcome.** `manuscript/` (`.Rbuildignore`d) holds an APA-7 Quarto doc (apaquarto 5.0.18 vendored)
that renders to PDF (13 pp) + docx. A `bfi25` walkthrough (Figure 1 + Table 1) and the Forbes-AMH
centerpiece (`ackwards(forbes2023, k_max=10, pairs="all")`, Figure 2) compute live; `references.bib`
resolves all citations; AMH fidelity anchored to `test-forbes-fidelity.R` (~1e-14, OSF pcwm8).
Abstract + Method + package sections drafted; Introduction + Discussion are marked `[AUTHOR TO DRAFT]`
stubs. `devtools::check()` 0 err/0 warn/0 note (manuscript excluded); no package R code changed.

**Key decisions.** Sole author (Girard, Univ. of Kansas); built outputs gitignored (sources only);
callout blocks replaced with plain blockquotes to drop the fontawesome5 LaTeX dep; manuscript is its
own deliverable, separate from the `forbes2023` vignette candidate row. Author-owned Intro/Discussion
prose deferred as a ROADMAP candidate.

**Review.** All 7 ACs verified with fresh evidence. Independent review fixed a real Waller (2007)
bib-title error (Crossref + roxygen) and a bfi25 n=1000-vs-875 prose inconsistency; one 54-vs-55
clarity nit logged, not actioned. Blame-history clean. The same Waller-title error found pre-existing
in two vignettes was fixed separately on master (a173f6f).

PR: https://github.com/jmgirard/ackwards/pull/56 · merged on local-green (docs-only, non-release).
