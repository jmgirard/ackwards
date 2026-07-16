# M62: Worked Forbes (2023) AMH example vignette — done 2026-07-16

**Goal.** Ship a precomputed vignette reproducing the Forbes (2023) AMH applied example
end-to-end on the bundled `forbes2023` dataset.

**Outcome.** `vignette("ackwards-forbes2023")` now walks the paper's 155-variable example:
10-level fit (`pairs = "all"`, 1320 edges/45 level-pairs), skip-level structure, the redundancy
chase (37/55 flagged under the default `direct` criterion, the worked d4→…→j4 chain, the
36-flag/`m3f3` divergence under the `adjacent` opt-in), and the pruned-factor diagram in her
publication style. Every headline number baked from a live run and pinned by
`test-forbes-fidelity.R`. Cross-links both ways with `ackwards-forbes`; `@seealso` on
`?forbes2023`; `R/data.R` Forbes-citation typo fixed (Crossref-verified "bass-ackward method");
pkgdown articles row; NEWS entry under a new "(development version)" heading.

**Key decisions.** New standalone vignette (bfi25 concept vignette untouched); reproduction
scope only — extended-tutorial extras declined at the plan gate. AC4 grep narrowed
mid-implement: Waller (2007)'s title legitimately uses "Bass-Ackwards"; the AC targets the
Forbes-citation pattern only (minor amendment, logged).

**Verification.** DoD gate passed (check 0 err/0 warn/0 note, coverage 100%, style/lint clean,
pkgdown complete); all 7 ACs evidence-verified; 3-lens independent review: 0 findings (diff-bug
reviewer re-knit the vignette — byte-identical output/figures). Two pre-existing observations
logged → candidate rows (Goldberg subtitle repo-wide; PCA `n_obs` message overpromise).

**PR.** https://github.com/jmgirard/ackwards/pull/63 · merged on local-green (non-release).
