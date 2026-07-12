# Lessons

Durable, append-only repo lessons (build quirks, testing tricks, gotchas).
One per line: `- YYYY-MM-DD (M<NN>): <lesson>`.

- 2026-07-12 (M54): `cairn_validate.py`'s ISO-date check flags any slash-separated `N/N/N` token
  (e.g. the zero-error/warning/note shorthand for a clean check) as a malformed date — spell it out
  ("0 err/0 warn/0 note") in tracking files instead.
- 2026-07-12 (M54): a review-complete milestone file with full AC evidence easily exceeds the
  150-line cap; tighten the Review prose (review-owned) rather than the append-only work-log, and
  remember it compresses to ≤25 lines at archive anyway.
- 2026-07-12 (M54): bundling a CC-BY dataset in an MIT package = data author as `Authors@R`
  `role = "cph"` scoped by `comment` to the file + a `LICENSE.note` (a WRE-recognized top-level
  filename that ships without tripping `R CMD check`). Gives a clean check with no license NOTE.
- 2026-07-12 (M54): to bundle a dataset that also backs a fidelity test, drive `data/<name>.rda`
  and the test fixture from one md5-pinned `data-raw/` download so they cannot drift; have the test
  read the matrix from the exported dataset (verifies the shipped object, not a fixture copy).
- 2026-07-12 (M55): CRAN's "cannot suppress console messages" complaint is fixed the R-idiomatic
  way by returning a small **classed** object and moving the rendering into its `print()` method —
  top-level display is unchanged while assignment/nesting go silent. Cheaper and more idiomatic than
  a `verbose=` arg. `cli_text()` writes to stderr (the message stream), so test its output with
  `capture.output(..., type = "message")` and the `cat()` literal with plain `capture.output()`.
- 2026-07-12 (M56): when local TinyTeX is a sync *ahead* of every CTAN mirror (`tlmgr`: "remote
  database … seems to be older than the local installation"), tlmgr refuses ALL installs — including
  quarto's render-time auto-install. Workaround: fetch the package tarball straight from the tlnet
  archive (`.../texlive/tlnet/archive/<pkg>.tar.xz`), `tar -xJf` into `~/Library/TinyTeX/texmf-dist`,
  `mktexlsr`. Map a missing `foo.sty` → its texlive package via `ctan.org/json/2.0/pkg/<name>`
  (`.texlive`/`.ctan.path`); e.g. `scalefnt.sty` lives in package `carlisle`.
- 2026-07-12 (M56): apaquarto (APA-7 Quarto) hard-requires `fontawesome5` (title-block ORCID) and its
  `man`/APA float mode uses `endfloat`, which emits `.fff`/`.ttt` aux files — gitignore them plus the
  usual LaTeX aux. Quarto callout blocks pull extra font deps; plain blockquotes render fine without.
- 2026-07-12 (M56): verify reference titles against Crossref *and* the package's own roxygen — the
  repo's intro/engines vignettes carried a wrong Waller (2007) title ("…Bass-Ackward factor analysis"
  vs the correct "…Goldberg's Bass-Ackwards method") that disagreed with `R/ackwards.R`.
- 2026-07-12 (M57): a fixture's *prose* provenance claim ("regenerated from OSF, set.seed(123)") is
  not reproducibility — the M44 sims matrices were provably unrecoverable (no fx/Phi/seed/n/method
  reproduced them; even same R/psych version strings) because no generator was committed. `set.seed`
  is deterministic on a fixed env but the specific draw is lost once the generation path is; pin a
  committed `data-raw/` generator, not a claim. `test-oracle-provenance.R` now enforces this.
- 2026-07-12 (M57): before swapping a test-input fixture, check downstream *hardcoded* test literals
  for topology-stability — the new Forbes sim matrices changed the chase labels but preserved the
  redundancy topology, so `test-forbes-fidelity.R`'s prune `{m3f1,m3f2,m3f3}`/`{m2f2,m4f1,m4f2}`
  literals held. Verify, don't assume. (data-raw is `.Rbuildignore`d, so a fixture guard can only
  check that provenance *names* a generator, not that the file exists at R CMD check time.)
