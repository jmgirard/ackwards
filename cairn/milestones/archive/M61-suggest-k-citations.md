# M61: Enrich suggest_k() docs with verified k-selection citations — done 2026-07-16

**Goal.** Source `suggest_k()`'s two uncited claims — the consensus-range stance and
"PA-PC overextracts" — with verified-at-ingest citations, in roxygen and the vignette.

**Outcome.** Help page + suggest-k vignette now cite Lim & Jahng (2019, PA estimate as a
±1 range resolved by interpretability) paired with Achim (2021, counterpoint — the
disagreement itself supports D-013's range-never-a-number stance), and Saucier (1997,
fn 14: PA suggesting ~30 factors in wide lexical sets) alongside Forbes (2023) for the
overextraction note. Three @references / vignette References entries; one NEWS 0.1.1
bullet. Doc-only; no behavior changes.

**Key decisions.**
- Vignette + NEWS scope added at the plan gate (candidate row was roxygen-only).
- precompute.R regenerates all vignettes; unrelated run-noise churn (cli timings, gt
  element IDs) restored rather than committed, keeping the PR diff to suggest-k files.

**Verification.** DOIs three-way-matched Rd ↔ rotation-and-k.md ↔ saucier1997.md;
dod-gate passed (check 0 err/0 warn/0 note, coverage 100%, style/lint/pkgdown clean);
suite 2300 pass / 0 fail. Independent review (3 lenses + scorer): 1 finding (85) — a
vignette clause overclaimed lim2019 (PA best among 13 PA *variants*, not all criteria);
fixed by dropping the clause + precompute re-run. Blame-history: no reintroduction of
the retired goldberg1990 citation. Prior-PR: no evidence (no-op).

**PR.** https://github.com/jmgirard/ackwards/pull/62
