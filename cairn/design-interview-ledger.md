# Design-interview banked-candidates ledger (Phase 1 → Phase 2 seam)

_Temporary working file for the /design-interview run started 2026-07-17. Owns nothing
durable — every candidate below gets a Phase-2 disposition (IP / GP / merged / skip /
ROADMAP candidate) and this file is then deleted at write-out. If found orphaned, resume
Phase 2 of /design-interview from it._

## A. Ingested CLAUDE.md Invariants (conservation: every #N needs a disposition)

| # | Invariant (short) | Proposed | Notes |
|---|---|---|---|
| 1 | One edge path — all between-level correlations via `compute_edges()`; standardize by real SDs `sqrt(diag(W'RW))` | IP | cited ~22 files deep in R/ + tests |
| 2 | Keep the cross-check — scores route retained; algebra-vs-scores agreement test | IP or GP | is it the capability or the current test idiom? |
| 3 | Light core, heavy opt-in — scores/fits/data NULL by default, recomputable | IP or GP | D-005 |
| 4 | Sign alignment anchors to primary parent, not "all positive" | IP | D-010 |
| 5 | Lineage lives in edges, never in IDs — `m{k}f{j}` stable labels | IP | D-009, D-029 lean on it |
| 6 | Loud defaults — announce consequential auto-choices; advise loudly, never switch silently | IP | D-006/D-008; most-cited in code |
| 7 | Convergence is data, not an error — warn + skip, never abort | IP | D-003; owner confirmed load-bearing (real-data ESEM wart) |
| 8 | Oracle-backed numerics — ≥2 independent oracle types; no unsourced reference value ships; provenance attrs | IP | interim home in CLAUDE.md; skill mandate: fold in as numbered principle |

## B. Banked from Phase-1 answers (2026-07-17)

- **Published-method fidelity is the capability bar** — features earn a place by implementing/
  serving a published method with verifiable fidelity; package-invented conventions stay out
  (D-020 precedent). Proposed: GP (or IP). Now recorded prose in DESIGN §2.
- **Forbes fidelity contract is permanent** — default output reproduces Forbes (2023) exactly,
  test-backed (M44/M53); reference-implementation status makes this non-negotiable. Proposed:
  IP (possibly merged with #8 or standing alone).
- **Contract boundary: tables maybe, prose never** — interpretation aids in-contract; gt table
  exports open future work; manuscript prose permanently out. Recorded in DESIGN §2; likely
  scope prose, not a numbered principle — decide.
- **Track CRAN-current upstream; floors best-effort** — recorded in DESIGN §2; posture prose,
  probably not a numbered principle.
- **Fluid until 1.0; 1.0 = BRM paper published** — recorded in DESIGN §3; posture prose.

## C. Phase-2 must add (not yet elicited)

- **Git/history-mined candidates:** report-first, flag-second (recurs: D-014 no pass/fail,
  D-017 chase reporting, D-028 trust tiering, D-030 advisory-only screens — a generalizable
  "values, never verdicts" principle); descriptive honesty (D-001 "series of linked solutions,
  not a fitted hierarchy"); wrap-don't-reimplement numerics (CLAUDE.md guardrail; D-011);
  lean-dependency posture (D-011: psych sole heavy Imports, no Rcpp); reproducibility bit-for-bit
  (D-023 serial ≡ parallel, seeds captured).
- **Domain-derived:** orthogonality as identity (D-002 varimax-only enables the whole algebra —
  essence vs. accident probe); zero-churn for unlabeled/default objects (D-029 byte-identical
  guarantee).
- **Wart follow-ups → ROADMAP candidates at write-out:** mechanize the precompute-staleness
  guard (currently prose-only in CLAUDE.md); upstream-churn watch (psych/lavaan release
  re-oracle run) — possibly a scheduled CI job (ties to "track CRAN-current").
