# ackwards departures from the canonical sources — Goldberg (2006) & Forbes (2023) (M72)

**Provenance.** Ingested 2026-07-23 by M72 from a first-hand audit of the repo's
own design records — `cairn/DESIGN.md` §9 and the IP/GP block, `cairn/DECISIONS.md`
(D-004, D-007, D-010, D-017, D-031), the source notes [[goldberg2006]],
[[forbes2023]], [[waller2007]], [[tong2025]], and `R/ackwards.R` — cross-read
against what the two canonical sources prescribe.
Pagination: —.
Extraction: derived — no external source of its own, only as current as its inputs (DESIGN §9, DECISIONS, and the cited reference notes), none re-read since 2026-07-23 — observed 2026-07-23.

**Scope.** This is the single auditable place answering "where, why, and on what
evidence does `ackwards` depart from Goldberg (2006) or Forbes (2023)?" It is
**not** a source summary (each source has its own note) and builds **no** new
policy — it *catalogues* choices recorded elsewhere. It covers **default-level
departures** (a default where `ackwards` does something different for a task the
source specifies) and the source-*matching* defaults; it deliberately does not
itemize additive capabilities *beyond* the sources (see the GP1 note below).
Tracking disclaimer: this is a reference, not an authority — status lives in
`ROADMAP.md`, decisions in `DECISIONS.md`, architecture and principles in
`DESIGN.md`. A synthesis note that starts asserting status is a second tracking
system.

**Evidence snapshot.**

- DESIGN §9 defaults table + IP1/IP2/IP4/IP9/GP1 — `cairn/DESIGN.md` — observed 2026-07-23.
- D-004, D-007, D-010, D-017, D-031 — `cairn/DECISIONS.md` — observed 2026-07-23.
- goldberg2006 / forbes2023 / waller2007 / tong2025 source notes — `cairn/references/` — observed 2026-07-23.
- `k_max` required + `cut_show = 0.3` defaults — `R/ackwards.R:34`, `R/ackwards.R:292` — observed 2026-07-23.

## What a "departure" is here

A departure is a **default** where `ackwards` does something different for a task
one of the canonical sources specifies — not a capability the sources never had.
Adding EFA/ESEM engines, FIML missing-data, out-of-sample scoring, `boot_edges()`,
and `comparability()` are **extensions beyond** the sources, governed by GP1
(published-method capability bar) with their own D-entries (D-003, D-020, D-021,
D-023, D-022) — they are "going beyond", not "going against", and are out of this
ledger by the Scope above.

## Departures ledger — tags: `match` · `depart-supported` · `depart-gap`

`match` = follows the source; `depart-supported` = differs, with a documented
reason **and** empirical/mathematical support; `depart-gap` = differs, reason
documented but empirical support missing (→ ROADMAP candidate). IDs are stable —
cite them, never renumber.

| # | Source behavior | `ackwards` behavior | Rationale (where) | Empirical / mathematical support | Tag |
|---|---|---|---|---|---|
| E1 | Goldberg: auto-stop the ladder when no variable's highest loading is on a new component (BA-maxLoading) | `k_max` **required**; extract and **retain all** levels `1..k_max` | DESIGN §9 `k_max` row ("don't silently pick"); IP6; [[goldberg2006]] | [[tong2025]] (p. 14): BA-maxLoading recovers the true count in 58% of cases vs PA 71% (BA-cutoff worse) — the auto-stop rule is unreliable, so we don't automate it | `depart-supported` |
| E2 | Goldberg: compute factor scores at each level, then correlate all scores | Between-level edges by **W′RW algebra** when linear; scores only otherwise | IP1 (one edge path); D-004; [[waller2007]] | Waller (2007) proof that the score step is unnecessary; IP2 "both routes must agree", test-enforced | `depart-supported` |
| E3 | Goldberg: PCA component scores (his "components over factors" stance) | PCA **default** → components (**matches** Goldberg); **factor** engines (EFA/ESEM) → **tenBerge**, not components | DESIGN §9 `scores` row; D-007 | Correlation-preserving + linear (algebra-eligible); trade-off formalized by [[grice2001]] / [[beauducel2024]]; Goldberg's own "virtually identical" for PCA-vs-factor | `depart-supported` |
| E4 | Forbes: `comp.corr = t(W_a) R W_b`, **no** sign alignment | Same algebra + **primary-parent sign alignment** (IP4) | IP4; D-010; [[forbes2023]] correspondence conventions | Fidelity suite: `|values|` identical entrywise to **1.3e-14** (M44/M53) — presentational only, numerically identical | `depart-supported` |
| E5 | Forbes: `cong = psych::factor.congruence`, rounded to 2 dp | **Exact** Tucker's φ | [[forbes2023]] correspondence conventions | Fidelity suite: agreement within **0.005** (M44) — a precision gain, not a divergence | `depart-supported` |
| M1 | Goldberg: display paths ≥ **.30** in figures | `cut_show = 0.3` (default show-cut) | `R/ackwards.R:292`; [[goldberg2006]] | — (matches the source) | `match` |
| M2 | Forbes: `ChaseCorrPaths` uses the **direct/skip-level** correlation to each ancestor | `redundancy_criterion = "direct"` (default) | D-017; [[forbes2023]]; legacy M53 | Fidelity: reproduced her redundancy chase 54/54 components (M53) | `match` |

## Disposition

- **E1–E5** are all `depart-supported`: each has a documented rationale **and**
  empirical or mathematical support (E1 tong2025; E2 Waller proof + IP2 test; E3
  grice2001/beauducel2024 + Goldberg; E4/E5 the fidelity suite). **No departure
  is currently a `depart-gap`, so no candidate row is spawned** — the "ideally
  empirical support" bar is met across the board.
- **M1, M2** are `match`: recorded so a future reader confirms they are deliberate
  alignments, not accidents.
- **The governing principle** is IP9 / D-031 (exact Forbes reproduction is a
  permanent *capability*, not a default lock-in; defaults may adopt a better
  method with loud docs + a D-entry). This ledger is the audit surface for that
  principle; it produces no new rule and locks no test.

## Open questions

- Additive extensions (engines, FIML, `boot_edges`, `comparability`) are summarized as GP1-governed rather than itemized; if a future reader wants them catalogued too, that is a separate ledger, not a gap here — observed 2026-07-23.
