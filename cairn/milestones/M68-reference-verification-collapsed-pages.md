# M68: Re-verify the 3 collapsed synthesis pages against their member sources

- **Status:** in-progress
- **Priority:** normal
- **Depends on:** M67
- **Branch:** `m68-reference-verification-collapsed-pages`
- **Principles touched:** IP8

## Goal

Read every standing fact on the three collapsed `cairn/references/` synthesis pages
back against the shelf PDFs of their 17 member sources, clearing the last of the
`references staleness` advisory.

## Scope

**In:** `rotation-and-k.md` (6 members: crawford1970, browne2001a, horn1965,
velicer1976, lim2019, achim2021), `applications.md` (6 members: kim2015, markon2005,
wright2014a, forbush2018, forbush2024, cowan2024), `background.md` (5 members:
kotov2017, schmid1957, yung1999, saucier1996, widiger2019). All 17 have shelf PDFs.
Same bar as M67 — standing facts only, corrections marked in place.

**Out:**
- The nine single-source pages → M67 (this milestone depends on it).
- Interpretive prose framing each collapsed set (the italic scope note, the
  "precedent for" clauses) — repo voice, not a source claim.
- User-facing propagation of any corrected fact → `/hotfix`, enumerated by T4.
- Splitting any collapsed page back into per-source notes. These pages were
  deliberately collapsed; if verification argues one member deserves its own page,
  that is a candidate row, not a scope grab.

## Acceptance criteria

- [ ] All three pages carry an extraction status naming a verification verb and its
      own `— observed YYYY-MM-DD`.
- [ ] `cairn_validate`'s `references staleness` advisory reaches 0.
- [ ] Every member entry's standing facts are checked against that member's shelf PDF;
      each page's work-log line names the members checked and any correction made.
- [ ] M61's already-verified `lim2019` / `achim2021` entries are absorbed by citation,
      not re-derived — including its finding that a vignette clause overclaimed
      `lim2019` (best among 13 PA *variants*, not all criteria).
- [ ] Each corrected fact is grepped across `R/`, `man/`, `vignettes/`, `NEWS.md`;
      any propagation is enumerated with `file:line` and routed to `/hotfix`.
- [ ] Docs-only — no files outside `cairn/` modified; `cairn_validate` all checks PASS.

## Coverage

- AC1 → T1, T2, T3
- AC2 → T5
- AC3 → T1, T2, T3
- AC4 → T1
- AC5 → T4
- AC6 → T5

## Tasks

- [x] T1: `rotation-and-k.md` — the varimax-default and k-advice supports; check the
      CF-varimax κ = 1/p claim, Horn's PA and Velicer's MAP descriptions, and absorb
      M61's verified DOIs for `lim2019` / `achim2021`.
- [ ] T2: `applications.md` — six published-application entries; verify each citation's
      standing facts and the "precedent for" claims that name a specific finding.
- [ ] T3: `background.md` — five framing entries; check the Schmid–Leiman and
      Yung et al. descriptions especially, since DESIGN cites them for what is
      deliberately out of scope, and the `saucier1996` "ruled out as comparability
      source" note against what that paper actually says.
- [ ] T4: grep every corrected fact across `R/`, `man/`, `vignettes/`, `NEWS.md`;
      enumerate propagation with `file:line`; open `/hotfix` if any is user-facing.
- [ ] T5: re-run `cairn_validate`; confirm staleness 0 and all checks PASS; commit.

## Work log

- 2026-07-19: created by /milestone-plan.
- 2026-07-19: status → in-progress; branch `m68-reference-verification-collapsed-pages` cut from master.
- 2026-07-19: T1 `rotation-and-k.md` verified — crawford1970/browne2001a/horn1965/velicer1976 accurate (anchors added: CF Eq. 7 p. 323, orthomax equivalence pp. 324–326, Browne Table 1 p. 119 for κ = 1/p, MAP Eq. 9 p. 323); lim2019 abstract-verbatim, M61's variant-scoping absorbed as a ⚠ note; **achim2021 corrected** — the ~17% is the share of noise-factor simulations (16.8%, p. 70), not of "correct" retentions. Added crawford1970 p. 331's explicit varimax-when-k-unknown recommendation. Noted lim2019's shelf PDF is online-first (no journal pagination). No propagation into `R/`, `man/`, `vignettes/`, `NEWS.md`.

## Decisions

## Review
