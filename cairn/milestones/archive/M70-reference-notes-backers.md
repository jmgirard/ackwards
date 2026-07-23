# M70: Author + verify the 5 default-rationale backer notes and wire them into DESIGN §9 + roxygen

**Status:** done (2026-07-23, PR #74 https://github.com/jmgirard/ackwards/pull/74)

**Goal:** Author a verified `cairn/references/` source note for each of 5 default-rationale backer PDFs (grice2001, beauducel2024, tong2025, williams2025, kaiser1958) and wire each citation into the DESIGN §9 default it supports and the relevant roxygen `@references`.

**Outcome:** 5 new source notes under a new INDEX section "Default-rationale backers", all `Extraction: verified 2026-07-23` against rendered shelf pages. DESIGN §9 rationale rows now cite their evidence — `rotation`→kaiser1958; `scores`→grice2001 + beauducel2024 + williams2025; `redundancy_phi`→grice2001; `k_max`→tong2025 — citations only, no default value changed. `@references` added to `R/ackwards.R` (kaiser1958, grice2001, beauducel2024, williams2025) + `R/suggest_k.R` (tong2025); `man/*.Rd` regenerated, NAMESPACE unchanged. Shelf `williams2025a.pdf → williams2025.pdf`. DOIs Crossref-confirmed; beauducel third author is "Kuhl" (not "Kühl"); kaiser1958 prints no DOI (registered BF02289233); **tong2025 shelf copy is the preprint** (page anchors are preprint pages, not the published 17–44).

**Decisions:** none (Principles touched: —; no NEWS — documentation-only refinement).

**Review:** 6/6 ACs verified with fresh evidence; DoD gate check 0/0/0, coverage 100%, style/lint/pkgdown clean; `cairn_validate` exit 0. Three fresh-context lenses (diff-bug/Opus, blame-history/Sonnet, prior-review/Sonnet) all **0 findings** — the Opus and prior-review lenses independently rendered every cited PDF page and re-checked all 5 DOIs against Crossref; scorer no-op. One in-review fix: dropped a stale `williams2025a` token from williams2025.md provenance so AC1's grep resolved cleanly.
