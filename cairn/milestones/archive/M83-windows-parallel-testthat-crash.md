# M83: Windows parallel-testthat crash — measure, diagnose, mitigate

**Status:** done (2026-07-27, PR #91 https://github.com/jmgirard/ackwards/pull/91)

**Goal:** Make the Windows CI job a trustworthy signal again by measuring, diagnosing, and mitigating the intermittent access violation that kills a testthat parallel worker.

**Outcome:** `.github/workflows/windows-stress.yaml` — a dispatch-only sharded harness that runs the suite N times on `windows-latest` and tallies `-1073741819` access violations, scoring by crash signature (the parent exits 1 and only reports the worker's code) and counting ordinary test failures separately. 530 iterations measured: 1 worker 15/120 (12.5%), 2 workers 11/230 (4.8%), 4 workers 6/60 (10.0%), 4 workers minus EFAtools 2/120 (1.7%). **No setting eliminates the crash**; the relationship is non-monotonic, serial being worst. Memory excluded (a crash landed at its shard's highest free-memory reading). Nine distinct test files affected, `test-suggest_k.R` over-represented. EFAtools implicated (p≈0.02) but not isolated — its removal also removes computation. `mnormt` untestable this way (`psych` requires it). `R-CMD-check.yaml` now runs 2 testthat workers on Windows only, via a `matrix.config.testthat-cpus` key, with the measurements and the operating policy recorded at that site. Confined to `.github/`; no tarball content changed.

**Decisions:** AC1 reworded (detection keys on the crash signature, not the iteration's exit code). AC4 amended mid-implementation from "zero crashes over ≥60 iterations" — unreachable on any axis — to a measured characterisation plus a stated re-run policy.

**Review:** Diff-bug lens 6 findings, all scored below the 80 threshold; blame-history and prior-review lenses none (PR-comment probe returned empty, so archived `## Review` sections were the evidence). F4 (45) and F5 (68) fixed despite being sub-threshold: an overclaimed memory comment, and shell fallbacks that disagreed with documented input defaults. F1/F2/F3/F6 rejected as scored. One review return: AC4's policy existed only inside its own wording.
