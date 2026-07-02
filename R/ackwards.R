#' Bass-ackwards hierarchical structural analysis
#'
#' Extracts factor/component solutions at levels 1 through `k`, then
#' characterises the hierarchy by computing between-level factor-score
#' correlations. The "hierarchy" is descriptive: edges are score correlations,
#' not a fitted higher-order SEM.
#'
#' @section Defaults and why:
#' * **`engine = "pca"`** -- the original Goldberg (2006) method; fastest; never
#'   fails to converge; the Waller (2007) algebra is exact for components.
#' * **`rotation = "varimax"`** -- the `T'=T^-1` property of orthogonal
#'   rotation enables the closed-form `W'RW` edge algebra and keeps
#'   within-level factors uncorrelated so cross-level edges reflect only the
#'   hierarchical signal. Matches Goldberg (2006), Kim & Eaton (2015), and
#'   Forbush et al. (2024). Varimax is the only supported rotation; oblique
#'   rotation would confound the between-level signal that is the method's
#'   core output.
#' * **`cor = "pearson"`** -- no silent basis switching. If your items look
#'   ordinal (<= 7 distinct integer values), a cli warning will suggest
#'   `cor = "polychoric"`, which is available for all three engines.
#' * **`align_signs = TRUE`** -- unaligned signs make the output unreadable.
#'   Anchor: m1f1 is oriented toward the positive manifold; each subsequent
#'   factor is flipped so its edge to its primary parent is positive.
#' * **`keep_scores = FALSE` / `keep_fits = FALSE`** -- memory and privacy.
#'   Scores are O(n x Sigmak) and often sensitive; raw engine fits can be large.
#'   Both are recomputable from the stored `r` matrix.
#'
#' @param data A data frame or numeric matrix of observed variables (items in
#'   columns, observations in rows). Alternatively, a pre-computed **correlation
#'   matrix** may be supplied (a square, symmetric, numeric matrix with unit
#'   diagonal). When a correlation matrix is supplied, `engine` must be `"pca"`
#'   or `"efa"` (ESEM requires raw data), and the `missing` and `cor` arguments
#'   are ignored. See the *Correlation-matrix input* section below.
#' @param k_max Maximum number of factors/components to extract. Required; use
#'   [suggest_k()] if uncertain. Sets the *depth* of the hierarchy: levels
#'   1 through `k_max` are all extracted and retained. (Note: `suggest_k()`'s
#'   `k_max` means something related but distinct -- the max number of
#'   factors/components it *evaluates* when recommending a depth, not a depth
#'   itself; see that function's docs.)
#' @param n_obs Number of observations, or (on the raw-data FIML path) a string
#'   selecting which N feeds the fit indices.
#'
#'   * **Correlation-matrix input:** a positive integer. Required for
#'     `engine = "efa"` ([psych::fa()] needs N for chi-square / RMSEA / TLI);
#'     optional for `"pca"` (stored as `NA_integer_` if omitted, disabling
#'     N-dependent fit statistics).
#'   * **Raw data:** N is normally taken from `nrow(data)` and a numeric `n_obs`
#'     is ignored (with a warning). The exception is `missing = "fiml"` with
#'     `engine = "pca"`/`"efa"` (M38): because [psych::corFiml()] estimates the
#'     correlation matrix from incomplete rows, `n_obs` may be `"total"`
#'     (default -- every row contributing to the FIML likelihood, matching the
#'     FIML convention; Enders, 2010) or `"complete"` (complete-case N, a
#'     conservative lower bound). Point estimates (loadings, edges) do not
#'     depend on this choice; only the EFA fit indices do, and those are
#'     *approximate* under this two-step (FIML matrix into normal-theory EFA)
#'     route regardless of N (Zhang & Savalei, 2020). A string `n_obs` is
#'     accepted only on this path.
#' @param engine Extraction engine: `"pca"` (default), `"efa"`, or `"esem"`.
#'   `"esem"` uses [lavaan::efa()] with rotation-aware SEs and per-level fit
#'   indices; recommended for the clinical/HiTOP workflow (Kim & Eaton, 2015;
#'   Forbush et al., 2024). Requires lavaan >= 0.6-13.
#' @param fm Factor extraction method passed to [psych::fa()]; only used when
#'   `engine = "efa"`. One of `"minres"` (default, robust OLS), `"ml"`
#'   (maximum likelihood, gives chi-square fit but converges less reliably at
#'   deep levels), or `"pa"` (principal axis). Ignored for `engine = "pca"`.
#' @param cor Correlation basis: `"pearson"` (default), `"spearman"`, or
#'   `"polychoric"`. For PCA/EFA, `"polychoric"` computes a polychoric
#'   correlation matrix via `psych::polychoric()` (requires psych). For ESEM,
#'   it triggers WLSMV estimation via lavaan. A cli warning is emitted when
#'   ordinal-looking columns are detected and `cor` is not `"polychoric"`.
#' @param estimator Estimation method for the ESEM engine. `NULL` (default)
#'   auto-selects: `"WLSMV"` when `cor = "polychoric"`, `"ML"` otherwise.
#'   Pass explicitly to override: `"ULSMV"` (unweighted WLS), `"MLR"`
#'   (robust ML). `cor = "polychoric"` with `estimator = "ML"`/`"MLR"` errors
#'   (lavaan itself does not support ML/MLR on ordered indicators); `"WLSMV"`/
#'   `"ULSMV"` with a continuous `cor` is allowed (a valid, if atypical,
#'   continuous WLS/ADF estimator). Ignored for PCA and EFA engines. The
#'   effective value (after auto-selection) is recorded in `x$meta$estimator`
#'   (`NA` for PCA/EFA).
#' @param missing How to handle missing item responses. One of:
#'   * `"pairwise"` (default) -- use all available observations pairwise.
#'     For PCA/EFA this feeds `stats::cor(use = "pairwise.complete.obs")`.
#'     For ESEM with WLSMV/ULSMV (ordinal), lavaan uses `available.cases`,
#'     which computes polychoric thresholds and correlations from all rows that
#'     contribute to each pair -- MCAR-valid and uses the full N.
#'     For ESEM with ML/MLR (continuous), lavaan uses listwise deletion
#'     internally while edges are computed from a pairwise correlation matrix;
#'     this minor inconsistency is documented in `$meta`. A warning is emitted
#'     when incomplete rows are detected.
#'   * `"listwise"` -- only complete rows are used. Reduces data to
#'     `stats::complete.cases()` before fitting, so the correlation matrix,
#'     the engine fit, and the edges are all consistent. `n_obs` in the result
#'     reflects the reduced N.
#'   * `"fiml"` -- Full Information Maximum Likelihood. For `engine = "esem"`
#'     (with `estimator = "ML"`/`"MLR"`), passes `missing = "fiml"` to
#'     [lavaan::efa()] and derives edge correlations from lavaan's
#'     FIML-estimated saturated model. For `engine = "pca"`/`"efa"` (M38), the
#'     correlation matrix is estimated via [psych::corFiml()] and fed to the
#'     usual `W'RW` algebra; this requires `cor = "pearson"` (corFiml estimates
#'     a multivariate-normal matrix) and the route is announced via a cli
#'     message. Errors for WLSMV/ULSMV and for a non-Pearson PCA/EFA basis.
#'     Note: FIML improves estimation under missingness but does not impute item
#'     responses; score materialisation (`keep_scores = TRUE`) still produces
#'     `NA` rows for incomplete observations. See `n_obs` for the fit-index
#'     sample size on the PCA/EFA path.
#' @param align_signs Logical; sign-align factors to primary-parent lineage?
#'   Default `TRUE`.
#' @param keep_scores Logical; store factor scores in the result? Default
#'   `FALSE` (recomputable via [augment.ackwards()]). When `TRUE`,
#'   per-observation scores are stored in `x$scores` as a named list of
#'   `n x k_j` matrices, one per level, standardized by real score SDs
#'   (see [augment.ackwards()]).
#' @param keep_fits Logical; store raw engine fit objects? Default `FALSE`.
#'   When `TRUE`, the per-level fit objects (psych or lavaan) are stored in
#'   `x$fits` as a named list indexed by level.
#' @param seed Integer seed for stochastic engines (not used by PCA but
#'   captured for reproducibility metadata). Default `NULL`.
#' @param pairs Which level pairs to compute edges for: `"adjacent"` (default,
#'   classic Goldberg -- only consecutive levels) or `"all"` (Forbes extension --
#'   every pair of levels, including skip-level correlations). [prune()]
#'   recomputes its own all-pairs edges on demand regardless of this setting,
#'   so pruning does not require `pairs = "all"` here.
#' @param cut_show Edges with `|r| >= cut_show` are flagged `above_cut` in
#'   `tidy()` output. Default `0.3`.
#' @param ... Reserved for future arguments.
#'
#' @return An object of class `"ackwards"`. See [print.ackwards()],
#'   [tidy.ackwards()], [glance.ackwards()], and [augment.ackwards()] for
#'   output methods.
#'
#' @seealso [print.ackwards()], [tidy.ackwards()], [glance.ackwards()],
#'   [prune()] (Forbes-extension redundancy/artifact flagging, piped off the
#'   result of this function)
#'
#' @references
#' Goldberg, L. R. (2006). Doing it all bass-ackwards. *Journal of Research in
#'   Personality*, 40(4), 347--358. \doi{10.1016/j.jrp.2006.01.001}
#'
#' Waller, N. G. (2007). A general method for computing hierarchical component
#'   structures by Goldberg's bass-ackwards method. *Journal of Research in
#'   Personality*, 41(4), 745--752. \doi{10.1016/j.jrp.2006.08.005}
#'
#' Forbes, M. K. (2023). Improving hierarchical models of individual
#'   differences: An extension of Goldberg's bass-ackward method.
#'   *Psychological Methods*. \doi{10.1037/met0000546}
#'
#' Enders, C. K. (2010). *Applied Missing Data Analysis*. Guilford Press.
#'
#' Zhang, X., & Savalei, V. (2020). Examining the effect of missing data on
#'   RMSEA and CFI under normal theory full-information maximum likelihood.
#'   *Structural Equation Modeling*, 27(2), 219--239.
#'   \doi{10.1080/10705511.2019.1642111}
#'
#' @section Performance (ESEM, large item sets):
#' The ESEM engine fits a separate `lavaan` model at every level 1..`k_max`.
#' For ordinal data (`cor = "polychoric"`, WLSMV) the costly sample statistics
#' lavaan derives from the raw data -- thresholds, the polychoric correlation
#' matrix, and the asymptotic weight matrix -- depend only on the data, not on
#' the number of factors, so they are **computed once** at the first level and
#' **reused** for every deeper level (identical solutions, much less work). This
#' matters most when you have many items (hundreds), where recomputing those
#' statistics at each level dominated the run time.
#'
#' The per-level model fits are mutually independent and are dispatched through
#' the \pkg{future} framework when \pkg{future.apply} is installed. By default
#' the plan is sequential (no behaviour change). To run the levels in parallel,
#' set a plan once before calling `ackwards()`:
#' \preformatted{
#'   future::plan(future::multisession, workers = 4)  # or multicore on Unix
#'   x <- ackwards(items, k_max = 8, engine = "esem", cor = "polychoric")
#' }
#' Parallelism pays off when the per-level fits are heavy (large `p`, several
#' levels); for small problems the worker startup cost can outweigh it. Results
#' are reproducible across plans when `seed` is supplied. PCA and EFA already
#' compute their correlation matrix once and are unaffected.
#'
#' @section Correlation-matrix input:
#' When `data` is a pre-computed correlation matrix (square, symmetric, unit
#' diagonal), `ackwards()` runs entirely from that matrix using the `W'RW`
#' algebra -- no raw item responses are needed. This is useful when you have
#' a published correlation table or a polychoric matrix computed externally.
#'
#' Constraints and behaviour when a correlation matrix is supplied:
#' * **Engine:** only `"pca"` and `"efa"` are supported. `"esem"` requires raw
#'   data (for lavaan's own polychoric computation, WLSMV estimation, and
#'   per-level fit indices) and will error clearly.
#' * **`n_obs`:** required for `"efa"` (psych needs N for chi-square / RMSEA /
#'   TLI); optional for `"pca"` (stored as `NA` if omitted).
#' * **`cor` argument:** ignored -- the basis is already determined by the
#'   matrix you supply. A warning is emitted if you set `cor` explicitly.
#' * **`missing` argument:** ignored -- missingness was handled when computing
#'   the matrix. A warning is emitted if you set `missing` explicitly.
#' * **Factor scores:** `keep_scores = TRUE` will error; `augment()` and
#'   `tidy(what = "scores")` will also error because individual-level scores
#'   require row-level item responses.
#' * **`$cor` field:** stored as `NA_character_`; printed as
#'   `"(user-supplied matrix)"`.
#'
#' @examples
#' x <- ackwards(bfi25, k_max = 5)
#' print(x)
#' tidy(x)
#' glance(x)
#'
#' # Correlation-matrix input (PCA engine; n_obs optional)
#' R <- cor(bfi25, use = "pairwise.complete.obs")
#' x_R <- ackwards(R, k_max = 5)
#' print(x_R)
#'
#' @export
ackwards <- function(
  data,
  k_max,
  engine = "pca",
  cor = "pearson",
  fm = "minres",
  estimator = NULL,
  missing = "pairwise",
  n_obs = NULL,
  align_signs = TRUE,
  keep_scores = FALSE,
  keep_fits = FALSE,
  seed = NULL,
  pairs = "adjacent",
  cut_show = 0.3,
  ...
) {
  cl <- match.call()

  # M34: pruning moved to a standalone prune() verb; the five args below no
  # longer exist here. Without this guard they would be silently absorbed by
  # `...` (a masked-argument footgun, not a clean break) -- error loudly
  # instead (Invariant 6).
  dots <- list(...)
  moved_args <- intersect(
    names(dots),
    c("prune", "redundancy_r", "redundancy_phi", "min_items", "orphan_r")
  )
  if (length(moved_args) > 0L) {
    cli::cli_abort(c(
      "!" = "{.arg {moved_args}} {?is/are} no longer {?an argument/arguments} \\
             to {.fn ackwards}.",
      "i" = "Pruning moved to a standalone verb (M34): use \\
             {.code ackwards(...) |> prune(...)} instead."
    ))
  }

  # --- Detect correlation-matrix vs. raw-data input ---------------------------
  # Guard first: a square symmetric matrix with non-unit diagonal is almost
  # certainly a covariance matrix -- give a targeted error before it silently
  # routes to the raw-data branch and produces a confusing message.
  if (is.matrix(data)) .check_maybe_cov_matrix(data)
  input_type <- if (.is_cor_matrix(data)) "cor_matrix" else "data"
  # Effective ESEM estimator, resolved in BRANCH B below; stays NULL for
  # PCA/EFA and for BRANCH A (cor_matrix input, which errors for engine =
  # "esem" before reaching this point). Recorded in $meta$estimator (M32).
  estimator_eff <- NULL

  # Variable labels (M36): captured from a data.frame's per-column "label"
  # attribute (the attribute labelled/haven write) before as.matrix() strips it,
  # so top_items() can display "How much..." instead of the bare colname. NULL
  # for matrix / correlation-matrix input, or when no column carries a label.
  item_labels <- .capture_item_labels(data)

  # --- Shared argument validation ---------------------------------------------
  engine <- rlang::arg_match(engine, c("pca", "efa", "esem"))
  fm <- rlang::arg_match(fm, c("minres", "ml", "pa"))
  pairs <- rlang::arg_match(pairs, c("adjacent", "all"))

  if (!is.numeric(k_max) || length(k_max) != 1L || k_max < 2L || k_max != as.integer(k_max)) {
    cli::cli_abort("{.arg k_max} must be an integer >= 2 (need at least two levels for a hierarchy).")
  }
  k_max <- as.integer(k_max)

  if (!is.numeric(cut_show) || length(cut_show) != 1L || is.na(cut_show) ||
    cut_show < 0 || cut_show > 1) {
    cli::cli_abort(
      "{.arg cut_show} must be a single number in [0, 1] (it thresholds |r|)."
    )
  }

  # ============================================================
  # BRANCH A -- correlation-matrix input
  # ============================================================
  if (input_type == "cor_matrix") {
    # Gate ESEM: lavaan must fit on raw data
    if (engine == "esem") {
      cli::cli_abort(c(
        "!" = "{.code engine = \"esem\"} requires raw item data.",
        "i" = "ESEM (lavaan) needs individual responses for WLSMV estimation, \\
               polychoric computation, and per-level fit indices.",
        "i" = "Use {.code engine = \"pca\"} or {.code engine = \"efa\"} \\
               when supplying a correlation matrix."
      ))
    }

    # Warn if cor/missing are set explicitly (they're ignored for R input).
    cor_was_set <- !missing(cor) && cor != "pearson"
    if (cor_was_set) {
      cli::cli_warn(
        c(
          "!" = "{.arg cor} is ignored when a correlation matrix is supplied \\
               (the basis is already determined by your matrix).",
          "i" = "Remove {.arg cor} from your call to suppress this warning."
        ),
        .frequency = "once",
        .frequency_id = "ackwards_cor_ignored"
      )
    }
    # missing() is a base R function that tests whether the caller supplied the
    # argument; it coincidentally has the same name as the formal `missing`.
    missing_was_set <- !missing(missing) && missing != "pairwise"
    if (missing_was_set) {
      cli::cli_warn(
        c(
          "!" = "{.arg missing} is ignored when a correlation matrix is supplied \\
               (missingness was handled before computing the matrix).",
          "i" = "Remove {.arg missing} from your call to suppress this warning."
        ),
        .frequency = "once",
        .frequency_id = "ackwards_missing_ignored"
      )
    }

    # n_obs: required for EFA, optional for PCA
    if (!is.null(n_obs)) {
      if (!is.numeric(n_obs) || length(n_obs) != 1L ||
        n_obs < 1L || n_obs != as.integer(n_obs)) {
        cli::cli_abort("{.arg n_obs} must be a positive integer when supplied.")
      }
      n_obs <- as.integer(n_obs)
    }
    if (engine == "efa" && is.null(n_obs)) {
      cli::cli_abort(c(
        "!" = "{.arg n_obs} is required when {.code engine = \"efa\"} and \\
               a correlation matrix is supplied.",
        "i" = "{.fn psych::fa} needs N to compute chi-square / RMSEA / TLI.",
        "i" = "Supply the number of observations: {.code n_obs = <N>}."
      ))
    }

    # Warn if n_obs supplied alongside raw data (handled below)
    # -- handled in the data branch instead.

    # Block keep_scores=TRUE: no row-level data to score
    if (isTRUE(keep_scores)) {
      cli::cli_abort(c(
        "!" = "{.code keep_scores = TRUE} requires raw item data.",
        "i" = "Factor scores are individual-level (n x k) projections and \\
               cannot be computed from a correlation matrix alone.",
        "i" = "Either supply raw data or set {.code keep_scores = FALSE}."
      ))
    }

    # Validate and normalise R; synthesise dimnames if absent
    R <- .validate_cor_matrix(as.matrix(data))
    p <- nrow(R)

    if (k_max > p) {
      cli::cli_abort(
        "{.arg k_max} ({k_max}) cannot exceed the number of variables ({p})."
      )
    }

    # Store effective n_obs (NA when not supplied for PCA)
    n_obs_eff <- if (is.null(n_obs)) NA_integer_ else n_obs
    if (engine == "pca" && is.null(n_obs)) {
      cli::cli_inform(c(
        "i" = "{.arg n_obs} not supplied; stored as {.code NA}.",
        "i" = "Fit statistics requiring N (chi-square, RMSEA, TLI) \\
               are unavailable. Pass {.code n_obs = <N>} to enable them."
      ))
    }

    cor_eff <- NA_character_
    missing_eff <- NA_character_
    n_complete <- NA_integer_
    is_ordinal <- FALSE

    if (!is.null(seed)) set.seed(seed)

    engine_out <- switch(engine,
      pca = pca_levels(R,
        k_max = k_max, cor = "pearson",
        keep_fits = keep_fits
      ),
      efa = efa_levels(R,
        k_max = k_max, fm = fm, n_obs = n_obs_eff,
        cor = "pearson", keep_fits = keep_fits
      )
    )
    levels_list <- engine_out$levels
    fits_stored <- engine_out$fits
    data_mat <- NULL
  } else {
    # ============================================================
    # BRANCH B -- raw data input (existing behaviour)
    # ============================================================

    # n_obs on the raw-data path is resolved after `missing`/`engine` below
    # (its meaning depends on them); see the n_obs_mode block.

    cor <- rlang::arg_match(cor, c("pearson", "spearman", "polychoric"))
    if (!is.null(estimator)) {
      estimator <- rlang::arg_match(estimator, c("ML", "MLR", "WLSMV", "ULSMV"))
    }
    missing <- rlang::arg_match(missing, c("pairwise", "listwise", "fiml"))

    # --- n_obs on the raw-data path (M38) ------------------------------------
    # Numeric n_obs is for correlation-matrix input only; with raw data N comes
    # from the data. Under missing = "fiml" (pca/efa) a *string* n_obs selects
    # which N feeds the fit indices: "total" (all rows contributing to the FIML
    # likelihood -- the FIML convention, Enders 2010) or "complete" (complete
    # cases; a conservative lower bound). Default "total".
    n_obs_mode <- "total"
    fiml_pcaefa <- missing == "fiml" && engine %in% c("pca", "efa")
    if (!is.null(n_obs)) {
      if (is.character(n_obs)) {
        n_obs <- rlang::arg_match(n_obs, c("total", "complete"))
        if (!fiml_pcaefa) {
          cli::cli_abort(c(
            "!" = "{.code n_obs = \"{n_obs}\"} is only valid with \\
                   {.code missing = \"fiml\"} and {.code engine = \"pca\"} or \\
                   {.code \"efa\"}.",
            "i" = "For correlation-matrix input pass a numeric N; with raw data \\
                   and no FIML, N is taken from {.code nrow(data)}."
          ))
        }
        n_obs_mode <- n_obs
      } else {
        cli::cli_warn(
          c(
            "!" = "{.arg n_obs} is ignored when raw data are supplied \\
                 (N is determined from {.code nrow(data)}).",
            if (fiml_pcaefa) {
              c("i" = "Under {.code missing = \"fiml\"}, pass \\
                       {.code n_obs = \"total\"} or {.code \"complete\"} to \\
                       choose which N feeds the fit indices.")
            } else {
              c("i" = "Remove {.arg n_obs} from your call to suppress this warning.")
            }
          ),
          .frequency = "once",
          .frequency_id = "ackwards_nobs_ignored"
        )
      }
    }

    # cor = "polychoric" marks every item `ordered` for lavaan; ML/MLR flatly
    # do not support ordered indicators (lavaan itself errors with "estimator
    # ML for ordered data is not supported yet"). Without this guard the
    # failure surfaces many calls deep as a per-level ESEM warning, then a
    # generic "failed to build at least 2 converged levels" abort that
    # misdiagnoses the cause as multicollinearity (Invariant 6: loud, not
    # buried). WLSMV/ULSMV with a continuous `cor` is left alone -- lavaan
    # runs it as a valid (if atypical) continuous WLS/ADF estimator.
    if (engine == "esem" && cor == "polychoric" && !is.null(estimator) &&
      estimator %in% c("ML", "MLR")) {
      cli::cli_abort(c(
        "!" = "{.code cor = \"polychoric\"} is incompatible with \\
               {.code estimator = \"{estimator}\"}.",
        "i" = "Polychoric correlations require ordinal (WLS-family) estimation; \\
               {.val {estimator}} assumes continuous indicators and lavaan \\
               itself errors on the combination.",
        "i" = "Use {.code estimator = \"WLSMV\"} (default) or \\
               {.code \"ULSMV\"}, or switch to {.code cor = \"pearson\"} \\
               to use {.val {estimator}}."
      ))
    }

    # Resolve effective estimator early (needed for missing= validation before
    # the engine dispatch block; also avoids repeating the logic below).
    estimator_eff <- if (engine == "esem") {
      if (!is.null(estimator)) {
        estimator
      } else if (cor == "polychoric") {
        "WLSMV"
      } else {
        "ML"
      }
    } else {
      NULL
    }
    .resolve_missing(missing, engine, estimator_eff, cor)

    # cor = "spearman" + engine = "esem" is semantically inconsistent: lavaan fits
    # Pearson ML on raw data while compute_edges() uses Spearman R for scoring
    # (DESIGN.md s.14 known limitations). Warn loudly rather than silently mixing bases.
    if (engine == "esem" && cor == "spearman") {
      cli::cli_warn(
        c(
          "!" = "{.code cor = \"spearman\"} with {.code engine = \"esem\"} \\
                 uses inconsistent bases: lavaan fits a Pearson-ML model on raw \\
                 data while edges are computed from the Spearman correlation matrix.",
          "i" = "Consider {.code cor = \"polychoric\"} for ordinal data or \\
                 {.code cor = \"pearson\"} for a consistent continuous path."
        ),
        .frequency = "once",
        .frequency_id = "ackwards_esem_spearman"
      )
    }

    if (!is.data.frame(data) && !is.matrix(data)) {
      cli::cli_abort("{.arg data} must be a data frame or numeric matrix.")
    }
    data_mat <- as.matrix(data)
    if (!is.numeric(data_mat)) {
      cli::cli_abort("{.arg data} must contain only numeric columns.")
    }

    p <- ncol(data_mat)

    if (k_max > p) {
      cli::cli_abort("{.arg k_max} ({k_max}) cannot exceed the number of variables ({p}).")
    }

    # --- Missing-data preparation (M16, Invariant 6) --------------------------
    # n_complete is computed from the original data regardless of deletion mode.
    n_complete <- sum(stats::complete.cases(data_mat))

    # Listwise: reduce to complete rows now; all downstream steps (R computation,
    # engine calls, score materialisation) will use the same consistent set.
    if (missing == "listwise" && n_complete < nrow(data_mat)) {
      data_mat <- data_mat[stats::complete.cases(data_mat), , drop = FALSE]
    }

    # n_obs reflects listwise reduction (if applied); pairwise/fiml keep total N.
    n_obs_eff <- nrow(data_mat)

    # Pairwise advisory when data has NAs (Invariant 6 -- loud defaults).
    # For ESEM ML/MLR this is especially relevant: lavaan defaults to listwise for
    # model fitting while edges use pairwise correlations, which can make fit
    # statistics slightly anti-conservative.
    # No .frequency = "once": unlike the ordinal-detection warning (a fixed
    # property of the data), this is a per-call data-state advisory -- the user
    # may be iterating k_max or engine and should be reminded each time NAs are
    # present with the pairwise default.
    if (missing == "pairwise" && n_complete < nrow(data_mat)) {
      n_incomplete <- nrow(data_mat) - n_complete
      supports_fiml <- !is.null(estimator_eff) && !estimator_eff %in% c("WLSMV", "ULSMV")
      cli::cli_warn(
        c(
          "!" = "{n_incomplete} row{?s} have missing values; correlations are \\
                 computed pairwise.",
          if (supports_fiml) {
            c("i" = "Use {.code missing = \"listwise\"} for consistent \\
                     complete-case analysis, or {.code missing = \"fiml\"} \\
                     for full-information ML.")
          } else {
            c("i" = "Use {.code missing = \"listwise\"} for consistent \\
                     complete-case analysis.")
          }
        )
      )
    }

    # --- Ordinal-detection warning (DESIGN.md s.9, Invariant 6) ---------------
    # Compute once; reused in meta$ordinal_warned below.
    # Only warn when the user has NOT already opted into the polychoric basis.
    ordinal_cols <- detect_ordinal(as.data.frame(data_mat))
    is_ordinal <- length(ordinal_cols) > 0L
    if (is_ordinal && cor != "polychoric") {
      # Name the flagged columns so the advice is actionable for mixed data
      # (M42/e3); truncate long lists to keep the warning readable.
      cols_show <- cli::cli_vec(ordinal_cols, style = list("vec-trunc" = 8))
      cli::cli_warn(
        c(
          "!" = "{length(ordinal_cols)} column{?s} look{?s/} like \\
                 ordinal/Likert items (<= 7 distinct integer values): \\
                 {.val {cols_show}}.",
          "i" = "Results use a {.val {cor}} basis. Consider \\
                 {.code cor = \"polychoric\"} for ordinal data."
        ),
        .frequency = "once",
        .frequency_id = "ackwards_ordinal_warning"
      )
    }

    cor_eff <- cor
    missing_eff <- missing

    if (!is.null(seed)) set.seed(seed)

    # --- Compute correlation matrix + extract levels --------------------------
    # ESEM: lavaan owns R computation for polychoric (WLSMV uses its own polychoric
    # matrix internally); for continuous paths we pre-compute R and pass it through.
    # PCA/EFA: always compute R here and pass to the engine.
    if (engine == "esem") {
      # Pre-compute edge R for non-polychoric, non-FIML paths. For FIML, R is
      # extracted from lavaan's FIML saturated model inside esem_levels().
      # For listwise, data_mat is already complete so pairwise.complete.obs ==
      # complete.obs -- no special handling needed.
      R_ext <- if (cor != "polychoric" && missing != "fiml") {
        stats::cor(data_mat, method = cor, use = "pairwise.complete.obs")
      } else {
        NULL
      }
      # estimator_eff already computed above for .resolve_missing() validation.
      esem_out <- esem_levels(data_mat,
        k_max = k_max, estimator = estimator_eff, cor = cor,
        R_external = R_ext, keep_fits = keep_fits,
        missing = missing
      )
      levels_list <- esem_out$levels
      fits_stored <- esem_out$fits
      R <- esem_out$r_lv
      if (is.null(R)) R <- R_ext # nocov
      if (is.null(R)) { # nocov start
        R <- stats::cor(data_mat, method = "pearson", use = "pairwise.complete.obs")
      } # nocov end
    } else {
      # PCA / EFA: compute R then dispatch
      if (cor == "polychoric") {
        poly_out <- tryCatch(
          psych::polychoric(data_mat),
          error = function(e) { # nocov start
            cli::cli_abort(
              c(
                "!" = "{.fn psych::polychoric} failed: {conditionMessage(e)}",
                "i" = "Check that your data contains integer-like columns with few \\
                     distinct values, or use {.code cor = \"pearson\"}."
              )
            )
          } # nocov end
        )
        R <- poly_out$rho
        # Guard against non-positive-definite polychoric matrices (DESIGN.md s.14 remaining)
        min_eig <- min(eigen(R, symmetric = TRUE, only.values = TRUE)$values)
        if (min_eig <= 0) { # nocov start
          cli::cli_warn(
            c(
              "!" = "Polychoric correlation matrix is not positive definite \\
                     (min eigenvalue = {round(min_eig, 4)}).",
              "i" = "Applying smoothing via {.fn psych::cor.smooth}."
            )
          )
          R <- psych::cor.smooth(R)
        } # nocov end
      } else if (missing == "fiml") {
        # M38: FIML correlation via psych::corFiml (cor is guaranteed "pearson"
        # here by .resolve_missing). The estimated R feeds the normal W'RW
        # algebra unchanged; only the fit-index N choice is FIML-specific.
        R <- .corfiml_R(data_mat)
        n_obs_eff <- switch(n_obs_mode,
          total    = nrow(data_mat),
          complete = n_complete
        )
        cli::cli_inform(c(
          "i" = "{.code missing = \"fiml\"}: correlation matrix estimated via \\
                 {.fn psych::corFiml} (full-information ML).",
          if (engine == "efa") {
            c(
              "i" = "Fit indices use N = {n_obs_eff} \\
                     ({.code n_obs = \"{n_obs_mode}\"}); point estimates \\
                     (loadings, edges) are unaffected by this choice.",
              "!" = "Fit indices are approximate: a FIML correlation matrix is \\
                     fed into a normal-theory EFA (a two-step procedure). See \\
                     {.code ?ackwards} ({.arg n_obs})."
            )
          }
        ))
      } else {
        R <- stats::cor(data_mat, method = cor, use = "pairwise.complete.obs")
      }
      engine_out <- switch(engine,
        pca = pca_levels(R,
          k_max = k_max, cor = cor,
          keep_fits = keep_fits
        ),
        efa = efa_levels(R,
          k_max = k_max, fm = fm, n_obs = n_obs_eff,
          cor = cor, keep_fits = keep_fits
        )
      )
      levels_list <- engine_out$levels
      fits_stored <- engine_out$fits
    }
  } # end BRANCH B

  # --- Handle convergence truncation (PCA never truncates; EFA/ESEM may) -----
  k_eff <- length(levels_list)
  if (k_eff < 2L) {
    cli::cli_abort(
      c(
        "!" = "{toupper(engine)} failed to build at least 2 converged levels \\
               (k_eff = {k_eff}; at least 2 are required for a hierarchy).",
        "i" = "Try fewer factors, or check your data for \\
               perfect multicollinearity or near-singular correlation."
      )
    )
  }
  if (k_eff < k_max) {
    cli::cli_warn(
      c(
        "!" = "Hierarchy truncated: {k_max - k_eff} level{?s} did not converge \\
               (requested k_max = {k_max}, built k = {k_eff}).",
        "i" = "Set {.arg k_max = {k_eff}} to suppress this message."
      )
    )
  }

  # --- Lineage matching -------------------------------------------------------
  # Always adjacent: lineage and sign alignment are defined by parent-child
  # relationships between consecutive levels.
  raw_edges_adj <- compute_edges(
    levels      = levels_list,
    R           = R,
    edge_method = "auto",
    pairs       = "adjacent",
    align       = FALSE,
    cut_show    = cut_show
  )

  lineage <- vector("list", k_eff)
  names(lineage) <- as.character(seq_len(k_eff))
  for (ki in seq_len(k_eff - 1L) + 1L) {
    key <- paste0(ki - 1L, ":", ki)
    lineage[[as.character(ki)]] <- match_parents(raw_edges_adj$matrices[[key]])
  }

  # --- Sign alignment (DESIGN.md s.7) -----------------------------------------
  if (align_signs && k_eff >= 1L) {
    loadings_list <- lapply(levels_list, `[[`, "loadings")
    aligned <- .align_signs(loadings_list, raw_edges_adj$matrices, lineage)
    for (ki in seq_len(k_eff)) {
      levels_list[[as.character(ki)]]$loadings <- aligned$loadings[[ki]]
      # Flip weights consistently with loadings
      levels_list[[as.character(ki)]]$scoring$weights <-
        flip_weights(
          levels_list[[as.character(ki)]]$scoring$weights,
          aligned$signs[[ki]]
        )
    }
  }

  # --- Final edge set (adjacent or all-levels Forbes extension) ---------------
  # Recompute using aligned levels_list so skip-level edges inherit correct signs.
  # Adjacent-pair matrices match aligned$edges exactly; recomputing is safe.
  final_edges <- compute_edges(
    levels      = levels_list,
    R           = R,
    edge_method = "auto",
    pairs       = pairs,
    align       = FALSE,
    cut_show    = cut_show
  )

  # --- Assemble aligned edge object -------------------------------------------
  edges_obj <- list(
    matrices = final_edges$matrices,
    tidy     = fill_primary(final_edges$tidy, lineage, levels_list)
  )

  # --- Meta -------------------------------------------------------------------
  conv <- vapply(levels_list, `[[`, logical(1L), "converged")
  meta <- list(
    k_requested       = k_max, # what the user asked for (may exceed k_eff)
    converged_levels  = conv,
    deepest_converged = max(which(conv)),
    pairs             = pairs,
    cut_show          = cut_show,
    ordinal_warned    = is_ordinal,
    missing           = missing_eff,
    n_complete        = n_complete,
    input_type        = input_type,
    estimator         = estimator_eff %||% NA_character_,
    item_labels       = item_labels
  )

  # --- Assemble result --------------------------------------------------------
  x <- new_ackwards(
    call        = cl,
    engine      = engine,
    cor         = cor_eff,
    n_obs       = n_obs_eff,
    k_max       = k_eff, # effective depth (may be < k_max if truncated)
    seed        = seed,
    pkg_version = utils::packageVersion("ackwards"),
    levels      = levels_list,
    edges       = edges_obj,
    lineage     = lineage,
    scores      = NULL,
    fits        = NULL,
    r           = R,
    data        = NULL,
    meta        = meta
  )

  # --- Opt-in storage (scores and raw fits) -----------------------------------
  # Scores use the sign-aligned levels_list so column orientations match loadings.
  # keep_scores=TRUE is blocked at entry for cor_matrix input (already errored above).
  if (isTRUE(keep_scores)) {
    x$scores <- .compute_scores(levels_list, data_mat)
  }
  if (isTRUE(keep_fits)) {
    x$fits <- fits_stored
  }

  x
}

# S3 constructor -- validates structure and attaches class
new_ackwards <- function(
  call, engine, cor, n_obs, k_max, seed, pkg_version,
  levels, edges, lineage, scores, fits, r, data, meta
) {
  structure(
    list(
      call        = call,
      engine      = engine,
      rotation    = "varimax",
      cor         = cor,
      n_obs       = n_obs,
      k_max       = k_max,
      seed        = seed,
      pkg_version = pkg_version,
      levels      = levels,
      edges       = edges,
      lineage     = lineage,
      scores      = scores,
      fits        = fits,
      r           = r,
      data        = data,
      meta        = meta,
      prune       = NULL # populated by prune() (see R/prune.R); NULL until called
    ),
    class = "ackwards"
  )
}
