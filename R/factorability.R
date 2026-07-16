#' Screen a dataset for factorability and sampling adequacy
#'
#' @description
#' Before extracting a hierarchy, it is worth asking two prior questions:
#' *is this correlation matrix even worth factoring*, and *is the sample large
#' enough to trust the answer*. `factorability()` reports the standard
#' diagnostics for both, one object you can print or pull values from:
#'
#' * **KMO** -- the Kaiser-Meyer-Olkin measure of sampling adequacy, overall
#'   and per variable (via [psych::KMO()]). It contrasts correlations with
#'   partial correlations: a low value means the variables share little common
#'   variance once other variables are partialled out, so a factor model has
#'   little to recover.
#' * **Bartlett's test of sphericity** -- tests whether the correlation matrix
#'   differs from the identity (via [psych::cortest.bartlett()]). Needs `N`; a
#'   non-significant result says the correlations are too weak to factor.
#' * **Sample size** -- `N`, the number of variables `p`, and the `N:p` ratio.
#' * **Ledermann bound** -- the largest number of *common factors* identifiable
#'   from `p` variables (a hard limit for EFA/ESEM, though not for PCA).
#'
#' **Read these as conventions, not verdicts.** Every cutoff here -- Kaiser's
#' KMO bands, the 5:1 / 10:1 `N:p` rules of thumb, Bartlett significance at
#' `.05` -- is a widely repeated *rule of thumb*, not a settled threshold, and
#' each is contested in the methodological literature (the required `N` in
#' particular depends on communalities and factor overdetermination far more
#' than on any fixed ratio; MacCallum et al. 1999). `factorability()`
#' deliberately reports the numbers and their conventional bands rather than
#' returning a pass/fail flag. [ackwards()] runs the same screen internally and
#' warns only at the genuinely consequential extreme (KMO `< 0.5`, `N:p < 5`).
#'
#' @param data A data frame or numeric matrix of item responses, **or** a
#'   correlation matrix (detected automatically). KMO and Bartlett are computed
#'   on the correlation matrix in the basis you name via `cor`.
#' @param cor The correlation basis to screen on: `"pearson"` (default),
#'   `"spearman"`, or `"polychoric"`. Ignored when `data` is already a
#'   correlation matrix. Use the same basis you plan to pass to [ackwards()].
#' @param n_obs Number of observations. Taken from `nrow(data)` for raw data
#'   (a supplied value is ignored with a warning). **Required** for the N-based
#'   diagnostics (Bartlett, `N:p`) when `data` is a correlation matrix; those
#'   are reported as `NA` when it is absent.
#'
#' @return An object of class `factorability` (a named list): `cor`, `n_obs`,
#'   `n_vars`, `np_ratio`, `kmo_overall`, `kmo_items` (a data frame of per-item
#'   MSA), `bartlett` (a list of `chisq`/`df`/`p_value`, or `NULL` when `N` is
#'   unknown), and `ledermann` (the bound). Print it for a banded summary; index
#'   the list for the raw values.
#'
#' @references
#' Kaiser, H. F. (1974). An index of factorial simplicity. *Psychometrika*,
#' 39(1), 31--36.
#'
#' MacCallum, R. C., Widaman, K. F., Zhang, S., & Hong, S. (1999). Sample size
#' in factor analysis. *Psychological Methods*, 4(1), 84--99.
#'
#' @seealso [check_items()] (per-item screening), [suggest_k()] (how many
#'   factors), and [comparability()] (split-half replicability) -- the other
#'   pre-analysis diagnostics; [ackwards()], which runs this screen internally.
#'
#' @examples
#' # Continuous data, Pearson basis
#' factorability(sim16)
#'
#' # From a correlation matrix: supply n_obs for the N-based diagnostics
#' R <- cor(sim16)
#' factorability(R, n_obs = nrow(sim16))
#'
#' @export
factorability <- function(data, cor = c("pearson", "spearman", "polychoric"),
                          n_obs = NULL) {
  cor <- rlang::arg_match(cor)

  if (!is.null(n_obs)) n_obs <- .check_count(n_obs, "n_obs")

  if (.is_cor_matrix(data)) {
    R <- .validate_cor_matrix(as.matrix(data))$R
    n_obs_eff <- if (is.null(n_obs)) NA_integer_ else as.integer(n_obs)
    basis <- cor
  } else {
    if (is.matrix(data)) .check_maybe_cov_matrix(data)
    data_mat <- .as_numeric_matrix(data)
    if (ncol(data_mat) < 2L) {
      cli::cli_abort("{.arg data} must have at least two columns to screen.")
    }
    # A constant item makes cor()/KMO undefined; error early and point at the
    # per-item screen (mirrors ackwards() and check_items()).
    scr <- .screen_items(as.data.frame(data_mat), cor)
    constant_items <- scr$item[scr$flag == "constant"]
    if (length(constant_items) > 0L) {
      cli::cli_abort(c(
        "!" = "{length(constant_items)} item{?s} ha{?s/ve} no variance: \\
               {.val {constant_items}}.",
        "i" = "Drop {cli::qty(constant_items)}{?it/them} first; \\
               {.fn check_items} screens the rest."
      ))
    }
    if (!is.null(n_obs)) {
      cli::cli_warn(c(
        "!" = "{.arg n_obs} is ignored for raw data (N is {.code nrow(data)}).",
        "i" = "Remove {.arg n_obs} to suppress this warning."
      ))
    }
    R <- .factorability_R(data_mat, cor)
    n_obs_eff <- nrow(data_mat)
    basis <- cor
  }

  out <- .compute_factorability(R, n_obs_eff)
  out$cor <- basis
  structure(out, class = "factorability")
}

# Correlation matrix on the requested basis for raw data. Pearson/Spearman go
# through stats::cor (pairwise); polychoric through psych, with the same
# targeted failure message ackwards() uses.
.factorability_R <- function(data_mat, cor) {
  if (cor == "polychoric") {
    poly <- tryCatch(
      suppressWarnings(psych::polychoric(data_mat)),
      error = function(e) { # nocov start
        cli::cli_abort(c(
          "!" = "{.fn psych::polychoric} failed: {conditionMessage(e)}",
          "i" = "Screen the items with {.fn check_items}; a sparse or \\
                 near-constant category is the usual cause."
        ))
      } # nocov end
    )
    return(poly$rho)
  }
  stats::cor(data_mat, method = cor, use = "pairwise.complete.obs")
}

# Shared factorability core used by factorability() and the ackwards() screen.
# Returns the diagnostic values for a correlation matrix R and (possibly NA)
# sample size. KMO is wrapped in tryCatch: on a near-singular R the underlying
# solve() fails, in which case MSA is reported as NA rather than aborting.
.compute_factorability <- function(R, n_obs) {
  p <- nrow(R)
  # psych::KMO handles a singular R internally (it reports MSA anyway); the NULL
  # branch is purely defensive. Suppress its own notes -- factorability() and
  # the near-singular check own the messaging.
  kmo <- tryCatch(
    suppressWarnings(suppressMessages(psych::KMO(R))),
    error = function(e) NULL
  )
  if (is.null(kmo)) {
    kmo_overall <- NA_real_ # nocov
    msai <- stats::setNames(rep(NA_real_, p), rownames(R)) # nocov
  } else {
    kmo_overall <- unname(kmo$MSA)
    msai <- kmo$MSAi
    if (is.null(msai)) msai <- stats::setNames(rep(NA_real_, p), rownames(R)) # nocov
  }
  kmo_items <- data.frame(
    item = rownames(R) %||% paste0("V", seq_len(p)),
    msa = unname(msai),
    stringsAsFactors = FALSE
  )

  bartlett <- if (!is.na(n_obs)) {
    bt <- suppressWarnings(psych::cortest.bartlett(R, n = n_obs))
    list(chisq = unname(bt$chisq), df = unname(bt$df), p_value = unname(bt$p.value))
  } else {
    NULL
  }

  list(
    cor = NA_character_,
    n_obs = if (is.na(n_obs)) NA_integer_ else as.integer(n_obs),
    n_vars = p,
    np_ratio = if (is.na(n_obs)) NA_real_ else n_obs / p,
    kmo_overall = kmo_overall,
    kmo_items = kmo_items,
    bartlett = bartlett,
    ledermann = .ledermann_bound(p)
  )
}

# Conventional (Kaiser 1974) band label for an overall/per-item KMO value.
# These bands are a rule of thumb, not a settled threshold -- see the printout
# and ?factorability.
.kmo_band <- function(msa) {
  if (is.na(msa)) {
    return("uncomputable")
  }
  if (msa >= 0.90) {
    "marvellous"
  } else if (msa >= 0.80) {
    "meritorious"
  } else if (msa >= 0.70) {
    "middling"
  } else if (msa >= 0.60) {
    "mediocre"
  } else if (msa >= 0.50) {
    "miserable"
  } else {
    "unacceptable"
  }
}

# Consequential-threshold constants, shared by the printout and the internal
# ackwards() screen so the two never drift apart.
.factorability_thresholds <- function() {
  list(kmo_warn = 0.50, np_warn = 5, n_floor = 100)
}

# Internal factorability screen run by ackwards(). Emits at most two loud
# warnings (Invariant 6): a Ledermann under-identification warning for
# EFA/ESEM when k_max exceeds the bound (PCA is exempt -- components extract
# exactly up to p, no latent-model df constraint), and a single consolidated
# adequacy warning when KMO or N:p breach the consequential thresholds. Never
# aborts and never alters the fit -- it only advises.
.factorability_screen <- function(R, n_obs, p, k_max, engine) {
  thr <- .factorability_thresholds()

  # (a) Ledermann identification ceiling -- EFA/ESEM only.
  if (engine %in% c("efa", "esem")) {
    bound <- .ledermann_bound(p)
    if (k_max > bound) {
      first_bad <- bound + 1L
      cli::cli_warn(
        c(
          "!" = "{.arg k_max} = {k_max} exceeds the Ledermann bound for \\
                 {p} variable{?s}: at most {bound} common factor{?s} \\
                 {?is/are} identifiable ({.code engine = \"{engine}\"}).",
          "i" = "Level{?s} {first_bad}+ {?is/are} under-identified \\
                 (model df < 0) and may fail to converge or yield a \\
                 degenerate solution.",
          "i" = "Reduce {.arg k_max} to {bound} or fewer, or use \\
                 {.code engine = \"pca\"} (components have no such limit)."
        ),
        .frequency = "once",
        .frequency_id = "ackwards_ledermann"
      )
    }
  }

  # (b) Consolidated sampling-adequacy warning. Reuse the shared core so the
  # KMO call inherits its NA-safe tryCatch and its suppress wrappers -- an
  # unwrapped psych::KMO() here would leak psych's "matrix is not invertible"
  # message on a near-singular R, on top of ackwards()' own near-singular
  # warning.
  kmo_overall <- .compute_factorability(R, n_obs)$kmo_overall
  np_ratio <- if (is.na(n_obs)) NA_real_ else n_obs / p

  kmo_bad <- !is.na(kmo_overall) && kmo_overall < thr$kmo_warn
  np_bad <- !is.na(np_ratio) && np_ratio < thr$np_warn
  n_bad <- !is.na(n_obs) && n_obs < thr$n_floor

  if (kmo_bad || np_bad || n_bad) {
    bullets <- c("!" = "The data may be poorly suited to factor analysis:")
    if (kmo_bad) {
      bullets <- c(bullets, "*" = "KMO = {.val {round(kmo_overall, 2)}} \\
                    (below {thr$kmo_warn}: the variables share little common variance).")
    }
    if (np_bad || n_bad) {
      bullets <- c(bullets, "*" = "N = {n_obs} for {p} variable{?s} \\
                    (N:p = {.val {round(np_ratio, 1)}}; conventional floors are \\
                    N:p >= 5 and N >= {thr$n_floor}).")
    }
    bullets <- c(
      bullets,
      "i" = "These are contested rules of thumb, not hard cutoffs -- \\
             run {.fn factorability} for the full report."
    )
    cli::cli_warn(bullets,
      .frequency = "once", .frequency_id = "ackwards_factorability"
    )
  }
  invisible(NULL)
}

#' Print a factorability screen
#'
#' @param x A `factorability` object (produced by [factorability()]).
#' @param ... Ignored.
#' @return `x` invisibly.
#' @export
print.factorability <- function(x, ...) {
  cli::cli_h1("Factorability screen ({.pkg ackwards})")

  n_txt <- if (is.na(x$n_obs)) {
    cli::col_grey("not supplied")
  } else {
    format(x$n_obs, big.mark = ",")
  }
  np_txt <- if (is.na(x$np_ratio)) {
    cli::col_grey("-- (N not supplied)")
  } else {
    paste0(round(x$np_ratio, 1), ":1")
  }
  cli::cli_dl(c(
    "Basis" = cli::style_bold(x$cor),
    "Observations" = n_txt,
    "Variables" = as.character(x$n_vars),
    "N:p ratio" = np_txt
  ))

  cli::cli_h2("Sampling adequacy")
  if (is.na(x$kmo_overall)) {
    cli::cli_text(cli::col_yellow( # nocov start
      "! KMO could not be computed (correlation matrix near-singular)."
    )) # nocov end
  } else {
    band <- .kmo_band(x$kmo_overall)
    cli::cli_text(
      "Overall KMO: {.strong {round(x$kmo_overall, 2)}} ",
      "({.emph {band}})"
    )
    # Flag items below the same KMO threshold the internal screen warns at, so
    # the printout and the ackwards() screen never drift (single source below).
    low_cut <- .factorability_thresholds()$kmo_warn
    low_cut_txt <- .format_r(low_cut)
    low <- x$kmo_items[!is.na(x$kmo_items$msa) & x$kmo_items$msa < low_cut, , drop = FALSE]
    if (nrow(low) > 0L) {
      cli::cli_text(cli::col_yellow(
        "! {nrow(low)} item{?s} with MSA < {low_cut_txt}: {.val {low$item}}"
      ))
    }
  }

  cli::cli_text()
  if (is.null(x$bartlett)) {
    cli::cli_text(cli::col_grey(
      "Bartlett's test: not computed (N not supplied)."
    ))
  } else {
    b <- x$bartlett
    sig <- if (!is.na(b$p_value) && b$p_value < .05) "" else " -- NOT significant"
    p_txt <- if (!is.na(b$p_value) && b$p_value < .001) {
      "p < .001"
    } else {
      paste0("p = ", .format_r(b$p_value, 3L))
    }
    cli::cli_text(
      "Bartlett's test of sphericity: ",
      "chi-square({b$df}) = {round(b$chisq, 1)}, ",
      "{p_txt}{sig}"
    )
  }

  cli::cli_h2("Identifiability")
  cli::cli_text(
    "Ledermann bound: at most {.strong {x$ledermann}} common factor{?s} ",
    "{?is/are} identifiable from {x$n_vars} variable{?s} (EFA/ESEM; PCA is unbounded)."
  )

  cli::cli_rule()
  cli::cli_text(cli::col_grey(
    "KMO bands (Kaiser 1974), the N:p >= 5/10 rules, and Bartlett at .05 are \\
     widely used *rules of thumb*, not settled thresholds -- required N depends \\
     on communalities and factor overdetermination (MacCallum et al. 1999). \\
     Read the numbers, not a pass/fail."
  ))
  invisible(x)
}
