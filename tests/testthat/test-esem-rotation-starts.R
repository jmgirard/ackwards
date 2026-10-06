# ESEM counts the lavaan rotation starts that did not converge, under every
# rotation (M098). lavaan rotates each level from rotation.args$rstarts random
# starts and warns once per start that did not converge. Some failed starts:
# the level is kept with a counted warning. Every start failed: the kept
# rotation did not converge, so the hierarchy ends at the level before.

# Plants a rotation failure at level `k_plant` of the real lavaan::efa call.
# `max_iter` caps that level's rotation iterations inside lavaan, so real starts
# fail; `n_warn` adds that many warnings in lavaan's text after the real fit, a
# deterministic partial failure. The version check in .esem_rotation_args()
# reads lavaan::efa's formals, so it stays on the real ones while efa is
# replaced (M90 lesson), and the levels run serially in this process.
.local_rotation_plant <- function(k_plant, max_iter = NULL, n_warn = 0L, env = parent.frame()) {
  real_efa <- lavaan::efa
  real_formals <- names(formals(real_efa))
  list_form <- "rotation_args" %in% real_formals
  real_rot <- .esem_rotation_args
  local_mocked_bindings(
    .esem_rotation_args = function(rotation) real_rot(rotation, real_formals),
    .esem_lapply = function(X, FUN) lapply(X, FUN),
    .env = env
  )
  local_mocked_bindings(
    efa = function(...) {
      args <- list(...)
      planted <- args$nfactors == k_plant
      if (planted && !is.null(max_iter)) {
        if (list_form) {
          args$rotation <- c(as.list(args$rotation), max_iter = max_iter)
        } else {
          args$rotation.args <- c(args$rotation.args, max_iter = max_iter)
        }
      }
      fit <- do.call(real_efa, args)
      if (planted) {
        # Every second planted warning carries a line break inside the phrase,
        # as message formatting can.
        for (i in seq_len(n_warn)) {
          msg <- if (i %% 2L == 0L) {
            "GP rotation algorithm did not\n  converge after 10000 iterations"
          } else {
            "GP rotation algorithm did not converge after 10000 iterations"
          }
          warning(msg, call. = FALSE)
        }
      }
      fit
    },
    .package = "lavaan",
    .env = env
  )
}

# The start total lavaan used, from the options of a kept fit.
.rstarts_of <- function(x, k) {
  lavaan::lavInspect(x$fits[[as.character(k)]], "options")$rotation.args$rstarts
}

test_that("esem: a level where some rotation starts fail is kept, with a counted warning", {
  skip_if_not_installed("lavaan")
  # local() scopes each plant to one fit, so the next plant wraps the real
  # functions, not the previous mock.
  for (rot in c("varimax", "oblimin")) {
    local({
      .local_rotation_plant(3L, n_warn = 7L)
      w <- .warnings_of(x <- suppressMessages(
        ackwards(sim16, k_max = 4, engine = "esem", rotation = rot, seed = 1, keep_fits = TRUE)
      ))
      n_starts <- .rstarts_of(x, 3)
      expect_gt(n_starts, 7)
      expect_identical(names(x$levels), c("1", "2", "3", "4"), label = rot)
      hits <- grep("random starts", w, fixed = TRUE, value = TRUE)
      expect_length(grep("k = 3", hits, fixed = TRUE), 1L)
      expect_match(hits, "at k = 3:", fixed = TRUE, all = TRUE)
      expect_match(hits, sprintf("7 of %d random starts did not converge", n_starts), fixed = TRUE)
      expect_match(hits, "The level is kept.", fixed = TRUE)
    })
  }
})

test_that("esem: a level where every rotation start fails ends the hierarchy", {
  skip_if_not_installed("lavaan")
  fits <- list(
    list(data = sim16, rotation = "varimax", cor = "pearson", estimator = "ML"),
    list(data = sim16, rotation = "oblimin", cor = "pearson", estimator = "ML"),
    list(data = sim16, rotation = "geomin", cor = "pearson", estimator = "ML"),
    list(data = bfi25, rotation = "varimax", cor = "polychoric", estimator = "WLSMV")
  )
  for (f in fits) {
    local({
      lab <- paste(f$rotation, f$estimator)
      .local_rotation_plant(3L, max_iter = 2L)
      w <- .warnings_of(x <- suppressMessages(ackwards(
        f$data,
        k_max = 4, engine = "esem", cor = f$cor, rotation = f$rotation,
        seed = 1, keep_fits = TRUE
      )))
      expect_identical(x$meta$estimator, f$estimator, label = lab)
      expect_identical(names(x$levels), c("1", "2"), label = lab)
      n_starts <- .rstarts_of(x, 2)
      hit <- grep("rotation did not converge at k = 3", w, fixed = TRUE, value = TRUE)
      expect_false(any(grepl("at k = 2", w, fixed = TRUE)), label = lab)
      expect_length(hit, 1L)
      expect_match(hit, sprintf("%d of %d starts did not converge", n_starts, n_starts), fixed = TRUE)
      expect_match(hit, "the kept rotation did not converge", fixed = TRUE)
      expect_match(hit, "Truncating hierarchy at level 2.", fixed = TRUE)
    })
  }
})

test_that("esem: default varimax fits whose rotation converges raise no rotation warning", {
  skip_if_not_installed("lavaan")
  w <- c(
    .warnings_of(suppressMessages(ackwards(sim16, k_max = 4, engine = "esem", seed = 1))),
    .warnings_of(suppressMessages(ackwards(bfi25, k_max = 4, engine = "esem", cor = "polychoric", seed = 1)))
  )
  expect_false(any(grepl("random starts", w, fixed = TRUE)))
  expect_false(any(grepl("rotation algorithm", w, fixed = TRUE)))
  expect_false(any(grepl("did not converge", w, fixed = TRUE)))
})

test_that(".esem_rotation_starts() reads rstarts and counts rstarts = 0 or an unreadable total as one start", {
  skip_if_not_installed("lavaan")
  local_mocked_bindings(lavInspect = function(object, what) object, .package = "lavaan")
  expect_identical(.esem_rotation_starts(list(rotation.args = list(rstarts = 30L))), 30L)
  expect_identical(.esem_rotation_starts(list(rotation.args = list(rstarts = 0L))), 1L)
  # A total that cannot be read counts as one start.
  expect_identical(.esem_rotation_starts(list()), 1L)
  expect_identical(.esem_rotation_starts(list(rotation.args = list(rstarts = "30"))), 1L)
  local_mocked_bindings(lavInspect = function(object, what) stop("no options"), .package = "lavaan")
  expect_identical(.esem_rotation_starts(NULL), 1L)
})
