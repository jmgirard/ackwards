# The rotation argument: per-engine validation, the object stamp, and the
# pass-through to boot_edges() refits.

# Stated here, independently of R/utils.R's .supported_rotations, so the test
# fails if the package table drops one of these four rotations or moves one to
# another engine.
rotation_support <- list(
  pca  = c("varimax", "oblimin", "promax"),
  efa  = c("varimax", "oblimin", "promax"),
  esem = c("varimax", "oblimin", "geomin")
)
all_rotations <- c("varimax", "oblimin", "promax", "geomin")

test_that("every engine-by-rotation pair fits or aborts naming both", {
  skip_if_not_installed("GPArotation")
  skip_if_not_installed("lavaan")

  n_fit <- 0L
  n_abort <- 0L
  for (engine in names(rotation_support)) {
    for (rotation in all_rotations) {
      if (rotation %in% rotation_support[[engine]]) {
        x <- cached(ackwards(sim16,
          k_max = 3, engine = engine,
          rotation = rotation
        ))
        expect_s3_class(x, "ackwards")
        expect_identical(x$rotation, rotation)
        n_fit <- n_fit + 1L
      } else {
        err <- expect_error(
          ackwards(sim16, k_max = 3, engine = engine, rotation = rotation),
          class = "rlang_error"
        )
        msg <- cli::ansi_strip(conditionMessage(err))
        expect_match(msg, "is not available with", fixed = TRUE)
        expect_match(msg, paste0('rotation = "', rotation, '"'), fixed = TRUE)
        expect_match(msg, paste0('engine = "', engine, '"'), fixed = TRUE)
        n_abort <- n_abort + 1L
      }
    }
  }
  # The loop covered the whole cross-product: 9 supported, 3 unsupported.
  expect_identical(c(n_fit, n_abort), c(9L, 3L))
})

test_that("an unknown rotation name fails as an arg_match error", {
  err <- expect_error(
    ackwards(sim16, k_max = 3, rotation = "quartimax"),
    class = "rlang_error"
  )
  expect_match(conditionMessage(err), "must be one of", fixed = TRUE)
  err <- expect_error(
    ackwards(sim16, k_max = 3, rotation = 1),
    class = "rlang_error"
  )
  expect_match(conditionMessage(err), "must be a string", fixed = TRUE)
})

test_that("a psych rotation that loads GPArotation checks for it before fitting", {
  # Record each rlang::check_installed() call: the guard is the routing under
  # test, and the install prompt itself is rlang's behavior. psych loads
  # GPArotation for oblimin on both engines and for promax on EFA only
  # (psych::pca() uses stats::promax()).
  asked <- character()
  local_mocked_bindings(
    check_installed = function(pkg, ...) {
      asked <<- c(asked, pkg)
      invisible(NULL)
    },
    .package = "rlang"
  )
  pairs <- list(c("oblimin", "pca"), c("oblimin", "efa"), c("promax", "efa"))
  for (pr in pairs) {
    expect_identical(.check_rotation(pr[[1]], pr[[2]]), pr[[1]])
  }
  expect_identical(asked, rep("GPArotation", 3L))

  # Control: PCA promax, varimax, and every ESEM rotation never ask for it.
  asked <- character()
  expect_identical(.check_rotation("promax", "pca"), "promax")
  expect_identical(.check_rotation("varimax", "pca"), "varimax")
  expect_identical(.check_rotation("varimax", "efa"), "varimax")
  expect_identical(.check_rotation("geomin", "esem"), "geomin")
  expect_identical(.check_rotation("oblimin", "esem"), "oblimin")
  expect_identical(asked, character())
})

test_that("print() and summary() show the chosen rotation", {
  skip_if_not_installed("GPArotation")
  x <- cached(ackwards(sim16, k_max = 3, rotation = "promax"))
  # The header is a cli definition list, written to the message stream.
  header_text <- function(obj) {
    out <- capture.output(print(obj), type = "message")
    cli::ansi_strip(paste(out, collapse = "\n"))
  }
  expect_match(header_text(x), "Rotation: promax", fixed = TRUE)
  expect_match(header_text(summary(x)), "Rotation: promax", fixed = TRUE)
})

test_that("boot_edges() refits each replicate with the object's rotation", {
  skip_if_not_installed("GPArotation")
  seen <- character()
  real <- .fit_levels_muffled
  local_mocked_bindings(
    .fit_levels_muffled = function(..., rotation = "varimax") {
      seen <<- c(seen, rotation)
      real(..., rotation = rotation)
    }
  )
  # Serial dispatch, so the mocked binding is the one each replicate calls.
  local_mocked_bindings(.boot_lapply = function(X, FUN) lapply(X, FUN))

  x_obl <- cached(ackwards(sim16, k_max = 3, rotation = "oblimin"))
  suppressMessages(suppressWarnings(
    boot_edges(x_obl, sim16, n_boot = 2, seed = 1)
  ))
  expect_identical(seen, c("oblimin", "oblimin"))

  # Control: a varimax object refits varimax.
  seen <- character()
  x_var <- cached(ackwards(sim16, k_max = 3))
  suppressMessages(suppressWarnings(
    boot_edges(x_var, sim16, n_boot = 2, seed = 1)
  ))
  expect_identical(seen, c("varimax", "varimax"))
})

test_that("boot_edges() checks for GPArotation before refitting a rotation that loads it", {
  skip_if_not_installed("GPArotation")
  x_obl <- cached(ackwards(sim16, k_max = 3, rotation = "oblimin"))
  x_pro <- cached(ackwards(sim16, k_max = 3, rotation = "promax"))
  local_mocked_bindings(
    check_installed = function(pkg, ...) stop("planted: ", pkg, " is missing"),
    .package = "rlang"
  )
  expect_error(
    boot_edges(x_obl, sim16, n_boot = 2, seed = 1),
    "planted: GPArotation is missing",
    fixed = TRUE
  )
  # Control: PCA promax runs without GPArotation, so nothing asks for it.
  expect_no_error(suppressMessages(suppressWarnings(
    boot_edges(x_pro, sim16, n_boot = 2, seed = 1)
  )))
})
