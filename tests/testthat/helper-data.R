# Shared simulated datasets loaded by testthat before every test file.

# Memoized fits for tests that need a fitted object but assert nothing about
# fit-time conditions. Keyed on the *evaluated* arguments (rlang::hash), so the
# same call on different data never collides. Warnings/messages are muffled on
# the first evaluation so cache misses and hits behave identically; only use
# cached() where the surrounding test does not expect_warning()/expect_message()
# on the fit itself. Do not use for calls that consume RNG the test relies on
# afterwards (e.g. unseeded suggest_k()); ackwards() fits are deterministic.
# Treat returned objects as read-only: rebinding/modifying the returned value
# is copy-on-modify safe, but environment-bearing components (e.g. lavaan fits
# under keep_fits = TRUE) are shared by reference across cache hits.
.fit_cache <- new.env(parent = emptyenv())

cached <- function(call) {
  expr <- substitute(call)
  env <- parent.frame()
  # Key on the deparsed call text plus the *values* of every symbol it
  # references -- variables AND functions (all.names, not all.vars) -- so
  # `ackwards(d, ...)` on different `d` never collides, and neither do two
  # same-named local helpers called with identical text and data. Functions
  # are keyed by formals + body (hashing a closure's value would serialize
  # its environment, which is unstable whenever anything in the enclosing
  # test frame mutates). Nothing in the expression is evaluated for the key
  # (a nested call evaluated here would run twice and outside the suppress
  # wrappers below).
  syms <- unique(all.names(expr))
  syms <- syms[vapply(syms, exists, logical(1), envir = env)]
  vals <- lapply(syms, function(s) {
    v <- get(s, envir = env)
    if (is.function(v)) list(formals(v), body(v)) else v
  })
  key <- rlang::hash(list(deparse(expr), syms, vals))
  if (is.null(.fit_cache[[key]])) {
    .fit_cache[[key]] <- suppressMessages(suppressWarnings(eval(expr, env)))
  }
  .fit_cache[[key]]
}

.make_esem_data <- function(seed = 42, n = 200) {
  set.seed(seed)
  f1 <- rnorm(n)
  f2 <- rnorm(n)
  data.frame(
    x1 = f1 + 0.4 * rnorm(n), x2 = f1 + 0.4 * rnorm(n),
    x3 = f1 + 0.4 * rnorm(n), x4 = f2 + 0.4 * rnorm(n),
    x5 = f2 + 0.4 * rnorm(n), x6 = f2 + 0.4 * rnorm(n)
  )
}

.make_ordinal_data <- function(seed = 42, n = 300) {
  set.seed(seed)
  f1 <- rnorm(n)
  f2 <- rnorm(n)
  d <- data.frame(
    x1 = f1 + 0.4 * rnorm(n), x2 = f1 + 0.4 * rnorm(n),
    x3 = f1 + 0.4 * rnorm(n), x4 = f2 + 0.4 * rnorm(n),
    x5 = f2 + 0.4 * rnorm(n), x6 = f2 + 0.4 * rnorm(n)
  )
  breaks <- c(-Inf, -1, -0.5, 0.5, 1, Inf)
  as.data.frame(lapply(d, function(x) as.integer(cut(x, breaks))))
}
