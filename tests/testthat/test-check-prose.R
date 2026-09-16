# Guards the plain-English prose checker (tools/check-prose.R). The checker
# and the lists it reads live under tools/, which is .Rbuildignore'd, so this
# test SKIPS in the built package and runs only in the source checkout, as
# test-vignette-freshness.R does.

.prose_env <- function() {
  root <- normalizePath(test_path("..", ".."), mustWork = FALSE)
  checker <- file.path(root, "tools", "check-prose.R")
  skip_if_not(file.exists(checker), "tools/check-prose.R absent (built package)")
  env <- new.env()
  sys.source(checker, envir = env)
  env$tools_dir <- file.path(root, "tools")
  env
}

.write_fixture <- function(lines, ext) {
  path <- tempfile(fileext = ext)
  writeLines(lines, path)
  path
}

.run <- function(env, path) {
  env$check_prose(
    path,
    banned = env$read_prose_list("prose-banned.txt", env$tools_dir),
    abbrev = env$read_prose_list("prose-abbrev.txt", env$tools_dir)
  )
}

long_sentence <- paste(c("Word", rep("word", 30)), collapse = " ")

test_that("each report class fires at each location, for each dash form", {
  env <- .prose_env()
  dashes <- c("—", "–", " -- ")
  classes <- c("em dash", "en dash", "double hyphen")

  for (d in seq_along(dashes)) {
    dash <- dashes[[d]]

    # Markdown: YAML-adjacent prose, heading, bullet, table cell, link text.
    md <- .write_fixture(c(
      "---", "title: x", "---",
      paste0("Right after the yaml", dash, "prose."),
      paste0("# A heading", dash, "with a dash; and a semicolon"),
      paste0("- a bullet", dash, "item; that is robust."),
      "| a | b |", "|---|---|",
      paste0("| a cell", dash, "here; | x |"),
      paste0("[link text", dash, "with dash](https://example.org)."),
      paste0(long_sentence, ".")
    ), ".Rmd")
    res <- .run(env, md)
    expect_setequal(
      res$line[res$class == classes[[d]]], c(4L, 5L, 6L, 9L, 10L)
    )
    expect_setequal(res$line[res$class == "semicolon"], c(5L, 6L, 9L))
    expect_true(6L %in% res$line[res$class == "banned phrase: robust"])
    expect_true(any(grepl("^sentence over 30 words", res$class) & res$line == 11L))
    # No line reports twice for the same class.
    expect_false(any(duplicated(res[c("line", "class")])))

    # Roxygen: @param and @details.
    rf <- .write_fixture(c(
      "#' Title",
      "#'",
      paste0("#' @param x a param", dash, "with a dash; and more"),
      "#' @details",
      paste0("#' Some details", dash, "here, robust ones."),
      paste0("#' ", long_sentence, "."),
      "#' @export",
      "f <- function(x) x"
    ), ".R")
    res <- .run(env, rf)
    expect_setequal(res$line[res$class == classes[[d]]], c(3L, 5L))
    expect_equal(res$line[res$class == "semicolon"], 3L)
    expect_equal(res$line[res$class == "banned phrase: robust"], 5L)
    expect_true(any(grepl("^sentence over", res$class) & res$line == 6L))

    # NEWS entry: only the development section is swept.
    nf <- .write_fixture(c(
      "# pkg (development version)",
      "",
      paste0("* An entry", dash, "with a dash; robust."),
      paste0("* ", long_sentence, "."),
      "",
      "# pkg 1.0.0",
      "",
      paste0("* Old entry", dash, "left alone; robust.")
    ), ".md")
    file.rename(nf, file.path(dirname(nf), "NEWS.md"))
    nf <- file.path(dirname(nf), "NEWS.md")
    res <- .run(env, nf)
    expect_equal(res$line[res$class == classes[[d]]], 3L)
    expect_equal(res$line[res$class == "semicolon"], 3L)
    expect_equal(res$line[res$class == "banned phrase: robust"], 3L)
    expect_true(any(grepl("^sentence over", res$class) & res$line == 4L))
    expect_false(8L %in% res$line)

    # DESCRIPTION: only the Description field is swept.
    dd <- tempfile()
    dir.create(dd)
    df <- file.path(dd, "DESCRIPTION")
    writeLines(c(
      "Package: x",
      paste0("Title: A title", dash, "with a dash; robust"),
      paste0("Description: A description", dash, "with a dash;"),
      paste0("    continued, robust. ", long_sentence, "."),
      "License: MIT"
    ), df)
    res <- .run(env, df)
    expect_equal(res$line[res$class == classes[[d]]], 3L)
    expect_equal(res$line[res$class == "semicolon"], 3L)
    expect_equal(res$line[res$class == "banned phrase: robust"], 4L)
    expect_true(any(grepl("^sentence over", res$class) & res$line == 4L))
    expect_false(2L %in% res$line)
  }
})

test_that("the checker is silent on every protected location", {
  env <- .prose_env()
  all_classes <- paste0("— – -- ; robust in order to ", long_sentence, ".")

  md <- .write_fixture(c(
    "---",
    paste0("title: ", all_classes),
    "---",
    paste0("<!-- precompute-stamp: source=x.Rmd.orig md5=abc ", all_classes, " -->"),
    "",
    paste0("Inline `", all_classes, "` span."),
    paste0("Double ``", all_classes, "`` span."),
    "A span crossing `a line",
    paste0(all_classes, "` break."),
    "",
    "```{r}",
    all_classes,
    "```",
    "",
    "---",
    "",
    "Abbreviations e.g. this, i.e. that, et al. those, vs. these, p. 3, pp. 4, cf. it, Fig. 2, No. 1 and",
    "some more words to pass thirty if the periods above split this sentence apart."
  ), ".Rmd")
  res <- .run(env, md)
  expect_equal(nrow(res), 0L, info = paste(res$class, res$line, collapse = "; "))

  rf <- .write_fixture(c(
    "#' Title",
    "#'",
    paste0("#' Inline `", all_classes, "` span."),
    "#'",
    "#' ```",
    paste0("#' ", all_classes),
    "#' ```",
    "#'",
    "#' \\preformatted{",
    paste0("#' ", all_classes),
    "#' }",
    "#'",
    "#' @examples",
    paste0("#' # ", all_classes),
    paste0("#' f(1) # ", all_classes),
    "#' @export",
    paste0("f <- function(x) x # ", all_classes)
  ), ".R")
  res <- .run(env, rf)
  expect_equal(nrow(res), 0L, info = paste(res$class, res$line, collapse = "; "))
})

test_that("a multi-word banned phrase is found when it wraps across a line", {
  env <- .prose_env()
  md <- .write_fixture(c(
    "We did this in order",
    "to make it work. It is worth",
    "noting that robustness matters.",
    "",
    "in order",
    "",
    "to: a paragraph break never joins a phrase."
  ), ".Rmd")
  res <- .run(env, md)
  expect_equal(res$line[res$class == "banned phrase: in order to"], 1L)
  expect_equal(res$line[res$class == "banned phrase: it is worth noting"], 2L)
  expect_equal(res$line[res$class == "banned phrase: robustness"], 3L)
  expect_false(any(res$line %in% c(5L, 7L)))
})

test_that("an unclosed span, fence, or YAML header is an error, not silence", {
  env <- .prose_env()
  stray <- .write_fixture(c(
    "A stray backtick ` here opens nothing.",
    "",
    "em dash — and robust.",
    "",
    "closes it ` finally."
  ), ".Rmd")
  expect_error(.run(env, stray), "unmatched backtick", fixed = TRUE)
  expect_error(.run(env, stray), ":1: ", fixed = TRUE)

  fence <- .write_fixture(c("Prose.", "```{r}", "code", "em dash — and robust."), ".Rmd")
  expect_error(.run(env, fence), "fenced code block is never closed", fixed = TRUE)
  expect_error(.run(env, fence), ":2: ", fixed = TRUE)

  rox <- .write_fixture(c(
    "#' Title", "#'", "#' ```", "#' code", "#' em dash — and robust.", "f <- 1"
  ), ".R")
  expect_error(.run(env, rox), "roxygen fenced code is never closed", fixed = TRUE)
  expect_error(.run(env, rox), ":3: ", fixed = TRUE)

  yaml <- .write_fixture(c("---", "title: x", "", "em dash — and robust."), ".Rmd")
  expect_error(.run(env, yaml), "YAML header", fixed = TRUE)
})

test_that("a paragraph never spans a code line, a one-line preformatted block, or a leading rule", {
  env <- .prose_env()
  # Two roxygen blocks separated by a code line are two paragraphs: a phrase
  # split across them is not matched, and a backtick in each is not a span.
  rf <- .write_fixture(c(
    "#' Alpha beta not",
    "f <- 1",
    "#' just gamma with a stray ` mark",
    "g <- 2",
    "#' and another ` mark here"
  ), ".R")
  expect_error(.run(env, rf), "unmatched backtick", fixed = TRUE)
  rf2 <- .write_fixture(c("#' Alpha beta not", "f <- 1", "#' just gamma"), ".R")
  expect_equal(nrow(.run(env, rf2)), 0L)

  one_line <- .write_fixture(c(
    "#' Title",
    "#'",
    "#' \\preformatted{x <- 1}",
    "#' Prose with an em dash — after it.",
    "f <- function(x) x"
  ), ".R")
  res <- .run(env, one_line)
  expect_equal(res$line[res$class == "em dash"], 4L)

  rule_first <- .write_fixture(c("---", "", "Prose; with a semicolon."), ".Rmd")
  res <- .run(env, rule_first)
  expect_equal(res$line[res$class == "semicolon"], 3L)

  yaml_first <- .write_fixture(c("---", "title: x", "---", "Prose; here."), ".Rmd")
  res <- .run(env, yaml_first)
  expect_equal(res$line[res$class == "semicolon"], 4L)
})

test_that("a sentence is counted between terminators, not per line", {
  env <- .prose_env()
  md <- .write_fixture(c(
    "Short one. Short two? Short three! Then a fourth that wraps",
    "onto the next line and stays short.",
    paste(rep("w", 16), collapse = " "),
    paste0(paste(rep("w", 15), collapse = " "), ".")
  ), ".Rmd")
  res <- .run(env, md)
  expect_equal(nrow(res), 1L)
  expect_equal(res$line, 3L)
  expect_match(res$class, "^sentence over 30 words \\(31\\)")
})

test_that("the default domain and the directory expansion are non-empty", {
  env <- .prose_env()
  root <- normalizePath(test_path("..", ".."), mustWork = FALSE)
  dom <- env$default_prose_domain(root)
  expect_true(all(file.exists(dom)))
  expect_gt(length(dom), 10L)
  rs <- env$expand_prose_paths(file.path(root, "R"))
  expect_gt(length(rs), 10L)
  expect_true(all(grepl("\\.R$", rs)))
})

test_that("check_code_unchanged sees only code and reports a changed code line", {
  env <- .prose_env()
  skip_if_not(nzchar(Sys.which("git")), "git not available")
  root <- tempfile("prose-git-")
  dir.create(file.path(root, "R"), recursive = TRUE)
  old <- setwd(root)
  on.exit(setwd(old), add = TRUE)
  git <- function(...) system2("git", c(...), stdout = TRUE, stderr = TRUE)
  git("init", "-q", "-b", "master")
  git("config", "user.email", "t@t")
  git("config", "user.name", "t")
  writeLines(c(
    "#' Title -- old prose",
    "#' @examples",
    "#' f(1)",
    "#' @export",
    "f <- function(x) x"
  ), "R/f.R")
  writeLines(c(
    "Package: x",
    "Title: A title",
    "Description: Old prose.",
    "License: MIT"
  ), "DESCRIPTION")
  writeLines(c(
    "---", "output: github_document", "---",
    "Old prose `r 1 + 1` here.",
    "```{r}", "x <- 1", "```"
  ), "README.Rmd")
  dir.create("vignettes")
  vig <- c(
    "---", "title: v", "---",
    "Old vignette prose `r 2 + 2` here.",
    "```{r chunk, eval = FALSE}", "y <- 1", "```"
  )
  writeLines(vig, "vignettes/v.Rmd.orig")
  git("add", "-A")
  git("commit", "-q", "-m", "base")
  git("checkout", "-q", "-b", "work")

  # Prose-only edits: silent.
  writeLines(c(
    "#' Title, new prose.",
    "#' @examples",
    "#' f(1)",
    "#' @export",
    "f <- function(x) x"
  ), "R/f.R")
  writeLines(c(
    "Package: x",
    "Title: A title",
    "Description: New prose.",
    "License: MIT"
  ), "DESCRIPTION")
  writeLines(c(
    "---", "output: github_document", "---",
    "New prose `r 1 + 1` here.",
    "```{r}", "x <- 1", "```"
  ), "README.Rmd")
  vig[4] <- "New vignette prose `r 2 + 2` here."
  writeLines(vig, "vignettes/v.Rmd.orig")
  git("commit", "-q", "-am", "prose")
  expect_equal(env$check_code_unchanged("master", root), character(0L))

  # Each code edit is reported by name.
  writeLines(c(
    "#' Title, new prose.",
    "#' @examples",
    "#' f(2)",
    "#' @export",
    "f <- function(x) x"
  ), "R/f.R")
  expect_match(env$check_code_unchanged("master", root), "R/f.R: code line 2")
  git("checkout", "-q", "--", "R/f.R")

  writeLines(c(
    "#' Title, new prose.",
    "#' @examples",
    "#' f(1)",
    "#' @export",
    "f <- function(x) x + 1"
  ), "R/f.R")
  expect_match(env$check_code_unchanged("master", root), "R/f.R: code line 3")
  git("checkout", "-q", "--", "R/f.R")

  writeLines(c(
    "Package: x",
    "Title: Another title",
    "Description: New prose.",
    "License: MIT"
  ), "DESCRIPTION")
  expect_match(env$check_code_unchanged("master", root), "DESCRIPTION: field Title")
  git("checkout", "-q", "--", "DESCRIPTION")

  writeLines(c(
    "---", "output: github_document", "---",
    "New prose `r 1 + 2` here.",
    "```{r}", "x <- 1", "```"
  ), "README.Rmd")
  expect_match(env$check_code_unchanged("master", root), "README.Rmd: chunk or inline-code item 4")
  git("checkout", "-q", "--", "README.Rmd")

  writeLines(c(
    "---", "output: github_document", "---",
    "New prose `r 1 + 1` here.",
    "```{r}", "x <- 2", "```"
  ), "README.Rmd")
  expect_match(env$check_code_unchanged("master", root), "README.Rmd: chunk or inline-code item 2")
  git("checkout", "-q", "--", "README.Rmd")

  # A vignette source is guarded the same way: a chunk option, a chunk body
  # line, and an inline span each report by file and item.
  v2 <- vig
  v2[5] <- "```{r chunk, eval = TRUE}"
  writeLines(v2, "vignettes/v.Rmd.orig")
  expect_match(
    env$check_code_unchanged("master", root),
    "vignettes/v.Rmd.orig: chunk or inline-code item 1"
  )
  v2 <- vig
  v2[6] <- "y <- 2"
  writeLines(v2, "vignettes/v.Rmd.orig")
  expect_match(
    env$check_code_unchanged("master", root),
    "vignettes/v.Rmd.orig: chunk or inline-code item 2"
  )
  v2 <- vig
  v2[4] <- "New vignette prose `r 2 + 3` here."
  writeLines(v2, "vignettes/v.Rmd.orig")
  expect_match(
    env$check_code_unchanged("master", root),
    "vignettes/v.Rmd.orig: chunk or inline-code item 4"
  )
})
