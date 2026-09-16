# Plain-English prose checker (M85).
#
# Sweeps the user-facing prose of the package and reports, with file and line,
# every em dash, en dash, ` -- `, semicolon, banned phrase
# (tools/prose-banned.txt), and sentence over `max_words` words. It never
# touches code: fenced chunks, roxygen `@examples` blocks, roxygen fenced code,
# YAML headers, and inline code spans are removed before any pattern runs.
# It errors, rather than going quiet, on a YAML header or fenced block that is
# never closed and on a backtick left unmatched within its paragraph.
#
# `check_code_unchanged(ref)` is the companion guard: it proves a prose
# rewrite left the code alone by comparing, against the merge base with `ref`,
# the non-roxygen lines of R/*.R, the `#'` lines inside `@examples`, the
# fenced chunks and inline `r` spans of README.Rmd, and every DESCRIPTION
# field other than `Description:`.
#
# Base R only, so it runs before any dependency install.
#
# USAGE (from the package root):
#   Rscript tools/check-prose.R [paths]        # prose sweep; no path = default domain
#   Rscript tools/check-prose.R --code-unchanged [ref]   # code guard (ref default: master)
#
# A directory argument expands to its *.R, *.Rmd, and *.Rmd.orig files. The
# default domain is README.Rmd, vignettes/*.Rmd.orig, vignettes/ackwards-interpret.Rmd,
# DESCRIPTION, NEWS.md, and R/*.R.
#
# It is also sourced by tests/testthat/test-check-prose.R, which reuses
# `check_prose()` directly; sourcing defines the functions without running the
# script body (guarded by `sys.nframe()`).

# ---- lists -------------------------------------------------------------------

.prose_tools_dir <- function() {
  # When sourced, sys.frame carries the file; when run, commandArgs does.
  f <- tryCatch(sys.frame(1L)$ofile, error = function(e) NULL)
  if (is.null(f)) {
    arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(arg) > 0L) f <- sub("^--file=", "", arg[[1L]])
  }
  if (is.null(f) || !nzchar(f)) {
    return("tools")
  }
  dirname(normalizePath(f, mustWork = FALSE))
}

read_prose_list <- function(name, dir = .prose_tools_dir()) {
  path <- file.path(dir, name)
  if (!file.exists(path)) {
    stop("prose list not found: ", path, call. = FALSE)
  }
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  lines <- trimws(lines)
  lines[nzchar(lines) & !startsWith(lines, "#")]
}

# ---- extraction: which lines are prose ---------------------------------------

# Each extractor returns a list of blocks; a block is a data.frame(line, text)
# of consecutive prose lines (paragraph breaks are kept as empty text so the
# sentence splitter sees them).

.is_fence <- function(x) grepl("^\\s*(```|~~~)", x)

.extract_markdown <- function(lines, skip_yaml = TRUE, file = "<text>") {
  keep <- rep(TRUE, length(lines))
  i <- 1L
  # YAML header: a `---` on line 1 through the next `---`. An unterminated
  # header is an error, never a silently empty file.
  if (skip_yaml && length(lines) > 0L && grepl("^---\\s*$", lines[[1L]])) {
    end <- which(grepl("^(---|\\.\\.\\.)\\s*$", lines))
    end <- end[end > 1L]
    if (length(end) == 0L) {
      stop(file, ":1: YAML header opened by `---` is never closed.", call. = FALSE)
    }
    stop_at <- end[[1L]]
    keep[1L:stop_at] <- FALSE
    i <- stop_at + 1L
  }
  in_fence <- FALSE
  fence_at <- NA_integer_
  while (i <= length(lines)) {
    if (.is_fence(lines[[i]])) {
      in_fence <- !in_fence
      fence_at <- i
      keep[[i]] <- FALSE
    } else if (in_fence) {
      keep[[i]] <- FALSE
    }
    i <- i + 1L
  }
  if (in_fence) {
    stop(file, ":", fence_at, ": fenced code block is never closed.", call. = FALSE)
  }
  # A markdown horizontal rule and the precompute stamp comment carry no prose.
  keep[grepl("^\\s*[-*_]{3,}\\s*$", lines)] <- FALSE
  keep[grepl("precompute-stamp:", lines, fixed = TRUE)] <- FALSE
  data.frame(line = seq_along(lines)[keep], text = lines[keep], stringsAsFactors = FALSE)
}

.extract_news <- function(lines, file = "<text>") {
  heads <- which(grepl("^# ", lines))
  if (length(heads) == 0L) {
    return(data.frame(line = integer(0L), text = character(0L)))
  }
  from <- heads[[1L]]
  to <- if (length(heads) >= 2L) heads[[2L]] - 1L else length(lines)
  md <- .extract_markdown(lines[from:to], skip_yaml = FALSE, file = file)
  md$line <- md$line + from - 1L
  md
}

.extract_description <- function(lines) {
  start <- which(grepl("^Description:", lines))
  if (length(start) == 0L) {
    return(data.frame(line = integer(0L), text = character(0L)))
  }
  start <- start[[1L]]
  end <- start
  while (end < length(lines) && grepl("^\\s", lines[[end + 1L]])) end <- end + 1L
  text <- lines[start:end]
  text[[1L]] <- sub("^Description:\\s*", "", text[[1L]])
  data.frame(line = start:end, text = text, stringsAsFactors = FALSE)
}

.extract_roxygen <- function(lines, file = "<text>") {
  is_rox <- grepl("^#'", lines)
  text <- sub("^#' ?", "", lines)
  keep <- is_rox
  in_examples <- FALSE
  in_fence <- FALSE
  fence_at <- NA_integer_
  for (i in c(seq_along(lines), length(lines) + 1L)) {
    if (i > length(lines) || !is_rox[[i]]) {
      # Leaving a roxygen block resets both states. An open fence at that
      # point is an error, never a silently dropped block.
      if (in_fence) {
        stop(file, ":", fence_at, ": roxygen fenced code is never closed.", call. = FALSE)
      }
      in_examples <- FALSE
      in_fence <- FALSE
      next
    }
    t <- text[[i]]
    if (grepl("^@examples", t)) {
      in_examples <- TRUE
      keep[[i]] <- FALSE
      next
    }
    if (in_examples) {
      if (grepl("^@[A-Za-z]", t)) {
        in_examples <- FALSE
      } else {
        keep[[i]] <- FALSE
        next
      }
    }
    if (.is_fence(t)) {
      in_fence <- !in_fence
      fence_at <- i
      keep[[i]] <- FALSE
      next
    }
    # Rd's own code block, `\preformatted{ ... }`, closed by a lone `}`.
    if (grepl("\\\\preformatted\\{", t)) {
      in_fence <- TRUE
      fence_at <- i
      keep[[i]] <- FALSE
      next
    }
    if (in_fence && grepl("^\\s*\\}\\s*$", t)) {
      in_fence <- FALSE
      keep[[i]] <- FALSE
      next
    }
    if (in_fence) keep[[i]] <- FALSE
  }
  data.frame(line = seq_along(lines)[keep], text = text[keep], stringsAsFactors = FALSE)
}

.prose_kind <- function(path) {
  base <- basename(path)
  if (identical(base, "DESCRIPTION")) {
    return("description")
  }
  if (identical(base, "NEWS.md")) {
    return("news")
  }
  if (grepl("\\.(Rmd|Rmd\\.orig|md)$", base, ignore.case = TRUE)) {
    return("markdown")
  }
  if (grepl("\\.[Rr]$", base)) {
    return("roxygen")
  }
  "markdown"
}

extract_prose <- function(path, kind = .prose_kind(path)) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  switch(kind,
    description = .extract_description(lines),
    news = .extract_news(lines, file = path),
    roxygen = .extract_roxygen(lines, file = path),
    .extract_markdown(lines, file = path)
  )
}

# ---- code-span removal -------------------------------------------------------

# Paragraph ids for a vector of prose lines: consecutive non-empty lines share
# an id, and an empty line ends the paragraph. Spans and wrapped phrases are
# matched within a paragraph, never across one.
.paragraph_ids <- function(text) {
  blank <- !nzchar(trimws(text))
  cumsum(c(TRUE, blank[-length(blank)] & !blank[-1L])) * !blank
}

.strip_one_paragraph <- function(text) {
  joined <- paste(text, collapse = "\n")
  # Double backticks first so a single backtick inside them is not a span start.
  pat <- "``[\\s\\S]*?``|`[^`\\n][\\s\\S]*?`|``"
  m <- gregexpr(pat, joined, perl = TRUE)[[1L]]
  if (m[[1L]] != -1L) {
    starts <- as.integer(m)
    lens <- attr(m, "match.length")
    out <- character(0L)
    pos <- 1L
    for (j in seq_along(starts)) {
      out <- c(out, substr(joined, pos, starts[[j]] - 1L))
      span <- substr(joined, starts[[j]], starts[[j]] + lens[[j]] - 1L)
      nl <- gsub("[^\n]", "", span)
      out <- c(out, paste0(" ", nl))
      pos <- starts[[j]] + lens[[j]]
    }
    out <- c(out, substr(joined, pos, nchar(joined)))
    joined <- paste(out, collapse = "")
  }
  strsplit(joined, "\n", fixed = TRUE)[[1L]] |>
    (\(x) if (length(x) < length(text)) c(x, rep("", length(text) - length(x))) else x)()
}

# Remove single- and double-backtick spans (including spans that cross a line
# break, never a paragraph break) from a vector of consecutive prose lines,
# preserving line count. A span's content is replaced by a single space so
# neighbouring words never merge; the newlines inside a span are kept so line
# numbers stay aligned. A backtick left unmatched within its paragraph is an
# error: a stray backtick would otherwise swallow the prose up to the next one.
strip_code_spans <- function(text, file = "<text>", line = seq_along(text)) {
  ids <- .paragraph_ids(text)
  out <- text
  for (id in setdiff(unique(ids), 0L)) {
    at <- which(ids == id)
    stripped <- .strip_one_paragraph(text[at])
    left <- grepl("`", stripped, fixed = TRUE)
    if (any(left)) {
      stop(
        file, ":", line[at][left][[1L]], ": unmatched backtick in this paragraph.",
        call. = FALSE
      )
    }
    out[at] <- stripped
  }
  out
}

# ---- reports -----------------------------------------------------------------

.report <- function(file, line, class, text) {
  data.frame(
    file = file, line = as.integer(line), class = class, text = text,
    stringsAsFactors = FALSE
  )
}

.empty_report <- function() .report(character(0L), integer(0L), character(0L), character(0L))

.char_reports <- function(file, prose, banned) {
  out <- list()
  flag <- function(pat, class, fixed = FALSE, ignore.case = FALSE, perl = FALSE) {
    hit <- grepl(pat, prose$text, fixed = fixed, ignore.case = ignore.case, perl = perl)
    if (any(hit)) {
      out[[length(out) + 1L]] <<- .report(file, prose$line[hit], class, trimws(prose$text[hit]))
    }
  }
  flag("—", "em dash", fixed = TRUE)
  flag("–", "en dash", fixed = TRUE)
  flag(" -- ", "double hyphen", fixed = TRUE)
  flag(";", "semicolon", fixed = TRUE)
  # Banned phrases are matched over each paragraph joined with newlines, so a
  # multi-word phrase wrapped across a line break is still found. The report
  # names the line where the match starts.
  ids <- .paragraph_ids(prose$text)
  for (b in banned) {
    esc <- gsub("([][{}()+*^$|\\\\.?])", "\\\\\\1", b)
    esc <- gsub("\\s+", "\\\\s+", esc)
    pat <- paste0("\\b", esc, "\\b")
    for (id in setdiff(unique(ids), 0L)) {
      at <- which(ids == id)
      joined <- paste(prose$text[at], collapse = "\n")
      m <- gregexpr(pat, joined, ignore.case = TRUE, perl = TRUE)[[1L]]
      if (m[[1L]] == -1L) next
      # Line of each match: count the newlines before its start.
      nl_before <- vapply(
        as.integer(m),
        function(s) nchar(gsub("[^\n]", "", substr(joined, 1L, s - 1L))),
        integer(1L)
      )
      rows <- unique(at[nl_before + 1L])
      out[[length(out) + 1L]] <- .report(
        file, prose$line[rows], paste0("banned phrase: ", b), trimws(prose$text[rows])
      )
    }
  }
  if (length(out) == 0L) {
    return(.empty_report())
  }
  do.call(rbind, out)
}

# Sentences. A sentence is the text between terminators: `.`, `?`, or `!`
# (optionally followed by a closing quote, bracket, or emphasis mark) followed
# by whitespace and an uppercase letter, or standing at a line end. Headings,
# bullet markers, table rows, roxygen tags, and the listed abbreviations are
# removed first. A paragraph break (empty prose line, new bullet, new roxygen
# tag, or block end) also ends a sentence.
.sentence_reports <- function(file, prose, abbrev, max_words) {
  text <- prose$text
  line <- prose$line
  # Headings and table rows carry no sentence.
  drop <- grepl("^\\s*#{1,6}\\s", text) | grepl("^\\s*\\|", text) | grepl("^@section\\b", text)
  starts <- grepl("^\\s*([-*+]|\\d+[.)])\\s", text) | grepl("^@[A-Za-z]", text)
  # Roxygen `@section Title:` is a heading; other tags are markers whose tag
  # word (and, for @param, the argument name) is not prose.
  text <- sub("^@param\\s+\\S+\\s*", "", text)
  text <- sub("^@[A-Za-z]+\\s*", "", text)
  text <- sub("^\\s*([-*+]|\\d+[.)])\\s+", "", text)
  text <- gsub("\\\\item\\{[^}]*\\}", "", text)
  text[drop] <- ""
  # Abbreviations: remove so their period never terminates a sentence.
  for (a in abbrev) {
    esc <- gsub("([][{}()+*^$|\\\\.?])", "\\\\\\1", a)
    text <- gsub(paste0("(?<![A-Za-z])", esc, "(?![A-Za-z])"), "", text, perl = TRUE)
  }
  term <- "[.?!][)\"'*_\\]]*"
  # PCRE lookbehind must be fixed-width, so mark the split point with a
  # sentinel and split on that instead.
  mark_re <- paste0("(", term, ")\\s+(?=[(\\[\"'*_]*[A-Z])")
  sentinel <- ""

  reports <- list()
  cur_words <- 0L
  cur_line <- NA_integer_
  cur_text <- character(0L)
  flush <- function() {
    if (cur_words > max_words) {
      reports[[length(reports) + 1L]] <<- .report(
        file, cur_line, sprintf("sentence over %d words (%d)", max_words, cur_words),
        paste(cur_text, collapse = " ")
      )
    }
    cur_words <<- 0L
    cur_line <<- NA_integer_
    cur_text <<- character(0L)
  }
  for (i in seq_along(text)) {
    t <- trimws(text[[i]])
    if (!nzchar(t) || starts[[i]]) flush()
    if (!nzchar(t)) next
    marked <- gsub(mark_re, paste0("\\1", sentinel), t, perl = TRUE)
    pieces <- strsplit(marked, sentinel, fixed = TRUE)[[1L]]
    for (j in seq_along(pieces)) {
      p <- trimws(pieces[[j]])
      if (!nzchar(p)) next
      w <- length(strsplit(p, "\\s+")[[1L]])
      if (is.na(cur_line)) cur_line <- line[[i]]
      cur_words <- cur_words + w
      cur_text <- c(cur_text, p)
      ends <- j < length(pieces) || grepl(paste0(term, "$"), p, perl = TRUE)
      if (ends) flush()
    }
  }
  flush()
  if (length(reports) == 0L) {
    return(.empty_report())
  }
  do.call(rbind, reports)
}

# ---- public: check_prose -----------------------------------------------------

default_prose_domain <- function(root = ".") {
  c(
    file.path(root, "README.Rmd"),
    list.files(file.path(root, "vignettes"), pattern = "\\.Rmd\\.orig$", full.names = TRUE),
    file.path(root, "vignettes", "ackwards-interpret.Rmd"),
    file.path(root, "DESCRIPTION"),
    file.path(root, "NEWS.md"),
    list.files(file.path(root, "R"), pattern = "\\.[Rr]$", full.names = TRUE)
  )
}

expand_prose_paths <- function(paths) {
  out <- character(0L)
  for (p in paths) {
    if (dir.exists(p)) {
      p <- sub("/+$", "", p)
      out <- c(out, list.files(p, pattern = "\\.([Rr]|Rmd|Rmd\\.orig)$", full.names = TRUE))
    } else {
      out <- c(out, p)
    }
  }
  unique(out)
}

# Returns a data.frame(file, line, class, text); zero rows means clean.
check_prose <- function(paths = NULL, banned = read_prose_list("prose-banned.txt"),
                        abbrev = read_prose_list("prose-abbrev.txt"), max_words = 30L) {
  if (is.null(paths)) paths <- default_prose_domain()
  paths <- expand_prose_paths(paths)
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("file not found: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  out <- list()
  for (p in paths) {
    prose <- extract_prose(p)
    if (nrow(prose) == 0L) next
    prose$text <- strip_code_spans(prose$text, file = p, line = prose$line)
    out[[length(out) + 1L]] <- .char_reports(p, prose, banned)
    out[[length(out) + 1L]] <- .sentence_reports(p, prose, abbrev, max_words)
  }
  if (length(out) == 0L) {
    return(.empty_report())
  }
  res <- do.call(rbind, out)
  res <- res[order(res$file, res$line), , drop = FALSE]
  rownames(res) <- NULL
  res
}

# ---- public: check_code_unchanged --------------------------------------------

.git <- function(...) {
  out <- suppressWarnings(system2("git", c(...), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0L) {
    stop("git ", paste(c(...), collapse = " "), " failed: ", paste(out, collapse = "\n"),
      call. = FALSE
    )
  }
  out
}

.git_show <- function(ref, path) {
  out <- suppressWarnings(system2("git", c("show", paste0(ref, ":", path)), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0L) {
    return(NULL)
  }
  out
}

# The code lines of an R file: every non-roxygen line, plus the roxygen lines
# inside `@examples` blocks.
.code_lines_r <- function(lines) {
  is_rox <- grepl("^#'", lines)
  text <- sub("^#' ?", "", lines)
  keep <- !is_rox
  in_examples <- FALSE
  for (i in seq_along(lines)) {
    if (!is_rox[[i]]) {
      in_examples <- FALSE
      next
    }
    if (grepl("^@examples", text[[i]])) {
      in_examples <- TRUE
      keep[[i]] <- TRUE
      next
    }
    if (in_examples && grepl("^@[A-Za-z]", text[[i]])) in_examples <- FALSE
    if (in_examples) keep[[i]] <- TRUE
  }
  lines[keep]
}

# The code of README.Rmd: fenced chunk lines (fences included) and inline `r`
# spans, in order.
.code_lines_rmd <- function(lines) {
  in_fence <- FALSE
  keep <- rep(FALSE, length(lines))
  for (i in seq_along(lines)) {
    if (.is_fence(lines[[i]])) {
      in_fence <- !in_fence
      keep[[i]] <- TRUE
    } else if (in_fence) {
      keep[[i]] <- TRUE
    }
  }
  inline <- regmatches(lines[!keep], gregexpr("`r [^`]*`", lines[!keep]))
  c(lines[keep], unlist(inline))
}

.fields_description <- function(lines) {
  tmp <- tempfile()
  on.exit(unlink(tmp), add = TRUE)
  writeLines(lines, tmp)
  d <- read.dcf(tmp)
  d <- d[1L, , drop = TRUE]
  d[setdiff(names(d), "Description")]
}

.first_diff <- function(a, b) {
  n <- max(length(a), length(b))
  length(a) <- n
  length(b) <- n
  same <- !is.na(a) & !is.na(b) & a == b
  which(!same)[1L]
}

# Returns a character vector of problems (empty = code unchanged since the
# merge base of `ref` and HEAD).
check_code_unchanged <- function(ref = "master", root = ".") {
  old <- setwd(root)
  on.exit(setwd(old), add = TRUE)
  base <- .git("merge-base", ref, "HEAD")[[1L]]
  problems <- character(0L)

  r_files <- sort(unique(c(
    list.files("R", pattern = "\\.[Rr]$", full.names = TRUE),
    grep("^R/.*\\.[Rr]$", .git("ls-tree", "--name-only", base, "R/"), value = TRUE)
  )))
  for (f in r_files) {
    before <- .git_show(base, f)
    after <- if (file.exists(f)) readLines(f, warn = FALSE) else NULL
    if (is.null(before) || is.null(after)) {
      problems <- c(problems, sprintf("%s exists on only one side of %s.", f, substr(base, 1, 7)))
      next
    }
    a <- .code_lines_r(before)
    b <- .code_lines_r(after)
    d <- .first_diff(a, b)
    if (!is.na(d)) {
      problems <- c(problems, sprintf(
        "%s: code line %d differs from the merge base (%s).", f, d,
        if (is.na(b[d])) "line removed" else b[d]
      ))
    }
  }

  before <- .git_show(base, "README.Rmd")
  if (!is.null(before) && file.exists("README.Rmd")) {
    a <- .code_lines_rmd(before)
    b <- .code_lines_rmd(readLines("README.Rmd", warn = FALSE))
    d <- .first_diff(a, b)
    if (!is.na(d)) {
      problems <- c(problems, sprintf(
        "README.Rmd: chunk or inline-code item %d differs from the merge base (%s).", d,
        if (is.na(b[d])) "item removed" else b[d]
      ))
    }
  }

  before <- .git_show(base, "DESCRIPTION")
  if (!is.null(before) && file.exists("DESCRIPTION")) {
    a <- .fields_description(before)
    b <- .fields_description(readLines("DESCRIPTION", warn = FALSE))
    for (nm in union(names(a), names(b))) {
      if (!identical(unname(a[nm]), unname(b[nm]))) {
        problems <- c(problems, sprintf("DESCRIPTION: field %s differs from the merge base.", nm))
      }
    }
  }

  problems
}

# ---- script body -------------------------------------------------------------

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 0L && identical(args[[1L]], "--code-unchanged")) {
    ref <- if (length(args) >= 2L) args[[2L]] else "master"
    problems <- check_code_unchanged(ref)
    if (length(problems) > 0L) {
      message("Code-unchanged check FAILED:")
      for (p in problems) message("  - ", p)
      quit(status = 1L, save = "no")
    }
    message("Code unchanged OK: no code line differs from the merge base with ", ref, ".")
  } else {
    res <- check_prose(if (length(args) > 0L) args else NULL)
    if (nrow(res) > 0L) {
      message("Prose check FAILED (", nrow(res), " report(s)):")
      for (i in seq_len(nrow(res))) {
        message(sprintf("  %s:%d  [%s]  %s", res$file[i], res$line[i], res$class[i], res$text[i]))
      }
      quit(status = 1L, save = "no")
    }
    message("Prose OK: no dash, semicolon, banned phrase, or over-long sentence found.")
  }
}
