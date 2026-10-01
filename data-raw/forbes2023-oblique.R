# Generate the oblique-rotation fidelity fixture
#   tests/testthat/fixtures/forbes2023_oblique.rds
#
# Forbes's (2023) reference implementation is "equipped to do orthogonal or
# oblique PCA or EFA" (her fn. 1): ExtendedBassAckwards() passes its `rotate`
# argument straight to psych::pca() / psych::fa(). This script runs that
# oblique branch on the three simulation matrices of the paper and records the
# values ackwards must reproduce under `rotation = "oblimin"` / `"promax"`.
#
# INPUTS
# -------------------------------------------------------------------------
# The three Spearman correlation matrices in
# tests/testthat/fixtures/forbes2023_sims.rds, the deterministic realizations
# of Forbes's simulation recipe (generator: data-raw/oracle-forbes-sims.R).
# Each is run with her functions file (OSF guid 7jfkw, md5-pinned below) at
# num.comp = 4 for:
#   fm     = "pca" (psych::pca) and "minres" (psych::fa, her default
#            scores = "tenBerge"), matching ackwards engine = "pca" and "efa";
#   rotate = "oblimin" and "promax".
#
# WHAT IS RECORDED (per sim x fm x rotate)
# -------------------------------------------------------------------------
#   comp_corr  her comp.corr list, t(W_a) %*% R %*% W_b with her weights and no
#              sign alignment, pairs enumerated as for c in 2..4, i in 1..c-1;
#   D          per level, diag(t(W) %*% R %*% W) for her weights. ackwards
#              standardizes every edge by the real score SDs, so the test
#              compares against D_a^{-1/2} comp_corr D_b^{-1/2}. Under oblique
#              rotation D = I holds only for correlation-preserving weights,
#              so it is stored rather than assumed (RR02 Q5 item 13);
#   loadings   per level, her psych pattern loadings (the test reads the sign
#              flips ackwards applied from them);
#   Phi        per level, psych's factor correlation ($Phi);
#   corr_chase her FindRedundantComp(comp.corr, cong, "d4")$corr.chase, which
#              reads the unstandardized comp.corr;
#   corr_chase_std  the same function on the D-standardized comp.corr, that
#              is, her chase on the correlations ackwards reports. The two
#              differ where D != I: psych::fa's one-factor weights are not
#              unit-variance, so her raw level-1 EFA products sit below the
#              correlations, and a b-to-a link can fall under .9 raw only;
#   chase_unbroken  components at level 3 or deeper whose chase on the
#              D-standardized comp.corr is unbroken to level a under her rule
#              (signed column max >= .9 at every ancestor level). For these
#              her ChaseCorrPaths() returns "X--null": with no FALSE in its
#              vector, which.min() returns 1, so it counts zero consecutive
#              links. (Level-2 components are special-cased to "a1" in her
#              code.) The test expects a1 for them.
#
# Re-run after any change to the inputs or the recorded quantities:
#   Rscript data-raw/forbes2023-oblique.R
# Requires network access to osf.io, psych, and GPArotation (psych's oblimin
# and promax both load it).

# Forbes's functions call fa.sort()/pca() unqualified, so psych must be attached.
suppressPackageStartupMessages(library(psych))
stopifnot(requireNamespace("GPArotation", quietly = TRUE))

osf_functions <- list(
  guid = "7jfkw", name = "ExtendedBassAckwards functions with annotation.R",
  url = "https://osf.io/download/7jfkw/", md5 = "a3e85df897d2a4a4310b9c45dcb068d9"
)

fun_path <- tempfile(fileext = "_ExtendedBassAckwards.R")
utils::download.file(osf_functions$url, fun_path, mode = "wb", quiet = TRUE)
stopifnot(unname(tools::md5sum(fun_path)) == osf_functions$md5)
source(fun_path, local = TRUE) # defines ExtendedBassAckwards(), FindRedundantComp()

sims_path <- file.path("tests", "testthat", "fixtures", "forbes2023_sims.rds")
sims <- readRDS(sims_path)

K <- 4L
one_run <- function(R, fm, rotate) {
  # psych::fa() draws 20 random starts for an oblimin rotation (n.rotations),
  # so the minres + oblimin values depend on the RNG state. A fixed seed per
  # run makes the fixture bit-reproducible; PCA and promax use no RNG.
  set.seed(2023)
  fb <- suppressWarnings(suppressMessages(
    ExtendedBassAckwards(R, num.comp = K, fm = fm, rotate = rotate)
  ))
  chase <- unlist(FindRedundantComp(fb$comp.corr, fb$cong, "d4")$corr.chase)
  D <- lapply(seq_len(K), function(c) {
    W <- unclass(fb$pcas[[c]]$weights)
    diag(crossprod(W, R %*% W))
  })
  loadings <- lapply(seq_len(K), function(c) unclass(fb$pcas[[c]]$loadings))
  Phi <- lapply(seq_len(K), function(c) {
    P <- fb$pcas[[c]]$Phi
    if (is.null(P)) diag(c) else unname(P)
  })
  # Her chase on D-standardized comp.corr (dimnames kept: her chase reads
  # them). Pairs run as in her list: for c in 2..K, for i in 1..c-1.
  pairs <- do.call(rbind, lapply(2:K, function(c) cbind(i = seq_len(c - 1L), c = c)))
  comp_corr_std <- lapply(seq_len(nrow(pairs)), function(j) {
    fb$comp.corr[[j]] /
      sqrt(outer(unname(D[[pairs[j, "i"]]]), unname(D[[pairs[j, "c"]]])))
  })
  chase_std <- unlist(FindRedundantComp(comp_corr_std, fb$cong, "d4")$corr.chase)
  # Level >= 3 components whose every ancestor level passes her signed rule.
  deep <- unlist(lapply(3:K, function(c) paste0(letters[c], seq_len(c))))
  unbroken <- deep[vapply(deep, function(comp) {
    cols <- Filter(function(m) comp %in% colnames(m), comp_corr_std)
    all(vapply(cols, function(m) max(m[, comp]) >= 0.9, logical(1)))
  }, logical(1))]
  stopifnot(
    length(fb$comp.corr) == choose(K, 2),
    length(chase) == sum(2:K),
    length(chase_std) == sum(2:K),
    all(vapply(Phi[2:K], function(P) max(abs(P[upper.tri(P)])), numeric(1)) > 0.05)
  )
  list(
    comp_corr = lapply(fb$comp.corr, unname),
    D = D, loadings = lapply(loadings, unname), Phi = Phi,
    corr_chase = chase, corr_chase_std = chase_std,
    chase_unbroken = unbroken
  )
}

runs <- list()
for (nm in names(sims)) {
  for (fm in c("pca", "minres")) {
    for (rotate in c("oblimin", "promax")) {
      runs[[paste(nm, fm, rotate, sep = "_")]] <- c(
        list(sim = nm, fm = fm, rotate = rotate),
        one_run(sims[[nm]]$R, fm, rotate)
      )
    }
  }
}

attr(runs, "provenance") <- list(
  source = paste0(
    "Forbes, M. K. (2023). Improving hierarchical models of individual ",
    "differences: An extension of Goldberg's bass-ackward method. ",
    "Psychological Methods. doi:10.1037/met0000546 (fn. 1: the reference ",
    "implementation does orthogonal or oblique PCA or EFA)"
  ),
  osf = "https://osf.io/pcwm8/",
  license = "CC-BY 4.0 International (https://creativecommons.org/licenses/by/4.0/)",
  files = sprintf(
    "%s (guid %s, md5 %s) <%s>",
    osf_functions$name, osf_functions$guid, osf_functions$md5, osf_functions$url
  ),
  inputs = paste0(
    sims_path, " (the three simulation Spearman matrices; generator ",
    "data-raw/oracle-forbes-sims.R)"
  ),
  recipe = paste0(
    "Per simulation, fm in {pca, minres}, rotate in {oblimin, promax}: Forbes's ",
    "set.seed(2023); ExtendedBassAckwards(R, num.comp = 4, fm, rotate) with her default ",
    "scores = 'tenBerge'; comp_corr, per-level D = diag(W'RW), loadings, Phi, ",
    "and FindRedundantComp(..., 'd4')$corr.chase on the raw and on the ",
    "D-standardized comp.corr recorded."
  ),
  generator = "data-raw/forbes2023-oblique.R",
  generated = as.character(Sys.Date()),
  R_version = R.version.string,
  psych = as.character(utils::packageVersion("psych")),
  GPArotation = as.character(utils::packageVersion("GPArotation"))
)

out <- file.path("tests", "testthat", "fixtures", "forbes2023_oblique.rds")
saveRDS(runs, out, compress = "xz")
message(sprintf("Wrote %s (%.0f KB)", out, file.size(out) / 1024))
