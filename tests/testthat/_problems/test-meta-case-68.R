# Extracted from test-meta-case.R:68

# prequel ----------------------------------------------------------------------
.caps_meta <- function(case = "Sample") {
  m <- data.frame(x = paste0("s", 1:8),
                  region = rep(c("north", "south"), each = 4),
                  stringsAsFactors = FALSE)
  names(m)[1] <- case
  m
}
.small_geno <- function() {
  set.seed(4)
  G <- matrix(stats::rbinom(8 * 30, 2, 0.4), 8, 30)
  dimnames(G) <- list(paste0("s", 1:8),
                      paste0("Pf3D7_07_v3:", seq(1000, by = 900, length.out = 30)))
  G
}

# test -------------------------------------------------------------------------
ns <- asNamespace("plasgenomicsutilsR")
takes_meta <- Filter(function(nm) {
    f <- get(nm, envir = ns)
    is.function(f) && "meta" %in% names(formals(f))
  }, getNamespaceExports(ns))
missing <- Filter(function(nm) {
    txt <- paste(deparse(body(get(nm, envir = ns))), collapse = " ")
    !grepl(".normalise_meta(meta)", txt, fixed = TRUE)
  }, takes_meta)
expect_equal(missing, character(0))
