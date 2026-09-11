# `scan_hh()` orders a two-allele marker's columns by FREQUENCY -- IHH_MAJ then IHH_MIN --
# and says nothing about which allele is which. The contrast rows are labelled by allele
# INDEX ("0>1"). Wherever allele 0 is the minor allele the two disagree, and the row is an
# exact sign flip of what its own label claims: log(iHH_1/iHH_0) reported as "0>1".
#
# This is not a rounding difference and it is not rare. On a 60-haplotype panel simulated at
# p(alt) = 0.62, 59 of 60 markers were inverted, |unihs difference| up to 0.66. iHS is read
# by sign -- a negative score means the *second* allele carries the longer haplotypes -- so
# an inverted row points selection at the wrong allele.
#
# The truth these check against is `calc_ehh()`, which returns one integral per allele in
# allele-index order and therefore cannot be ambiguous.

.orient_hap <- function(p_alt, n = 60, m = 61, seed = 7) {
  set.seed(seed)
  G <- matrix(stats::rbinom(n * m, 1, p_alt), n, m)
  rownames(G) <- sprintf("s%02d", seq_len(n))
  colnames(G) <- paste0("Pf3D7_01_v3:", seq(1000, by = 500, length.out = m))
  parasite_haplotypes(G, maf = 0.02, alleles = "index", impute = FALSE)
}

.truth_01 <- function(hap, rows, snp_id, polarized = FALSE) {
  a <- plasgenomicsutilsR:::.marker_allele_ihh(hap, rows, snp_id, polarized)
  if (is.null(a)) NA_real_ else log(a$ihh[["0"]] / a$ihh[["1"]])
}

.contrast_rows <- function(hap, polarized = FALSE) {
  rows <- list(all = seq_len(nrow(hap$hap)))
  plasgenomicsutilsR:::.run_ihs_contrasts(hap, rows, "ref", polarized, 1, 0.05, NULL,
                                          NA, NA, NULL, 1, FALSE)
}

test_that("a 0>1 row is log(iHH_0/iHH_1) even when allele 0 is the MINOR allele", {
  skip_if_not_installed("rehh")
  hap <- .orient_hap(0.62)                      # allele 1 is the major allele nearly everywhere
  rows <- seq_len(nrow(hap$hap))
  res <- .contrast_rows(hap)
  minor0 <- vapply(res$snp_id, function(s) mean(hap$hap[, s] == 0L, na.rm = TRUE) < 0.5,
                   logical(1))
  expect_gt(sum(minor0), 40)                    # the case actually exercised, not a no-op
  truth <- vapply(res$snp_id, function(s) .truth_01(hap, rows, s), numeric(1))
  expect_equal(res$unihs, unname(truth), tolerance = 1e-9)
})

test_that("the marker where allele 0 IS the major allele is unchanged", {
  skip_if_not_installed("rehh")
  hap <- .orient_hap(0.35)                      # allele 0 major nearly everywhere
  rows <- seq_len(nrow(hap$hap))
  res <- .contrast_rows(hap)
  truth <- vapply(res$snp_id, function(s) .truth_01(hap, rows, s), numeric(1))
  expect_equal(res$unihs, unname(truth), tolerance = 1e-9)
})

test_that("a 50/50 marker is scored from the per-allele integrals, not guessed", {
  skip_if_not_installed("rehh")
  # at an exact tie the scan's MAJ/MIN order is unknowable from its output, so the code has
  # to fall through to calc_ehh rather than pick one and hope
  hap <- .orient_hap(0.5, seed = 3)
  n <- nrow(hap$hap)
  tie <- "Pf3D7_01_v3:11000"
  hap$hap[, tie] <- as.integer(rep(0:1, each = n / 2))
  rows <- seq_len(n)
  expect_equal(mean(hap$hap[, tie] == 0L), 0.5)
  res <- .contrast_rows(hap)
  got <- res$unihs[res$snp_id == tie]
  expect_length(got, 1L)
  expect_equal(got, .truth_01(hap, rows, tie), tolerance = 1e-9)
})

test_that("a polarized scan agrees with the per-allele integrals too", {
  skip_if_not_installed("rehh")
  hap <- .orient_hap(0.62)
  rows <- seq_len(nrow(hap$hap))
  res <- suppressWarnings(.contrast_rows(hap, polarized = TRUE))
  truth <- vapply(res$snp_id, function(s) .truth_01(hap, rows, s, TRUE), numeric(1))
  expect_equal(res$unihs, unname(truth), tolerance = 1e-9)
})
