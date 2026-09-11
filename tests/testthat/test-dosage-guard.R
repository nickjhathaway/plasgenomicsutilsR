# A dosage matrix and an allele-index matrix are both integer matrices, and telling them
# apart by looking is only possible when an index happens to exceed 2. So the object carries
# its own `encoding`, and these guards read it -- falling back to the value check for a bare
# matrix, which catches the unambiguous half.
#
# Without them an index panel is not rejected, it is *mangled*: `.haploid_calls()` maps 0->0
# and 2->1 and leaves everything else NA, so allele 2 becomes allele 1, allele 1 becomes
# missing, and allele 3 becomes missing. The result is a plausible-looking matrix of the
# right shape and the wrong contents.

.idx_mat <- function() {
  m <- matrix(c(0L, 1L, 2L, 3L, 0L, 1L), nrow = 2,
              dimnames = list(c("s1", "s2"), c("c1:1", "c1:2", "c1:3")))
  m
}
.dos_mat <- function() {
  matrix(c(0L, 2L, 2L, 0L, 1L, 2L), nrow = 2,
         dimnames = list(c("s1", "s2"), c("c1:1", "c1:2", "c1:3")))
}

test_that("a matrix carrying an allele index above 2 is refused, not mangled", {
  expect_error(plasgenomicsutilsR:::.require_dosage(.idx_mat(), "pop_diversity()"),
               "allele index")
  expect_error(plasgenomicsutilsR:::.require_dosage(.idx_mat(), "pop_diversity()"),
               "encoding")
})

test_that("an ordinary dosage matrix passes", {
  expect_silent(plasgenomicsutilsR:::.require_dosage(.dos_mat(), "pop_diversity()"))
  expect_silent(plasgenomicsutilsR:::.require_dosage(
    matrix(c(0L, NA_integer_, 2L, 1L), 2), "x()"))
})

test_that("a declared non-dosage encoding is refused even when its values look like dosage", {
  # the case the value check cannot catch: a triallelic marker coded 0/1/2 is a valid
  # allele-index column and an equally valid dosage column
  g <- list(genotype = .dos_mat(), encoding = "allele_index")
  expect_error(plasgenomicsutilsR:::.require_dosage(g, "pop_diff()"), "allele_index")
  g$encoding <- "dosage"
  expect_silent(plasgenomicsutilsR:::.require_dosage(g, "pop_diff()"))
})

test_that("the error names the caller, so it says which analysis refused", {
  expect_error(plasgenomicsutilsR:::.require_dosage(.idx_mat(), "run_snmf()"), "run_snmf")
})

test_that("pop_diversity refuses an allele-index matrix rather than silently voiding it", {
  # `.haploid_calls()` would turn allele 1 into NA and allele 2 into allele 1
  expect_error(pop_diversity(.idx_mat()), "allele index")
})

test_that("pop_structure refuses one too, since imputation would invent an allele", {
  # the mean of allele indices {0, 2, 3} is 1.67, which is not an allele
  expect_error(pop_structure(.idx_mat()), "allele index")
})

test_that("run_snmf refuses one, since .geno has no symbol for a 3", {
  skip_if_not_installed("LEA")
  expect_error(run_snmf(.idx_mat(), K = 2), "allele index")
})

test_that("a biallelic panel still reaches all three unchanged", {
  d <- .dos_mat()
  expect_no_error(plasgenomicsutilsR:::.haploid_calls(d))
  expect_no_error(plasgenomicsutilsR:::.impute_geno(d))
})
