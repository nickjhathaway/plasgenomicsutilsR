# One object, not two.
#
# A dosage cannot carry a multiallelic site and an allele index cannot be fed to PCA, so the
# obvious workaround is to keep two PopStructures and remember which analysis takes which.
# That is the thing this avoids: the object holds the richest panel it was given, and each
# analysis says what it *needs*. A dosage view of the biallelic sites is derivable from an
# allele-index panel, so it is derived on demand and cached.
#
# The reverse is not derivable -- which is the whole point of the phase -- so an object built
# from a dosage panel cannot serve an allele-index analysis, and says so.

.multi_ps <- function(umap = FALSE) example_pop_structure("multiallelic", umap = umap)

test_that("the multiallelic example object holds an allele-index panel", {
  skip_if_not_installed("SeqArray")
  ps <- .multi_ps()
  expect_equal(ps$encoding(), "allele_index")
  s <- ps$sites()
  expect_gte(sum(s$n_alt_real > 1), 20L)
  expect_equal(nrow(s), ncol(ps$genotype()))
})

test_that("a dosage analysis gets a biallelic dosage view without being asked to make one", {
  skip_if_not_installed("SeqArray")
  ps <- .multi_ps()
  d <- ps$genotype(needs = "dosage")
  expect_true(all(d %in% c(0L, 2L, NA_integer_)))
  # only the biallelic sites are in it, since a dosage cannot carry the others
  expect_equal(ncol(d), sum(ps$sites()$n_alt_real == 1L))
  expect_lt(ncol(d), ncol(ps$genotype()))
})

test_that("the derived view is cached as a panel rather than rebuilt each time", {
  skip_if_not_installed("SeqArray")
  ps <- .multi_ps()
  expect_false("biallelic_dosage" %in% ps$panels())   # nothing has asked for it yet
  d <- ps$genotype(needs = "dosage")
  expect_true("biallelic_dosage" %in% ps$panels())
  expect_identical(d, ps$genotype(needs = "dosage"))
  # and it is a panel like any other, so it can be asked for by name
  expect_identical(ps$genotype(panel = "biallelic_dosage"), ps$genotype(needs = "dosage"))
})

test_that("the object says once where the derived panel came from", {
  skip_if_not_installed("SeqArray")
  # construction derives the ONE-HOT view, because that is what the PCA now runs on: a
  # multiallelic site reaches the ordination as one indicator column per ALT instead of
  # being dropped for not fitting in a dosage
  msgs <- paste(testthat::capture_messages(.multi_ps()), collapse = " ")
  expect_match(msgs, "one-hot panel")
  # and the biallelic dosage view still announces itself when something asks for it
  ps <- suppressMessages(.multi_ps())
  d <- paste(testthat::capture_messages(ps$genotype(needs = "dosage")), collapse = " ")
  expect_match(d, "derived a biallelic dosage panel")
  expect_match(d, "multiallelic ones left out")
})

test_that("the dosage statistics run straight off a multiallelic object", {
  skip_if_not_installed("SeqArray")
  ps <- .multi_ps()
  # no filtering, no second object, no error
  expect_no_error(div <- pop_diversity(ps))
  expect_true(nrow(div) > 0)
  expect_no_error(pd <- pop_diff(ps, group = "country"))
  expect_true(!is.null(pd))
})

test_that("an allele-index analysis gets the index panel", {
  skip_if_not_installed("SeqArray")
  ps <- .multi_ps()
  g <- ps$genotype(needs = "allele_index")
  expect_true(any(g > 1L, na.rm = TRUE), info = "a third allele must be present")
  expect_equal(ncol(g), ncol(ps$genotype()))
})

test_that("asking a dosage-only object for allele indices says why it cannot", {
  skip_if_not_installed("SNPRelate")
  ps <- example_pop_structure(umap = FALSE)
  expect_equal(ps$encoding(), "dosage")
  expect_error(ps$genotype(needs = "allele_index"), "cannot be recovered")
})

test_that("the toggle back to biallelic is one argument, not a second object", {
  skip_if_not_installed("SeqArray")
  ps <- .multi_ps()
  auto <- ps$genotype(needs = "allele_index")
  bi <- ps$genotype(needs = "dosage")
  expect_gt(ncol(auto), ncol(bi))
  # and the biallelic columns are a subset, so the two describe the same samples and sites
  expect_true(all(colnames(bi) %in% colnames(auto)))
  expect_equal(rownames(bi), rownames(auto))
})

test_that("a dosage object is unaffected: needs = dosage is what it already had", {
  skip_if_not_installed("SNPRelate")
  ps <- example_pop_structure(umap = FALSE)
  expect_identical(ps$genotype(needs = "dosage"), ps$genotype())
})

test_that("parasite_haplotypes follows the object by default, and toggles back with one word", {
  skip_if_not_installed("SeqArray")
  ps <- .multi_ps()
  # the fixture is a 1-in-500 stride through a real callset, so most sites are patchy;
  # `max_snp_missing` is what decides the panel size here, not the encoding
  # `impute = FALSE` drops any site with a missing call, and the multiallelic ones are the
  # patchiest, so imputation is what lets them reach the haplotypes at all here
  args <- list(maf = 0.02, impute = TRUE, max_snp_missing = 0.6)
  auto <- do.call(parasite_haplotypes, c(list(ps), args))
  expect_equal(auto$alleles, "index")
  expect_true(any(auto$hap > 1L, na.rm = TRUE), info = "a third allele reaches the haplotypes")

  bi <- do.call(parasite_haplotypes, c(list(ps), args, list(alleles = "dosage")))
  expect_equal(bi$alleles, "dosage")
  expect_true(all(bi$hap %in% c(0, 1, NA)))
  expect_lt(ncol(bi$hap), ncol(auto$hap))
})

test_that("a dosage object keeps reading as dosage under auto", {
  skip_if_not_installed("SNPRelate")
  ps <- example_pop_structure(umap = FALSE)
  h <- parasite_haplotypes(ps, maf = 0.02, impute = FALSE)
  expect_equal(h$alleles, "dosage")
})

# --- one-hot: how a multiallelic site can still reach PCA and UMAP -------------------

test_that("one-hot expansion reduces to the dosage matrix on a biallelic panel", {
  skip_if_not_installed("SeqArray")
  # This is why it is the right generalisation rather than an alternative: dropping the
  # reference column and keeping one indicator per ALT gives back exactly the alt-dosage
  # column when there is one ALT. So a biallelic panel's PCA does not move.
  ps <- .multi_ps()
  oh <- ps$genotype(needs = "onehot")
  dos <- ps$genotype(needs = "dosage")
  shared <- intersect(colnames(oh), colnames(dos))
  expect_gt(length(shared), 100)
  expect_equal(oh[, shared, drop = FALSE], dos[, shared, drop = FALSE])
})

test_that("a multiallelic site becomes one column per alternate, named by the allele", {
  skip_if_not_installed("SeqArray")
  ps <- .multi_ps()
  oh <- ps$genotype(needs = "onehot")
  s <- ps$sites()
  multi <- s[s$n_alt_real > 1, ][1, ]
  cols <- grep(paste0("^", multi$site_key, ":"), colnames(oh), value = TRUE)
  expect_equal(length(cols), multi$n_alt_real)
  # the allele is in the name, so a loading can be traced back to a base
  expect_true(all(sub("^.*:", "", cols) %in% multi$alt[[1]]))
  # a sample carries at most one of them, since it has one allele
  expect_true(all(rowSums(oh[, cols, drop = FALSE] > 0, na.rm = TRUE) <= 1))
})

test_that("the expanded panel is wider than the dosage one and reaches PCA", {
  skip_if_not_installed("SeqArray")
  ps <- .multi_ps()
  expect_gt(ncol(ps$genotype(needs = "onehot")), ncol(ps$genotype(needs = "dosage")))
  p <- pop_structure(list(genotype = ps$genotype(needs = "onehot"), encoding = "dosage"),
                     n_pcs = 3)
  expect_equal(nrow(p$pca), nrow(ps$genotype()))
})

test_that("a dosage-only object can still be asked for one-hot, and is unchanged", {
  skip_if_not_installed("SNPRelate")
  ps <- example_pop_structure(umap = FALSE)
  expect_identical(ps$genotype(needs = "onehot"), ps$genotype())
})

test_that("admixture takes the derived biallelic subset instead of refusing", {
  skip_if_not_installed("SeqArray")
  skip_if_not_installed("LEA")
  # sNMF's `.geno` alphabet is 0/1/2/9, so a multiallelic locus has no representation in it,
  # and one-hot would not help: sNMF would read the indicator columns as independent loci
  # when they are perfectly anti-correlated. The biallelic subset is the honest answer, and a
  # genome-wide ancestry summary is where dropping 50 of 577 sites costs least.
  ps <- .multi_ps()
  expect_no_error(ps$run_snmf(K = 2, rep = 1))
  q <- ps$q(K = 2)
  expect_equal(nrow(q), nrow(ps$genotype()))
})
