# Diversity and differentiation, written for k alleles.
#
# Every one of these is a Gini-Simpson heterozygosity, `H = 1 - sum(p_i^2)` -- the chance two
# randomly drawn gene copies differ. At two alleles that is `2p(1-p)`, which is what the code
# wrote, and the biallelic form is a special case rather than a different formula. So the
# generalisation is a substitution, and on a biallelic panel it must move nothing.
#
# The blocker was never the formula: it was that `.group_freqs()` produced one ALT frequency
# per SNP from a dosage matrix, so there was no `p_i` vector to sum over.

test_that("the site heterozygosity is the Gini-Simpson form", {
  h <- plasgenomicsutilsR:::.site_het_k
  # two alleles: 2p(1-p), with the n/(n-1) correction
  expect_equal(h(cbind(0.7, 0.3), 100), 2 * 0.7 * 0.3 * 100 / 99)
  # three: 1 - sum(p^2)
  expect_equal(h(cbind(0.5, 0.25, 0.25), 100), (1 - (0.25 + 0.0625 + 0.0625)) * 100 / 99)
  # monomorphic is zero, not NA
  expect_equal(h(cbind(1, 0, 0), 100), 0)
  # fewer than two gene copies says nothing
  expect_true(is.na(h(cbind(0.5, 0.5), 1)))
})

test_that("it reduces to the old biallelic formula exactly", {
  old <- plasgenomicsutilsR:::.site_het
  new <- plasgenomicsutilsR:::.site_het_k
  set.seed(1)
  p <- runif(200, 0.01, 0.99)
  n <- sample(4:200, 200, TRUE)
  expect_equal(new(cbind(1 - p, p), n), old(p, n), tolerance = 1e-12)
})

test_that("Jost's D and its relatives take a per-allele frequency matrix", {
  d <- plasgenomicsutilsR:::.pair_diff_k
  # two groups, one site, three alleles: A fixed for allele 1, B fixed for allele 2
  pa <- array(c(1, 0, 0), dim = c(1, 3))
  pb <- array(c(0, 1, 0), dim = c(1, 3))
  v <- d(pa, pb, 100, 100, "jost_d", clamp = TRUE)
  expect_equal(unname(v), 1, tolerance = 1e-6)   # complete differentiation
  # and no differentiation when the two are identical
  expect_equal(unname(d(pa, pa, 100, 100, "jost_d", clamp = TRUE)), 0, tolerance = 1e-6)
})

test_that("the k-allele differentiation reduces to the biallelic one", {
  set.seed(2)
  pa <- runif(150, 0.05, 0.95); pb <- runif(150, 0.05, 0.95)
  Na <- rep(80, 150); Nb <- rep(60, 150)
  for (st in c("jost_d", "gst_hedrick", "fst")) {
    old <- plasgenomicsutilsR:::.pair_diff(pa, pb, Na, Nb, st, clamp = TRUE)
    new <- plasgenomicsutilsR:::.pair_diff_k(cbind(1 - pa, pa), cbind(1 - pb, pb),
                                             Na, Nb, st, clamp = TRUE)
    expect_equal(unname(new), unname(old), tolerance = 1e-12, info = st)
  }
})

test_that("a third allele changes the answer, which is the point", {
  # A is half allele-1 half allele-2; B is half allele-1 half allele-3. Collapsing the
  # alternates would make them look identical; they share only half their diversity.
  d <- plasgenomicsutilsR:::.pair_diff_k
  pa <- array(c(0.5, 0.5, 0.0), dim = c(1, 3))
  pb <- array(c(0.5, 0.0, 0.5), dim = c(1, 3))
  v <- d(pa, pb, 200, 200, "jost_d", clamp = TRUE)
  expect_gt(v, 0.1)
})

# --- end to end, off a real multiallelic object -------------------------------------

test_that("pop_diversity computes k-allele diversity when asked", {
  skip_if_not_installed("SeqArray")
  ps <- example_pop_structure("multiallelic", umap = FALSE)
  bi <- pop_diversity(ps)                       # the derived biallelic view
  k <- pop_diversity(ps, alleles = "index")     # every site, k-allele forms
  expect_gt(k$n_snps, bi$n_snps)
  # more sites and more alleles per site, so more diversity to find
  expect_gt(k$he, 0)
})

test_that("on a biallelic panel the two routes agree exactly", {
  skip_if_not_installed("SeqArray")
  # The k-allele forms must not move a biallelic number, so the two routes are run over the
  # same sites. The index view is *derived* from the dosage one here, which is only allowed
  # because the panel's `sites` table says every site has one alternate -- a dosage of 2 at a
  # collapsed multiallelic site could be any of the alternates.
  vcf <- tempfile(fileext = ".vcf")
  set.seed(4)
  n <- 12
  gts <- replicate(30, paste(sample(c("0/0", "1/1", "./."), n, TRUE,
                                    prob = c(.5, .4, .1)), collapse = "\t"))
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=Pf3D7_01_v3,length=640851>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
            sprintf("s%02d", seq_len(n))), collapse = "\t"),
    vapply(seq_along(gts), function(i)
      paste(c("Pf3D7_01_v3", i * 1000, ".", "A", "G", ".", "PASS", ".", "GT", gts[i]),
            collapse = "\t"), character(1))), vcf)
  g <- load_genotypes(vcf, gds = tempfile(fileext = ".gds"), prune = FALSE)
  ps <- PopStructure$new(g, meta = data.frame(sample = g$sample.id, country = "X"))

  a <- pop_diversity(ps)
  b <- pop_diversity(ps, alleles = "index")
  expect_equal(a$he, b$he, tolerance = 1e-12)
  expect_equal(a$pi, b$pi, tolerance = 1e-12)
  expect_equal(a$seg_sites, b$seg_sites)
})

test_that("a dosage panel with no sites table refuses to be read as indices", {
  skip_if_not_installed("SNPRelate")
  # nothing says its sites are biallelic, and a dosage of 2 at a collapsed multiallelic site
  # could be any alternate, so guessing is the one thing not to do
  ps <- example_pop_structure(umap = FALSE)
  expect_error(pop_diversity(ps, alleles = "index"), "cannot be recovered")
})

test_that("pop_diff scores every site under the index route, not just the biallelic ones", {
  skip_if_not_installed("SeqArray")
  ps <- example_pop_structure("multiallelic", umap = FALSE)
  a <- pop_diff(ps, group = "country")                        # derived biallelic view
  b <- pop_diff(ps, group = "country", alleles = "index")     # every site
  # `D` is sites x group-pairs
  expect_equal(ncol(a$D), 1L)
  expect_gt(nrow(b$D), nrow(a$D))
  expect_true(any(is.finite(b$D)))
  # the multiallelic sites are the ones the dosage route could not reach at all
  s <- ps$sites()
  multi <- s$site_key[s$n_alt_real > 1]
  expect_true(all(multi %in% rownames(b$D)))
  expect_false(any(multi %in% rownames(a$D)))
})

test_that("segregating sites count alleles, not the collapsed frequency", {
  skip_if_not_installed("SeqArray")
  # a site where every sample carries an alternate is polymorphic if the alternates differ,
  # and `0 < p < 1` on the collapsed frequency says it is not
  G <- matrix(c(1L, 1L, 2L, 2L), nrow = 4,
              dimnames = list(paste0("s", 1:4), "Pf3D7_01_v3:100"))
  d <- pop_diversity(list(genotype = G, encoding = "allele_index"), alleles = "index")
  expect_equal(d$seg_sites, 1L)
})

test_that("beta_score is biallelic by definition and says so rather than collapsing", {
  skip_if_not_installed("SeqArray")
  # Beta1 folds the SFS, `min(x, 1-x)`, which presupposes one allele and its complement. A
  # three-allele site has no single folded frequency, so there is nothing to generalise --
  # unlike Jost's D or pi, where the biallelic form was only a special case.
  G <- matrix(c(0L, 1L, 2L, 0L, 1L, 2L), nrow = 3,
              dimnames = list(paste0("s", 1:3), c("Pf3D7_01_v3:100", "Pf3D7_01_v3:200")))
  expect_error(beta_score(list(genotype = G, encoding = "allele_index")), "allele_index")
})

test_that("a multiallelic object still scores its biallelic sites", {
  skip_if_not_installed("SeqArray")
  # the derived dosage view is what it gets, so the score exists over the sites Beta is
  # defined for and the multiallelic ones are absent rather than mangled
  ps <- example_pop_structure("multiallelic", umap = FALSE)
  b <- beta_score(ps, window = 20000, min_window_snps = 2)
  expect_true(nrow(b) > 0)
  expect_true(all(b$snp_id %in% colnames(ps$genotype(needs = "dosage"))))
})
