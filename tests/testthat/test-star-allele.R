# `*` is not a base.
#
# It marks a position whose sequence is deleted on that haplotype -- a confident observation,
# and the reason `spanning_del_filter` is default-off, but not one of the alleles being
# compared. Every k-allele estimator in the package works off the genotype matrix alone
# (`he`, pi, Jost's D, the one-hot expansion, `allele_states()`), and each one counts
# whatever distinct values it finds there. Left in, a site where most haplotypes are deleted
# reads as a *highly diverse* site rather than a mostly-absent one, and the error is large:
# a two-base site with half its haplotypes deleted has a true H of 0.5 and reads as 0.67.
#
# So the encoding blanks them: a haplotype with no base here gets NA, which is what NA means
# everywhere else in the matrix. Allele indices are untouched, so `sites$alt` still lists `*`
# in its own slot and index 2 still names the same base. `star = "allele"` keeps them, for
# the deliberate case where presence/absence IS the state.

.star_vcf <- function() {
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=Pf3D7_01_v3,length=640851>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
    paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",
           paste(sprintf("s%d", 1:8), collapse = "\t")),
    # 4 x ref, 4 x `*`: two "alleles" that are really one base and an absence
    "Pf3D7_01_v3\t100\t.\tA\t*\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1",
    # ref, one real ALT, and `*`: index 2 must keep meaning G whatever happens to `*`
    "Pf3D7_01_v3\t200\t.\tA\t*,G\t.\tPASS\t.\tGT\t0/0\t0/0\t1/1\t1/1\t2/2\t2/2\t2/2\t2/2",
    "Pf3D7_01_v3\t300\t.\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t1/1"),
    vcf)
  vcf
}

.load_star <- function(star = "missing") {
  suppressMessages(load_genotypes(.star_vcf(), gds = tempfile(fileext = ".gds"),
                                  prune = FALSE, variants = "all",
                                  encoding = "allele_index", star = star))
}

test_that("a `*` call is missing by default, and the record is not dropped", {
  skip_if_not_installed("SeqArray")
  g <- .load_star()
  expect_equal(ncol(g$genotype), 3L)                 # the site is kept, only the calls go
  expect_equal(unname(g$genotype[, "Pf3D7_01_v3:99"]),
               c(0L, 0L, 0L, 0L, NA, NA, NA, NA))
  expect_equal(g$star, "missing")
})

test_that("blanking `*` does not renumber the alleles beside it", {
  skip_if_not_installed("SeqArray")
  # `*` is ALT index 1 here and G is index 2. Dropping the allele from the list would make
  # index 1 point at G and silently relabel every call; blanking the calls cannot.
  g <- .load_star()
  st <- g$sites[g$sites$site_key == "Pf3D7_01_v3:199", ]
  expect_equal(st$alt[[1]], c("*", "G"))
  expect_true(st$has_spanning_del)
  expect_equal(unname(g$genotype[, "Pf3D7_01_v3:199"]),
               c(0L, 0L, NA, NA, 2L, 2L, 2L, 2L))
})

test_that("star = \"allele\" keeps them, for presence/absence as a state", {
  skip_if_not_installed("SeqArray")
  g <- .load_star("allele")
  expect_equal(unname(g$genotype[, "Pf3D7_01_v3:99"]), c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L))
  expect_equal(g$star, "allele")
})

test_that("counting `*` as an allele inflates heterozygosity, and the default does not", {
  skip_if_not_installed("SeqArray")
  # site 199 truly has A x2 and G x4 among the 6 haplotypes with a base: H = 1 - (1/9+4/9)
  # = 4/9, times the 6/5 small-sample correction. Counting `*` makes it 3 alleles at
  # 2/2/4 of 8, H = 1 - (1/16+1/16+1/4) = 0.625 -- the site reads as more diverse for
  # being more deleted.
  ps_m <- PopStructure$new(.load_star(), meta = NULL)
  ps_a <- PopStructure$new(.load_star("allele"), meta = NULL)
  d_m <- pop_diversity(ps_m, alleles = "index", min_snps = 1)
  d_a <- pop_diversity(ps_a, alleles = "index", min_snps = 1)
  expect_lt(d_m$he, d_a$he)
  expect_gt(d_a$he / d_m$he - 1, 0.2)     # not a rounding difference
})

test_that("a deleted haplotype is not scored as reference in the one-hot panel", {
  skip_if_not_installed("SeqArray")
  # the one-hot expansion writes 2 for "carries this ALT" and 0 for "does not". A `*` call
  # is neither -- there is no base to carry or not carry -- and 0 there hands PCA a fake
  # reference haplotype for every deleted sample.
  ps <- PopStructure$new(.load_star(), meta = NULL)
  oh <- ps$genotype(needs = "onehot")
  col <- grep("Pf3D7_01_v3:199", colnames(oh), value = TRUE)
  expect_length(col, 1L)
  expect_equal(unname(oh[, col]), c(0L, 0L, NA, NA, 2L, 2L, 2L, 2L))
})
