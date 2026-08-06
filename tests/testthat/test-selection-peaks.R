# Collapsing a per-SNP scan into loci. The point of the function is that a count of
# significant SNPs is not a count of findings, so most of these check the merging.

.scan <- function(pos, value, chr = "Pf3D7_07_v3", group = "a") {
  data.frame(group = group, chr = chr, pos = pos, neg_log10_p = value,
             stringsAsFactors = FALSE)
}

test_that("neighbouring hits become one peak and distant ones do not", {
  s <- .scan(c(1000, 1500, 2000, 90000, 95000), c(6, 8, 7, 6, 6))
  pk <- selection_peaks(s, criterion = "value", cutoff = 5, gap = 20000)
  expect_equal(nrow(pk), 2)
  wide <- pk[which.max(pk$n_snps), ]
  expect_equal(wide$n_snps, 3L)
  expect_equal(wide$start, 1000)
  expect_equal(wide$end, 2001)              # half-open, past the last hit
  expect_equal(wide$peak_pos, 1500)         # the strongest SNP, not the first
  expect_equal(wide$peak_value, 8)
})

test_that("gap controls how aggressively hits merge", {
  s <- .scan(c(1000, 30000), c(6, 6))
  expect_equal(nrow(selection_peaks(s, criterion = "value", cutoff = 5, gap = 20000)), 2)
  expect_equal(nrow(selection_peaks(s, criterion = "value", cutoff = 5, gap = 50000)), 1)
})

test_that("peaks never span chromosomes or groups", {
  s <- rbind(.scan(1000, 6, chr = "Pf3D7_07_v3"),
             .scan(1100, 6, chr = "Pf3D7_08_v3"),
             .scan(1000, 6, chr = "Pf3D7_07_v3", group = "b"))
  pk <- selection_peaks(s, criterion = "value", cutoff = 5)
  expect_equal(nrow(pk), 3)
  expect_setequal(as.character(pk$group), c("a", "a", "b"))
})

test_that("min_snps drops isolated single-SNP hits", {
  s <- .scan(c(1000, 1200, 50000), c(6, 6, 9))
  expect_equal(nrow(selection_peaks(s, criterion = "value", cutoff = 5, min_snps = 1)), 2)
  lone <- selection_peaks(s, criterion = "value", cutoff = 5, min_snps = 2)
  expect_equal(nrow(lone), 1)
  expect_equal(lone$n_snps, 2L)             # the tall singleton is gone, despite being tallest
})

test_that("the criteria select different SNPs", {
  s <- .scan(seq(1000, by = 100000, length.out = 100), seq(1, 10, length.out = 100))
  s$significant <- s$neg_log10_p > 9
  s$significant_fdr <- s$neg_log10_p > 8
  expect_lt(sum(selection_peaks(s, criterion = "bonferroni")$n_snps),
            sum(selection_peaks(s, criterion = "fdr")$n_snps))
  top <- selection_peaks(s, criterion = "top", top = 0.1)
  expect_equal(sum(top$n_snps), 10L)
})

test_that("a threshold table is used when the scan carries no flag", {
  s <- .scan(c(1000, 2000), c(4, 7))
  thr <- data.frame(group = "a", threshold = 5, neg_log10_p_fdr_threshold = 3)
  expect_equal(sum(selection_peaks(s, "bonferroni", thresholds = thr)$n_snps), 1L)
  expect_equal(sum(selection_peaks(s, "fdr", thresholds = thr)$n_snps), 2L)
  expect_error(selection_peaks(s, "bonferroni"), "no `significant` column")
})

test_that("genes are named at the peak SNP, not at the interval's midpoint", {
  # the top SNP is in pfcrt at the far left; the interval's midpoint falls in `other`
  s <- .scan(c(403500, 404000, 500000), c(6, 9, 6))
  genes <- data.frame(name = c("pfcrt", "other"), chrom = "Pf3D7_07_v3",
                      start = c(403221, 480000), end = c(406317, 520000))
  pk <- selection_peaks(s, criterion = "value", cutoff = 5, gap = 200000, genes = genes)
  expect_equal(nrow(pk), 1)
  expect_equal(pk$peak_pos, 404000)         # the strongest SNP
  expect_equal(pk$peak_genes, "pfcrt")      # ...and the gene containing it
  expect_equal(pk$distance_to_gene, 0)
  expect_equal(pk$n_genes, 2L)              # the interval reaches both
  # the midpoint (~451,750) is inside neither, which is why it is not the anchor
  expect_true(mean(c(pk$start, pk$end)) > 406317)
})

test_that("an intergenic peak names a gene it covers elsewhere", {
  # the peak spans pfcrt, but its top SNP sits past the gene's end
  s <- .scan(c(404000, 410500), c(6, 8))
  genes <- data.frame(name = "pfcrt", chrom = "Pf3D7_07_v3",
                      start = 403221, end = 406317)
  pk <- selection_peaks(s, criterion = "value", cutoff = 5, genes = genes)
  expect_equal(pk$peak_genes, "")           # the top SNP is not inside it
  expect_equal(pk$nearest_gene, "pfcrt")    # ...but the peak does cover it
  expect_equal(pk$distance_to_gene, 410500 - (406317 - 1))
  expect_equal(pk$n_genes, 1L)
})

test_that("a gene outside the peak is never named, however close it is", {
  # the whole peak sits past pfcrt: covering nothing means naming nothing
  s <- .scan(c(410000, 410500), c(6, 8))
  genes <- data.frame(name = "pfcrt", chrom = "Pf3D7_07_v3",
                      start = 403221, end = 406317)
  pk <- selection_peaks(s, criterion = "value", cutoff = 5, genes = genes)
  expect_equal(pk$n_genes, 0L)
  expect_equal(pk$peak_genes, "")
  expect_equal(pk$nearest_gene, "")         # not "pfcrt, 3.7 kb away"
  expect_true(is.na(pk$distance_to_gene))

  # padding the peak far enough to reach it brings it back
  padded <- selection_peaks(s, criterion = "value", cutoff = 5, genes = genes, pad = 5000)
  expect_equal(padded$n_genes, 1L)
  expect_equal(padded$nearest_gene, "pfcrt")
})

test_that("the three gene columns cannot contradict each other", {
  # a short gene list, which is when this bites: most peaks cover none of it
  s <- rbind(.scan(c(404000, 405000), c(6, 9)),                    # inside pfcrt
             .scan(c(900000, 901000), c(6, 7)),                    # nowhere near
             .scan(c(404000, 405000), c(6, 9), group = "b"))
  genes <- data.frame(name = "pfcrt", chrom = "Pf3D7_07_v3",
                      start = 403221, end = 406317)
  pk <- selection_peaks(s, criterion = "value", cutoff = 5, genes = genes)
  expect_true(all(pk$nearest_gene[pk$n_genes == 0] == ""))
  expect_true(all(is.na(pk$distance_to_gene[pk$n_genes == 0])))
  expect_true(all(pk$peak_genes[pk$n_genes == 0] == ""))
  expect_true(all(pk$distance_to_gene[pk$peak_genes != ""] == 0))
  expect_true(all(pk$nearest_gene[pk$peak_genes != ""] == pk$peak_genes[pk$peak_genes != ""]))
})

test_that("peak_genes holds every gene covering the peak SNP when spans overlap", {
  s <- .scan(1000, 9)
  genes <- data.frame(name = c("a", "b", "c"), chrom = "Pf3D7_07_v3",
                      start = c(500, 900, 5000), end = c(1500, 1200, 6000))
  pk <- selection_peaks(s, criterion = "value", cutoff = 5, genes = genes)
  expect_equal(pk$peak_genes, "a,b")        # comma-separated, sorted
  expect_equal(pk$distance_to_gene, 0)
})

test_that("distance is measured to the closer edge, half-open end included", {
  genes <- data.frame(name = "g", chrom = "Pf3D7_07_v3", start = 1000, end = 2000)
  # peaks wide enough to cover the gene, with the top SNP on either side of it
  before <- selection_peaks(.scan(c(900, 1500), c(9, 6)), criterion = "value",
                            cutoff = 5, genes = genes)
  after <- selection_peaks(.scan(c(1500, 2100), c(6, 9)), criterion = "value",
                           cutoff = 5, genes = genes)
  expect_equal(before$distance_to_gene, 100)
  expect_equal(after$distance_to_gene, 2100 - 1999)   # `end` is exclusive
  # the last base inside the gene is 1999, and it counts as inside
  expect_equal(selection_peaks(.scan(c(1999, 2500), c(9, 6)), criterion = "value",
                               cutoff = 5, genes = genes)$distance_to_gene, 0)
})

test_that("nothing qualifying returns an empty frame with a warning", {
  s <- .scan(c(1000, 2000), c(1, 2))
  expect_warning(pk <- selection_peaks(s, criterion = "value", cutoff = 5), "no SNP met")
  expect_equal(nrow(pk), 0)
})

test_that("it reads any of the package's scans", {
  skip_if_not_installed("ggplot2")
  ps <- example_pop_structure(umap = FALSE)
  b <- beta_score(ps, group = "country", window = 300000, min_window_snps = 1)
  pk <- selection_peaks(b, criterion = "top", top = 0.3)
  expect_true(nrow(pk) > 0)
  expect_equal(attr(pk, "peak_metric"), "beta")   # picked the right column unprompted
  expect_true(all(pk$end > pk$start))
})

test_that("pad widens the reported interval", {
  s <- .scan(c(10000, 10500), c(6, 6))
  a <- selection_peaks(s, criterion = "value", cutoff = 5)
  b <- selection_peaks(s, criterion = "value", cutoff = 5, pad = 5000)
  expect_equal(b$start, a$start - 5000)
  expect_equal(b$end, a$end + 5000)
  # never negative
  expect_gte(selection_peaks(.scan(100, 6), criterion = "value", cutoff = 5,
                             pad = 5000)$start, 0)
})
