test_that("annotate_snps finds the interval each SNP sits in", {
  iv <- data.frame(name = c("a", "b"), chr = "Pf3D7_07_v3",
                   start = c(1000, 5000), end = c(2000, 6000), stringsAsFactors = FALSE)
  s <- data.frame(snp_id = c("Pf3D7_07_v3:1500", "Pf3D7_07_v3:3000", "Pf3D7_07_v3:5500"),
                  v = 1:3, stringsAsFactors = FALSE)
  a <- annotate_snps(s, iv)
  expect_equal(a$name, c("a", NA, "b"))
  expect_equal(nrow(a), 3)                       # keep = "all" keeps the miss
  expect_equal(a$v, 1:3)                         # original columns and order preserved
  expect_equal(nrow(annotate_snps(s, iv, keep = "hits")), 2)

  # half-open: a SNP at `end` is outside, at `start` is inside
  edge <- data.frame(snp_id = c("Pf3D7_07_v3:1000", "Pf3D7_07_v3:2000"), stringsAsFactors = FALSE)
  expect_equal(annotate_snps(edge, iv)$name, c("a", NA))
  # ...unless widened
  expect_equal(annotate_snps(edge, iv, within = 1)$name, c("a", "a"))
})

test_that("chromosome naming and the chr/pos form both work", {
  iv <- data.frame(name = "a", chrom = "7", start = 1000, end = 2000, stringsAsFactors = FALSE)
  # Pf3D7_07_v3 in the scan vs "7" in the intervals, and `chrom` instead of `chr`
  expect_equal(annotate_snps(data.frame(snp_id = "Pf3D7_07_v3:1500"), iv)$name, "a")
  expect_equal(annotate_snps(data.frame(chr = "chr7", pos = 1500), iv)$name, "a")
  expect_error(annotate_snps(data.frame(v = 1), iv), "needs a `snp_id`")
  expect_error(annotate_snps(data.frame(snp_id = "nope"), iv), "could not read a position")
})

test_that("a SNP in several intervals gives a row each, or one collapsed row", {
  iv <- data.frame(name = c("outer", "inner"), chr = "Pf3D7_07_v3",
                   start = c(1000, 1400), end = c(2000, 1600), stringsAsFactors = FALSE)
  s <- data.frame(snp_id = "Pf3D7_07_v3:1500", v = 9, stringsAsFactors = FALSE)
  a <- annotate_snps(s, iv)
  expect_equal(nrow(a), 2)                                    # not silently collapsed
  expect_setequal(a$name, c("outer", "inner"))

  c1 <- annotate_snps(s, iv, collapse = TRUE)
  expect_equal(nrow(c1), 1)
  expect_equal(c1$name, "outer,inner")                        # genomic order
  expect_equal(c1$n_intervals, 2L)
  # prefix keeps a second annotation alongside
  p <- annotate_snps(s, iv, collapse = TRUE, prefix = "core_")
  expect_true(all(c("core_name", "core_n_intervals") %in% names(p)))
})

test_that("gene_id rides along when the interval table has one", {
  a <- annotate_snps(data.frame(snp_id = "Pf3D7_07_v3:403500"), PF_EXAMPLE_DRUG_GENES)
  expect_equal(a$name, "pfcrt")
  expect_equal(a$gene_id, "PF3D7_0709000")
  expect_true(a$interval_start <= 403500 && 403500 < a$interval_end)
})

test_that("1-based SNP positions can be annotated without shifting them by hand", {
  # a 3 bp interval, 0-based half-open: covers 1-based positions 101, 102, 103
  iv <- data.frame(name = "codon", chr = "1", start = 100, end = 103)
  snps <- data.frame(snp_id = c("1:100", "1:101", "1:103", "1:104"))

  # the package convention: positions are 0-based, so 100/101/102 are inside
  zero <- annotate_snps(snps, iv, keep = "hits")
  expect_setequal(zero$snp_id, c("1:100", "1:101"))

  # 1-based (a genotype-matrix id from load_genotypes()): 101/102/103 are inside, and 100 is
  # not -- without saying so, every SNP shifts a base and a near-miss reads as a hit
  one <- annotate_snps(snps, iv, keep = "hits", one_based_snps = TRUE)
  expect_setequal(one$snp_id, c("1:101", "1:103"))
  expect_false("1:100" %in% one$snp_id)
  # the reported interval bounds are untouched; only the test position moves
  expect_equal(unique(one$interval_start), 100)
  expect_equal(unique(one$interval_end), 103)
})
