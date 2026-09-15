# iHS at a marker with more than two alleles.
#
# `log(iHH_A / iHH_B)` is a ratio, so it needs exactly two terms and there is no k-allele iHS.
# rehh's `scan_hh` resolves that by keeping the two commonest alleles and dropping the rest,
# which is a two-allele contrast on a subset of the haplotypes reported as though it were the
# whole marker. The right answer is k-1 *pairwise contrasts*, each naming which two alleles it
# compared.
#
# The thing that makes this cheap: per-allele iHH is **invariant** to which other alleles are
# present, because EHH for an allele class only ever involves that class's haplotypes. So a
# contrast is a choice of which two of `calc_ehh()`'s integrals to divide -- no subsetting, no
# re-scanning, and the numbers are identical to what a hand-built subset would give.

.chap <- function(counts, m = 41, seed = 11) {
  set.seed(seed)
  n <- sum(counts)
  G <- matrix(stats::rbinom(n * m, 1, 0.35), n, m)
  G[, (m + 1) %/% 2] <- rep(seq_along(counts) - 1L, times = counts)
  rownames(G) <- sprintf("s%02d", seq_len(n))
  colnames(G) <- paste0("Pf3D7_01_v3:", seq(1000, by = 500, length.out = m))
  parasite_haplotypes(G, maf = 0.02, alleles = "index", impute = FALSE)
}

test_that("per-allele iHH does not depend on which other alleles are at the marker", {
  skip_if_not_installed("rehh")
  # the property the whole design rests on; if it ever stopped holding, contrasts would have
  # to be computed on explicit subsets again
  hap <- .chap(c(26, 23, 11))
  full <- plasgenomicsutilsR:::.marker_allele_ihh(hap, seq_len(nrow(hap$hap)), "Pf3D7_01_v3:11000")
  expect_equal(length(full$ihh), 3L)

  keep <- which(hap$hap[, "Pf3D7_01_v3:11000"] %in% c(0L, 1L))
  sub <- hap
  sub$hap <- sub$hap[keep, , drop = FALSE]
  sub$hap[, "Pf3D7_01_v3:11000"] <- as.integer(sub$hap[, "Pf3D7_01_v3:11000"] == 1L)
  s <- plasgenomicsutilsR:::.marker_allele_ihh(sub, seq_len(nrow(sub$hap)), "Pf3D7_01_v3:11000")
  expect_equal(unname(s$ihh[1]), unname(full$ihh[1]), tolerance = 1e-9)
  expect_equal(unname(s$ihh[2]), unname(full$ihh[2]), tolerance = 1e-9)
})

test_that("a biallelic scan is unchanged by asking for contrasts", {
  skip_if_not_installed("rehh")
  # the safety property: `contrast = "ref"` on a biallelic panel has exactly one contrast per
  # marker, so every number must be the one the old path gave
  hap <- .chap(c(35, 25))
  a <- run_ihs(hap, min_maf = 0.02, contrast = "none")
  b <- run_ihs(hap, min_maf = 0.02, contrast = "ref")
  expect_equal(nrow(a), nrow(b))
  expect_equal(a$unihs, b$unihs, tolerance = 1e-9)
  expect_equal(a$ihs, b$ihs, tolerance = 1e-9)
  expect_equal(a$neg_log10_p, b$neg_log10_p, tolerance = 1e-9)
})

test_that("a three-allele marker gives one row per alternate, each naming its contrast", {
  skip_if_not_installed("rehh")
  hap <- .chap(c(26, 23, 11))
  r <- run_ihs(hap, min_maf = 0.02, contrast = "ref")
  focal <- r[r$snp_id == "Pf3D7_01_v3:11000", ]
  expect_equal(nrow(focal), 2L)
  expect_setequal(focal$contrast, c("0>1", "0>2"))
  # and the two are different measurements, not the same number twice
  expect_false(isTRUE(all.equal(focal$unihs[1], focal$unihs[2])))
})

test_that("the reported frequency is the allele's own, and the marker's sum to one", {
  skip_if_not_installed("rehh")
  hap <- .chap(c(26, 23, 11))
  r <- run_ihs(hap, min_maf = 0.02, contrast = "ref")
  focal <- r[r$snp_id == "Pf3D7_01_v3:11000", ]
  # 23/60 and 11/60 -- the alternates' real frequencies, not "the second commonest allele"
  expect_equal(sort(round(focal$freq_minor, 4)), sort(round(c(23, 11) / 60, 4)))
  expect_equal(round(sum(focal$freq_minor) + 26 / 60, 6), 1)
})

test_that("a biallelic marker keeps one row even when the panel has multiallelic ones", {
  skip_if_not_installed("rehh")
  hap <- .chap(c(26, 23, 11))
  r <- run_ihs(hap, min_maf = 0.02, contrast = "ref")
  other <- r[r$snp_id != "Pf3D7_01_v3:11000", ]
  expect_true(all(table(other$snp_id) == 1L))
  expect_true(all(other$contrast == "0>1"))
})

test_that("pairwise adds the alternate-against-alternate comparison", {
  skip_if_not_installed("rehh")
  hap <- .chap(c(26, 23, 11))
  r <- run_ihs(hap, min_maf = 0.02, contrast = "pairwise")
  focal <- r[r$snp_id == "Pf3D7_01_v3:11000", ]
  expect_setequal(focal$contrast, c("0>1", "0>2", "1>2"))
})

test_that("contrast = none still warns, since it is still dropping alleles", {
  skip_if_not_installed("rehh")
  hap <- .chap(c(26, 23, 11))
  expect_warning(run_ihs(hap, min_maf = 0.02, contrast = "none"), "more than two alleles")
  # and the contrast paths do not, because they are not dropping any
  expect_no_warning(run_ihs(hap, min_maf = 0.02, contrast = "ref"))
})

test_that("the aggregators keep two contrasts at one position apart", {
  skip_if_not_installed("rehh")
  hap <- .chap(c(26, 23, 11))
  r <- run_ihs(hap, min_maf = 0.02, contrast = "ref")
  w <- ihs_windows(r, window = 50000, step = 50000, min_snps = 1)
  expect_true("contrast" %in% names(w))
  # without the split, the focal marker's two rows would be counted as two SNPs in one window
  expect_true(any(w$contrast == "0>2"))
  expect_true(all(w$n_snps[w$contrast == "0>2"] == 1L))

  pk <- selection_peaks(r, criterion = "top", top = 0.5)
  expect_true("contrast" %in% names(pk) || any(grepl("0>", pk$group)))

  gn <- ihs_genes(r, genes = data.frame(name = "focal", chr = "Pf3D7_01_v3",
                                        start = 10500, end = 11500))
  expect_true("contrast" %in% names(gn))
  # the gene is summarised once per contrast, not once with whichever allele scored higher
  expect_setequal(gn$contrast, c("0>1", "0>2"))
})

test_that("the Manhattan plot gives each contrast its own facet", {
  skip_if_not_installed("rehh")
  skip_if_not_installed("ggplot2")
  hap <- .chap(c(26, 23, 11))
  r <- run_ihs(hap, min_maf = 0.02, contrast = "ref")
  # the synthetic panel is too small for the p-values to be finite, so plot the score
  p <- plot_ihs(r, metric = "ihs")
  b <- ggplot2::ggplot_build(p)
  # two contrasts, so two panels: overplotted they would be two points on one x with nothing
  # to tell them apart
  expect_gte(length(unique(b$data[[1]]$PANEL)), 2L)
})

test_that("a biallelic scan is drawn exactly as before, with no extra facet", {
  skip_if_not_installed("rehh")
  skip_if_not_installed("ggplot2")
  hap <- .chap(c(35, 25))
  r <- run_ihs(hap, min_maf = 0.02, contrast = "ref")
  expect_equal(length(unique(r$contrast)), 1L)
  b <- ggplot2::ggplot_build(plot_ihs(r, metric = "ihs"))
  expect_equal(length(unique(b$data[[1]]$PANEL)), 1L)
})


# `parasite_haplotypes()` given a raw `load_genotypes()` list must honour the list's own
# encoding. The raw-list branch used to hardcode the dosage requirement, so a list built with
# `encoding = "allele_index"` -- whose whole purpose is to keep the alternates apart -- was
# refused here as "not dosages", and the only way through was to wrap it in a PopStructure
# first. "auto" now reads the list's `encoding`, and `alleles = "index"` is honoured.

.index_list <- function() {
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=Pf3D7_01_v3,length=640851>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
            sprintf("s%02d", 1:12)), collapse = "\t"),
    paste(c("Pf3D7_01_v3", "1000", ".", "A", "C,G", ".", "PASS", ".", "GT",
            rep(c("0/0", "1/1", "2/2"), length.out = 12)), collapse = "\t"),
    paste(c("Pf3D7_01_v3", "2000", ".", "A", "T", ".", "PASS", ".", "GT",
            rep(c("0/0", "1/1"), length.out = 12)), collapse = "\t"),
    paste(c("Pf3D7_01_v3", "3000", ".", "C", "G", ".", "PASS", ".", "GT",
            rep(c("0/0", "1/1"), length.out = 12)), collapse = "\t")), vcf)
  suppressMessages(load_genotypes(vcf, gds = tempfile(fileext = ".gds"), prune = FALSE,
                                  variants = "all", encoding = "allele_index"))
}

test_that("a raw allele-index list reaches the haplotype path via auto", {
  skip_if_not_installed("SeqArray")
  g <- .index_list()
  expect_equal(g$encoding, "allele_index")
  hap <- suppressMessages(suppressWarnings(
    parasite_haplotypes(g, maf = 0, max_snp_missing = 0.9, max_sample_missing = 0.9,
                        impute = FALSE)))
  expect_equal(hap$alleles, "index")
  # the triallelic site keeps its three states rather than being collapsed to a dosage
  j <- match("Pf3D7_01_v3:999", hap$map$snp_id)
  expect_length(unique(hap$hap[!is.na(hap$hap[, j]), j]), 3L)
})

test_that("alleles = \"index\" on a raw list is honoured, not refused as non-dosage", {
  skip_if_not_installed("SeqArray")
  g <- .index_list()
  expect_no_error(suppressMessages(suppressWarnings(
    parasite_haplotypes(g, alleles = "index", maf = 0,
                        max_snp_missing = 0.9, max_sample_missing = 0.9, impute = FALSE))))
})
