# Haplotype preparation and the rehh-backed scans. The preparation is where the
# judgement calls live, so most of these test what it drops and what it reports.

.hap_geno <- function(n = 30, m = 30, seed = 1) {
  set.seed(seed)
  p <- stats::runif(m, 0.15, 0.85)
  G <- matrix(2 * stats::rbinom(n * m, 1, rep(p, each = n)), n, m)
  dimnames(G) <- list(paste0("s", seq_len(n)),
                      paste0("Pf3D7_01_v3:", seq(1000, by = 2000, length.out = m)))
  G
}

test_that("haplotypes come back complete, 0/1 and integer", {
  G <- .hap_geno()
  G[1, 1] <- NA
  h <- parasite_haplotypes(G, maf = 0)
  expect_s3_class(h, "parasite_haplotypes")
  expect_false(anyNA(h$hap))
  expect_setequal(unique(as.vector(h$hap)), c(0L, 1L))
  expect_type(h$hap[1, 1], "integer")            # rehh rejects a double matrix
  expect_equal(nrow(h$map), ncol(h$hap))
})

test_that("the Fws gate keeps only monoclonal infections", {
  G <- .hap_geno(n = 10)
  fws <- stats::setNames(c(rep(0.99, 6), rep(0.5, 4)), rownames(G))
  h <- parasite_haplotypes(G, fws = fws, min_fws = 0.95, maf = 0)
  expect_equal(nrow(h$hap), 6)
  expect_equal(h$filtering$n_dropped_polyclonal, 4)
  # a data frame works the same way
  df <- data.frame(sample = names(fws), fws = as.numeric(fws))
  expect_equal(nrow(parasite_haplotypes(G, fws = df, maf = 0)$hap), 6)
})

test_that("a sample with no Fws is an error, not a silent drop", {
  G <- .hap_geno(n = 5)
  fws <- stats::setNames(rep(0.99, 3), rownames(G)[1:3])
  expect_error(parasite_haplotypes(G, fws = fws), "no Fws for 2 sample")
})

test_that("mixed calls are resolved by an allele draw, or left to imputation", {
  G <- .hap_geno(n = 20, m = 10)
  G[1:5, 1] <- 1                                  # five mixed calls at one SNP
  drawn <- parasite_haplotypes(G, het = "sample", maf = 0, seed = 3)
  expect_equal(drawn$filtering$n_het_calls, 5)
  expect_equal(drawn$filtering$n_imputed, 0)      # the draw already filled them
  # left as missing they are 25% of that SNP, so it survives only a laxer ceiling
  missing <- parasite_haplotypes(G, het = "missing", maf = 0, max_snp_missing = 0.3,
                                 seed = 3)
  expect_equal(missing$filtering$n_imputed, 5)
  expect_equal(parasite_haplotypes(G, het = "missing", maf = 0,
                                   seed = 3)$filtering$n_dropped_snp_missing, 1)
})

test_that("the same seed gives the same haplotypes and a different one does not", {
  G <- .hap_geno(n = 20, m = 12)
  G[1:8, 2] <- 1
  a <- parasite_haplotypes(G, maf = 0, seed = 11)
  b <- parasite_haplotypes(G, maf = 0, seed = 11)
  c2 <- parasite_haplotypes(G, maf = 0, seed = 12)
  expect_identical(a$hap, b$hap)
  expect_false(identical(a$hap, c2$hap))
})

test_that("MAF and missingness filters are applied and counted", {
  # a deterministic panel, so only the two SNPs planted below can fail a filter
  G <- matrix(rep(c(0, 2), each = 10), 20, 10,
              dimnames = list(paste0("s", 1:20), paste0("Pf3D7_01_v3:", 1:10 * 1000)))
  G[, 1] <- 0; G[1, 1] <- 2                        # one rare-allele SNP (MAF 0.05)
  G[2:19, 2] <- NA                                 # one mostly-missing SNP
  h <- parasite_haplotypes(G, maf = 0.1, max_snp_missing = 0.1, max_sample_missing = 1)
  expect_equal(h$filtering$n_dropped_snp_missing, 1)
  expect_equal(h$filtering$n_dropped_maf, 1)
  expect_equal(ncol(h$hap), 8)
})

test_that("impute = FALSE drops incomplete SNPs instead of filling them", {
  G <- .hap_geno(n = 20, m = 10)
  G[1, 3] <- NA
  kept <- parasite_haplotypes(G, maf = 0, impute = TRUE)
  dropped <- parasite_haplotypes(G, maf = 0, impute = FALSE)
  expect_equal(ncol(kept$hap), 10)
  expect_equal(ncol(dropped$hap), 9)
  expect_equal(dropped$filtering$n_imputed, 0)
  expect_false(anyNA(dropped$hap))
})

test_that("filtering everything away is an error with a reason", {
  G <- .hap_geno(n = 20, m = 10)
  expect_error(parasite_haplotypes(G, maf = 0.6), "no SNP passed maf")
  fws <- stats::setNames(rep(0.1, 20), rownames(G))
  expect_error(parasite_haplotypes(G, fws = fws), "no sample has Fws")
})

test_that("print reports how the haplotypes were made", {
  G <- .hap_geno(n = 20, m = 10)
  G[1:3, 1] <- 1
  out <- utils::capture.output(print(parasite_haplotypes(G, maf = 0)))
  expect_match(paste(out, collapse = " "), "parasite_haplotypes")
  expect_match(paste(out, collapse = " "), "mixed calls")
  expect_match(paste(out, collapse = " "), "seed")
})


# --------------------------------------------------------------------------- #
#  Scans                                                                       #
# --------------------------------------------------------------------------- #

test_that("run_ihs returns one row per scored SNP per group", {
  skip_if_not_installed("rehh")
  ps <- example_pop_structure(umap = FALSE)
  hap <- parasite_haplotypes(ps, maf = 0.05)
  ihs <- suppressWarnings(run_ihs(hap, group = "country"))
  expect_s3_class(ihs, "data.frame")
  expect_true(all(c("group", "chr", "pos", "snp_id", "ihs", "neg_log10_p") %in% names(ihs)))
  expect_setequal(as.character(unique(ihs$group)), c("Ghana", "Cambodia"))
  expect_true(all(ihs$pos > 0))
})

test_that("the cross-population scans cover every group pair", {
  skip_if_not_installed("rehh")
  ps <- example_pop_structure(umap = FALSE)
  hap <- parasite_haplotypes(ps, maf = 0.05)
  rsb <- suppressWarnings(run_rsb(hap, group = "country"))
  xp <- suppressWarnings(run_xpehh(hap, group = "country"))
  expect_equal(unique(rsb$pair), "Cambodia vs Ghana")
  expect_true(all(c("pop1", "pop2", "value", "neg_log10_p") %in% names(rsb)))
  expect_equal(nrow(rsb), nrow(xp))            # same SNPs, different statistic
  expect_false(isTRUE(all.equal(rsb$value, xp$value)))
})

test_that("a cross-population scan needs two groups", {
  skip_if_not_installed("rehh")
  ps <- example_pop_structure(umap = FALSE)
  hap <- parasite_haplotypes(ps, maf = 0.05)
  meta <- ps$get_meta()
  meta$one <- "all"
  expect_error(run_rsb(hap, group = "one", meta = meta), "at least two groups")
})

test_that("ihs_genes keeps the strongest SNP in each gene", {
  scan <- data.frame(
    group = factor("a"), chr = "Pf3D7_07_v3",
    pos = c(403500, 404000, 900000),
    snp_id = c("x", "y", "z"), ihs = c(1, -3, 0.5),
    neg_log10_p = c(1.0, 4.2, 0.3), stringsAsFactors = FALSE)
  genes <- data.frame(name = "pfcrt", chrom = "Pf3D7_07_v3", start = 403221, end = 406317)
  out <- ihs_genes(scan, genes = genes, min_snps = 1)
  expect_equal(nrow(out), 1)
  expect_equal(out$n_snps, 2L)                 # the third SNP is outside the CDS
  expect_equal(out$max_neg_log10_p, 4.2)
  expect_equal(out$max_abs_value, 3)           # takes the magnitude, sign is unpolarised
  expect_equal(out$peak_pos, 404000)
})

test_that("within widens the gene window", {
  scan <- data.frame(group = factor("a"), chr = "Pf3D7_07_v3", pos = 402000,
                     snp_id = "x", ihs = 2, neg_log10_p = 3, stringsAsFactors = FALSE)
  genes <- data.frame(name = "pfcrt", chrom = "Pf3D7_07_v3", start = 403221, end = 406317)
  expect_equal(nrow(ihs_genes(scan, genes = genes, min_snps = 1)), 0)
  expect_equal(nrow(ihs_genes(scan, genes = genes, within = 5000, min_snps = 1)), 1)
})


# --------------------------------------------------------------------------- #
#  Plots                                                                       #
# --------------------------------------------------------------------------- #

test_that("the scan plots build", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("rehh")
  ps <- example_pop_structure(umap = FALSE)
  hap <- parasite_haplotypes(ps, maf = 0.05)
  ihs <- suppressWarnings(run_ihs(hap, group = "country"))
  p <- plot_ihs(ihs, genes = PF_EXAMPLE_DRUG_GENES)
  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(p))
  expect_length(attr(p, "plasgenomics_dims"), 2)
  # a pair-keyed scan facets on `pair` instead of `group`
  rsb <- suppressWarnings(run_rsb(hap, group = "country"))
  expect_s3_class(plot_ihs(rsb), "ggplot")
})

test_that("plot_diversity insists on a windowed result", {
  skip_if_not_installed("ggplot2")
  ps <- example_pop_structure(umap = FALSE)
  expect_error(plot_diversity(pop_diversity(ps, group = "country")), "by = \"window\"")
  d <- pop_diversity(ps, group = "country", by = "window", window = 500000)
  expect_s3_class(plot_diversity(d, metric = "he"), "ggplot")
})

test_that("read_ld_decay round-trips the Python output and plots", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("readr")
  f <- tempfile(fileext = ".tsv")
  writeLines(c("#max_dist=50000\t#max_snps=3000\t#maf=0.05",
    "group\tbin_start\tbin_end\tbin_mid\tn_pairs\tmean_r2\tmedian_r2",
    "a\t0\t5000\t2500\t100\t0.40\t0.35",
    "a\t5000\t10000\t7500\t90\t0.30\t0.25",
    "a\t10000\t15000\t12500\t80\t0.15\t0.10",
    "b\t0\t5000\t2500\t50\t0.20\t0.18",
    "b\t5000\t10000\t7500\t40\t0.19\t0.17",
    "b\t10000\t15000\t12500\t30\t0.18\t0.16"), f)
  ld <- read_ld_decay(f)
  expect_equal(nrow(ld), 6)
  expect_equal(attr(ld, "ld_max_dist"), 50000)   # the scan settings ride along
  expect_equal(attr(ld, "ld_maf"), 0.05)

  hd <- attr(ld, "ld_half_decay")
  # group a halves between the 7500 and 12500 bins (0.40 -> target 0.20)
  expect_true(is.finite(hd$half_decay_bp[hd$group == "a"]))
  expect_gt(hd$half_decay_bp[hd$group == "a"], 7500)
  expect_lt(hd$half_decay_bp[hd$group == "a"], 12500)
  # group b never falls to half, so there is no half-decay to report
  expect_true(is.na(hd$half_decay_bp[hd$group == "b"]))

  p <- plot_ld_decay(ld)
  expect_s3_class(p, "ggplot")
  expect_true(any(vapply(p$layers, function(l) inherits(l$geom, "GeomVline"), logical(1))))
})

test_that("read_ld_decay rejects a table missing its columns", {
  skip_if_not_installed("readr")
  f <- tempfile(fileext = ".tsv")
  writeLines(c("group\tbin_mid", "a\t100"), f)
  expect_error(read_ld_decay(f), "missing column")
})
