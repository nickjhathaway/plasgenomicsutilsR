# Diversity, LD and selection statistics. Where an independent implementation exists
# (pegas, poppr) the tests check against it rather than against a remembered number.

# a small haploid panel with known structure
.fake_geno <- function(n = 30, m = 40, chr = "Pf3D7_01_v3", seed = 1, spacing = 1000) {
  set.seed(seed)
  p <- stats::runif(m, 0.1, 0.9)
  G <- matrix(2 * stats::rbinom(n * m, 1, rep(p, each = n)), n, m)
  dimnames(G) <- list(paste0("s", seq_len(n)),
                      paste0(chr, ":", seq(1000, by = spacing, length.out = m)))
  G
}

test_that("pi divides by accessible base pairs, not by SNP count", {
  G <- .fake_geno(m = 20, spacing = 100)
  # 20 SNPs spread over a 100 kb window of a chromosome
  acc <- data.frame(chrom = "Pf3D7_01_v3", start = 0, end = 100000)
  d <- pop_diversity(G, accessible = acc)
  expect_equal(d$n_sites, 100000)
  # he is the per-SNP average; pi is the same numerator over the bp denominator
  expect_equal(d$pi, d$he * d$n_snps / d$n_sites, tolerance = 1e-12)
  expect_lt(d$pi, d$he)                       # and so pi is far smaller than he

  # doubling the accessible span halves pi but leaves he alone -- the whole point
  acc2 <- data.frame(chrom = "Pf3D7_01_v3", start = 0, end = 200000)
  d2 <- pop_diversity(G, accessible = acc2)
  expect_equal(d2$pi, d$pi / 2, tolerance = 1e-12)
  expect_equal(d2$he, d$he)
})

test_that("Tajima's D, haplotype diversity and pi agree with pegas", {
  skip_if_not_installed("pegas")
  skip_if_not_installed("ape")
  set.seed(7)
  n <- 24; m <- 35
  H <- matrix(stats::rbinom(n * m, 1, rep(stats::runif(m, 0.1, 0.9), each = n)), n, m)
  d <- ape::as.DNAbin(t(apply(H, 1, function(r) ifelse(r == 1, "t", "a"))))

  f <- plasgenomicsutilsR:::.snp_freqs(H)
  h <- plasgenomicsutilsR:::.site_het(f$p, f$n)
  expect_equal(tajima_d(h, n), pegas::tajima.test(d)$D, tolerance = 1e-8)
  expect_equal(sum(h) / m, as.numeric(pegas::nuc.div(d)), tolerance = 1e-10)

  hs <- plasgenomicsutilsR:::.haplotype_stats(H)
  expect_equal(hs$hap_div, as.numeric(pegas::hap.div(d)[1]), tolerance = 1e-10)
  expect_equal(hs$n_hap, nrow(pegas::haplotype(d)))
})

test_that("a monomorphic panel has no diversity and no Tajima's D", {
  G <- matrix(0, 10, 8, dimnames = list(paste0("s", 1:10),
                                        paste0("Pf3D7_01_v3:", 1:8 * 100)))
  d <- pop_diversity(G)
  expect_equal(d$seg_sites, 0L)
  expect_equal(d$he, 0)
  expect_equal(d$pi, 0)
  expect_true(is.na(d$tajima_d))              # no segregating sites -> undefined
  expect_equal(d$n_hap, 1L)                   # everyone shares the one haplotype
  expect_equal(d$hap_div, 0)
})

test_that("heterozygous calls are dropped as mixed infections by default", {
  G <- matrix(c(0, 2, 1, 1), 4, 1,
              dimnames = list(paste0("s", 1:4), "Pf3D7_01_v3:100"))
  drop <- pop_diversity(G)                    # het = "missing": 2 calls, p = 0.5
  expect_equal(drop$n_samples, 4L)
  expect_equal(drop$he, 1)                    # 2*0.5*0.5 * 2/1 with n = 2
  dose <- pop_diversity(G, het = "dosage")    # hets contribute 0.5 each
  expect_equal(dose$he, 2 * 0.5 * 0.5 * 4 / 3)
})

test_that("group ordering follows the metadata's factor levels", {
  G <- .fake_geno(n = 20)
  meta <- data.frame(sample = rownames(G),
                     site = factor(rep(c("b", "a"), each = 10), levels = c("b", "a")))
  d <- pop_diversity(G, group = "site", meta = meta)
  expect_equal(levels(d$group), c("b", "a"))
  expect_equal(as.character(d$group), c("b", "a"))
})

test_that("small groups are skipped rather than silently reported", {
  G <- .fake_geno(n = 20)
  meta <- data.frame(sample = rownames(G), site = c(rep("big", 18), "tiny", "tiny"))
  d <- pop_diversity(G, group = "site", meta = meta, min_samples = 4)
  expect_equal(as.character(d$group), "big")
})

test_that("per-gene units pick up only the SNPs inside the CDS", {
  G <- .fake_geno(m = 30, spacing = 100)      # SNPs at 1000, 1100, ... 3900
  genes <- data.frame(name = c("inside", "elsewhere"), chrom = c("Pf3D7_01_v3", "Pf3D7_02_v3"),
                      start = c(1000, 0), end = c(1500, 5000))
  d <- pop_diversity(G, by = "gene", genes = genes)
  expect_equal(d$unit, c("inside", "elsewhere"))
  expect_equal(d$n_snps[d$unit == "inside"], 5L)      # 1000..1400, end exclusive
  expect_equal(d$n_snps[d$unit == "elsewhere"], 0L)
  expect_equal(d$n_sites[d$unit == "inside"], 500)    # gene length, no accessible BED
})

test_that("windows tile the chromosome and report their own span", {
  G <- .fake_geno(m = 30, spacing = 100)
  d <- pop_diversity(G, by = "window", window = 1000, step = 1000)
  expect_true(all(d$end - d$start == 1000))
  expect_equal(sum(d$n_snps), 30L)                    # every SNP lands in exactly one
  expect_true(all(diff(d$start) > 0))
})

test_that("pop_diversity accepts a PopStructure and its grouping", {
  ps <- example_pop_structure(umap = FALSE)
  d <- pop_diversity(ps, group = "country")
  expect_equal(nrow(d), 2)
  expect_setequal(as.character(d$group), c("Ghana", "Cambodia"))
  expect_true(all(d$he > 0 & d$he < 1))
  # these read the full panel, not the thinned one PCA uses
  expect_true(all(d$n_snps == ncol(ps$genotype(prefer = "full"))))
})


# --------------------------------------------------------------------------- #
#  Linkage disequilibrium                                                      #
# --------------------------------------------------------------------------- #

test_that("Ia and rbarD agree with poppr", {
  skip_if_not_installed("poppr")
  skip_if_not_installed("adegenet")
  set.seed(4)
  H <- matrix(stats::rbinom(30 * 20, 1, 0.4), 30, 20)
  gi <- adegenet::df2genind(apply(H, 2, as.character), ploidy = 1, type = "codom")
  ref <- poppr::ia(gi, quiet = TRUE)
  ours <- plasgenomicsutilsR:::.ia_rbard(H)
  expect_equal(ours$ia, unname(ref[["Ia"]]), tolerance = 1e-10)
  expect_equal(ours$rbar_d, unname(ref[["rbarD"]]), tolerance = 1e-10)
})

test_that("perfectly linked loci give the maximum index of association", {
  H <- matrix(rep(stats::rbinom(20, 1, 0.5), 6), 20, 6)   # six copies of one locus
  s <- plasgenomicsutilsR:::.ia_rbard(H)
  expect_equal(s$rbar_d, 1, tolerance = 1e-8)             # rbarD is 1 when fully linked
  expect_equal(s$ia, ncol(H) - 1, tolerance = 1e-8)       # Ia grows with the locus count
})

test_that("beta matches the folded Beta1 definition", {
  # independent transcription of Siewert & Voight's folded statistic
  d_ref <- function(freq, x, p) {
    xf <- min(x, 1 - x); f <- pmin(freq, 1 - freq); md <- max(xf, 0.5 - xf)
    ((md - abs(xf - f)) / md)^p
  }
  ref <- function(wf, x, n, p = 2) {
    i <- seq_len(n - 1)
    sum(d_ref(wf, x, p)) / sum((1 / i) * d_ref(i / n, x, p)) - length(wf) / sum(1 / i)
  }
  set.seed(2)
  wf <- stats::runif(15, 0.05, 0.95)
  expect_equal(plasgenomicsutilsR:::.beta_folded(wf, 0.4, 50, 2), ref(wf, 0.4, 50),
               tolerance = 1e-12)
})

test_that("beta is larger when neighbouring frequencies cluster on the core", {
  clustered <- plasgenomicsutilsR:::.beta_folded(rep(0.45, 12), 0.45, 40, 2)
  dispersed <- plasgenomicsutilsR:::.beta_folded(seq(0.02, 0.5, length.out = 12), 0.45, 40, 2)
  expect_gt(clustered, dispersed)
})

test_that("beta_score excludes the core SNP from its own window", {
  G <- .fake_geno(m = 12, spacing = 500)
  b <- beta_score(G, window = 100000, min_window_snps = 1)
  expect_equal(unique(b$n_window_snps + 1L), 12L)   # every other SNP, never itself
})

test_that("beta_genes summarises to one row per gene", {
  G <- .fake_geno(m = 20, spacing = 100)
  b <- beta_score(G, window = 100000, min_window_snps = 1)
  genes <- data.frame(name = "g1", chrom = "Pf3D7_01_v3", start = 1000, end = 2000)
  bg <- beta_genes(b, genes = genes, min_snps = 1)
  expect_equal(nrow(bg), 1)
  expect_equal(bg$gene, "g1")
  expect_true(bg$beta_max >= bg$beta_mean)
})

test_that("Tajima's D p-values agree with pegas", {
  skip_if_not_installed("pegas")
  skip_if_not_installed("ape")
  set.seed(21)
  for (k in 1:4) {
    n <- sample(15:45, 1); m <- sample(25:60, 1)
    H <- matrix(stats::rbinom(n * m, 1, rep(stats::runif(m, 0.1, 0.9), each = n)), n, m)
    d <- ape::as.DNAbin(t(apply(H, 1, function(r) ifelse(r == 1, "t", "a"))))
    ref <- pegas::tajima.test(d)
    f <- plasgenomicsutilsR:::.snp_freqs(H)
    h <- plasgenomicsutilsR:::.site_het(f$p, f$n)
    S <- sum(f$p > 0 & f$p < 1)
    D <- tajima_d(h, n)
    expect_equal(tajima_d_pvalue(D, n, S, "normal"), ref$Pval.normal, tolerance = 1e-10)
    expect_equal(tajima_d_pvalue(D, n, S, "beta"), ref$Pval.beta, tolerance = 1e-9)
  }
})

test_that("the p-value is two-sided and undefined where D is", {
  expect_equal(tajima_d_pvalue(0, n = 40, S = 20, method = "normal"), 1)
  # symmetric under the normal approximation
  expect_equal(tajima_d_pvalue(-2, 40, 20, "normal"), tajima_d_pvalue(2, 40, 20, "normal"))
  expect_lt(tajima_d_pvalue(-2.5, 40, 20), tajima_d_pvalue(-0.5, 40, 20))
  expect_true(is.na(tajima_d_pvalue(NA, 40, 20)))
  expect_true(is.na(tajima_d_pvalue(-1, n = 3, S = 20)))     # too few samples
  expect_true(is.na(tajima_d_pvalue(-1, n = 40, S = 0)))     # nothing segregating
})

test_that("pop_diversity carries the p-value and an empirical percentile", {
  G <- .fake_geno(n = 30, m = 60, spacing = 100)
  genes <- data.frame(name = paste0("g", 1:5), chrom = "Pf3D7_01_v3",
                      start = seq(1000, by = 1200, length.out = 5),
                      end = seq(2200, by = 1200, length.out = 5))
  d <- pop_diversity(G, by = "gene", genes = genes)
  expect_true(all(c("tajima_p", "tajima_percentile") %in% names(d)))
  ok <- is.finite(d$tajima_d)
  if (sum(ok) > 2) {
    expect_true(all(d$tajima_percentile[ok] >= 0 & d$tajima_percentile[ok] <= 1))
    # the percentile ranks D within the group, so the smallest D is the smallest percentile
    expect_equal(which.min(d$tajima_d[ok]), which.min(d$tajima_percentile[ok]))
  }
  expect_true(all(is.na(d$tajima_p[!ok])))
})

test_that("step slides the window and defaults to abutting", {
  G <- .fake_geno(m = 60, spacing = 250)          # SNPs every 250 bp from 1000
  abut <- pop_diversity(G, by = "window", window = 5000)
  slide <- pop_diversity(G, by = "window", window = 5000, step = 2500)

  expect_equal(abut$start, pop_diversity(G, by = "window", window = 5000,
                                         step = 5000)$start)   # step defaults to window
  expect_true(all(diff(abut$start) == 5000))
  expect_true(all(diff(slide$start) == 2500))
  expect_true(all(slide$end - slide$start == 5000))            # width is still `window`
  expect_gt(nrow(slide), nrow(abut))

  # consecutive sliding windows overlap by window - step, so SNPs are counted twice
  expect_gt(sum(slide$n_snps), sum(abut$n_snps))
  expect_lt(slide$start[2], abut$end[1])
})
