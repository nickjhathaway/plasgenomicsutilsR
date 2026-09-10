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

test_that("freqbin defaults to one bin for an unpolarized scan", {
  rf <- plasgenomicsutilsR:::.resolve_freqbin
  # unpolarized has no derived allele to bin by, so a single bin is the honest choice
  expect_equal(rf(NULL, polarized = FALSE), 1)
  expect_equal(rf(NULL, polarized = TRUE), 0.05)
  # an explicit choice is respected, but binning an unpolarized scan is called out
  expect_warning(rf(0.05, polarized = FALSE), "not the derived-allele frequency")
  expect_silent(rf(0.05, polarized = TRUE))
  expect_silent(rf(1, polarized = FALSE))
})

test_that("subset_haplotypes keeps samples, metadata groups, or both", {
  ps <- example_pop_structure(umap = FALSE)
  hap <- parasite_haplotypes(ps, maf = 0.05)
  n_all <- nrow(hap$hap)

  one <- subset_haplotypes(hap, country = "Ghana")
  expect_lt(nrow(one$hap), n_all)
  expect_setequal(unique(as.character(one$meta$country)), "Ghana")
  # the SNP panel is deliberately untouched, so two subsets stay comparable
  expect_identical(colnames(one$hap), colnames(hap$hap))
  expect_identical(one$map, hap$map)

  # several values for one column
  both <- subset_haplotypes(hap, country = c("Ghana", "Cambodia"))
  expect_equal(nrow(both$hap), n_all)

  keep <- head(rownames(hap$hap), 8)
  expect_setequal(rownames(subset_haplotypes(hap, samples = keep)$hap), keep)
  # samples and metadata together intersect
  gh <- rownames(one$hap)
  mix <- subset_haplotypes(hap, samples = c(gh[1:3], setdiff(rownames(hap$hap), gh)[1:3]),
                           country = "Ghana")
  expect_setequal(rownames(mix$hap), gh[1:3])

  # the restriction is recorded, since it changes what every scan off the object means
  expect_match(paste(capture.output(print(one)), collapse = " "), "subset")
  expect_equal(one$subset$from, n_all)
})

test_that("subset_haplotypes refuses filters that cannot mean anything", {
  ps <- example_pop_structure(umap = FALSE)
  hap <- parasite_haplotypes(ps, maf = 0.05)
  # a typo'd level would otherwise come back as "no data" rather than as a mistake
  expect_error(subset_haplotypes(hap, country = "Narnia"), "no sample has country = Narnia")
  expect_error(subset_haplotypes(hap, nope = "x"), "not a metadata column: nope")
  expect_error(subset_haplotypes(hap, samples = "nobody"), "none of `samples`")
  expect_warning(subset_haplotypes(hap, samples = c(rownames(hap$hap)[1], "nobody")),
                 "not in these haplotypes")
  expect_error(subset_haplotypes(hap$hap), "must be a parasite_haplotypes")
  # metadata-free haplotypes cannot match `...`
  bare <- hap; bare$meta <- NULL
  expect_error(subset_haplotypes(bare, country = "Ghana"), "no metadata")
})

test_that("a subset feeds the scans and the EHH curve", {
  skip_if_not_installed("rehh")
  skip_if_not_installed("ggplot2")
  ps <- example_pop_structure(umap = FALSE)
  hap <- parasite_haplotypes(ps, maf = 0.05)
  gh <- subset_haplotypes(hap, country = "Ghana")
  # a scan on the subset sees only those haplotypes
  scan <- run_ihs(gh)
  expect_true(nrow(scan) > 0)
  expect_s3_class(plot_ehh(gh, hap$map$snp_id[10], span = 5e5), "ggplot")
})

test_that("maxgap stops the integration at a hole instead of crossing it", {
  skip_if_not_installed("rehh")
  ps <- example_pop_structure(umap = FALSE)
  hap <- parasite_haplotypes(ps, maf = 0.05)
  # a SNP-free stretch has nothing to break the haplotype, so EHH runs flat across it and
  # the integral accrues the width of the hole rather than anything about the haplotypes
  gaps <- unlist(lapply(split(hap$map$pos, hap$map$chr), function(p) diff(sort(p))))
  free <- suppressWarnings(run_ihs(hap, group = "country"))

  # hold the border rule fixed so this measures `maxgap` alone: a gap rule that no gap in
  # the data exceeds has to leave the scan exactly where it was. It has to clear the gaps
  # *within a group*, which are wider than the map's -- a SNP monomorphic in one group is
  # dropped from that group's scan, and the hole it leaves is real for that scan.
  wide <- suppressWarnings(run_ihs(hap, group = "country", maxgap = 1e9,
                                   discard_at_border = FALSE))
  expect_equal(wide$ihs, free$ihs)
  # and one under most of the spacing has to change it, by cutting the integrals short
  capped <- suppressWarnings(run_ihs(hap, group = "country",
                                     maxgap = stats::quantile(gaps, 0.25),
                                     discard_at_border = FALSE))
  expect_false(isTRUE(all.equal(capped$ihs, free$ihs)))
})

test_that("a maxgap that leaves nothing to score says so rather than returning NAs", {
  skip_if_not_installed("rehh")
  ps <- example_pop_structure(umap = FALSE)
  hap <- parasite_haplotypes(ps, maf = 0.05)
  gaps <- unlist(lapply(split(hap$map$pos, hap$map$chr), function(p) diff(sort(p))))
  # setting `maxgap` turns the border rule on, and on markers this sparse EHH never decays
  # before the data runs out -- every marker is then at a border and scores NA
  expect_warning(out <- run_ihs(hap, group = "country",
                                maxgap = stats::quantile(gaps, 0.25)),
                 "reached a border")
  expect_true(all(is.na(out$ihs)))
})

test_that("discard_at_border follows maxgap unless it is set outright", {
  expect_false(.resolve_border(NULL, NA))       # no gap rule: keep the telomeric markers
  expect_true(.resolve_border(NULL, 20000))     # a gap rule makes the border a data hole
  expect_true(.resolve_border(TRUE, NA))
  expect_false(.resolve_border(FALSE, 20000))
})

test_that("ihs_windows counts the extreme SNPs in each window", {
  scan <- data.frame(
    group = factor(rep(c("b", "a"), each = 8), levels = c("b", "a")),
    chr = "Pf3D7_01_v3",
    pos = rep(c(10, 20, 30, 40, 1010, 1020, 1030, 1040), 2),
    ihs = c(3, 3, -3, 0.1,   0.1, 0.2, 0.3, 0.4,      # group b: 3/4 then 0/4
            0, 0, 0, 0,      3, -3, 0.1, 0.2))        # group a: 0/4 then 2/4
  w <- ihs_windows(scan, window = 1000, min_snps = 4)
  expect_equal(nrow(w), 4)
  expect_equal(w$n_snps, rep(4L, 4))
  expect_equal(w$frac_extreme[w$group == "b"], c(0.75, 0))
  expect_equal(w$frac_extreme[w$group == "a"], c(0, 0.5))
  # the magnitude is what counts, so a negative iHS is extreme too
  expect_equal(w$n_extreme[w$group == "a"], c(0L, 2L))
  expect_equal(w$max_abs[w$group == "b"], c(3, 0.4))
  # the window midpoint plots inside the data, never past the last SNP
  expect_true(all(w$pos <= max(scan$pos)))
  expect_equal(levels(w$group), c("b", "a"))   # the scan's order, not the alphabet
})

test_that("ihs_windows drops thin windows rather than letting them read 100%", {
  scan <- data.frame(chr = "Pf3D7_01_v3", pos = c(10, 20, 5010),
                     ihs = c(0.1, 0.2, 9))
  # the lone SNP in the second window is extreme, which would be a 100% window
  expect_equal(nrow(ihs_windows(scan, window = 1000, min_snps = 1)), 2)
  w <- ihs_windows(scan, window = 1000, min_snps = 2)
  expect_equal(nrow(w), 1)
  expect_equal(w$frac_extreme, 0)
  expect_error(ihs_windows(scan, window = 1000, min_snps = 5), "no window holds")
})

test_that("ihs_windows slides when given a step, and reports what it summarised", {
  scan <- data.frame(chr = "Pf3D7_01_v3", pos = seq(0, 900, by = 100),
                     ihs = c(rep(0.1, 5), rep(3, 5)))
  tiled <- ihs_windows(scan, window = 500, min_snps = 1)
  slid <- ihs_windows(scan, window = 500, step = 250, min_snps = 1)
  expect_gt(nrow(slid), nrow(tiled))
  expect_true(all(slid$end - slid$start == 500))
  # a window straddling the switch has to land between the two flat halves
  expect_true(any(slid$frac_extreme > 0 & slid$frac_extreme < 1))
  expect_equal(sum(tiled$n_extreme), 5L)
})

test_that("ihs_windows says which column it cannot find", {
  scan <- data.frame(chr = "Pf3D7_01_v3", pos = 1:10, ihs = 0)
  expect_error(ihs_windows(scan, metric = "rsb"), "no 'rsb' column")
  expect_error(ihs_windows(scan[, c("chr", "ihs")]), "no 'pos' column")
  expect_error(ihs_windows(scan, window = -1), "positive width")
  expect_error(ihs_windows(transform(scan, ihs = NA_real_)), "no finite")
})

test_that("maf_bands standardises within frequency bands rather than over all of them", {
  skip_if_not_installed("rehh")
  ps <- example_pop_structure(umap = FALSE)
  hap <- parasite_haplotypes(ps, maf = 0.05)
  one <- suppressWarnings(run_ihs(hap, group = "country"))
  banded <- suppressWarnings(run_ihs(hap, group = "country", maf_bands = 4))

  # same SNPs, different scores
  expect_equal(nrow(banded), nrow(one))
  expect_equal(banded$pos, one$pos)
  expect_false(isTRUE(all.equal(banded$ihs, one$ihs)))

  # a standardised score is centred and scaled -- within each band now, not just overall
  ok <- is.finite(banded$ihs)
  expect_lt(abs(mean(banded$ihs[ok])), 0.15)
  expect_lt(abs(stats::sd(banded$ihs[ok]) - 1), 0.15)

  # and the point of it: the spread should no longer track the minor-allele frequency
  spread <- function(d) {
    d <- d[is.finite(d$ihs) & is.finite(d$freq_minor), ]
    b <- cut(d$freq_minor, stats::quantile(d$freq_minor, c(0, 0.5, 1)),
             include.lowest = TRUE)
    s <- tapply(abs(d$ihs), b, mean)
    unname(s[1] / s[2])            # rare-allele half over common-allele half
  }
  expect_lt(abs(spread(banded) - 1), abs(spread(one) - 1))
})

test_that("a band too thin to standardise in says so", {
  raw <- data.frame(CHR = "Pf3D7_01_v3", POSITION = seq_len(40), FREQ_MAJ = 0.8,
                    FREQ_MIN = c(rep(0.2, 36), rep(0.45, 4)),
                    UNIHS = stats::rnorm(40))
  expect_warning(.band_standardise(raw, raw$FREQ_MIN, 10), "fewer than 10 markers")
})


test_that("the unstandardised log ratio comes back, so the ratio is recoverable", {
  set.seed(11)
  raw <- data.frame(CHR = "Pf3D7_01_v3", POSITION = seq_len(60), FREQ_MAJ = 0.75,
                    FREQ_MIN = runif(60, 0.05, 0.5), UNIHS = stats::rnorm(60, -0.4, 0.9))
  b <- .band_standardise(raw, raw$FREQ_MIN, 3)

  expect_true("UNIHS" %in% names(b))
  expect_equal(b$UNIHS, raw$UNIHS)                       # passed through untouched
  # and it is the thing `ihs` was standardised from: undo the z-score band by band and the
  # raw values come back, which is exactly what makes exp(unihs) the EHH ratio
  br <- unique(stats::quantile(raw$FREQ_MIN, seq(0, 1, length.out = 4)))
  bands <- cut(raw$FREQ_MIN, br, include.lowest = TRUE)
  m <- tapply(raw$UNIHS, bands, mean)
  s <- tapply(raw$UNIHS, bands, stats::sd)
  expect_equal(b$IHS, as.vector((raw$UNIHS - m[bands]) / s[bands]))
  expect_equal(unname(exp(b$UNIHS)), exp(raw$UNIHS))

  # a z-score of 0 is the band mean, which is not a ratio of 1 unless the mean happens to be
  expect_false(isTRUE(all.equal(exp(mean(raw$UNIHS)), 1)))
})

test_that("an undefined log ratio stays NA through standardisation", {
  # rehh returns NA for UNIHS when it cannot integrate EHH for one of the two alleles; that
  # must travel through as NA rather than becoming a score
  raw <- data.frame(CHR = "Pf3D7_01_v3", POSITION = seq_len(40), FREQ_MAJ = 0.8,
                    FREQ_MIN = runif(40, 0.05, 0.5),
                    UNIHS = c(NA_real_, stats::rnorm(39)))
  b <- .band_standardise(raw, raw$FREQ_MIN, 2)
  expect_true(is.na(b$UNIHS[1]))
  expect_true(is.na(b$IHS[1]))
  expect_true(is.na(b$LOGPVALUE[1]))
  expect_equal(sum(is.na(b$IHS)), 1L)                    # only that one
})

test_that("run_ihs says when a marker had more than two alleles", {
  skip_if_not_installed("rehh")
  mk <- function(codes, probs, seed = 3) {
    set.seed(seed); n <- 70; m <- 30
    H <- matrix(stats::rbinom(n * m, 1, 0.4), n, m)
    H[, 15] <- sample(codes, n, TRUE, prob = probs)
    map <- data.frame(chr = "Pf3D7_01_v3", pos = seq(1000, by = 500, length.out = m),
                      snp_id = paste0("Pf3D7_01_v3:",
                                      seq(1000, by = 500, length.out = m)),
                      stringsAsFactors = FALSE)
    colnames(H) <- map$snp_id; rownames(H) <- paste0("s", seq_len(n))
    structure(list(hap = H, map = map,
                   meta = data.frame(sample = rownames(H), grp = "all",
                                     stringsAsFactors = FALSE),
                   filtering = list()), class = "parasite_haplotypes")
  }
  # a biallelic scan is silent, so the warning cannot be background noise
  expect_no_warning(run_ihs(mk(0:1, c(.6, .4)), min_maf = 0.02))

  # rehh keeps the two commonest alleles and drops the rest without saying so; this is the
  # saying so, and it reports how much of the group went missing
  expect_warning(run_ihs(mk(0:2, c(.5, .35, .15)), min_maf = 0.02),
                 "more than two alleles")
  expect_warning(run_ihs(mk(0:2, c(.5, .35, .15)), min_maf = 0.02),
                 "Pf3D7_01_v3:8000 \\(20% of haplotypes excluded\\)")
  # four alleles drop more than three do
  expect_warning(run_ihs(mk(0:3, c(.4, .3, .2, .1)), min_maf = 0.02),
                 "39% of haplotypes excluded")

  # the count is of markers, not of marker-by-group pairs
  h <- mk(0:2, c(.5, .35, .15))
  h$meta$grp <- rep(c("a", "b"), length.out = nrow(h$hap))
  w <- tryCatch(run_ihs(h, group = "grp", min_maf = 0.02),
                warning = function(e) conditionMessage(e))
  expect_match(w, "^1 marker\\(s\\)")

  # the detector itself, away from rehh
  d <- plasgenomicsutilsR:::.multiallelic_drop(mk(0:2, c(.5, .35, .15)),
                                               seq_len(nrow(mk(0:2, c(.5, .35, .15))$hap)))
  expect_equal(nrow(d), 1L)
  expect_equal(d$n_alleles, 3L)
  expect_null(plasgenomicsutilsR:::.multiallelic_drop(
    mk(0:1, c(.6, .4)), seq_len(70)))
})
