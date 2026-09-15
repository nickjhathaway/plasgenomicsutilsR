# The lossless `allele_set` encoding: one panel that keeps a mixed call AND a multiallelic
# identity per cell, and derives the allele_index / dosage views every analysis reads.

# A small VCF written to a temp file, one record per (pos, alt, gts) triple.
.allele_set_vcf <- function(recs, samp, chrom = "Pf3D7_07_v3") {
  hdr <- c("##fileformat=VCFv4.2",
           sprintf("##contig=<ID=%s,length=1445207>", chrom),
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
           paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",
                  paste(samp, collapse = "\t")))
  row <- function(r) paste(c(chrom, format(r$pos, scientific = FALSE, trim = TRUE), ".", "A",
                             r$alt, ".", "PASS", ".", "GT", r$gts), collapse = "\t")
  v <- tempfile(fileext = ".vcf")
  writeLines(c(hdr, vapply(recs, row, character(1))), v)
  v
}

test_that(".encode_allele_sets names a diploid het and a multiallelic mix as their allele sets", {
  # one triallelic column (A > C,G): the states span a pure ref, both pure alternates, a
  # ref+alt1 het and an alt1+alt2 mix -- exactly what a dosage or an allele index cannot keep.
  # ploidy x sample x variant, allele indices (0 = ref, 1 = C, 2 = G), NA = missing
  gt <- array(NA_integer_, dim = c(2, 6, 1))
  gt[, 1, 1] <- c(0L, 0L)   # {0}         reference
  gt[, 2, 1] <- c(1L, 1L)   # {1}         alternate 1
  gt[, 3, 1] <- c(2L, 2L)   # {2}         alternate 2
  gt[, 4, 1] <- c(0L, 1L)   # {0,1}       reference + alternate 1
  gt[, 5, 1] <- c(1L, 2L)   # {1,2}       alternate 1 + alternate 2
  gt[, 6, 1] <- c(NA, NA)   # missing
  es <- plasgenomicsutilsR:::.encode_allele_sets(gt, list(c("C", "G")), star = "missing")

  # the code is a 0-based index into that column's own state list
  states <- es$levels[[1]][es$codes[, 1] + 1L]
  expect_equal(states[1:5],
               c("reference", "alternate 1", "alternate 2",
                 "reference + alternate 1", "alternate 1 + alternate 2"))
  expect_true(is.na(states[6]))                       # missing stays missing
  # only the states that occur are listed (a triallelic site has seven possible)
  expect_length(es$levels[[1]], 5L)
  # the integer decomposition lines up with the names, in the same order
  expect_equal(es$sets[[1]][[match("reference + alternate 1", es$levels[[1]])]], c(0L, 1L))
  expect_equal(es$sets[[1]][[match("alternate 1 + alternate 2", es$levels[[1]])]], c(1L, 2L))
})

test_that(".encode_allele_sets drops the `*` allele under star = 'missing'", {
  # A > *,T : `*` is index 1, T is index 2. A cell that is only `*` becomes missing; a T/*
  # mix drops the `*` and reads as the singleton T -- which is more than an allele_index panel
  # keeps (there a het is missing outright).
  gt <- array(NA_integer_, dim = c(2, 4, 1))
  gt[, 1, 1] <- c(0L, 0L)   # {0}          reference
  gt[, 2, 1] <- c(1L, 1L)   # {*}   -> {}  missing
  gt[, 3, 1] <- c(2L, 2L)   # {T}          alternate 2 (index kept, matching the index panel)
  gt[, 4, 1] <- c(1L, 2L)   # {*,T} -> {T} alternate 2
  es <- plasgenomicsutilsR:::.encode_allele_sets(gt, list(c("*", "T")), star = "missing")
  states <- es$levels[[1]][es$codes[, 1] + 1L]
  expect_equal(states, c("reference", NA, "alternate 2", "alternate 2"))
  expect_equal(es$n_star, 2L)                          # two cells touched a `*` (pure, and mixed)
  expect_equal(es$star_col, 1L)

  # star = "allele" keeps `*` as a state instead
  es2 <- plasgenomicsutilsR:::.encode_allele_sets(gt, list(c("*", "T")), star = "allele")
  st2 <- es2$levels[[1]][es2$codes[, 1] + 1L]
  expect_equal(st2[2], "alternate 1")                  # the pure `*` is its own state now
  expect_equal(es2$n_star, 0L)
})

test_that("load_genotypes(encoding='allele_set') carries state levels and derives the index", {
  skip_if_not_installed("SeqArray")
  samp <- sprintf("s%02d", 1:12)
  recs <- list(
    list(pos = 429000, alt = "G",   gts = rep(c("0/0", "1/1", "0/1"), length.out = 12)),
    list(pos = 429500, alt = "C,G", gts = rep(c("0/0", "1/1", "2/2", "1/2"), length.out = 12)),
    list(pos = 430000, alt = "T",   gts = rep(c("0/0", "1/1"), length.out = 12)))
  vcf <- .allele_set_vcf(recs, samp)

  set <- suppressMessages(load_genotypes(vcf, gds = tempfile(fileext = ".gds"),
                                         variants = "all", encoding = "allele_set"))
  idx <- suppressMessages(load_genotypes(vcf, gds = tempfile(fileext = ".gds"),
                                         variants = "all", encoding = "allele_index"))

  expect_equal(set$encoding, "allele_set")
  # the list carries the per-column states and their integer decomposition
  expect_named(set$state_levels, colnames(set$genotype))
  expect_true("mixed" %in% set$state_levels[["Pf3D7_07_v3:428999"]])                  # biallelic het
  expect_true("alternate 1 + alternate 2" %in% set$state_levels[["Pf3D7_07_v3:429499"]])

  # a PopStructure DERIVES an allele-index view from the set panel, and it is column-for-column
  # identical to a directly loaded allele_index panel: singletons keep their allele, mixes -> NA
  ps <- suppressMessages(PopStructure$new(set,
          meta = data.frame(sample = set$sample.id, region = "X")))
  di <- suppressMessages(ps$genotype(needs = "allele_index"))
  direct <- idx$genotype[rownames(di), colnames(di), drop = FALSE]
  expect_equal(di, direct)

  # a dosage view derives too: the multiallelic column is dropped (a dosage cannot hold it), and
  # it equals a DIRECT dosage load column for column -- a biallelic mixed call keeps its
  # intermediate 1 rather than being forced to missing the way an allele index would.
  dose <- suppressMessages(ps$genotype(needs = "dosage"))
  expect_false("Pf3D7_07_v3:429499" %in% colnames(dose))
  expect_setequal(colnames(dose), c("Pf3D7_07_v3:428999", "Pf3D7_07_v3:429999"))
  direct_dose <- suppressWarnings(suppressMessages(
    load_genotypes(vcf, gds = tempfile(fileext = ".gds"),
                   variants = "all", encoding = "dosage", prune = FALSE)))
  expect_equal(dose,
               direct_dose$genotype[rownames(dose), colnames(dose), drop = FALSE])
  # the mixed (0/1) samples at the biallelic het site are dosage 1, not NA
  mixed_rows <- samp[seq(3, 12, by = 3)]                 # every third sample was 0/1
  expect_true(all(dose[mixed_rows, "Pf3D7_07_v3:428999"] == 1L))
})

test_that("state_levels() and the derived panels are reachable through the accessors", {
  skip_if_not_installed("SeqArray")
  samp <- sprintf("s%02d", 1:9)
  recs <- list(list(pos = 500000, alt = "C,G",
                    gts = rep(c("0/0", "1/1", "2/2", "0/1", "1/2"), length.out = 9)))
  set <- suppressMessages(load_genotypes(.allele_set_vcf(recs, samp),
          gds = tempfile(fileext = ".gds"), variants = "all", encoding = "allele_set"))
  ps <- suppressMessages(PopStructure$new(set, meta = data.frame(sample = samp, region = "X")))

  sl <- ps$state_levels()
  expect_named(sl, colnames(ps$genotype()))
  expect_true(all(c("reference + alternate 1", "alternate 1 + alternate 2") %in% sl[[1]]))
  # a dosage panel does not carry sets, so the accessor is NULL there
  bare <- PopStructure$new(matrix(0L, 3, 2, dimnames = list(c("a","b","c"), c("1:1","1:2"))))
  expect_null(bare$state_levels())
})

test_that("plot_region_haplotypes auto-prefers the set panel and draws mixed + multiallelic cells", {
  skip_if_not_installed("SeqArray")
  skip_if_not_installed("ggplot2")
  samp <- sprintf("s%02d", 1:12)
  recs <- list(
    list(pos = 429000, alt = "G",   gts = rep(c("0/0", "1/1", "0/1"), length.out = 12)),  # het
    list(pos = 429500, alt = "C,G", gts = rep(c("0/0", "1/1", "2/2", "1/2"), length.out = 12)),
    list(pos = 430000, alt = "T",   gts = rep(c("0/0", "1/1"), length.out = 12)))
  set <- suppressMessages(load_genotypes(.allele_set_vcf(recs, samp),
          gds = tempfile(fileext = ".gds"), variants = "all", encoding = "allele_set"))
  ps <- suppressMessages(PopStructure$new(set, meta = data.frame(sample = samp, region = "X")))

  # the set panel is the one the plot reaches for, above index/full
  expect_equal(plasgenomicsutilsR:::.set_panel_of(ps), ps$panels()[1])

  p <- suppressMessages(plot_region_haplotypes(ps, "7", genes = PF_EXAMPLE_DRUG_GENES,
                                               pad = 5000))
  hm <- NULL
  for (i in seq_len(8)) {
    q <- tryCatch(p[[i]], error = function(e) NULL)
    if (!is.null(q) && is.data.frame(q$data) && "call" %in% names(q$data)) { hm <- q; break }
  }
  expect_false(is.null(hm))
  # the biallelic het lands on the shared green "mixed"
  bi <- unique(as.character(hm$data$call[hm$data$snp_id == "Pf3D7_07_v3:428999"]))
  expect_true("mixed" %in% bi)
  # the multiallelic site keeps every allele apart AND shows the polyclonal mix as its own state
  mu <- unique(as.character(hm$data$call[hm$data$snp_id == "Pf3D7_07_v3:429499"]))
  expect_true(all(c("reference", "alternate 1", "alternate 2",
                    "alternate 1 + alternate 2") %in% mu))
})

test_that("parasite_haplotypes takes an allele_set input and gives the same haplotypes as allele_index", {
  skip_if_not_installed("SeqArray")
  samp <- sprintf("s%02d", 1:12)
  recs <- list(
    list(pos = 429000, alt = "G",   gts = rep(c("0/0", "1/1", "0/1"), length.out = 12)),  # het
    list(pos = 429500, alt = "C,G", gts = rep(c("0/0", "1/1", "2/2", "1/2"), length.out = 12)),  # triallelic
    list(pos = 430000, alt = "T",   gts = rep(c("0/0", "1/1"), length.out = 12)))
  vcf <- .allele_set_vcf(recs, samp)
  L <- function(enc) suppressWarnings(suppressMessages(
    load_genotypes(vcf, gds = tempfile(fileext = ".gds"), prune = FALSE,
                   variants = "all", encoding = enc)))
  idx <- L("allele_index")
  set <- L("allele_set")
  meta <- data.frame(sample = samp, region = "X", stringsAsFactors = FALSE)

  # generous missingness thresholds so the mixed-heavy triallelic site is not filtered out --
  # this test is about representation, not the missingness filter (which drops it identically
  # under either encoding)
  ph <- function(g) suppressMessages(parasite_haplotypes(
    g, meta = if (inherits(g, "PopStructure")) NULL else meta, maf = 0,
    max_snp_missing = 0.9, max_sample_missing = 0.9))

  # the raw-list path the qmd uses: an index list resolves alleles = "index" today
  h_idx <- ph(idx)
  expect_equal(h_idx$alleles, "index")
  # a bare allele_set list now reaches the same path (auto -> index, index derived in
  # .coerce_geno) instead of erroring, and yields an identical haplotype matrix
  h_set <- ph(set)
  expect_equal(h_set$alleles, "index")
  expect_equal(h_set$hap[rownames(h_idx$hap), colnames(h_idx$hap)], h_idx$hap)

  # and the triallelic codon site is carried as distinct alleles (0/1/2), not collapsed
  tri <- "Pf3D7_07_v3:429499"
  expect_true(tri %in% colnames(h_set$hap))
  expect_setequal(sort(unique(h_set$hap[, tri])), c(0L, 1L, 2L))

  # a PopStructure holding an allele_set panel also auto-resolves to the index reading now
  ps <- suppressMessages(PopStructure$new(set, meta = meta))
  h_ps <- ph(ps)
  expect_equal(h_ps$alleles, "index")
  expect_equal(h_ps$hap[rownames(h_idx$hap), colnames(h_idx$hap)], h_idx$hap)
})

# A multi-contig VCF: a focal chromosome carrying a triallelic codon plus biallelic SNPs, and
# background contigs for the null windows. Returns the vcf path.
.ibs_multi_vcf <- function(samp, seed = 1) {
  set.seed(seed)
  contigs <- c("Pf3D7_07_v3", "Pf3D7_01_v3", "Pf3D7_02_v3")
  hdr <- c("##fileformat=VCFv4.2",
           sprintf("##contig=<ID=%s,length=1500000>", contigs),
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
           paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",
                  paste(samp, collapse = "\t")))
  ns <- length(samp)
  gt_bi <- function() {
    v <- sample(c("0/0", "1/1", "0/1", "./."), ns, replace = TRUE, prob = c(.55, .3, .1, .05))
    v
  }
  gt_tri <- function() sample(c("0/0", "1/1", "2/2", "1/2", "./."), ns, replace = TRUE,
                              prob = c(.4, .3, .15, .1, .05))
  rows <- character(0)
  for (ci in seq_along(contigs)) {
    ch <- contigs[ci]
    for (k in seq_len(15)) {
      pos <- 1000 * k
      # one triallelic codon on the focal chromosome, everything else biallelic
      if (ci == 1L && k == 8L) { alt <- "C,G"; gts <- gt_tri() }
      else                     { alt <- "T";   gts <- gt_bi() }
      rows <- c(rows, paste(c(ch, pos, ".", "A", alt, ".", "PASS", ".", "GT", gts),
                            collapse = "\t"))
    }
  }
  v <- tempfile(fileext = ".vcf"); writeLines(c(hdr, rows), v); v
}

test_that(".as_ibs_panel derives the same index matrix from allele_set as from allele_index", {
  skip_if_not_installed("SeqArray")
  samp <- sprintf("s%02d", 1:24)
  vcf <- .ibs_multi_vcf(samp)
  idx <- suppressWarnings(suppressMessages(load_genotypes(vcf, gds = tempfile(fileext = ".gds"),
           prune = FALSE, variants = "all", encoding = "allele_index")))
  set <- suppressWarnings(suppressMessages(load_genotypes(vcf, gds = tempfile(fileext = ".gds"),
           prune = FALSE, variants = "all", encoding = "allele_set")))

  pi <- plasgenomicsutilsR:::.as_ibs_panel(idx)
  ps <- plasgenomicsutilsR:::.as_ibs_panel(set)
  expect_identical(pi$Gt, ps$Gt)          # SNPs x samples, integer -- the matrix IBS reads
  expect_identical(pi$chr, ps$chr)
  expect_identical(pi$pos, ps$pos)
  # the triallelic codon is present and carries distinct alleles, not collapsed
  tri <- "Pf3D7_07_v3:7999"
  expect_true(tri %in% rownames(ps$Gt))
  expect_true(any(ps$Gt[tri, ] == 2L, na.rm = TRUE))
})

test_that("locus_ibs_by_allele gives identical results from an allele_set and an allele_index list", {
  skip_if_not_installed("SeqArray")
  samp <- sprintf("s%02d", 1:24)
  vcf <- .ibs_multi_vcf(samp)
  meta <- data.frame(sample = samp,
                     region = rep(c("A", "B"), each = 12),
                     car = rep(c("Present", "Absent"), length.out = 24),
                     stringsAsFactors = FALSE)
  locus <- data.frame(name = "codon", chr = "Pf3D7_07_v3", start = 7500, end = 8500)
  idx <- suppressWarnings(suppressMessages(load_genotypes(vcf, gds = tempfile(fileext = ".gds"),
           prune = FALSE, variants = "all", encoding = "allele_index")))
  set <- suppressWarnings(suppressMessages(load_genotypes(vcf, gds = tempfile(fileext = ".gds"),
           prune = FALSE, variants = "all", encoding = "allele_set")))

  run <- function(g) suppressWarnings(locus_ibs_by_allele(
    g, locus, allele = "car", carrier = "Present", reference = "Absent",
    group_a = "A", group_b = "B", group = "region", meta = meta,
    n_snps = 4L, min_sites = 2L, span_tol = NULL))
  ri <- run(idx); rs <- run(set)
  expect_equal(as.data.frame(rs), as.data.frame(ri))

  # an allele_set list must NOT trip the dosage warning (it is not dosage-encoded), and it must
  # not be read as codes -- reading it correctly is what the equality above proves
  expect_no_warning(plasgenomicsutilsR:::.as_ibs_panel(set))
  # a set panel with the decomposition stripped is refused rather than silently misread
  set_broken <- set; set_broken$state_sets <- NULL
  expect_error(plasgenomicsutilsR:::.as_ibs_panel(set_broken), "state-set decomposition")
})

test_that(".encode_allele_sets is vectorised, so a genome-wide diploid build stays fast", {
  # the naive form is a variants x samples double loop; the shipped one vectorises over samples
  # within each column, so 400 samples x 5000 triallelic records must not take seconds.
  set.seed(1)
  ns <- 400L; nv <- 5000L
  gt <- array(sample(c(0:2, NA), 2 * ns * nv, replace = TRUE,
                     prob = c(0.5, 0.25, 0.2, 0.05)),
              dim = c(2, ns, nv))
  alt <- rep(list(c("C", "G")), nv)
  t <- system.time(es <- plasgenomicsutilsR:::.encode_allele_sets(gt, alt, star = "missing"))
  expect_equal(dim(es$codes), c(ns, nv))
  expect_length(es$levels, nv)
  expect_lt(t[["elapsed"]], 20)                      # generous, but a per-cell loop blows past it
})
