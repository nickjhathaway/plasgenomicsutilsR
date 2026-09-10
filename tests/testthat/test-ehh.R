hap_for_ehh <- function() parasite_haplotypes(example_pop_structure(umap = FALSE), maf = 0.05)

ehh_panel <- function(p) {
  if (!inherits(p, "patchwork")) return(p)
  for (i in seq_len(6)) {
    q <- tryCatch(p[[i]], error = function(e) NULL)
    if (!is.null(q) && is.data.frame(q$data) && "ehh" %in% names(q$data)) return(q)
  }
  stop("no EHH panel")
}

test_that("plot_ehh draws one decay curve per allele at the focal SNP", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("rehh")
  hap <- hap_for_ehh()
  p <- plot_ehh(hap, hap$map$snp_id[10], span = 5e5)
  d <- ehh_panel(p)$data
  # the two states at the focal SNP, which is the mutant-vs-reference split without needing
  # the SNPs annotated
  expect_setequal(as.character(unique(d$allele)), c("reference", "alternate"))
  expect_true(all(d$ehh >= 0 & d$ehh <= 1))
  # EHH is 1 at the focal SNP itself and decays away from it
  focal <- hap$map$pos[10]
  expect_equal(max(d$ehh[d$pos == focal]), 1)
  far <- d$ehh[abs(d$pos - focal) > 3e5]
  if (length(far)) expect_lte(max(far), 1)
  # the gene track is opt-in, since `genes` is usually only there to resolve `focal`
  expect_s3_class(p, "ggplot")
  expect_false(inherits(p, "patchwork"))
})

test_that("the focal SNP can be named three ways", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("rehh")
  hap <- hap_for_ehh()
  by_id <- ehh_panel(plot_ehh(hap, hap$map$snp_id[10], span = 5e5))$data
  by_pos <- ehh_panel(plot_ehh(hap, hap$map$pos[10], span = 5e5))$data
  expect_equal(by_id$ehh, by_pos$ehh)

  # a gene holding several SNPs picks the most balanced one, saying which and how many it
  # chose between; the shortlist itself comes from ehh_candidates()
  g <- data.frame(name = "wide", chr = normalise_chr(hap$map$chr[10]),
                  start = min(hap$map$pos) - 1, end = max(hap$map$pos) + 1)
  expect_message(plot_ehh(hap, "wide", genes = g, span = 5e5), "holds [0-9]+ SNPs")
  expect_message(plot_ehh(hap, "wide", genes = g, span = 5e5), "ehh_candidates")

  expect_error(plot_ehh(hap, "nowhere"), "not a SNP in the haplotypes")
  expect_error(plot_ehh(hap, 1), "no SNP at position 1")
  expect_error(plot_ehh(hap, c("a", "b")), "must be a `chr:pos` id")
  # a gene with no genotyped SNP in it
  empty <- data.frame(name = "empty", chr = "1", start = 1, end = 2)
  expect_error(plot_ehh(hap, "empty", genes = empty), "no genotyped SNP inside")
})

test_that("group facets in the metadata's order and skips groups too small to mean anything", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("rehh")
  hap <- hap_for_ehh()
  p <- plot_ehh(hap, hap$map$snp_id[10], group = "country", span = 5e5, min_haplotypes = 5)
  d <- ehh_panel(p)$data
  if (is.factor(d$group)) {
    expect_identical(levels(d$group),
                     intersect(levels(.as_group_factor(hap$meta$country)),
                               as.character(unique(d$group))))
  }
  # a group with fewer haplotypes than asked for is dropped, with a reason
  expect_message(plot_ehh(hap, hap$map$snp_id[10], group = "country", span = 5e5,
                          min_haplotypes = 25),
                 "skipping|fewer than|not variable")
})

test_that("span crops the curves and takes two sides", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("rehh")
  hap <- hap_for_ehh()
  focal <- hap$map$pos[10]
  wide <- ehh_panel(plot_ehh(hap, hap$map$snp_id[10], span = 6e5))$data
  narrow <- ehh_panel(plot_ehh(hap, hap$map$snp_id[10], span = 2e5))$data
  expect_lt(diff(range(narrow$pos)), diff(range(wide$pos)))
  expect_true(all(abs(narrow$pos - focal) <= 2e5))
  # Asymmetric, as everywhere else in the package. Assert on the axis, not the data: a curve
  # stops where EHH falls below `limehh`, which can be well inside the window on either side.
  asym <- plot_ehh(hap, hap$map$snp_id[10], span = c(left = 5e4, right = 5e5))
  rng <- ggplot2::ggplot_build(ehh_panel(asym))$layout$panel_params[[1]]$x.range
  # more room on the right than the left -- the window is still clamped to the chromosome,
  # so assert the asymmetry rather than an absolute reach
  expect_gt(rng[2] - focal, focal - rng[1])
  expect_lt(focal - rng[1], 1e5)
  expect_true(all(ehh_panel(asym)$data$pos >= focal - 5e4))
})

test_that("the gene track can be asked for", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("rehh")
  skip_if_not_installed("patchwork")
  hap <- hap_for_ehh()
  g <- data.frame(name = "amark", chr = normalise_chr(hap$map$chr[10]),
                  start = hap$map$pos[10] - 500, end = hap$map$pos[10] + 500)
  p <- plot_ehh(hap, hap$map$snp_id[10], span = 5e5, genes = g, gene_track = TRUE)
  expect_s3_class(p, "patchwork")
  expect_true("amark" %in% p[[2]]$data$name)
})

test_that("the frequency note can be moved or turned off", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("rehh")
  hap <- hap_for_ehh()
  note <- function(...) {
    b <- ggplot2::ggplot_build(ehh_panel(plot_ehh(hap, hap$map$snp_id[10], span = 5e5, ...)))
    d <- b$data[vapply(b$data, function(z) "label" %in% names(z), logical(1))]
    if (!length(d)) return(NULL)
    d[[1]]
  }
  top <- note()
  # default: the top of the panel, where the curves have not decayed to and cannot cross it
  expect_false(is.null(top))
  expect_gt(top$y[1], 0.9)
  expect_equal(top$hjust[1], 0)
  expect_equal(top$vjust[1], 1)
  expect_match(top$label[1],
               "^n = [0-9]+; reference [0-9]+ \\([0-9]+%\\), alternate [0-9]+ \\([0-9]+%\\)$")
  # the counts are counted, not recovered from the share, so they have to reconstruct n --
  # and matching the percentages is what says the two alleles were not labelled the wrong
  # way round, since the count is derived independently of the frequency rehh reports
  n <- as.integer(sub("^n = ([0-9]+);.*", "\\1", top$label[1]))
  k <- as.integer(regmatches(top$label[1],
                             gregexpr("[0-9]+(?= \\()", top$label[1], perl = TRUE))[[1]])
  pct <- as.integer(regmatches(top$label[1],
                               gregexpr("[0-9]+(?=%)", top$label[1], perl = TRUE))[[1]])
  expect_length(k, 2)
  expect_equal(sum(k), n)
  expect_equal(round(100 * k / n), pct)

  bottom <- note(freq_position = "bottomleft")
  expect_lt(bottom$y[1], 0.1)
  right <- note(freq_position = "topright")
  expect_gt(right$x[1], top$x[1])
  expect_equal(right$hjust[1], 1)

  expect_null(note(show_freq = FALSE))
  expect_error(plot_ehh(hap, hap$map$snp_id[10], freq_position = "middle"),
               "should be one of")
})

test_that("ehh_candidates lists the shortlist with the chosen SNP on top", {
  testthat::skip_if_not_installed("rehh")
  ps <- example_pop_structure(umap = FALSE)
  hap <- parasite_haplotypes(ps, maf = 0.05)

  cand <- ehh_candidates(hap, "pfcrt", genes = PF_EXAMPLE_DRUG_GENES)
  expect_true(nrow(cand) > 1)
  expect_true(all(c("snp_id", "chr", "pos", "maf", "n_hap", "chosen") %in% names(cand)))
  expect_equal(sum(cand$chosen), 1L)
  expect_true(cand$chosen[1])                              # chosen first
  expect_equal(cand$maf[-1], sort(cand$maf[-1], decreasing = TRUE))   # then by maf
  expect_equal(cand$maf[1], max(cand$maf))                 # and it is the most balanced
})

test_that("the SNP ehh_candidates marks chosen is the one plot_ehh measures from", {
  testthat::skip_if_not_installed("rehh")
  testthat::skip_if_not_installed("ggplot2")
  ps <- example_pop_structure(umap = FALSE)
  hap <- parasite_haplotypes(ps, maf = 0.05)

  chosen <- ehh_candidates(hap, "pfcrt", genes = PF_EXAMPLE_DRUG_GENES)$snp_id[1]
  msg <- testthat::capture_messages(
    plot_ehh(hap, "pfcrt", genes = PF_EXAMPLE_DRUG_GENES, span = 30000))
  expect_true(any(grepl(chosen, msg, fixed = TRUE)))
})

test_that("naming a chr:pos leaves a single candidate, already chosen", {
  testthat::skip_if_not_installed("rehh")
  ps <- example_pop_structure(umap = FALSE)
  hap <- parasite_haplotypes(ps, maf = 0.05)

  one <- ehh_candidates(hap, hap$map$snp_id[1], genes = PF_EXAMPLE_DRUG_GENES)
  expect_equal(nrow(one), 1L)
  expect_true(one$chosen)
  expect_equal(one$snp_id, hap$map$snp_id[1])
})

test_that("per-group columns show a SNP that is balanced overall but flat within a group", {
  testthat::skip_if_not_installed("rehh")
  ps <- example_pop_structure(umap = FALSE)
  hap <- parasite_haplotypes(ps, maf = 0.05)

  cand <- ehh_candidates(hap, "pfcrt", group = "country", genes = PF_EXAMPLE_DRUG_GENES)
  expect_true(all(c("maf_Ghana", "maf_Cambodia", "n_groups_variable") %in% names(cand)))
  expect_true(all(cand$n_groups_variable <= 2))
  # the fixture's most balanced pfcrt SNP separates the two countries, so it is monomorphic
  # inside each -- which is the case that leaves plot_ehh(group=) with no curve to draw
  expect_equal(cand$n_groups_variable[cand$chosen], 0L)
  expect_true(any(cand$n_groups_variable > 0))     # ... and better choices do exist
  # a group MAF can never beat the pooled one
  expect_true(all(cand$maf_Ghana <= cand$maf + 1e-9, na.rm = TRUE))
})

test_that("ehh_candidates accepts a PopStructure and rejects anything else", {
  testthat::skip_if_not_installed("rehh")
  ps <- example_pop_structure(umap = FALSE)
  expect_s3_class(ehh_candidates(ps, "pfcrt", genes = PF_EXAMPLE_DRUG_GENES), "tbl_df")
  expect_error(ehh_candidates(list(a = 1), "pfcrt"), "parasite_haplotypes")
})

test_that("the title is set on the curves, not on the gene track underneath", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("rehh")
  skip_if_not_installed("patchwork")
  hap <- hap_for_ehh()
  snp <- hap$map$snp_id[10]
  args <- list(hap, snp, span = 5e5, genes = PF_EXAMPLE_DRUG_GENES, gene_track = TRUE)

  # NULL keeps the written title, a string replaces it, NA/FALSE drops it
  expect_equal(ehh_panel(do.call(plot_ehh, args))$labels$title, paste0("EHH around ", snp))
  custom <- do.call(plot_ehh, c(args, list(title = "EHH around A675V")))
  expect_equal(ehh_panel(custom)$labels$title, "EHH around A675V")
  expect_null(ehh_panel(do.call(plot_ehh, c(args, list(title = NA))))$labels$title)
  expect_null(ehh_panel(do.call(plot_ehh, c(args, list(title = FALSE))))$labels$title)
  expect_equal(ehh_panel(do.call(plot_ehh, c(args, list(subtitle = "north vs south"))))$labels$subtitle,
               "north vs south")

  # the point of the argument: `+ labs()` on the returned patchwork would land on the gene
  # track, so the title has to be on the panel that sits on top
  track <- custom[[2]]
  expect_null(track$labels$title)
  outer <- patchwork::patchworkGrob(custom)
  cells <- outer$layout[outer$layout$name %in% c("panel-1", "panel-2"), ]
  expect_lt(cells$t[cells$name == "panel-1"], cells$t[cells$name == "panel-2"])
  g <- ggplot2::ggplotGrob(custom[[1]])
  expect_lt(g$layout$t[g$layout$name == "title"],
            min(g$layout$t[grepl("^panel", g$layout$name)]))

  # and without a track it is a plain ggplot carrying the same title
  plain <- plot_ehh(hap, snp, span = 5e5, title = "plain")
  expect_false(inherits(plain, "patchwork"))
  expect_equal(plain$labels$title, "plain")
})


# A focal marker with more than two alleles, which ALT dosage cannot express.
.multi_hap <- function(counts, m = 41, seed = 7) {
  set.seed(seed)
  n <- sum(counts)
  G <- matrix(stats::rbinom(n * m, 1, 0.35), n, m)
  G[, (m + 1) %/% 2] <- rep(seq_along(counts) - 1L, times = counts)
  rownames(G) <- sprintf("s%02d", seq_len(n))
  colnames(G) <- paste0("c1:", seq(1000, by = 500, length.out = m))
  parasite_haplotypes(G, maf = 0.02, alleles = "index")
}

test_that("a focal marker with three alleles gets three curves, correctly labelled", {
  skip_if_not_installed("rehh")
  skip_if_not_installed("ggplot2")
  hap <- .multi_hap(c(30, 10, 20))
  expect_setequal(unique(as.vector(hap$hap)), c(0L, 1L, 2L))

  p <- plot_ehh(hap, "c1:11000", span = 12000)
  expect_equal(levels(p$data$allele), c("reference", "alternate 1", "alternate 2"))

  b <- ggplot2::ggplot_build(p)
  lab <- unique(unlist(lapply(b$data, function(d)
    if ("label" %in% names(d)) as.character(d$label))))
  # rehh names the columns EHH_MAJ / EHH_MIN1 / EHH_MIN2 positionally, not by frequency, so
  # allele 0 stays "reference" even though it is not the rarest or the commonest by design
  expect_match(lab[1], "reference 30 \\(50%\\)")
  expect_match(lab[1], "alternate 1 10 \\(17%\\)")
  expect_match(lab[1], "alternate 2 20 \\(33%\\)")
  k <- as.integer(regmatches(lab[1], gregexpr("[0-9]+(?= \\()", lab[1], perl = TRUE))[[1]])
  expect_equal(sum(k), 60)                      # every haplotype accounted for, none dropped

  # three curves need three colours; two keep the pair the plot has always used
  curve_colours <- function(pp) {
    bb <- ggplot2::ggplot_build(pp)
    d <- bb$data[[which(vapply(bb$data, function(z)
      "colour" %in% names(z) && nrow(z) > 50, logical(1)))[1]]]
    sort(unique(d$colour))
  }
  expect_length(curve_colours(p), 3)
  bi <- .multi_hap(c(35, 25))
  expect_equal(curve_colours(plot_ehh(bi, "c1:11000", span = 12000)),
               sort(unname(plasgenomicsutilsR:::.EHH_FILL)))
})

test_that("nothing about the curves is fixed at three alleles", {
  skip_if_not_installed("rehh")
  skip_if_not_installed("ggplot2")
  # 1 ref + 3 alts is rare but legal, and rehh answers it with FREQ_MIN3. The column
  # selector, the labels and the palette are all sized from the data, so the only way to
  # know they stay in step is to ask at more than one arity.
  for (counts in list(c(24, 6, 18, 12), c(20, 5, 15, 10, 10))) {
    hap <- .multi_hap(counts)
    p <- plot_ehh(hap, "c1:11000", span = 12000)
    expect_equal(levels(p$data$allele),
                 c("reference", paste("alternate", seq_len(length(counts) - 1))))
    b <- ggplot2::ggplot_build(p)
    lab <- unique(unlist(lapply(b$data, function(d)
      if ("label" %in% names(d)) as.character(d$label))))
    k <- as.integer(regmatches(lab[1], gregexpr("[0-9]+(?= \\()", lab[1], perl = TRUE))[[1]])
    expect_equal(k, counts)                     # in allele order, not frequency order
    expect_equal(sum(k), sum(counts))           # no allele quietly left out
  }
})

test_that("the curve builder refuses a marker it cannot label rather than guessing", {
  # the guard: whatever rehh returns, the frequencies and the curves have to correspond, or
  # the labels would be attached to the wrong lines
  fake <- list(ehh = data.frame(POSITION = 1:3, EHH_MAJ = 1, EHH_MIN1 = 1, EHH_MIN2 = 1),
               freq = c(FREQ_MAJ = 0.5, FREQ_MIN = 0.5))
  local_mocked_bindings(calc_ehh = function(...) fake, .package = "rehh")
  hap <- .multi_hap(c(30, 10, 20))
  msg <- plasgenomicsutilsR:::.ehh_curve(hap, seq_len(nrow(hap$hap)), 11000, "c1", FALSE, 0.05)
  expect_type(msg, "character")
  expect_match(msg, "3 curves and 2 frequencies for the 3 allele\\(s\\)")
})

test_that("allele = 'index' refuses what is not an allele index", {
  G <- matrix(c(0, 1, 2, -1), 2, 2,
              dimnames = list(c("a", "b"), c("c1:100", "c1:200")))
  expect_error(parasite_haplotypes(G, maf = 0, alleles = "index"), "non-negative whole-number")
  G2 <- matrix(c(0, 1, 2, 0.5), 2, 2, dimnames = dimnames(G))
  expect_error(parasite_haplotypes(G2, maf = 0, alleles = "index"), "non-negative whole-number")
})


# A marker called on its own because the main callset dropped it for being multiallelic --
# the situation add_haplotype_markers() exists for.
.tri_bcf <- function(dir, chrom = "Pf3D7_13_v3", pos = 1725592,
                     samps = sprintf("s%02d", 1:12)) {
  gts <- c(rep("0/0", 6), rep("1/1", 3), rep("2/2", 2), "0/1")
  hdr <- c("##fileformat=VCFv4.2", sprintf("##contig=<ID=%s,length=2000000>", chrom),
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
           paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
                   "FORMAT", samps), collapse = "\t"),
           paste(c(chrom, pos, ".", "C", "T,A", ".", ".", ".", "GT", gts), collapse = "\t"))
  v <- file.path(dir, "tri.vcf")
  writeLines(hdr, v)
  v
}

.hap_on <- function(chrom, samps = sprintf("s%02d", 1:12), n_snp = 30, seed = 4) {
  set.seed(seed)
  G <- matrix(stats::rbinom(length(samps) * n_snp, 1, 0.4), length(samps), n_snp)
  rownames(G) <- samps
  colnames(G) <- paste0(chrom, ":", seq(1700000, by = 2000, length.out = n_snp))
  parasite_haplotypes(G, maf = 0.02, alleles = "index")
}

test_that("a triallelic marker called on its own joins the haplotypes with its alleles intact", {
  skip_if_not(nzchar(Sys.which("bcftools")))
  skip_if_not_installed("rehh")
  d <- tempfile(); dir.create(d)
  hap <- .hap_on("Pf3D7_13_v3")
  expect_setequal(unique(as.vector(hap$hap)), c(0L, 1L))     # biallelic to start with

  out <- add_haplotype_markers(hap, .tri_bcf(d), het = "draw")
  expect_setequal(unique(as.vector(out$hap)), c(0L, 1L, 2L))
  expect_equal(ncol(out$hap), ncol(hap$hap) + 1L)

  # a dosage matrix cannot hold this: copy.num.of.ref makes 1/1 and 2/2 the same number
  i <- which(out$map$snp_id == "Pf3D7_13_v3:1725591")       # 0-based, as the package counts
  expect_length(i, 1L)
  counts <- as.integer(table(out$hap[, i])[c("0", "1", "2")])
  expect_equal(counts, c(6L, 4L, 2L))            # the one het drew allele 1

  # inserted in coordinate order, not appended
  expect_false(is.unsorted(out$map$pos[out$map$chr == "Pf3D7_13_v3"]))
  # and it reaches plot_ehh as three curves
  p <- plot_ehh(out, "Pf3D7_13_v3:1725591", span = 30000)
  expect_equal(levels(p$data$allele), c("reference", "alternate 1", "alternate 2"))
})

test_that("the added marker adopts the chromosome spelling the haplotypes already use", {
  skip_if_not(nzchar(Sys.which("bcftools")))
  # SNPRelate keeps `Pf3D7_13_v3` from one file and reduces a recognised name to `13` in
  # another. A marker that keeps its own spelling lands on a chromosome of its own, with no
  # neighbours to decay against -- which surfaces as "too few polymorphic SNPs on that
  # chromosome", nothing like a naming problem.
  d <- tempfile(); dir.create(d)
  hap <- .hap_on("13")                                   # haplotypes say "13"
  out <- add_haplotype_markers(hap, .tri_bcf(d, chrom = "Pf3D7_13_v3"), het = "draw")
  expect_true(all(out$map$chr == "13"))
  expect_true("13:1725591" %in% out$map$snp_id)
  i <- which(out$map$snp_id == "13:1725591")
  expect_gt(i, 1L)                                       # sorted among its neighbours
  expect_lt(i, nrow(out$map))
})

test_that("add_haplotype_markers refuses what it cannot merge", {
  skip_if_not(nzchar(Sys.which("bcftools")))
  d <- tempfile(); dir.create(d)
  hap <- .hap_on("Pf3D7_13_v3")

  # a position already in the haplotypes is two answers, not a merge
  same <- .tri_bcf(d, pos = as.integer(sub(".*:", "", hap$map$snp_id[1])) + 1L)
  expect_error(add_haplotype_markers(hap, same), "already in the haplotypes")

  # a callset missing some of the haplotypes' samples cannot fill the column
  d2 <- tempfile(); dir.create(d2)
  few <- .tri_bcf(d2, samps = sprintf("s%02d", 1:8))
  expect_error(add_haplotype_markers(hap, few), "not in")
})


test_that("an allele keeps one name and one colour in every facet", {
  skip_if_not_installed("rehh")
  skip_if_not_installed("ggplot2")
  hap <- .multi_hap(c(30, 20, 10))
  i <- which(hap$map$pos == 11000)
  a <- hap$hap[, i]
  # one group without allele 2, which is what a region lacking a variant looks like. rehh
  # numbers its columns densely over the alleles it is shown, so that group's allele 1 would
  # come back as plain "alternate" -- one allele under two names in one plot, and a stray
  # factor level with no colour in the scale.
  grp <- ifelse(a == 2L, "has_all", rep(c("has_all", "no_alt2"), length.out = length(a)))
  hap$meta <- data.frame(sample = rownames(hap$hap), grp = grp)

  curves <- lapply(unique(grp), function(g)
    plasgenomicsutilsR:::.ehh_curve(hap, which(grp == g), hap$map$pos[i], hap$map$chr[i],
                                    FALSE, 0.05))
  lv <- lapply(curves, attr, "levels")
  expect_equal(lv[[1]], lv[[2]])                       # the same names in both groups
  expect_equal(lv[[1]], c("reference", "alternate 1", "alternate 2"))

  p <- plot_ehh(hap, "c1:11000", group = "grp", span = 12000)
  expect_equal(levels(p$data$allele), c("reference", "alternate 1", "alternate 2"))
  b <- ggplot2::ggplot_build(p)
  d <- b$data[[which(vapply(b$data, function(z)
    "colour" %in% names(z) && nrow(z) > 20, logical(1)))[1]]]
  expect_length(unique(d$colour), 3)                   # no fourth, uncoloured level
})

test_that("a group carrying no reference allele keeps both of its alternates", {
  skip_if_not_installed("rehh")
  hap <- .multi_hap(c(30, 20, 10))
  i <- which(hap$map$pos == 11000)
  a <- hap$hap[, i]
  rows <- which(a != 0L)
  # shown alleles 1 and 2, rehh reports FREQ_MAJ = 0 for the absent allele 0 and drops one
  # of the two real ones. The dense recoding is what stops that.
  cur <- plasgenomicsutilsR:::.ehh_curve(hap, rows, hap$map$pos[i], hap$map$chr[i],
                                         FALSE, 0.05)
  expect_false(is.character(cur))
  expect_setequal(levels(droplevels(cur$allele)), c("alternate 1", "alternate 2"))
  expect_equal(sum(attr(cur, "count")), length(rows))  # every haplotype still accounted for
  expect_equal(names(attr(cur, "count")), c("alternate 1", "alternate 2"))
})

test_that("add_ihs puts the focal SNP's score in the corner note", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("rehh")
  hap <- hap_for_ehh()
  # a pfcrt SNP that is variable in both countries, so both panels draw a curve
  foc <- "Pf3D7_07_v3:405361"
  G <- PF_EXAMPLE_DRUG_GENES
  sc <- suppressWarnings(run_ihs(hap, group = "country", min_maf = 0.02))
  sc_pooled <- suppressWarnings(run_ihs(hap, min_maf = 0.02))
  ann <- function(p) {
    i <- vapply(p$layers, function(l) is.data.frame(l$data) && "label" %in% names(l$data),
                logical(1))
    if (!any(i)) NULL else p$layers[[which(i)[1]]]$data
  }
  ehh <- function(...) plot_ehh(hap, foc, genes = G, span = 30000, ...)
  at <- function(d, g) d$label[as.character(d$group) == g]

  # off by default: the note is exactly what it was before
  base <- ann(ehh(group = "country"))
  expect_false(any(grepl("iHS", base$label)))

  # the score joins the counts on its own line, as the magnitude, matching the scan
  got <- suppressMessages(ann(ehh(group = "country", add_ihs = sc)))
  want <- abs(sc$ihs[sc$pos == 405361 & sc$group == "Cambodia"])
  expect_match(at(got, "Cambodia"), "^n = 30;.*\nabs\\(iHS\\) ")
  expect_equal(as.numeric(sub(".*abs\\(iHS\\) ", "", at(got, "Cambodia"))), round(want, 2))
  # the counts line itself is untouched by adding one
  expect_equal(sub("\n.*", "", at(got, "Cambodia")), at(base, "Cambodia"))
  # a panel the scan has no score for keeps its counts rather than borrowing the other's
  expect_equal(at(got, "Ghana"), at(base, "Ghana"))
  expect_message(ehh(group = "country", add_ihs = sc), "no iHS for Ghana")

  # an ungrouped scan labels a pooled plot
  pooled <- suppressMessages(ann(ehh(add_ihs = sc_pooled)))
  expect_match(pooled$label, "abs\\(iHS\\) ")

  # the bars of the conventional |iHS| render as letters at this text size, so the label
  # spells the magnitude out instead; pin that so it cannot drift back to bars
  expect_false(any(grepl("|", got$label, fixed = TRUE)))

  # polarized, the sign is meaningful and is kept
  signed <- sc; signed$ihs <- -abs(signed$ihs)
  pol <- suppressMessages(ann(ehh(group = "country", polarized = TRUE, add_ihs = signed)))
  expect_match(at(pol, "Cambodia"), "\niHS -")

  # `show_freq = FALSE` leaves the score alone rather than dropping the note entirely
  only <- suppressMessages(ann(ehh(group = "country", add_ihs = sc, show_freq = FALSE)))
  expect_equal(as.character(only$group), "Cambodia")
  expect_match(only$label, "^abs\\(iHS\\) ")

  # running the scan here agrees with handing one in
  ran <- suppressWarnings(suppressMessages(
    ann(ehh(group = "country", add_ihs = TRUE, ihs_args = list(min_maf = 0.02)))))
  expect_equal(at(ran, "Cambodia"), at(got, "Cambodia"))

  # a grouped scan cannot label a pooled plot, and says so instead of picking a group
  expect_message(ehh(add_ihs = sc), "pools every haplotype")
  # nor can a scan that does not carry the focal SNP
  expect_message(ehh(group = "country", add_ihs = sc[sc$pos != 405361, ]),
                 "no iHS for Pf3D7_07_v3:405361 in the scan")

  # a scan run on other haplotypes is caught by the frequency it reports at this SNP
  moved <- sc; moved$freq_minor[moved$pos == 405361] <- 0.3
  expect_message(ehh(group = "country", add_ihs = moved), "different set of haplotypes")
  # ...and a matching one does not trip it, which is what makes the check worth having
  msgs <- character(0)
  withCallingHandlers(ehh(group = "country", add_ihs = sc),
                      message = function(m) { msgs <<- c(msgs, conditionMessage(m))
                                              invokeRestart("muffleMessage") })
  expect_false(any(grepl("different set of haplotypes", msgs)))

  # the two arguments that would let the note describe a different scan than the curves
  expect_error(ehh(add_ihs = TRUE, ihs_args = list(group = "country")), "must not set")
  expect_error(ehh(add_ihs = TRUE, ihs_args = list(polarized = TRUE)), "must not set")
  expect_error(ehh(add_ihs = TRUE, ihs_args = list(1, 2)), "named list")
  expect_error(ehh(add_ihs = "yes"), "TRUE, FALSE, or a run_ihs")
  expect_error(ehh(add_ihs = sc[, c("group", "chr", "pos")]), "no 'ihs' column")
  expect_warning(ehh(group = "country", add_ihs = sc, ihs_args = list(min_maf = 0.02)),
                 "ignored")
})

test_that("a contrast column gives one corner line per contrast", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("rehh")
  hap <- hap_for_ehh()
  foc <- "Pf3D7_07_v3:405361"
  ann <- function(p) {
    i <- vapply(p$layers, function(l) is.data.frame(l$data) && "label" %in% names(l$data),
                logical(1))
    if (!any(i)) NULL else p$layers[[which(i)[1]]]$data
  }
  base <- data.frame(chr = "Pf3D7_07_v3", pos = 405361, group = c("Cambodia", "Ghana"),
                     ihs = c(-3.1, 2.4), stringsAsFactors = FALSE)
  # two named contrasts, the second missing for one panel
  two <- rbind(cbind(base, contrast = "A"),
               cbind(base[1, ], contrast = "B"))
  two$ihs[3] <- 1.75
  p <- plot_ehh(hap, foc, genes = PF_EXAMPLE_DRUG_GENES, span = 30000,
                group = "country", add_ihs = two)
  d <- ann(p)
  cam <- d$label[as.character(d$group) == "Cambodia"]
  gha <- d$label[as.character(d$group) == "Ghana"]
  # the panel with both gets both lines, each named
  expect_match(cam, "A abs\\(iHS\\) 3\\.10")
  expect_match(cam, "B abs\\(iHS\\) 1\\.75")
  # the panel with only one gets only that line, still named
  expect_match(gha, "A abs\\(iHS\\) 2\\.40")
  expect_false(grepl("B abs", gha))
  expect_equal(lengths(regmatches(cam, gregexpr("abs\\(iHS\\)", cam)))[[1]], 2L)

  # and a table without the column behaves exactly as it did before
  one <- ann(plot_ehh(hap, foc, genes = PF_EXAMPLE_DRUG_GENES, span = 30000,
                      group = "country", add_ihs = base))
  expect_match(one$label[as.character(one$group) == "Cambodia"], "^n = .*\nabs\\(iHS\\) 3\\.10$")
})
