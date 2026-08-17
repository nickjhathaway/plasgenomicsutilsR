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
  expect_match(top$label[1], "^n = [0-9]+; reference [0-9]+%, alternate [0-9]+%$")

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
