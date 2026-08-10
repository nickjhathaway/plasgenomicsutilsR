scan_on_ibd_groups <- function(x) {
  ps <- example_pop_structure(umap = FALSE)
  scan <- as.data.frame(run_ihs(parasite_haplotypes(ps, maf = 0.05), group = "country"))
  # the two example datasets are different cohorts, so relabel the scan onto the IBD
  # groups to have something the mirror can face
  lv <- unique(x$get_per_snp_group()$group)
  scan$group <- factor(lv[as.integer(factor(scan$group))], levels = lv)
  scan
}

test_that(".combine_groups keeps the declared order across column types", {
  a <- factor(c("Kenya", "DRC"), levels = c("DRC", "Kenya"))
  b <- c("Uganda", "DRC")
  out <- .combine_groups(a, b)
  expect_s3_class(out, "factor")
  expect_identical(levels(out), c("DRC", "Kenya", "Uganda"))
  # the failure this guards against: c(factor, character) yields the integer codes
  expect_false(any(levels(out) %in% c("1", "2")))
  expect_null(.combine_groups(NULL, NULL))
})

test_that("the tug-of-war can hang an external scan from the top", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("rehh")
  x <- example_ibd_results()
  scan <- scan_on_ibd_groups(x)
  p <- plot_ibd_tugofwar(x, top = scan, top_label = "iHS", draw_threshold = 1)
  expect_s3_class(p, "ggplot")
  # both tracks must carry one factor, or the panels come back alphabetised
  built <- ggplot2::ggplot_build(p)
  expect_s3_class(built$layout$layout$group, "factor")
  expect_identical(levels(built$layout$layout$group),
                   unique(as.character(x$get_per_snp_group()$group)))
})

test_that("a top track that cannot face the IBD groups is refused", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("rehh")
  x <- example_ibd_results()
  ps <- example_pop_structure(umap = FALSE)
  hap <- parasite_haplotypes(ps, maf = 0.05)
  expect_error(plot_ibd_tugofwar(x, top = run_ihs(hap, group = "country")),
               "share no group")
  expect_error(plot_ibd_tugofwar(x, top = run_rsb(hap, group = "country")),
               "population pair")
  expect_error(plot_ibd_tugofwar(x, top = data.frame(pos = 1)), "needs chr and pos")
  expect_error(plot_ibd_tugofwar(x, top = scan_on_ibd_groups(x), metric = "nope"),
               "no 'nope' column")
})

test_that("a threshold above every value in the top track is not drawn in the IBD half", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("rehh")
  x <- example_ibd_results()
  scan <- scan_on_ibd_groups(x)
  expect_message(plot_ibd_tugofwar(x, top = scan, draw_threshold = 99),
                 "above every value in the top track")
  expect_silent(plot_ibd_tugofwar(x, top = scan, draw_threshold = FALSE))
})

test_that("plot_ibd_locus draws both tracks on two axes over one window", {
  skip_if_not_installed("ggplot2")
  x <- example_ibd_results()
  p <- plot_ibd_locus(x, "pfcrt", pad = 50000)
  expect_s3_class(p, "patchwork")
  panel <- ggplot2::ggplot_build(p[[1]])
  # the window, not the genome
  rng <- panel$layout$panel_params[[1]]$x.range
  expect_lt(diff(rng), 2e5)
  # a secondary axis carries the scan
  expect_true(!is.null(panel$layout$panel_params[[1]]$y.sec))
  expect_true("pfcrt" %in% p[[2]]$data$name)
  expect_s3_class(plot_ibd_locus(x, "pfcrt", pad = 50000, gene_track = FALSE), "ggplot")
})

test_that("plot_ibd_locus widens a gene-sized window rather than plotting nothing", {
  skip_if_not_installed("ggplot2")
  x <- example_ibd_results()
  expect_s3_class(plot_ibd_locus(x, "pfcrt"), "patchwork")
  expect_error(plot_ibd_locus(x, "pfcrt", min_span = 0, pad = 0),
               "no IBD SNPs in 7:")
})

test_that("plot_ibd_locus maps the optional point aesthetics only when asked", {
  skip_if_not_installed("ggplot2")
  x <- example_ibd_results()
  sel <- as.data.frame(x$get_selection())
  sel$mutation_type <- rep(c("synonymous", "nonsynonymous"), length.out = nrow(sel))
  sel$delta_af <- abs(sel$z_score) / max(abs(sel$z_score), na.rm = TRUE)

  # detected from the scan
  shaped <- plot_ibd_locus(x, "pfcrt", scan = sel, pad = 50000, gene_track = FALSE)
  expect_gt(length(unique(ggplot2::ggplot_build(shaped)$data[[2]]$shape)), 1)
  # and refused when told never to
  plain <- plot_ibd_locus(x, "pfcrt", scan = sel, shape_by = NA, pad = 50000,
                          gene_track = FALSE)
  # every built layer has shape/size columns, so it is the number of distinct values that
  # says whether anything was mapped
  expect_length(unique(ggplot2::ggplot_build(plain)$data[[2]]$shape), 1)
  expect_length(unique(ggplot2::ggplot_build(plain)$data[[2]]$size), 1)
  sized <- plot_ibd_locus(x, "pfcrt", scan = sel, size_by = "delta_af", pad = 50000,
                          gene_track = FALSE)
  expect_gt(length(unique(ggplot2::ggplot_build(sized)$data[[2]]$size)), 1)
  expect_error(plot_ibd_locus(x, "pfcrt", scan = sel, size_by = "nope"),
               "not a column of the scan")
})

test_that("plot_ibd_locus honours the group selection", {
  skip_if_not_installed("ggplot2")
  x <- example_ibd_results()
  p <- plot_ibd_locus(x, "7:380000-460000", groups = c("Ethiopia", "Kenya"))
  expect_identical(as.character(ggplot2::ggplot_build(p[[1]])$layout$layout$group),
                   c("Ethiopia", "Kenya"))
})

hlines <- function(p) {
  b <- ggplot2::ggplot_build(if (inherits(p, "patchwork")) p[[1]] else p)
  d <- b$data[vapply(b$data, function(z) !is.null(z$yintercept), logical(1))]
  list(n = length(d), y = unlist(lapply(d, function(z) z$yintercept)),
       range = b$layout$panel_params[[1]]$y.range)
}

test_that("plot_ibd_locus takes a threshold kind, not just a height", {
  skip_if_not_installed("ggplot2")
  x <- example_ibd_results()
  loc <- function(...) plot_ibd_locus(x, "pfcrt", pad = 30000, gene_track = FALSE, ...)

  expect_equal(hlines(loc(threshold = NA))$n, 0)
  expect_equal(hlines(loc(threshold = 5))$n, 1)
  expect_equal(hlines(loc())$n, 1)                     # the object's Bonferroni line
  expect_equal(hlines(loc(threshold = "bonferroni"))$n, 1)

  # a kind the run did not write is an error, not a silently missing line -- which is what a
  # bare string used to be, since it reached the plot as a non-finite "height"
  expect_error(loc(threshold = "fdr"), "no FDR threshold")
  expect_error(loc(threshold = "permutation"), "no permutation threshold")
  expect_error(loc(threshold = "nope"), "should be one of")
  # and a kind only describes the object's own statistic
  expect_error(loc(scan = as.data.frame(x$get_selection()), threshold = "bonferroni"),
               "names a kind in the object's own threshold table")
})

test_that("a threshold above everything in the window is still drawn", {
  skip_if_not_installed("ggplot2")
  x <- example_ibd_results()
  low <- hlines(plot_ibd_locus(x, "pfcrt", pad = 30000, gene_track = FALSE, threshold = NA))
  high <- hlines(plot_ibd_locus(x, "pfcrt", pad = 30000, gene_track = FALSE,
                                threshold = 500))
  expect_equal(high$n, 1)
  # the scale grew to hold it rather than ggplot dropping it for being out of limits
  expect_gt(high$range[2], low$range[2])
  expect_lte(max(high$y), high$range[2])
})

test_that("a per-group threshold draws one line per panel", {
  skip_if_not_installed("ggplot2")
  x <- example_ibd_results()
  thr <- x$get_thresholds()
  skip_if(is.null(thr) || !"group" %in% names(thr))
  skip_if(length(unique(thr$threshold[is.finite(thr$threshold)])) < 2)
  h <- hlines(plot_ibd_locus(x, "pfcrt", pad = 30000, gene_track = FALSE,
                             threshold = "bonferroni"))
  expect_gt(length(unique(h$y)), 1)
})

test_that("a signed metric is plotted as a magnitude rather than clipped at zero", {
  skip_if_not_installed("ggplot2")
  x <- example_ibd_results()
  sel <- as.data.frame(x$get_selection())
  expect_true(any(sel$z_score < 0, na.rm = TRUE))     # the case this is about
  loc <- function(...) plot_ibd_locus(x, "pfcrt", pad = 30000, gene_track = FALSE,
                                      scan = sel, metric = "z_score", ...)

  # default: magnitudes, and the axis says so
  b <- ggplot2::ggplot_build(loc())
  expect_equal(sum(b$data[[2]]$y < 0), 0)
  expect_match(b$layout$panel_params[[1]]$y.sec$name, "|z_score|", fixed = TRUE)

  # keeping the sign extends the panel below zero instead of hiding the negative half
  bs <- ggplot2::ggplot_build(loc(scan_abs = FALSE))
  expect_gt(sum(bs$data[[2]]$y < 0), 0)
  expect_lt(bs$layout$panel_params[[1]]$y.range[1], 0)
  expect_true(all(bs$data[[2]]$y >= bs$layout$panel_params[[1]]$y.range[1]))
  expect_match(bs$layout$panel_params[[1]]$y.sec$name, "z_score", fixed = TRUE)

  # an unsigned metric is left exactly as it was
  b0 <- ggplot2::ggplot_build(plot_ibd_locus(x, "pfcrt", pad = 30000, gene_track = FALSE))
  expect_match(b0$layout$panel_params[[1]]$y.sec$name, "neg_log10_p", fixed = TRUE)
})

test_that("the tug-of-war mirrors a signed metric as a magnitude, or refuses", {
  skip_if_not_installed("ggplot2")
  x <- example_ibd_results()
  sel <- as.data.frame(x$get_selection())
  lv <- unique(x$get_per_snp_group()$group)
  sel$group <- factor(lv[as.integer(factor(sel$group))], levels = lv)

  p <- plot_ibd_tugofwar(x, top = sel, metric = "z_score")
  b <- ggplot2::ggplot_build(p)
  # the top half spans [0, 1]; nothing may end up above it, which is where a negative value
  # would land and be silently clipped
  tips <- unlist(lapply(b$data, function(d) d$yend[is.finite(d$yend)]))
  expect_false(any(tips > 1.001))
  expect_match(b$layout$panel_scales_y[[1]]$name, "|z_score|", fixed = TRUE)

  # the mirror has no room for the sign, so asking for it is an error, not a clipped plot
  expect_error(plot_ibd_tugofwar(x, top = sel, metric = "z_score", scan_abs = FALSE),
               "no room below zero")
  expect_s3_class(plot_ibd_tugofwar(x), "ggplot")     # unsigned default untouched
})
