ps_for_hap <- function() example_pop_structure(umap = FALSE)

# The heatmap panel, whether or not a dendrogram / gene track was stacked around it.
# `Filter()` over a patchwork returns a plain list rather than plots, so walk the indices.
hap_panel <- function(p) {
  if (!inherits(p, "patchwork")) return(p)
  for (i in seq_len(8)) {
    q <- tryCatch(p[[i]], error = function(e) NULL)
    if (is.null(q) || !is.data.frame(q$data)) next
    if ("call" %in% names(q$data)) return(q)
  }
  stop("no heatmap panel in this patchwork")
}

# the sample order the heatmap drew, top to bottom
drawn_rows <- function(p) {
  d <- hap_panel(p)$data
  unique(d$sample[order(d$.row)])
}

test_that("plot_region_haplotypes draws a heatmap over the window", {
  skip_if_not_installed("ggplot2")
  ps <- ps_for_hap()
  p <- plot_region_haplotypes(ps, "7", genes = PF_EXAMPLE_DRUG_GENES)
  expect_s3_class(p, "patchwork")
  # one row per sample, one column per SNP in the window
  hm <- hap_panel(p)
  expect_length(unique(hm$data$sample), length(ps$get_samples()))
  # every SNP of the FULL panel on that chromosome -- this plot prefers it, since pruning
  # removes the correlated SNPs a haplotype block is made of
  loci <- .parse_snp_ids(colnames(ps$genotype(prefer = "full")))
  in_win <- sum(normalise_chr(loci$chr) == "7")
  expect_length(unique(hm$data$snp_id), in_win)
  expect_gt(in_win, ncol(ps$genotype()) / 2)   # the dense windows are actually being read
  expect_true(all(levels(hm$data$call) == c("reference", "mixed", "alternate")))
})

test_that("split blocks the rows in the metadata's level order and clusters inside each", {
  skip_if_not_installed("ggplot2")
  ps <- ps_for_hap()
  p <- plot_region_haplotypes(ps, "7", split = "country")
  hm <- hap_panel(p)
  by_row <- hm$data[order(hm$data$.row), c("sample", ".split")]
  by_row <- by_row[!duplicated(by_row$sample), ]

  meta <- ps$get_meta()
  expect_identical(levels(by_row$.split), levels(.as_group_factor(meta$country)))
  # every block is contiguous: the split is what fixes the blocks, clustering only reorders
  # samples inside them
  runs <- rle(as.character(by_row$.split))
  expect_length(runs$values, nlevels(by_row$.split))
  expect_identical(runs$values, levels(by_row$.split))

  # and the order within a block is learned, not the order the samples arrived in
  ordered <- plot_region_haplotypes(ps, "7", split = "country", cluster = FALSE)
  expect_false(identical(drawn_rows(ordered), drawn_rows(p)))
})

test_that("the dendrogram lines up with the rows it labels", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  ps <- ps_for_hap()
  p <- plot_region_haplotypes(ps, "7", split = "country")
  dend <- p[[1]]; hm <- hap_panel(p)
  db <- ggplot2::ggplot_build(dend); hb <- ggplot2::ggplot_build(hm)
  expect_length(db$layout$panel_params, length(hb$layout$panel_params))
  for (i in seq_along(hb$layout$panel_params))
    expect_equal(db$layout$panel_params[[i]]$y.range,
                 hb$layout$panel_params[[i]]$y.range, tolerance = 1e-6)
  # no dendrogram when there is nothing to cluster
  expect_false(inherits(plot_region_haplotypes(ps, "7", cluster = FALSE), "patchwork"))
})

test_that("spacing decides what the x axis means", {
  skip_if_not_installed("ggplot2")
  ps <- ps_for_hap()
  xr <- function(...) {
    d <- hap_panel(plot_region_haplotypes(ps, "7", ...))$data
    range(c(d$xmin, d$xmax))
  }
  even <- xr(spacing = "even")
  genomic <- xr(spacing = "genomic")
  # even counts SNP columns, so its axis tops out at the number of them; genomic is in base
  # pairs, so it spans the region itself
  n_snp <- length(unique(hap_panel(plot_region_haplotypes(ps, "7"))$data$snp_id))
  expect_equal(even[2], n_snp + 0.5)
  expect_gt(genomic[2], 1e5)
  expect_gt(genomic[2], even[2] * 100)
})

test_that("mark_snps takes an id, a position or a gene", {
  skip_if_not_installed("ggplot2")
  ps <- ps_for_hap()
  loci <- .parse_snp_ids(colnames(ps$genotype()))
  loci$chr <- normalise_chr(loci$chr)
  on7 <- loci[loci$chr == "7", ]
  id <- colnames(ps$genotype())[on7$idx[1]]

  n_marks <- function(...) {
    b <- ggplot2::ggplot_build(
      hap_panel(plot_region_haplotypes(ps, "7", spacing = "genomic", ...)))
    sum(vapply(b$data, function(z) sum(!is.null(z$xintercept)), integer(1)))
  }
  expect_gt(n_marks(mark_snps = id), 0)
  expect_gt(n_marks(mark_snps = on7$pos[1]), 0)
  # a gene with SNPs in it marks them all
  expect_gt(n_marks(mark_snps = "pfcrt", genes = PF_EXAMPLE_DRUG_GENES), 0)
  # one with none says so, rather than looking like the argument was ignored
  empty <- data.frame(name = "nosnps", chr = "7", start = 1, end = 2)
  expect_message(plot_region_haplotypes(ps, "7", mark_snps = "nosnps", genes = empty),
                 "no genotyped SNP inside nosnps")
  expect_error(plot_region_haplotypes(ps, "7", mark_snps = "not-a-thing"),
               "not a SNP in the window")
})

test_that("plot_region_haplotypes refuses windows it cannot draw", {
  skip_if_not_installed("ggplot2")
  ps <- ps_for_hap()
  expect_error(plot_region_haplotypes(ps, "7:1-1000"), "no genotyped SNPs")
  expect_error(plot_region_haplotypes(ps, "7", max_snps = 2), "more than `max_snps`")
  expect_error(plot_region_haplotypes(ps, "7", split = "nope"),
               "is not a metadata column")
  expect_error(plot_region_haplotypes(ps, "7", samples = "nobody"),
               "none of `samples`")
})

test_that("samples can be narrowed, and the R6 method is the same plot", {
  skip_if_not_installed("ggplot2")
  ps <- ps_for_hap()
  keep <- head(ps$get_samples(), 12)
  p <- plot_region_haplotypes(ps, "7", samples = keep)
  hm <- hap_panel(p)
  expect_setequal(unique(hm$data$sample), keep)
  # few enough rows to be worth labelling
  expect_true(all(keep %in% ggplot2::ggplot_build(hm)$layout$panel_params[[1]]$y$get_labels()))
  expect_s3_class(ps$plot_region_haplotypes("7", samples = keep), class(p)[1])
})

test_that("the call labels follow which allele the dosages count", {
  skip_if_not_installed("ggplot2")
  ps <- ps_for_hap()
  expect_identical(ps$allele(), "alt")           # the fixture records it

  alt <- hap_panel(plot_region_haplotypes(ps, "7"))$data
  ref <- hap_panel(plot_region_haplotypes(ps, "7", allele = "ref"))$data
  # 2 is homozygous alternate under alt dosage and homozygous reference under ref dosage, so
  # the two readings of the same matrix are mirror images -- getting it wrong silently
  # mislabels every call
  two <- alt$value == 2 & !is.na(alt$value)
  expect_true(all(as.character(alt$call[two]) == "alternate"))
  expect_true(all(as.character(ref$call[two]) == "reference"))
  zero <- alt$value == 0 & !is.na(alt$value)
  expect_true(all(as.character(alt$call[zero]) == "reference"))
  expect_true(all(as.character(ref$call[zero]) == "alternate"))

  # an object that cannot say assumes alt, and says so rather than guessing silently
  bare <- PopStructure$new(ps$genotype(), meta = ps$get_meta())
  expect_null(bare$allele())
  expect_message(plot_region_haplotypes(bare, "7"), "does not record which allele")
})

test_that("every call has its own colour and only missing data is grey", {
  ps <- ps_for_hap()
  expect_setequal(names(.GENO_FILL), c("reference", "mixed", "alternate"))
  expect_length(unique(unname(.GENO_FILL)), 3L)
  # none of them may be near-white, or a call is indistinguishable from the panel and from
  # the grey of a missing call
  rgb_of <- function(h) grDevices::col2rgb(h)[, 1]
  expect_true(all(vapply(.GENO_FILL, function(h) mean(rgb_of(h)) < 220, logical(1))))
})

test_that("annotations draw one coloured strip per column, sharing the object's colours", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggnewscale")
  skip_if_not_installed("patchwork")
  ps <- ps_for_hap()
  # a second annotation to prove each gets its own scale rather than sharing one palette
  meta <- ps$get_meta()
  meta$half <- ifelse(seq_len(nrow(meta)) %% 2 == 0, "even", "odd")
  ps$add_meta(meta)
  p <- plot_region_haplotypes(ps, "7", split = "country",
                              annotations = c("country", "half"))
  # the strips live in their own panel, one x position per annotation
  ann <- NULL
  for (i in seq_len(6)) {
    q <- tryCatch(p[[i]], error = function(e) NULL)
    if (is.null(q)) next
    xs <- tryCatch(ggplot2::ggplot_build(q)$layout$panel_params[[1]]$x$get_labels(),
                   error = function(e) NULL)
    if (!is.null(xs) && all(c("country", "half") %in% xs)) { ann <- q; break }
  }
  expect_false(is.null(ann))
  # two annotations -> two fill scales, so a level of one never borrows the other's colour
  expect_gte(length(ann$layers), 2L)
  expect_error(plot_region_haplotypes(ps, "7", annotations = "nope"),
               "not a metadata column")
})

test_that("genomic spacing gives every SNP the same width at its own position", {
  skip_if_not_installed("ggplot2")
  ps <- ps_for_hap()
  d <- hap_panel(plot_region_haplotypes(ps, "7", spacing = "genomic"))$data
  w <- unique(round(d$xmax - d$xmin, 6))
  # one width for all of them: equal marks are what make the distances between SNPs readable,
  # and stretching each tile to its neighbours would fill the gaps back in
  expect_length(w, 1L)
  cent <- unique(round((d$xmin + d$xmax) / 2))
  expect_setequal(cent, unique(d$pos))
  # a wider mark on request
  d2 <- hap_panel(plot_region_haplotypes(ps, "7", spacing = "genomic",
                                         snp_width = 5000))$data
  expect_equal(unique(round(d2$xmax - d2$xmin)), 5000)
})

test_that("borders are drawn by default and can be turned off", {
  skip_if_not_installed("ggplot2")
  ps <- ps_for_hap()
  col_of <- function(...) {
    b <- ggplot2::ggplot_build(hap_panel(plot_region_haplotypes(ps, "7", ...)))
    unique(b$data[[1]]$colour)
  }
  expect_false(any(is.na(col_of())))
  expect_true(all(is.na(col_of(border = FALSE))))
})

test_that("each block is named in exactly one place", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  ps <- ps_for_hap()
  blank <- function(q) {
    if (is.null(q$theme)) return(TRUE)
    any(vapply(c("strip.text", "strip.text.y", "strip.text.y.right"),
               function(k) inherits(q$theme[[k]], "element_blank"), logical(1)))
  }
  named <- function(p) {
    n <- 0L
    for (i in seq_len(6)) {
      q <- tryCatch(p[[i]], error = function(e) NULL)
      if (is.null(q) || is.null(q$facet) || inherits(q$facet, "FacetNull")) next
      if (!blank(q)) n <- n + 1L
    }
    n
  }
  # with annotations the strips belong to the annotation panel; the dendrogram never names
  # them, or every label appears twice and the dendrogram is pushed off the genotypes
  expect_equal(named(plot_region_haplotypes(ps, "7", split = "country",
                                            annotations = "country")), 1L)
  # without them the heatmap is the one place
  expect_equal(named(plot_region_haplotypes(ps, "7", split = "country")), 1L)
})
