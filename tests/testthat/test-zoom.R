test_that(".resolve_region reads every accepted spelling", {
  chrom <- .resolve_region("7")
  expect_equal(chrom$chr, "7")
  expect_equal(chrom$start, 0)
  expect_gt(chrom$end, 1e6)

  rng <- .resolve_region("Pf3D7_07_v3:728,081-988,719")
  expect_equal(rng, list(chr = "7", start = 728081, end = 988719))

  gene <- .resolve_region("pfcrt", genes = PF_EXAMPLE_DRUG_GENES)
  expect_equal(gene$chr, "7")
  expect_lt(gene$end - gene$start, 1e4)

  df <- .resolve_region(data.frame(chr = "chr4", start = 10, end = 20))
  expect_equal(df, list(chr = "4", start = 10, end = 20))

  expect_null(.resolve_region(NULL))
  expect_error(.resolve_region("not-a-place"), "could not read")
  expect_error(.resolve_region(data.frame(chr = "7")), "needs chr, start and end")
})

test_that(".pad_region takes a fraction or base pairs and clamps to the chromosome", {
  layout <- .chrom_layout(DEFAULT_REFERENCE)
  iv <- list(chr = "7", start = 400000, end = 500000)
  expect_equal(.pad_region(iv, 0.1, layout)$start, 390000)
  expect_equal(.pad_region(iv, 0.1, layout)$end, 510000)
  expect_equal(.pad_region(iv, 25000, layout)$start, 375000)
  # clamped at both ends
  len <- layout$len[layout$chr == "7"]
  wide <- .pad_region(list(chr = "7", start = 100, end = len - 100), 1e6, layout)
  expect_equal(wide$start, 0)
  expect_equal(wide$end, len)
  expect_identical(.pad_region(iv, 0, layout), iv)
})

test_that("zoom crops the axis without moving the coordinates", {
  skip_if_not_installed("ggplot2")
  x <- example_ibd_results()
  full <- ggplot2::ggplot_build(plot_selection_manhattan(x))
  zoomed <- plot_selection_manhattan(x, zoom = "7:400000-500000", zoom_pad = 0)
  # a patchwork: the scan on top, the gene track underneath
  built <- ggplot2::ggplot_build(zoomed[[1]])
  rng <- built$layout$panel_params[[1]]$x.range
  layout <- x$chrom_layout()
  offset <- layout$offset[layout$chr == "7"]
  expect_lt(abs((rng[1] - offset) - 400000), 2000)   # 1% scale expansion only
  expect_lt(abs((rng[2] - offset) - 500000), 2000)
  # the same locus sits at the same x in both plots
  expect_true(rng[1] > full$layout$panel_params[[1]]$x.range[1])
  expect_true(diff(rng) < diff(full$layout$panel_params[[1]]$x.range))
})

test_that("zoom names the genes in the window in a track underneath", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  x <- example_ibd_results()
  p <- plot_ibd_sharing_manhattan(x, zoom = "pfcrt", zoom_pad = 20000)
  expect_s3_class(p, "patchwork")
  expect_true("pfcrt" %in% p[[2]]$data$name)
  # labelling off -> no track, so a plain ggplot again
  expect_false(inherits(plot_ibd_sharing_manhattan(x, zoom = "pfcrt", zoom_pad = 20000,
                                                  label_genes = FALSE), "patchwork"))
  # a window with no gene in it also stays a plain ggplot
  expect_false(inherits(plot_ibd_sharing_manhattan(x, zoom = "1:200000-260000"),
                        "patchwork"))
})

test_that("zoom keeps the declared group order and sizes for one window", {
  skip_if_not_installed("ggplot2")
  x <- example_ibd_results()
  p <- plot_selection_manhattan(x, zoom = "7")
  built <- ggplot2::ggplot_build(p[[1]])
  expect_identical(levels(built$layout$layout$group), levels(x$get_selection()$group))
  # a window is one width regardless of how much genome the object covers
  expect_equal(attr(p, "plasgenomics_dims")[["width"]],
               attr(plot_selection_manhattan(x, zoom = "1"), "plasgenomics_dims")[["width"]])
})

test_that("zoom scales the y axis to the window, not to the genome", {
  skip_if_not_installed("ggplot2")
  x <- example_ibd_results()
  full <- ggplot2::ggplot_build(plot_selection_manhattan(x))
  zoomed <- ggplot2::ggplot_build(
    plot_selection_manhattan(x, zoom = "7:400000-460000", label_genes = FALSE))
  # the window's tallest SNP is well below the genome-wide peak, so leaving the genome-wide
  # rows in would flatten the panel
  expect_lt(zoomed$layout$panel_params[[1]]$y.range[2],
            full$layout$panel_params[[1]]$y.range[2])
})

test_that("an empty zoom window is an error, not a blank panel", {
  skip_if_not_installed("ggplot2")
  x <- example_ibd_results()
  expect_error(plot_ibd_sharing_manhattan(x, zoom = "pfcrt", zoom_pad = 0),
               "no IBD SNPs in 7:")
  expect_error(plot_ibd_tugofwar(x, zoom = "pfcrt", zoom_pad = 0), "no ")
})

test_that("zoom is rejected when the chromosome is not being drawn", {
  skip_if_not_installed("ggplot2")
  x <- example_ibd_results()
  expect_error(plot_selection_manhattan(x, chroms = "1", zoom = "7"),
               "which this plot is not showing")
})

test_that("the scan plots zoom the same way", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("rehh")
  ps <- example_pop_structure(umap = FALSE)
  scan <- run_ihs(parasite_haplotypes(ps, maf = 0.05), group = "country")
  p <- plot_ihs(scan, genes = PF_EXAMPLE_DRUG_GENES, zoom = "7")
  expect_s3_class(p, "patchwork")
  built <- ggplot2::ggplot_build(p[[1]])
  expect_identical(levels(built$layout$layout$.facet), levels(scan$group))
})

test_that("genes_for_track fills the track while the plot's own genes mark the panel", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  x <- example_ibd_results()          # its track is the eight drug genes
  p <- plot_ibd_sharing_manhattan(x, zoom = "pfpx1", zoom_pad = 20000,
                                  genes_for_track = PF3D7_GENES)
  track <- p[[2]]$data
  # the full annotation has many genes in this window, the drug track only pfpx1
  expect_gt(nrow(track), 3)
  expect_true("pfpx1" %in% track$name)
  # ... while the reference lines in the panel still come from the object's own genes
  marks <- ggplot2::ggplot_build(p[[1]])
  n_marks <- sum(vapply(marks$data, function(d) sum(!is.null(d$xintercept)), integer(1)))
  expect_lt(n_marks, nrow(track))
  # without it, both come from the same short track
  q <- plot_ibd_sharing_manhattan(x, zoom = "pfpx1", zoom_pad = 20000)
  expect_lt(nrow(q[[2]]$data), nrow(track))
})

test_that("a gene only in the full annotation still resolves as a zoom target", {
  x <- example_ibd_results()
  expect_error(.resolve_region("pfaif", PF_EXAMPLE_DRUG_GENES), "could not read")
  iv <- .resolve_region("pfaif", PF_EXAMPLE_DRUG_GENES, fallback = PF3D7_GENES)
  expect_equal(iv$chr, "7")
  # and by systematic id
  expect_equal(.resolve_region("PF3D7_0720900", NULL, fallback = PF3D7_GENES)$chr, "7")
})

test_that("rotating the gene names moves them to one baseline and keeps them in the window", {
  skip_if_not_installed("ggplot2")
  x <- example_ibd_results()
  lay <- x$chrom_layout()
  z <- .zoom_setup("pfpx1", x$get_genes(), lay, TRUE, 20000, "pf3d7", 0, PF3D7_GENES)
  boxes <- .zoom_gene_boxes(z$track, lay)

  flat <- .gene_track_panel(boxes, z$xlim, angle = 0)
  turned <- .gene_track_panel(boxes, z$xlim, angle = 45)
  # horizontal: a label sits under its own box row, so wide names force extra rows
  expect_gt(max(flat$data$.row), max(turned$data$.row))
  # rotated: boxes pack only on where the genes actually are, and every name hangs below
  # all of them rather than under its own row
  expect_lte(max(turned$data$.row), max(flat$data$.row))
  expect_true(all(turned$data$.label_y < min(turned$data$.y)))
  expect_true(all(turned$data$.label_y %in% unique(turned$data$.label_y)[1:(max(turned$data$.tier) + 1)]))
  # and the strip is told how many inches it needs, so nothing is clipped
  expect_gt(attr(turned, "track_in"), attr(flat, "track_in"))

  # no name runs outside the window on either setting
  for (tr in list(flat, turned)) {
    d <- tr$data
    span <- diff(z$xlim)
    fp <- nchar(d$name) * .GENE_LABEL_SIZE * .CHAR_IN * cos(0) * span / (9 - 1.2)
    lo <- d$.mid - fp * d$.hjust
    expect_true(all(lo >= z$xlim[1] - span * 0.02))
    expect_true(all(lo + fp <= z$xlim[2] + span * 0.02))
  }
})

test_that("the page grows for the gene track that was actually drawn", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  x <- example_ibd_results()
  h <- function(...) attr(plot_ibd_sharing_manhattan(x, zoom = "pfpx1", zoom_pad = 20000,
                                                    ...), "plasgenomics_dims")[["height"]]
  expect_gt(h(genes_for_track = PF3D7_GENES, gene_label_angle = 45),
            h(genes_for_track = PF3D7_GENES))
  expect_gt(h(genes_for_track = PF3D7_GENES), h(label_genes = FALSE))
})

test_that("padding can differ on each side", {
  layout <- .chrom_layout(DEFAULT_REFERENCE)
  iv <- list(chr = "7", start = 400000, end = 410000)
  side <- function(pad) {
    r <- .pad_region(iv, pad, layout)
    c(left = 400000 - r$start, right = r$end - 410000)
  }
  expect_equal(side(20000), c(left = 20000, right = 20000))
  expect_equal(side(c(5000, 40000)), c(left = 5000, right = 40000))
  # named, in either order, and naming one side pads only that side
  expect_equal(side(c(right = 40000, left = 5000)), c(left = 5000, right = 40000))
  expect_equal(side(c(right = 40000)), c(left = 0, right = 40000))
  expect_equal(side(c(left = 5000)), c(left = 5000, right = 0))
  # each side judged on its own: a fraction below 1, base pairs at or above it
  expect_equal(side(c(0.1, 20000)), c(left = 1000, right = 20000))
  expect_equal(side(NULL), c(left = 0, right = 0))

  # still clamped to the chromosome
  wide <- .pad_region(list(chr = "7", start = 100, end = 200), c(left = 1e6, right = 1e9),
                      layout)
  expect_equal(wide$start, 0)
  expect_equal(wide$end, layout$len[layout$chr == "7"])

  expect_error(.pad_region(iv, "x", layout), "must be numeric")
  expect_error(.pad_region(iv, c(1, 2, 3), layout), "at most two values")
  expect_error(.pad_region(iv, c(lo = 1, hi = 2), layout), "names must be")
  expect_error(.pad_region(iv, c(left = 1, 2), layout), "names must be")
})

test_that("asymmetric padding reaches the plots", {
  skip_if_not_installed("ggplot2")
  x <- example_ibd_results()
  rng <- function(p) {
    b <- ggplot2::ggplot_build(if (inherits(p, "patchwork")) p[[1]] else p)
    b$layout$panel_params[[1]]$x.range
  }
  a <- rng(plot_selection_manhattan(x, zoom = "pfcrt", zoom_pad = 20000))
  b <- rng(plot_selection_manhattan(x, zoom = "pfcrt",
                                    zoom_pad = c(left = 20000, right = 60000)))
  expect_equal(a[1], b[1], tolerance = 0.01)     # same left edge
  expect_gt(b[2], a[2])                          # more room on the right
  l <- rng(plot_ibd_locus(x, "pfcrt", pad = c(left = 60000, right = 20000), min_span = 0))
  r <- rng(plot_ibd_locus(x, "pfcrt", pad = c(left = 20000, right = 60000), min_span = 0))
  expect_lt(l[1], r[1])
  expect_lt(l[2], r[2])
})
