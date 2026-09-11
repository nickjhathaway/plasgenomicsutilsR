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

test_that("mark_snps takes an interval table, so a codon table needs no conversion", {
  skip_if_not_installed("ggplot2")
  ps <- ps_for_hap()
  loci <- .parse_snp_ids(colnames(ps$genotype(prefer = "full")))
  on7 <- loci$pos[normalise_chr(loci$chr) == "7"]
  # an interval covering three known SNPs, in the package's 0-based half-open terms
  iv <- data.frame(chr = "7", start = min(on7), end = sort(on7)[3] + 1)
  n_marks <- function(m) {
    b <- ggplot2::ggplot_build(hap_panel(
      plot_region_haplotypes(ps, "7", spacing = "genomic", mark_snps = m)))
    length(unlist(lapply(b$data, function(z) z$xintercept)))
  }
  expect_equal(n_marks(iv), 3)
  # the same three by bare position
  expect_equal(n_marks(sort(on7)[1:3]), 3)
  # an interval with nothing in it says so rather than drawing nothing silently
  expect_message(plot_region_haplotypes(ps, "7", mark_snps = data.frame(chr = "7", start = 1,
                                                                       end = 2)),
                 "no genotyped SNP inside any")
  expect_error(plot_region_haplotypes(ps, "7", mark_snps = data.frame(x = 1)),
               "needs start and end")
})

# every SNP tile and every gene box, from a drawn plot
hap_geometry <- function(p) {
  hm <- tr <- NULL
  for (i in seq_len(6)) {
    q <- tryCatch(p[[i]], error = function(e) NULL)
    if (is.null(q) || !is.data.frame(q$data)) next
    if ("call" %in% names(q$data)) hm <- q
    if (".gene_xmin" %in% names(q$data)) tr <- q
  }
  list(tiles = hm$data[!duplicated(hm$data$snp_id), c("snp_id", "pos", "xmin", "xmax")],
       genes = if (is.null(tr)) NULL else
         tr$data[, c("name", "start", "end", ".gene_xmin", ".gene_xmax")])
}

test_that("a SNP is only ever drawn over the genes it actually falls in", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  ps <- ps_for_hap()
  for (sp in c("even", "genomic", "gapped")) {
    g <- hap_geometry(suppressMessages(
      plot_region_haplotypes(ps, "pfcrt", pad = 30000, genes = PF3D7_GENES, spacing = sp)))
    skip_if(is.null(g$genes) || !nrow(g$genes))
    bad <- 0
    for (k in seq_len(nrow(g$tiles))) for (j in seq_len(nrow(g$genes))) {
      overlaps <- g$tiles$xmax[k] > g$genes$.gene_xmin[j] &&
                  g$tiles$xmin[k] < g$genes$.gene_xmax[j]
      inside <- g$tiles$pos[k] >= g$genes$start[j] && g$tiles$pos[k] < g$genes$end[j]
      if (overlaps && !inside) bad <- bad + 1
    }
    expect_equal(bad, 0, info = paste(sp, "spacing: tiles over a gene they are not in"))
  }
})

test_that("under even spacing a gene's box is exactly the columns it holds", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  ps <- ps_for_hap()
  g <- hap_geometry(suppressMessages(
    plot_region_haplotypes(ps, "pfcrt", pad = 30000, genes = PF3D7_GENES)))
  skip_if(is.null(g$genes) || !nrow(g$genes))
  tiles <- g$tiles[order(g$tiles$xmin), ]
  for (j in seq_len(nrow(g$genes))) {
    cols <- which(tiles$pos >= g$genes$start[j] & tiles$pos < g$genes$end[j])
    expect_gt(length(cols), 0)                       # empty genes are dropped, not drawn
    expect_equal(g$genes$.gene_xmin[j], min(cols) - 0.5)
    expect_equal(g$genes$.gene_xmax[j], max(cols) + 0.5)
  }
  # a gene in the window with no genotyped SNP has no width on a SNP-index axis, so it is
  # left off rather than drawn at an interpolated spot under someone else's SNPs
  expect_message(plot_region_haplotypes(ps, "pfcrt", pad = 30000, genes = PF3D7_GENES),
                 "hold no genotyped SNP")
})

test_that("the genotype legend sits above the annotations, which follow the order asked for", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("patchwork")
  testthat::skip_if_not_installed("ggnewscale")
  # ggplot sorts guides by an internal hash unless each carries an `order`, so without one the
  # stack rearranged itself between datasets and two figures stopped being comparable.
  fill_orders <- function(p) {
    subs <- c(list(p), if (!is.null(p$patches)) p$patches$plots)
    rows <- list()
    for (q in subs) {
      b <- tryCatch(ggplot2::ggplot_build(q), error = function(e) NULL)
      if (is.null(b)) next
      for (sc in b$plot$scales$scales) {
        # ggnewscale renames the stashed scales (fill_new, ...), so match on the family
        if (!grepl("fill", sc$aesthetics[1])) next
        nm <- sc$name
        if (!is.character(nm) || length(nm) != 1) next
        g <- sc$guide
        ord <- if (inherits(g, "Guide")) g$params$order else if (is.list(g)) g$order else NA
        rows[[length(rows) + 1L]] <- data.frame(name = nm, order = as.numeric(ord)[1])
      }
    }
    d <- do.call(rbind, rows)
    d$name[order(d$order)]
  }

  build <- function(seed, nlev) {
    obj <- example_pop_structure(umap = FALSE)
    mm <- obj$get_meta()
    set.seed(seed)
    mm$marker <- sample(letters[seq_len(nlev)], nrow(mm), replace = TRUE)
    mm$site <- sample(LETTERS[seq_len(nlev + 1)], nrow(mm), replace = TRUE)
    obj$add_meta(mm[, c("sample", "marker", "site")])
    obj
  }
  seen <- lapply(list(c(1, 2), c(7, 3), c(11, 4)), function(cfg)
    fill_orders(suppressMessages(plot_region_haplotypes(
      build(cfg[1], cfg[2]), "pfcrt", pad = 20000, genes = PF_EXAMPLE_DRUG_GENES,
      annotations = c("marker", "site")))))

  expect_equal(seen[[1]], c("call", "marker", "site"))
  expect_equal(seen[[2]], seen[[1]])       # more marker levels must not reshuffle them
  expect_equal(seen[[3]], seen[[1]])

  # and the annotations follow the order they were listed in, not alphabetical
  rev_ann <- fill_orders(suppressMessages(plot_region_haplotypes(
    build(1, 2), "pfcrt", pad = 20000, genes = PF_EXAMPLE_DRUG_GENES,
    annotations = c("site", "marker"))))
  expect_equal(rev_ann, c("call", "site", "marker"))
})

# the fill each annotation level was actually drawn with, one entry per strip
.strip_map <- function(p) {
  for (i in seq_len(8)) {
    q <- tryCatch(p[[i]], error = function(e) NULL)
    if (is.null(q) || !length(q$layers)) next
    d1 <- tryCatch(q$layers[[1]]$data, error = function(e) list())
    if (!"value" %in% names(d1)) next
    b <- ggplot2::ggplot_build(q)
    return(lapply(seq_along(q$layers), function(k) {
      d <- q$layers[[k]]$data
      m <- unique(data.frame(level = as.character(d$value), fill = b$data[[k]]$fill,
                             stringsAsFactors = FALSE))
      stats::setNames(m$fill, m$level)[sort(m$level)]
    }))
  }
  stop("no annotation strip")
}

test_that("annotation colours can be set per call without touching the object", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggnewscale")
  ps <- example_pop_structure("africa", umap = FALSE)
  regs <- sort(unique(as.character(ps$get_meta()$region)))
  before <- ps$get_colors()$region

  # a named, partial palette recolours those levels and leaves the rest of the shared map
  p <- plot_region_haplotypes(ps, "7", annotations = c("country", "region"),
                              annotation_colours = list(
                                region = stats::setNames("#FF00FF", regs[1])))
  got <- .strip_map(p)[[2]]
  expect_equal(unname(got[regs[1]]), "#FF00FF")
  expect_equal(unname(got[regs[2]]), unname(before[regs[2]]))

  # unnamed is positional in level order
  pos <- plot_region_haplotypes(ps, "7", annotations = "region",
                                annotation_colours = list(region = c("#111111", "#222222")))
  expect_equal(unname(.strip_map(pos)[[1]][regs]), c("#111111", "#222222"))

  # either spelling
  amer <- plot_region_haplotypes(ps, "7", annotations = "region",
                                 annotation_colors = list(region = c("#111111", "#222222")))
  expect_equal(.strip_map(pos), .strip_map(amer))

  # the object's shared map is unchanged, so other plots keep their colours
  expect_equal(ps$get_colors()$region, before)
})

test_that("annotation colours say when they are given something unusable", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggnewscale")
  ps <- example_pop_structure("africa", umap = FALSE)
  expect_error(plot_region_haplotypes(ps, "7", annotations = "region",
                                      annotation_colours = c("#111111", "#222222")),
               "named list")
  expect_error(plot_region_haplotypes(ps, "7", annotations = "region",
                                      annotation_colours = list(region = "#111111")),
               "colour\\(s\\) for 2 level\\(s\\)")
  expect_warning(plot_region_haplotypes(ps, "7", annotations = "region",
                                        annotation_colours = list(nope = c(a = "#111111"))),
                 "not an annotation")
  expect_warning(plot_region_haplotypes(ps, "7", annotations = "region",
                                        annotation_colours = list(
                                          region = c(nowhere = "#111111"))),
                 "not in the data")
})

# the legend keys of the annotation strips: level -> the colour actually mapped to it
.strip_keys <- function(p) {
  for (i in seq_len(8)) {
    q <- tryCatch(p[[i]], error = function(e) NULL)
    if (is.null(q) || !length(q$layers)) next
    d1 <- tryCatch(q$layers[[1]]$data, error = function(e) list())
    if (!"value" %in% names(d1)) next
    b <- ggplot2::ggplot_build(q)
    sc <- Filter(function(s) grepl("fill", s$aesthetics[1]), b$plot$scales$scales)
    return(lapply(sc, function(s) {
      br <- s$get_breaks(); br <- br[!is.na(br)]
      stats::setNames(unname(s$map(br)), as.character(br))
    }))
  }
  stop("no annotation strip")
}

test_that("an annotation level with no samples here is not given a legend key", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggnewscale")
  ps <- example_pop_structure("africa", umap = FALSE)
  m <- ps$get_meta()
  # a factor annotation carries every level it was built with, so subsetting the samples
  # does not by itself remove a level -- which is how empty keys reached the legend
  m$region <- factor(as.character(m$region),
                     levels = c(unique(as.character(m$region)), "Nowhere"))
  ps$add_meta(m)

  full <- .strip_keys(plot_region_haplotypes(ps, "7", annotations = "region"))[[1]]
  expect_setequal(names(full), unique(as.character(m$region)))
  expect_false("Nowhere" %in% names(full))

  sub <- ps$subset(region = "East Africa")
  keys <- .strip_keys(plot_region_haplotypes(sub, "7", annotations = "region"))[[1]]
  expect_equal(names(keys), "East Africa")
  expect_false("Central Africa" %in% names(keys))
})

test_that("a level the shared colour map does not cover still gets a colour", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggnewscale")
  ps <- example_pop_structure("africa", umap = FALSE)
  regs <- sort(unique(as.character(ps$get_meta()$region)))
  # a map covering one level only: the other used to be drawn as a key with no swatch
  ps$set_colors(list(region = stats::setNames("#E20134", regs[1])))
  keys <- .strip_keys(plot_region_haplotypes(ps, "7", annotations = "region"))[[1]]
  expect_setequal(names(keys), regs)
  expect_false(anyNA(keys))
  expect_equal(unname(keys[regs[1]]), "#E20134")   # the one that was set is honoured
})


# ---- gapped spacing: even's readable columns, genomic's sense of distance ------------

test_that("gapped spacing keeps one full column per SNP", {
  skip_if_not_installed("ggplot2")
  ps <- ps_for_hap()
  d <- hap_panel(plot_region_haplotypes(ps, "7", spacing = "gapped"))$data
  expect_equal(unique(round(d$xmax - d$xmin, 9)), 1)
})

test_that("a blank column is spent per gap_unit of empty genome, up to gap_max", {
  # 100 bp apart, then a 40 kb desert, then 100 bp apart again
  pos <- c(1000, 1100, 1200, 41200, 41300)
  x <- plasgenomicsutilsR:::.gapped_x(pos, gap_unit = 10000, gap_max = 10)
  expect_equal(diff(x), c(1, 1, 5, 1))          # 4 blank columns bought by the desert
  expect_true(all(diff(x) >= 1))                # never overlapping, always in order

  # the cap is what stops one desert taking the panel
  capped <- plasgenomicsutilsR:::.gapped_x(c(1, 10e6), gap_unit = 10000, gap_max = 10)
  expect_equal(diff(capped), 11)

  # a smaller unit exaggerates the same gap, a larger one plays it down
  expect_gt(diff(plasgenomicsutilsR:::.gapped_x(pos, 2000, 100))[3],
            diff(plasgenomicsutilsR:::.gapped_x(pos, 20000, 100))[3])
})

test_that("the default gap_unit scales with the window, so ordinary spacing costs nothing", {
  # SNPs evenly spread over the window: none of the gaps is unusual, so none buys a column
  even_spread <- seq(1, 50000, length.out = 60)
  expect_equal(diff(plasgenomicsutilsR:::.gapped_x(even_spread)), rep(1, 59))
  # the same window with one desert in it: only the desert opens up
  with_desert <- c(seq(1, 20000, length.out = 40), seq(45000, 50000, length.out = 20))
  d <- diff(plasgenomicsutilsR:::.gapped_x(with_desert))
  expect_equal(sum(d > 1), 1L)
  expect_gt(max(d), 5)
})

test_that("gapped sits between even and genomic on the axis", {
  skip_if_not_installed("ggplot2")
  ps <- ps_for_hap()
  xr <- function(...) {
    d <- hap_panel(plot_region_haplotypes(ps, "7", ...))$data
    diff(range(c(d$xmin, d$xmax)))
  }
  even <- xr(spacing = "even")
  gapped <- xr(spacing = "gapped")
  expect_gte(gapped, even)                      # gaps add columns, never remove them
  expect_lt(gapped, xr(spacing = "genomic"))    # but the axis is still columns, not bp
})

test_that("a gene with no SNP of its own is drawn in the gap it sits in", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  ps <- ps_for_hap()
  g <- ps$genotype("full")
  pos <- as.numeric(sub(".*:", "", colnames(g)))
  # carve out every SNP in pfcrt, so it holds none
  thinned <- g[, !(pos > 403000 & pos < 407000), drop = FALSE]
  args <- list(ps, "pfcrt", pad = 30000, genotypes = thinned, genes = PF3D7_GENES)

  # under even it has no columns, so it cannot be drawn at all
  even <- hap_geometry(suppressMessages(do.call(plot_region_haplotypes,
                                                c(args, spacing = "even"))))
  expect_false("pfcrt" %in% even$genes$name)

  # under gapped the blank columns are room, and it lands between its flanking SNPs
  gap <- hap_geometry(suppressMessages(do.call(plot_region_haplotypes,
                                               c(args, spacing = "gapped"))))
  expect_true("pfcrt" %in% gap$genes$name)
  box <- gap$genes[gap$genes$name == "pfcrt", ]
  expect_gt(box$.gene_xmax, box$.gene_xmin)
})

test_that("a marked position with no SNP still lands under gapped spacing", {
  skip_if_not_installed("ggplot2")
  ps <- ps_for_hap()
  d <- hap_panel(plot_region_haplotypes(ps, "7", spacing = "gapped"))$data
  snps <- sort(unique(d$pos))
  between <- floor((snps[1] + snps[2]) / 2)
  if (between %in% snps) skip("no room between the first two SNPs")
  x <- plasgenomicsutilsR:::.marks_to_x(between, data.frame(pos = snps), 
                                        list(x = seq_along(snps)), "gapped")
  expect_length(x, 1L)
  expect_gt(x, 1); expect_lt(x, 2)              # between the first two columns, not dropped
})


# ---- nested splits: several metadata columns, each in its own level order ------------

test_that("split takes several columns and nests them in the order given", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  ps <- example_pop_structure("africa", umap = FALSE)
  m <- ps$get_meta()
  # a geographic order that is neither alphabetical nor the order the rows arrive in, so
  # the test can tell factor levels from either fallback
  regs <- rev(sort(unique(as.character(m$region))))
  m$region <- factor(as.character(m$region), levels = regs)
  m$half <- ifelse(seq_len(nrow(m)) %% 2 == 0, "odd", "even")   # "odd" sorts after "even"
  ps$add_meta(m)

  p <- plot_region_haplotypes(ps, "7", split = c("region", "half"))
  hm <- hap_panel(p)
  d <- hm$data[order(hm$data$.row), ]
  d <- d[!duplicated(d$sample), ]

  # the outer column's factor order is the block order, and the inner divides each block
  expect_identical(levels(d$.split1), regs)
  expect_identical(as.character(unique(d$.split1)), regs)
  expect_identical(levels(d$.split2), c("even", "odd"))
  within <- lapply(split(as.character(d$.split2), d$.split1), unique)
  for (w in within) expect_identical(w, c("even", "odd"))
  # every combination is one contiguous block, none of them drawn twice
  runs <- rle(as.character(d$.split))
  expect_length(runs$values, nlevels(d$.split))
  expect_identical(runs$values, levels(d$.split))
  expect_identical(levels(d$.split),
                   as.vector(t(outer(regs, c("even", "odd"), paste, sep = " / "))))

  # one panel per combination, and the dendrogram facets identically so the leaves line up
  hb <- ggplot2::ggplot_build(hm); db <- ggplot2::ggplot_build(p[[1]])
  expect_length(hb$layout$panel_params, nlevels(d$.split))
  expect_length(db$layout$panel_params, nlevels(d$.split))
  for (i in seq_along(hb$layout$panel_params))
    expect_equal(db$layout$panel_params[[i]]$y.range,
                 hb$layout$panel_params[[i]]$y.range, tolerance = 1e-6)
  # the strips name each column separately rather than pasting the two together
  expect_length(hb$layout$facet$params$rows, 2L)

  # a combination with no samples is not drawn as an empty block
  m2 <- m; m2$half[m2$region == regs[1]] <- "even"
  ps$add_meta(m2)
  q <- hap_panel(plot_region_haplotypes(ps, "7", split = c("region", "half")))
  expect_false(paste(regs[1], "odd", sep = " / ") %in% levels(q$data$.split))
  expect_true(paste(regs[1], "even", sep = " / ") %in% levels(q$data$.split))
})

test_that("nested splits keep the annotations and the block names in step", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("ggnewscale")
  ps <- example_pop_structure("africa", umap = FALSE)
  p <- plot_region_haplotypes(ps, "7", split = c("region", "country"),
                              annotations = c("region", "country"))
  hm <- hap_panel(p)
  n <- length(ggplot2::ggplot_build(hm)$layout$panel_params)
  expect_gt(n, 1L)
  for (i in seq_len(6)) {
    q <- tryCatch(p[[i]], error = function(e) NULL)
    if (is.null(q) || is.null(q$facet) || inherits(q$facet, "FacetNull")) next
    expect_length(ggplot2::ggplot_build(q)$layout$panel_params, n)
  }
  # a sample missing either column is dropped, and the message names both
  m <- ps$get_meta(); m$country[1] <- NA; ps$add_meta(m)
  expect_message(plot_region_haplotypes(ps, "7", split = c("region", "country")),
                 "1 sample\\(s\\) with no region / country")
  expect_error(plot_region_haplotypes(ps, "7", split = c("region", "nope")),
               "`split = \"nope\"` is not a metadata column")
})


.set_bcf <- function(dir, gts, alt = "T,A", chrom = "Pf3D7_13_v3", pos = 1725592,
                     samps = NULL) {
  if (is.null(samps)) samps <- sprintf("s%02d", seq_along(gts))
  hdr <- c("##fileformat=VCFv4.2", sprintf("##contig=<ID=%s,length=2000000>", chrom),
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
           paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
                   "FORMAT", samps), collapse = "\t"),
           paste(c(chrom, pos, ".", "C", alt, ".", ".", ".", "GT", gts), collapse = "\t"))
  v <- file.path(dir, paste0("m", pos, ".vcf"))
  writeLines(hdr, v)
  v
}

test_that("a call is named for the alleles it carries, and a biallelic one keeps its words", {
  nm <- plasgenomicsutilsR:::.allele_set_name
  # two alleles: exactly the wording the plot has always used, so adding a multiallelic
  # marker beside biallelic ones does not rename the calls they were already showing
  expect_equal(nm(0L, 2L), "reference")
  expect_equal(nm(1L, 2L), "alternate")
  expect_equal(nm(0:1, 2L), "mixed")
  # three: the alternates are told apart, and so are the mixtures between them
  expect_equal(nm(0L, 3L), "reference")
  expect_equal(nm(1L, 3L), "alternate 1")
  expect_equal(nm(2L, 3L), "alternate 2")
  expect_equal(nm(0:1, 3L), "reference + alternate 1")
  expect_equal(nm(c(1L, 2L), 3L), "alternate 1 + alternate 2")
  expect_equal(nm(0:2, 3L), "reference + alternate 1 + alternate 2")
})

test_that("allele sets are read from the calls, and only the states that occur are kept", {
  skip_if_not(nzchar(Sys.which("bcftools")))
  d <- tempfile(); dir.create(d)
  v <- .set_bcf(d, c(rep("0/0", 5), rep("1/1", 3), rep("2/2", 2), "0/1", "1/2", "./."))
  got <- plasgenomicsutilsR:::.read_genotype_sets(v)

  expect_equal(colnames(got$codes), "Pf3D7_13_v3:1725591")     # 0-based, as ids are here
  states <- got$levels[[1]][got$codes[, 1] + 1L]
  expect_equal(as.integer(table(states)[c("reference", "alternate 1", "alternate 2")]),
               c(5L, 3L, 2L))
  expect_true("reference + alternate 1" %in% states)
  expect_true("alternate 1 + alternate 2" %in% states)
  expect_true(is.na(states[length(states)]))                   # ./. stays missing
  # a triallelic site has seven possible states; only the five seen are levels
  expect_length(got$levels[[1]], 5L)
  expect_false("reference + alternate 2" %in% got$levels[[1]])
})

test_that("additional_genotypes puts a multiallelic marker in the heatmap and the legend", {
  skip_if_not(nzchar(Sys.which("bcftools")))
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("SNPRelate")
  ps <- example_pop_structure(umap = FALSE)
  ids <- colnames(ps$genotype(prefer = "full"))
  loc <- plasgenomicsutilsR:::.parse_snp_ids(ids)
  chrom <- loc$chr[1]
  samps <- rownames(ps$genotype(prefer = "full"))
  gts <- rep(c("0/0", "1/1", "2/2"), length.out = length(samps))
  d <- tempfile(); dir.create(d)
  # a position inside the window but not already genotyped
  free <- setdiff(seq(min(loc$pos), min(loc$pos) + 200), loc$pos)[1]
  v <- .set_bcf(d, gts, chrom = chrom, pos = free + 1L, samps = samps)

  p <- plot_region_haplotypes(ps, sprintf("%s:%d-%d", chrom, min(loc$pos), max(loc$pos)),
                              additional_genotypes = v, cluster = FALSE,
                              gene_track = FALSE)
  hm <- if (inherits(p, "patchwork")) p[[1]] else p
  keys <- hm$scales$scales[[which(vapply(hm$scales$scales, function(z)
    "fill" %in% z$aesthetics, logical(1)))[1]]]$limits
  expect_true(all(c("reference", "mixed", "alternate") %in% keys))
  expect_true(all(c("alternate 1", "alternate 2") %in% keys))
  # no key for a state nobody has: this marker has no mixed calls at all
  expect_false(any(grepl("\\+", keys)))
})

test_that("the extra call colours stay clear of the three already in use", {
  base <- plasgenomicsutilsR:::.GENO_FILL
  got <- plasgenomicsutilsR:::.distinct_fills(base, 3)
  expect_length(got, 3L)
  expect_length(intersect(got, unname(base)), 0L)
  # taking the next entries off the palette hands back an orange that sits beside the
  # existing `alternate`; picked on worst-case CIEDE2000 they stay apart under every
  # dichromacy the package checks
  d <- colour_blind_distance(c(unname(base), got))
  expect_true(all(d > 10, na.rm = TRUE))
})

test_that("extra samples in the marker's callset are left out, and said so", {
  skip_if_not(nzchar(Sys.which("bcftools")))
  skip_if_not_installed("SNPRelate")
  # a left join on the genotypes being plotted. The marker is normally called on the whole
  # cohort while the figure shows a subset, so extras are ordinary -- but a name mismatch
  # that drops most of the callset looks exactly the same, hence the count.
  ps <- example_pop_structure(umap = FALSE)
  loc <- plasgenomicsutilsR:::.parse_snp_ids(colnames(ps$genotype(prefer = "full")))
  samps <- rownames(ps$genotype(prefer = "full"))
  free <- setdiff(seq(min(loc$pos), min(loc$pos) + 300), loc$pos)[1]
  d <- tempfile(); dir.create(d)
  v <- .set_bcf(d, rep(c("0/0", "1/1", "2/2"), length.out = length(samps) + 5L),
                chrom = loc$chr[1], pos = free + 1L, samps = c(samps, paste0("ghost", 1:5)))
  region <- sprintf("%s:%d-%d", loc$chr[1], min(loc$pos), max(loc$pos))
  expect_message(p <- plot_region_haplotypes(ps, region, additional_genotypes = v,
                                             cluster = FALSE, gene_track = FALSE),
                 "5 sample\\(s\\) not in the genotypes")
  hm <- if (inherits(p, "patchwork")) p[[1]] else p
  expect_setequal(unique(hm$data$sample), samps)      # and none of the extras got in
})

test_that("additional_genotypes refuses what it cannot place", {
  skip_if_not(nzchar(Sys.which("bcftools")))
  skip_if_not_installed("SNPRelate")
  ps <- example_pop_structure(umap = FALSE)
  ids <- colnames(ps$genotype(prefer = "full"))
  loc <- plasgenomicsutilsR:::.parse_snp_ids(ids)
  samps <- rownames(ps$genotype(prefer = "full"))
  d <- tempfile(); dir.create(d)

  # a position already genotyped is two answers for one column
  same <- .set_bcf(d, rep("0/0", length(samps)), chrom = loc$chr[1],
                   pos = loc$pos[1] + 1L, samps = samps)
  expect_error(plot_region_haplotypes(ps, ids[1], additional_genotypes = same),
               "already in the genotypes")
  # and a callset missing samples cannot fill the column
  d2 <- tempfile(); dir.create(d2)
  few <- .set_bcf(d2, rep("0/0", 3), chrom = loc$chr[1], pos = loc$pos[1] + 5L,
                  samps = samps[1:3])
  expect_error(plot_region_haplotypes(ps, ids[1], additional_genotypes = few), "are not in")
})

test_that(".geno_calls never lets dosage arithmetic touch a nominal state column", {
  # `idx <- 3L - v` is a *dosage* transform. An additional_genotypes column's `v` is a
  # nominal index into that marker's own states, so under `allele = "ref"` a 4-state column
  # produced a 0 subscript (R drops it silently: 3 labels for 4 inputs) and a 5-state one
  # produced a negative subscript, which errors outright. Five states is ordinary -- the
  # test above this one already asserts one occurs.
  lv4 <- c("reference", "alternate 1", "alternate 2", "alternate 1 + alternate 2")
  lv5 <- c(lv4, "reference + alternate 1")
  for (lv in list(lv4, lv5)) {
    v <- seq_along(lv) - 1L
    id <- rep("c1:100", length(v))
    sl <- stats::setNames(list(lv), "c1:100")
    for (al in c("alt", "ref")) {
      got <- plasgenomicsutilsR:::.geno_calls(v, al, id, sl)
      expect_length(got, length(v))
      expect_equal(as.character(got), lv)
    }
  }
})

test_that(".geno_calls still reads a plain dosage column both ways round", {
  # the state-level branch must not change what a normal column does
  v <- c(0L, 1L, 2L, NA_integer_)
  expect_equal(as.character(plasgenomicsutilsR:::.geno_calls(v, "alt")),
               c("reference", "mixed", "alternate", NA))
  expect_equal(as.character(plasgenomicsutilsR:::.geno_calls(v, "ref")),
               c("alternate", "mixed", "reference", NA))
})

test_that(".geno_calls handles a state column sitting beside dosage columns", {
  # the failure mode that made this silent: the shortened vector only misaligns the calls
  # that come *after* the state column, so a marker at the end of the window looked fine
  v  <- c(0L, 2L, 3L, 0L, 2L)
  id <- c("c1:100", "c1:100", "c1:100", "c1:200", "c1:200")
  sl <- stats::setNames(list(c("reference", "alternate 1", "alternate 2")), "c1:100")
  got <- plasgenomicsutilsR:::.geno_calls(v, "ref", id, sl)
  expect_length(got, 5L)
  # v = 3 is beyond that marker's three states, so it is missing, not a silent drop
  expect_equal(as.character(got),
               c("reference", "alternate 2", NA, "alternate", "reference"))
})

test_that("the fill palette is short rather than recycled when it runs out", {
  # `.distinct_fills`' own comment says "a palette that ran out is better short than
  # recycled into a duplicate", and the caller then recycled with `rep(length.out=)`. Two
  # distinct allele states sharing a fill is exactly the confusion the greedy CIEDE2000 pick
  # exists to avoid, and it fails silently: the plot looks fine.
  got <- plasgenomicsutilsR:::.distinct_fills(plasgenomicsutilsR:::.GENO_FILL, 20L)
  expect_lt(length(got), 20L)                    # the palette really does run out
  expect_equal(length(unique(got)), length(got))
})

test_that("a marker with more states than colours warns instead of duplicating a fill", {
  skip_if_not_installed("ggplot2")
  fills <- plasgenomicsutilsR:::.GENO_FILL
  extra <- paste("state", 1:20)
  expect_warning(
    out <- plasgenomicsutilsR:::.assign_extra_fills(fills, extra),
    "colours")
  used <- out[extra]
  used <- used[!is.na(used)]
  expect_equal(length(unique(used)), length(used), info = "no fill is used twice")
})

test_that("nominal state codes are not treated as an ordered scale when clustering", {
  # `alternate 2` is not "twice as far from reference as alternate 1" -- the codes are an
  # arbitrary sorted index. Euclidean distance on them also lets one triallelic marker carry
  # up to 4 units of distance where a biallelic SNP carries 2, so it outweighs several SNPs
  # in the Ward ordering that decides row order.
  G <- matrix(c(0L, 0L, 2L,
                0L, 0L, 2L,
                0L, 0L, 2L), nrow = 3, byrow = TRUE,
              dimnames = list(c("a", "b", "c"), c("c1:1", "c1:2", "c1:3")))
  G["a", "c1:3"] <- 0L; G["b", "c1:3"] <- 1L; G["c", "c1:3"] <- 4L
  nominal <- "c1:3"
  d_plain <- as.matrix(stats::dist(G))
  d_nom <- as.matrix(plasgenomicsutilsR:::.geno_dist(G, nominal))
  # on the raw scale c looks 4x further from a than b does; on state identity they are equal
  expect_gt(d_plain["a", "c"], d_plain["a", "b"])
  expect_equal(unname(d_nom["a", "c"]), unname(d_nom["a", "b"]))
})

test_that("with no nominal columns the distance is the ordinary one", {
  G <- matrix(c(0L, 2L, 0L, 2L, 2L, 0L), nrow = 3,
              dimnames = list(c("a", "b", "c"), c("c1:1", "c1:2")))
  expect_equal(as.matrix(plasgenomicsutilsR:::.geno_dist(G, character(0))),
               as.matrix(stats::dist(G)))
})
