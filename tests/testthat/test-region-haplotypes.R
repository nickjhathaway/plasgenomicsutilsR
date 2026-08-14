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
  for (sp in c("even", "genomic")) {
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
