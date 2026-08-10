# Population structure: the PopStructure R6 workspace (PCA/UMAP/admixture), the
# shared colour map, variance-based PC selection, sNMF caching, and the plots.
# The heavy analysis packages are optional, so each test skips when they're absent.

test_that("example_pop_structure returns a PopStructure over the public fixture", {
  ps <- example_pop_structure(umap = FALSE)
  expect_s3_class(ps, "PopStructure")
  expect_length(ps$get_samples(), 60)
  expect_true(nrow(ps$pca_variance()) > 0)
  expect_s3_class(ps$prcomp(), "prcomp")
  expect_true(all(c("sample", "country") %in% names(ps$get_meta())))
  expect_setequal(unique(ps$get_meta()$country), c("Ghana", "Cambodia"))
})

test_that("colour maps are shared and stable across metadata levels", {
  cols <- meta_colors(data.frame(sample = 1:4, region = c("B", "A", "A", "B")))
  expect_named(cols, "region")
  expect_named(cols$region, c("A", "B"))            # sorted level order
  ps <- example_pop_structure(umap = FALSE)
  cc <- ps$get_colors()$country
  expect_setequal(names(cc), c("Cambodia", "Ghana"))
  # overrides win
  ov <- meta_colors(data.frame(sample = 1:2, region = c("A", "B")),
                    overrides = list(region = c(A = "#123456")))
  expect_equal(unname(ov$region["A"]), "#123456")
})

test_that("n_pcs_for_variance accepts proportion and percent", {
  ps <- example_pop_structure(umap = FALSE)
  p50 <- n_pcs_for_variance(ps$prcomp(), 0.5)
  expect_identical(p50, n_pcs_for_variance(ps$prcomp(), 50))   # percent == proportion
  expect_true(p50 >= 1 && p50 <= nrow(ps$pca_variance()))
  # fewer PCs for a smaller target
  expect_true(n_pcs_for_variance(ps$prcomp(), 0.1) <= p50)
})

test_that("run_umap accepts a variance fraction and the scatter plots build", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("uwot")
  ps <- example_pop_structure(umap = FALSE)
  ps$run_umap(pca_components = 0.5)        # PCs covering 50% of variance
  expect_equal(nrow(ps$umap_df()), 60)
  expect_s3_class(ps$plot_pca(colour = "country"), "ggplot")
  expect_s3_class(ps$plot_umap(colour = "country"), "ggplot")
  g <- ggplot2::ggplot_build(ps$plot_pca(colour = "country"))
  expect_true(length(g$data) >= 1)
})

test_that("subset by samples and by metadata narrows the active set", {
  ps <- example_pop_structure(umap = FALSE)
  gh <- ps$subset(country = "Ghana")
  expect_length(gh$get_samples(), 30)
  expect_setequal(unique(gh$get_meta()$country), "Ghana")
  expect_length(ps$get_samples(), 60)     # original untouched
  two <- ps$subset(samples = ps$get_samples()[1:2])
  expect_length(two$get_samples(), 2)
  expect_equal(nrow(two$pca_scores()), 2)
})

test_that("plot_umap errors without a UMAP embedding", {
  testthat::skip_if_not_installed("ggplot2")
  mat <- matrix(stats::rbinom(20 * 40, 2, 0.4), nrow = 20)
  ps <- pop_structure(mat, umap = FALSE)
  expect_null(ps$umap)
  expect_error(plot_umap(ps), "no UMAP")
})

test_that("admixture_order is a stable within-group ordering", {
  q <- matrix(c(0.9, 0.1, 0.85, 0.15, 0.2, 0.8, 0.1, 0.9), ncol = 2, byrow = TRUE)
  rownames(q) <- paste0("s", 1:4)
  meta <- data.frame(sample = paste0("s", 1:4), pop = c("A", "A", "B", "B"))
  ord <- admixture_order(q, meta = meta, group = "pop")
  expect_setequal(ord, rownames(q))
  expect_identical(admixture_order(q, meta = meta, group = "pop"), ord)   # deterministic
})

test_that("plot_admixture builds from a bare Q, with grouping and a group strip", {
  testthat::skip_if_not_installed("ggplot2")
  q <- matrix(c(0.8, 0.2, 0.3, 0.7, 0.9, 0.1, 0.4, 0.6), ncol = 2, byrow = TRUE)
  rownames(q) <- paste0("s", 1:4)
  expect_s3_class(plot_admixture(q), "ggplot")
  meta <- data.frame(sample = paste0("s", 1:4), pop = c("A", "A", "B", "B"))
  expect_s3_class(plot_admixture(q, meta = meta, group = "pop"), "ggplot")
  testthat::skip_if_not_installed("ggnewscale")
  cols <- meta_colors(meta)$pop
  p <- plot_admixture(q, meta = meta, group = "pop", group_bar = TRUE,
                      group_colours = cols)
  expect_s3_class(p, "ggplot")
})

test_that("the group strip is one colour per facet (not repeated across facets)", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("ggnewscale")
  set.seed(1)
  q <- matrix(stats::runif(60), ncol = 3); rownames(q) <- paste0("s", 1:20)
  meta <- data.frame(sample = rownames(q), region = rep(c("A", "B"), each = 10))
  p <- plot_admixture(q, meta = meta, group = "region", group_bar = TRUE, border = FALSE)
  b <- ggplot2::ggplot_build(p)
  strip <- Filter(function(d) "fill" %in% names(d) && any(abs(d$y - 1.06) < 1e-6), b$data)[[1]]
  per_panel <- tapply(strip$fill, strip$PANEL, function(f) length(unique(f)))
  expect_true(all(per_panel == 1))                    # each facet shows only its own colour
})

test_that("admixture bars are bordered by default and the border toggles off", {
  testthat::skip_if_not_installed("ggplot2")
  q <- matrix(stats::runif(2 * 10), ncol = 2); rownames(q) <- paste0("s", 1:10)
  bar_colour <- function(p) {
    d <- ggplot2::ggplot_build(p)$data[[which(vapply(ggplot2::ggplot_build(p)$data,
      function(x) "ymin" %in% names(x), logical(1)))[1]]]
    unique(d$colour)
  }
  expect_equal(bar_colour(plot_admixture(q)), "black")               # default: black borders
  expect_equal(bar_colour(plot_admixture(q, border_colour = "grey30")), "grey30")
  expect_true(all(is.na(bar_colour(plot_admixture(q, border = FALSE)))))   # toggled off
})

test_that("run_snmf caches: a second identical call reuses the project", {
  testthat::skip_if_not_installed("LEA")
  cache <- tempfile("snmf_cache_")
  dir.create(cache)
  geno <- list(genotype = readRDS(system.file("extdata",
    "pop_structure_ghana_cambodia.rds", package = "plasgenomicsutilsR"))$genotype)
  geno$sample.id <- rownames(geno$genotype)
  fit1 <- run_snmf(geno, K = 1:3, rep = 2, cpu = 1, cache_dir = cache)
  proj_files <- list.files(cache, pattern = "\\.snmfProject$")
  fit2 <- run_snmf(geno, K = 1:3, rep = 2, cpu = 1, cache_dir = cache)
  # no new project created on the cached re-run
  expect_equal(list.files(cache, pattern = "\\.snmfProject$"), proj_files)
  expect_s3_class(fit1, "snmf_fit")
  expect_true(snmf_best_k(fit2) %in% 1:3)
})

test_that("the PopStructure sNMF path recovers structure across the two countries", {
  testthat::skip_if_not_installed("LEA")
  ps <- example_pop_structure(umap = FALSE)
  ps$run_snmf(K = 1:4, rep = 3, cpu = 1)
  bk <- ps$best_k()
  expect_true(bk %in% 1:4)
  q <- ps$q(bk)
  expect_equal(nrow(q), 60)
  expect_equal(rownames(q), ps$get_samples())
  testthat::skip_if_not_installed("ggplot2")
  ord <- admixture_order(q, meta = ps$get_meta(), group = "country")
  expect_s3_class(ps$plot_admixture(K = bk, group = "country", sample_order = ord),
                  "ggplot")
})

test_that("the richer 'africa' example loads with multi-region metadata", {
  ps <- example_pop_structure("africa", umap = FALSE)
  expect_s3_class(ps, "PopStructure")
  expect_length(ps$get_samples(), 258)
  m <- ps$get_meta()
  expect_true(all(c("sample", "country", "site", "region") %in% names(m)))
  expect_setequal(unique(m$country), c("DRC", "Kenya", "Tanzania", "Uganda"))
  expect_true(length(unique(m$site)) >= 6)
  expect_false(any(grepl("-Past", m$site)))     # the "-Past" tag was dropped
})

test_that("set_levels controls the order used across metadata and plots", {
  ps <- example_pop_structure("africa", umap = FALSE)
  ord <- c("Uganda_North", "Uganda_East", "Kenya_West", "Kenya_East",
           "Tanzania_West", "Tanzania_East", "Uganda_Northeast", "DRC")
  ps$set_levels("site", ord)
  expect_identical(levels(ps$get_meta()$site), ord)
  # colours are re-ordered to the new level order (same colour per level)
  expect_identical(names(ps$get_colors()$site), ord)
})

test_that("save() / load_pop_structure() round-trips the workspace", {
  ps <- example_pop_structure(umap = FALSE)
  f <- tempfile(fileext = ".rds")
  ps$save(f)
  ps2 <- load_pop_structure(f)
  expect_s3_class(ps2, "PopStructure")
  expect_equal(ps2$get_samples(), ps$get_samples())
  expect_equal(ps2$pca_scores(), ps$pca_scores())
  expect_equal(ps2$get_colors(), ps$get_colors())
  # loading a non-PopStructure object errors clearly
  g <- tempfile(fileext = ".rds"); saveRDS(list(1, 2), g)
  expect_error(load_pop_structure(g), "does not contain a PopStructure")
})

test_that("plot_structure_figure composes UMAP + admixture in one object", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("patchwork")
  testthat::skip_if_not_installed("uwot")
  testthat::skip_if_not_installed("LEA")
  ps <- example_pop_structure(umap = TRUE)          # ghana/cambodia, fast
  ps$run_snmf(K = 1:3, rep = 2, cpu = 1)
  fig <- plot_structure_figure(ps, group = "country", orientation = "vertical")
  expect_s3_class(fig, "patchwork")
  expect_true(is.numeric(attr(fig, "width")) && is.numeric(attr(fig, "height")))
  # horizontal orientation and a custom row layout both build
  figh <- ps$plot_figure(group = "country", orientation = "horizontal",
                         rows = list("Ghana", "Cambodia"))
  expect_s3_class(figh, "patchwork")
})

test_that("load_genotypes converts a VCF and returns a genotype matrix", {
  testthat::skip_if_not_installed("SNPRelate")
  testthat::skip_if_not_installed("gdsfmt")
  set.seed(9)
  n_s <- 8; n_v <- 60
  hdr <- c("##fileformat=VCFv4.2", "##contig=<ID=chr1,length=100000>",
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
           paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
                   "FORMAT", sprintf("s%d", seq_len(n_s))), collapse = "\t"))
  gts <- c("0/0", "0/1", "1/1")
  body <- vapply(seq_len(n_v), function(i) {
    g <- sample(gts, n_s, replace = TRUE)
    paste(c("chr1", i * 100, ".", "A", "T", ".", ".", ".", "GT", g), collapse = "\t")
  }, character(1))
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c(hdr, body), vcf)
  res <- load_genotypes(vcf, gds = tempfile(fileext = ".gds"))
  expect_true(is.matrix(res$genotype))
  expect_equal(nrow(res$genotype), n_s)
  expect_length(res$sample.id, n_s)
})

test_that("group factor levels drive the order of every output", {
  testthat::skip_if_not_installed("ggplot2")
  set.seed(3)
  n <- 24
  meta <- data.frame(sample = paste0("s", 1:n),
                     site = rep(c("Delta", "Alpha", "Charlie"), each = n / 3),
                     stringsAsFactors = FALSE)
  want <- c("Delta", "Charlie", "Alpha")            # deliberately not alphabetical
  meta$site <- factor(meta$site, levels = want)
  q <- matrix(runif(n * 3), nrow = n, dimnames = list(meta$sample, paste0("K", 1:3)))
  q <- q / rowSums(q)

  facets <- function(p) as.character(ggplot2::ggplot_build(p)$layout$layout$site)
  # the group colour bar used to force alphabetical facets by carrying a character column
  expect_equal(facets(plot_admixture(q, meta$sample, meta, "site")), want)
  expect_equal(facets(plot_admixture(q, meta$sample, meta, "site", group_bar = TRUE)), want)

  # sample sweep follows the levels too
  ord <- admixture_order(q, meta$sample, meta, "site")
  site_of <- stats::setNames(as.character(meta$site), meta$sample)
  expect_equal(unique(site_of[ord]), want)

  # samples outside the declared levels are dropped, but not silently
  meta2 <- meta
  meta2$site <- factor(as.character(meta2$site), levels = c("Delta", "Charlie"))
  expect_warning(admixture_order(q, meta2$sample, meta2, "site"), "no level")
})

test_that("diff heatmap follows group levels, exactly when unclustered", {
  testthat::skip_if_not_installed("ggplot2")
  set.seed(5)
  n <- 30
  meta <- data.frame(sample = paste0("s", 1:n),
                     site = factor(rep(c("Delta", "Alpha", "Charlie"), each = n / 3),
                                   levels = c("Delta", "Charlie", "Alpha")),
                     stringsAsFactors = FALSE)
  g <- matrix(rbinom(n * 40, 2, 0.4), nrow = n, dimnames = list(meta$sample, NULL))
  pd <- pop_diff(g, group = "site", meta = meta)
  expect_equal(rownames(pop_diff_matrix(pd)), levels(meta$site))
  xo <- function(p) ggplot2::ggplot_build(p)$layout$panel_params[[1]]$x$get_labels()
  # unclustered follows the levels; clustered follows the clustering, which is the point
  expect_equal(xo(plot_diff_heatmap(pd, cluster = FALSE)), levels(meta$site))
  expect_setequal(xo(plot_diff_heatmap(pd)), levels(meta$site))
})

test_that("annotation colours follow the annotation's levels, not the axis order", {
  testthat::skip_if_not_installed("ggplot2")
  meta <- data.frame(site = c("S_a", "S_b", "S_c", "S_d"),
                     country = factor(c("Zed", "Alpha", "Zed", "Alpha"),
                                      levels = c("Zed", "Alpha")),
                     stringsAsFactors = FALSE)
  ord <- c("S_b", "S_d", "S_a", "S_c")          # as a clustering would order them
  ann <- .resolve_annotations("country", ord, "site", meta)
  expect_true(is.factor(ann$country))           # the column's factor must survive
  p <- .annotation_panel("country", ann$country, ord, NULL, 11)
  expect_equal(levels(p$data$value), c("Zed", "Alpha"))   # not the axis order Alpha, Zed

  meta$country <- c("c10", "c2", "c10", "c2")   # no factor -> natural sort
  a2 <- .resolve_annotations("country", ord, "site", meta)
  expect_equal(levels(.annotation_panel("country", a2$country, ord, NULL, 11)$data$value),
               c("c2", "c10"))
})

test_that("without a factor, group order is a natural sort", {
  expect_equal(.natural_sort(c("site10", "site2", "site1")), c("site1", "site2", "site10"))
  expect_equal(.levels_of(c("s10", "s2", "s1")), c("s1", "s2", "s10"))
  # a factor keeps its own order untouched
  expect_equal(.levels_of(factor(c("b", "a"), levels = c("b", "a"))), c("b", "a"))
  expect_equal(.natural_sort(c("chr10", "chrX", "chr2")), c("chr2", "chr10", "chrX"))
})

test_that("plot_snmf_cross_entropy reads an elbow table and marks the best K", {
  testthat::skip_if_not_installed("ggplot2")
  ce <- tibble::tibble(K = 1:5, n_runs = 3L,
                       min  = c(0.70, 0.66, 0.64, 0.645, 0.65),
                       mean = c(0.71, 0.67, 0.65, 0.655, 0.66),
                       max  = c(0.72, 0.68, 0.66, 0.665, 0.67),
                       best_run = 1L)
  p <- plot_snmf_cross_entropy(ce)                       # a table is accepted directly
  b <- ggplot2::ggplot_build(p)
  # the line follows `min` by default -- the replicate snmf_q() actually returns
  expect_equal(p$data$.y, ce$min)
  expect_equal(p$data$K[p$data$.best], 3L)               # lowest min
  labels <- vapply(b$plot$layers,
                   function(l) paste0(l$aes_params$label %||% "", collapse = ""), character(1))
  expect_true(any(grepl("best K = 3", labels)))
  # stat = "mean" follows the other column and can pick a different K
  expect_equal(plot_snmf_cross_entropy(ce, stat = "mean")$data$.y, ce$mean)
  # best_k is overridable, and NA marks none
  expect_equal(plot_snmf_cross_entropy(ce, best_k = 5)$data$K[
    plot_snmf_cross_entropy(ce, best_k = 5)$data$.best], 5L)
  expect_false(any(plot_snmf_cross_entropy(ce, best_k = NA)$data$.best))
  # the min-max band is what shows a flat / unreproducible K
  ribbons <- vapply(b$plot$layers, function(l) inherits(l$geom, "GeomRibbon"), logical(1))
  expect_true(any(ribbons))
  expect_false(any(vapply(ggplot2::ggplot_build(plot_snmf_cross_entropy(ce, show_range = FALSE))$plot$layers,
                          function(l) inherits(l$geom, "GeomRibbon"), logical(1))))
})

test_that("snmf_cross_entropy summarises replicates per K", {
  testthat::skip_if_not_installed("LEA")
  set.seed(2)
  g <- matrix(rbinom(40 * 60, 2, rep(c(0.2, 0.8), each = 20 * 60)), nrow = 40,
              dimnames = list(paste0("s", 1:40), NULL))
  fit <- run_snmf(g, K = 1:3, rep = 3, cache = FALSE)
  ce <- snmf_cross_entropy(fit)
  expect_equal(ce$K, 1:3)
  expect_true(all(ce$n_runs == 3))
  expect_true(all(ce$min <= ce$mean & ce$mean <= ce$max))
  # best_run is the replicate index snmf_q() defaults to
  k <- ce$K[which.min(ce$min)]
  expect_equal(ce$best_run[ce$K == k],
               which.min(LEA::cross.entropy(fit$project, K = k)))
})

test_that("plot_admixture_multi_k pages every K, best K first-page-marked", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("LEA")
  set.seed(4)
  n <- 40
  g <- matrix(rbinom(n * 60, 2, rep(c(0.2, 0.8), each = (n / 2) * 60)), nrow = n,
              dimnames = list(paste0("s", 1:n), NULL))
  meta <- data.frame(sample = rownames(g), pop = rep(c("A", "B"), each = n / 2),
                     stringsAsFactors = FALSE)
  ps <- PopStructure$new(g, meta = meta)
  ps$run_snmf(K = 1:4, rep = 2, cpu = 1)

  # pin best_k so the assertion does not depend on which K this tiny fit prefers
  pages <- ps$plot_admixture_multi_k(group = "pop", best_k = 3)
  # elbow first, then one page per K >= 2 (K = 1 carries no structure)
  expect_equal(names(pages), c("cross_entropy", paste0("K=", 2:4)))
  expect_true(all(vapply(pages, function(p) inherits(p, "ggplot"), logical(1))))
  expect_match(pages[["K=3"]]$labels$title, "best by")
  expect_equal(pages[["K=2"]]$labels$title, "K = 2")          # the others are plain
  expect_equal(pages$cross_entropy$labels$title, "sNMF cross-entropy by K")

  # one shared sample order keeps a sample at the same x on every page
  ords <- lapply(pages[-1], function(p) levels(p$data$sample))
  expect_length(unique(ords), 1)

  # turning it off lets each page cluster its own samples
  free <- ps$plot_admixture_multi_k(group = "pop", sample_order_best_k = FALSE)
  expect_true(length(unique(lapply(free[-1], function(p) levels(p$data$sample)))) >= 1)

  # an explicit order wins, and the elbow page can be dropped
  want <- rev(rownames(g))
  ex <- ps$plot_admixture_multi_k(group = "pop", sample_order = want,
                                  cross_entropy_first = FALSE)
  expect_false("cross_entropy" %in% names(ex))
  expect_equal(levels(ex[[1]]$data$sample), want)

  # K can be narrowed, and best_k overridden
  two <- ps$plot_admixture_multi_k(K = 3:4, best_k = 4, cross_entropy_first = FALSE)
  expect_equal(names(two), c("K=3", "K=4"))
  expect_match(two[["K=4"]]$labels$title, "best by")

  # the pages go straight into a multi-page PDF
  f <- tempfile(fileext = ".pdf")
  save_plot(f, pages)
  expect_true(file.exists(f) && file.size(f) > 0)
})

test_that("the diff heatmap has no dendrogram tip dots and keeps legends inside", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("patchwork")
  testthat::skip_if_not_installed("ggnewscale")
  set.seed(7)
  n <- 40
  meta <- data.frame(sample = paste0("s", 1:n),
                     site = rep(c("a", "b", "c", "d", "e"), each = n / 5),
                     country = rep(c("X", "Y"), each = n / 2), stringsAsFactors = FALSE)
  g <- matrix(rbinom(n * 60, 2, 0.4), nrow = n, dimnames = list(meta$sample, NULL))
  pd <- pop_diff(g, group = "site", meta = meta)

  p <- plot_diff_heatmap(pd, annotate = "country", meta = meta)
  panels <- Filter(function(z) inherits(z, "ggplot"), p$patches$plots)
  geoms <- unlist(lapply(panels, function(pl) vapply(pl$layers,
                    function(l) class(l$geom)[1], character(1))))
  # the dendrogram used to carry a coloured dot per group at each leaf
  expect_false("GeomPoint" %in% geoms)
  expect_true("GeomSegment" %in% geoms)

  # legends live on the heatmap panel (inside the empty triangle), so nothing is collected
  scale_names <- function(pl) {
    out <- lapply(pl$scales$scales, function(s) s$name)
    unlist(Filter(function(z) is.character(z) && length(z) == 1, out))
  }
  hm <- p[[length(p)]]
  expect_true("country" %in% scale_names(hm))       # annotation legend moved onto the matrix
  expect_null(p$patches$layout$guides)              # not collected to a margin column

  # the fallback keeps working, and puts them back in the margin
  m <- plot_diff_heatmap(pd, annotate = "country", meta = meta, legend_inside = FALSE)
  expect_s3_class(m, "patchwork")
  expect_equal(m$patches$layout$guides, "collect")
  expect_false("country" %in% scale_names(m[[length(m)]]))
})

test_that("inside legends are anchored top-left so a stack grows downward", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if(utils::packageVersion("ggplot2") < "3.5.0")
  th <- .legend_upper_triangle()
  # centre-anchoring let a tall stack of legends ride up over the panel above
  expect_equal(th$legend.justification.inside, c(0, 1))
  expect_equal(th$legend.position, "inside")
})

test_that("admixture clusters read K1..K15, not K1 K10 K11 K2", {
  testthat::skip_if_not_installed("ggplot2")
  q <- matrix(stats::runif(12 * 12), 12, 12,
              dimnames = list(paste0("s", 1:12), paste0("K", 1:12)))
  q <- q / rowSums(q)
  p <- plot_admixture(q)
  # reshape() leaves `cluster` a character, which ggplot would sort as K1, K10, K11, K12, K2...
  expect_equal(levels(p$data$cluster), paste0("K", 1:12))
})

test_that("legends wrap and the suggested height clears them", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("ggnewscale")
  set.seed(2)
  n <- 60; K <- 15
  q <- matrix(stats::runif(n * K), n, K,
              dimnames = list(paste0("s", 1:n), paste0("K", 1:K)))
  q <- q / rowSums(q)
  meta <- data.frame(sample = rownames(q),
                     region = rep(c("a", "b", "c", "d", "e", "f"), length.out = n),
                     stringsAsFactors = FALSE)
  h <- function(p) unname(attr(p, "plasgenomics_dims")[["height"]])
  legend_h <- function(p) .guide_box_height(ggplot2::ggplotGrob(p))

  side <- plot_admixture(q, meta = meta, group = "region", group_bar = TRUE)
  # 15 keys wrap into columns of 10, and the canvas is tall enough for the whole stack
  expect_gt(h(side), 4)
  expect_gte(h(side), legend_h(side))

  # along the bottom the legend is a layout row, so it adds to the height
  bottom <- plot_admixture(q, meta = meta, group = "region", group_bar = TRUE,
                           legend_position = "bottom")
  expect_gt(h(bottom), 4)
  expect_equal(.guide_box_height(ggplot2::ggplotGrob(bottom)), 0)   # nothing on the side

  # no legend at all -> back to the bare panel height
  expect_equal(h(plot_admixture(q, meta = meta, group = "region", group_bar = TRUE,
                                legend_position = "none")), 4)

  # legend_rows controls the wrap: fewer per column -> more columns -> a shorter legend
  wide <- plot_admixture(q, meta = meta, group = "region", group_bar = TRUE, legend_rows = 4)
  expect_lt(legend_h(wide), legend_h(side))
})

test_that("guide order is fixed, so the legends do not swap between plots", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("ggnewscale")
  q <- matrix(stats::runif(8 * 3), 8, 3, dimnames = list(paste0("s", 1:8), paste0("K", 1:3)))
  q <- q / rowSums(q)
  meta <- data.frame(sample = rownames(q), region = rep(c("a", "b"), 4),
                     stringsAsFactors = FALSE)
  orders <- replicate(4, {
    p <- plot_admixture(q, meta = meta, group = "region", group_bar = TRUE)
    paste(vapply(p$scales$scales, function(s) {
      o <- s$guide$params$order %||% s$guide$order %||% NA
      paste0(s$aesthetics[1], ":", o)
    }, character(1)), collapse = "|")
  })
  expect_length(unique(orders), 1)
})

test_that("a saved workspace is re-bound to the installed class", {
  ps <- example_pop_structure(umap = FALSE)
  ps$set_levels("country", c("Ghana", "Cambodia"))
  f <- tempfile(fileext = ".rds")
  ps$save(f)
  back <- load_pop_structure(f)
  expect_equal(back$get_samples(), ps$get_samples())
  expect_equal(levels(back$get_meta()$country), c("Ghana", "Cambodia"))
  expect_equal(dim(back$pca_scores()), dim(ps$pca_scores()))
  # an .rds written before a method existed still gets it: an R6 object serialises its own
  # methods, so the loader has to take them from the class as it stands now
  expect_true(all(names(PopStructure$public_methods) %in% ls(back)))
  expect_identical(body(back$plot_admixture),
                   body(PopStructure$public_methods$plot_admixture))
})

test_that("a PopStructure can hold both a pruned and a full SNP panel", {
  m <- matrix(c(0L, 1L, 2L, 0L, 1L, 2L), 3,
              dimnames = list(c("a", "b", "c"), c("1:10", "1:20")))
  big <- cbind(m, matrix(1L, 3, 2, dimnames = list(c("a", "b", "c"), c("1:30", "1:40"))))

  ps <- PopStructure$new(m, allele = "alt", pruned = TRUE)
  # the primary panel is named for what it is, so `$genotype("full")` finds an unpruned object
  expect_identical(ps$panels(), "pruned")
  expect_identical(PopStructure$new(m, pruned = FALSE)$panels(), "full")
  expect_identical(PopStructure$new(m)$panels(), "genotypes")

  ps$add_panel("full", big, allele = "alt", pruned = FALSE)
  expect_setequal(ps$panels(), c("pruned", "full"))
  expect_equal(ncol(ps$genotype()), 2L)              # primary is still primary
  expect_equal(ncol(ps$genotype("full")), 4L)
  expect_equal(ncol(ps$genotype(prefer = "full")), 4L)
  # `prefer` is a wish, `panel` a requirement
  expect_equal(ncol(ps$genotype(prefer = "nope")), 2L)
  expect_error(ps$genotype("nope"), "no panel called")
  # facts are per panel
  expect_true(ps$pruned("pruned"))
  expect_false(ps$pruned("full"))

  expect_setequal(PopStructure$new(m, full = big, pruned = TRUE)$panels(),
                  c("pruned", "full"))
  # a panel that does not cover the samples is refused rather than silently recycled
  expect_error(ps$add_panel("short", big[1:2, , drop = FALSE]), "missing 1 of")
  expect_error(ps$add_panel("bad", unname(big)), "needs sample row names")
})

test_that("asking for the full panel says so when only a pruned one is there", {
  m <- matrix(c(0L, 1L, 2L, 0L, 1L, 2L), 3,
              dimnames = list(c("a", "b", "c"), c("1:10", "1:20")))
  ps <- PopStructure$new(m, pruned = TRUE)
  expect_message(ps$genotype(prefer = "full"), "only holds a pruned panel")
  # said once, not on every call
  expect_silent(ps$genotype(prefer = "full"))
  # and never when the object does not claim to be pruned
  expect_silent(PopStructure$new(m)$genotype(prefer = "full"))
})

test_that("the analyses whose signal is SNP correlation read the full panel", {
  skip_if_not_installed("ggplot2")
  ps <- example_pop_structure(umap = FALSE)
  G <- ps$genotype()
  # a "full" panel with an extra SNP, so which panel was read is visible in the output
  extra <- cbind(G, matrix(G[, 1], ncol = 1,
                           dimnames = list(rownames(G), "14:12345")))
  ps$add_panel("full", extra, allele = "alt", pruned = FALSE)
  expect_equal(ncol(.geno_for(ps)), ncol(extra))
  expect_equal(ncol(.geno_for(ps, prefer = "pruned")), ncol(G))
  # an explicit matrix still wins over either panel
  expect_equal(ncol(.geno_for(ps, genotype = G[, 1:3, drop = FALSE])), 3L)
})

test_that("a saved object without panels rebuilds them from its genotypes", {
  m <- matrix(c(0L, 1L, 2L, 0L, 1L, 2L), 3,
              dimnames = list(c("a", "b", "c"), c("1:10", "1:20")))
  ps <- PopStructure$new(m)
  # simulate an object saved before panels existed: the fields simply are not there
  e <- ps$.__enclos_env__$private
  e$panel_list <- NULL
  e$primary <- NULL
  f <- tempfile(fileext = ".rds")
  saveRDS(ps, f)

  back <- load_pop_structure(f)
  # .refresh_pop_structure() copies onto a two-sample placeholder and skips NULL fields, so a
  # stale panel list would still describe that placeholder rather than these genotypes
  expect_identical(back$panels(), "genotypes")
  expect_equal(dim(back$genotype()), dim(m))
  expect_setequal(rownames(back$genotype()), rownames(m))
  expect_equal(unname(back$genotype()), unname(m))
})
