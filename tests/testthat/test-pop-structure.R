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

test_that("run_ld_prune converts a VCF and returns a genotype matrix", {
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
  res <- run_ld_prune(vcf, gds = tempfile(fileext = ".gds"))
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
