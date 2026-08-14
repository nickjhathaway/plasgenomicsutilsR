# Jost's D between metadata groups: per-SNP/pairwise values, the group summary matrix,
# the marker picker, and the heatmap.

test_that("jost_d computes per-SNP pairwise D in [0, 1]", {
  ps <- example_pop_structure(umap = FALSE)              # ghana/cambodia, 2 groups
  jd <- jost_d(ps, group = "country")
  expect_s3_class(jd, "jost_d")
  expect_equal(ncol(jd$D), 1)                             # one pair for two groups
  # differentiation reads the full panel: pruning removes the differentiating SNPs
  expect_equal(nrow(jd$D), ncol(ps$genotype(prefer = "full")))
  v <- jd$D[is.finite(jd$D)]
  expect_true(all(v >= 0 & v <= 1))
  # strongly differentiated countries -> some clearly informative SNPs
  expect_true(max(v) > 0.3)
})

test_that("a fixed-difference SNP gives D near 1 and a shared SNP near 0", {
  # 20 samples, two groups; SNP 1 is fixed-different, SNP 2 is identical across groups
  G <- cbind(fixed = c(rep(0L, 10), rep(2L, 10)),
             shared = rep(c(0L, 2L), 10))
  rownames(G) <- paste0("s", 1:20)
  grp <- rep(c("A", "B"), each = 10)
  jd <- jost_d(G, group = grp)
  d <- jd$D[, 1]
  expect_gt(d[["fixed"]], 0.9)
  expect_lt(d[["shared"]], 0.1)
})

test_that("jost_d_matrix is symmetric with a zero diagonal", {
  ps <- example_pop_structure("africa", umap = FALSE)
  jd <- jost_d(ps, group = "site")
  M <- jost_d_matrix(jd, stat = "mean")
  expect_equal(dim(M), c(length(jd$groups), length(jd$groups)))
  expect_equal(M, t(M))
  expect_true(all(diag(M) == 0))
  expect_true(all(M[upper.tri(M)] >= 0))
})

test_that("top_differentiating_snps returns the requested count of real SNPs", {
  ps <- example_pop_structure("africa", umap = FALSE)
  jd <- jost_d(ps, group = "site")
  top <- top_differentiating_snps(jd, 200)
  expect_length(top, 200)
  expect_true(all(top %in% jd$snp))
  expect_false(any(duplicated(top)))
  expect_length(top_differentiating_snps(jd, 50, method = "max"), 50)
})

test_that("all three differentiation statistics compute in a sane range", {
  ps <- example_pop_structure("africa", umap = FALSE)
  for (s in c("jost_d", "gst_hedrick", "fst")) {
    pd <- ps$pop_diff(group = "site", statistic = s)
    expect_s3_class(pd, "pop_diff")
    v <- pd$D[is.finite(pd$D)]
    expect_true(all(v >= 0 & v <= 1))
    # the top-percentile summary lifts the near-zero genome-wide mean
    M_mean <- pop_diff_matrix(pd, "mean")
    M_top  <- pop_diff_matrix(pd, "top_mean", top = 0.05)
    expect_true(mean(M_top[upper.tri(M_top)]) > mean(M_mean[upper.tri(M_mean)]))
  }
})

test_that("plot_diff_heatmap builds plain and with dendrogram + annotation", {
  testthat::skip_if_not_installed("ggplot2")
  ps <- example_pop_structure("africa", umap = FALSE)
  plain <- ps$plot_diff_heatmap(group = "site", statistic = "fst",
                                dendrogram = FALSE, triangle = TRUE)
  expect_s3_class(plain, "ggplot")
  testthat::skip_if_not_installed("patchwork")
  # multiple annotations with custom colours for one of them
  full <- ps$plot_diff_heatmap(group = "site", stat = "top_mean",
                               annotate = c("country", "region"), dendrogram = TRUE,
                               annotate_colours = list(country = c(DRC = "#000000")))
  expect_s3_class(full, "patchwork")
})

test_that("pop_diff_table gathers all statistics per pair", {
  ps <- example_pop_structure("africa", umap = FALSE)
  tbl <- ps$pop_diff_table(group = "site")
  expect_s3_class(tbl, "data.frame")
  expect_equal(nrow(tbl), choose(length(ps$get_meta()$site |> unique()), 2))
  expect_true(all(c("a", "b", "n_snps",
                    "jost_d_mean", "jost_d_top_mean", "jost_d_max",
                    "gst_hedrick_mean", "gst_hedrick_max", "fst_top_mean") %in% names(tbl)))
  expect_false(any(grepl("^gst_mean$|^gst_top", names(tbl))))   # plain Nei's Gst dropped
  # top_mean >= mean for every pair/statistic
  expect_true(all(tbl$jost_d_top_mean >= tbl$jost_d_mean))
})

test_that("pop_diff_snps unpacks per-SNP values with parsed coordinates", {
  ps <- example_pop_structure("africa", umap = FALSE)
  pd <- ps$pop_diff(group = "site", statistic = "jost_d")
  sn <- pop_diff_snps(pd)
  expect_true(all(c("snp", "chr", "pos", "a", "b", "pair", "value") %in% names(sn)))
  expect_equal(nrow(sn), length(pd$snp) * ncol(pd$D))
  expect_true(all(is.finite(sn$pos[!is.na(sn$pos)])))          # chr:pos parsed to numbers
})

test_that("plot_diff_manhattan builds combined and per-pair", {
  testthat::skip_if_not_installed("ggplot2")
  ps <- example_pop_structure("africa", umap = FALSE)
  pd <- ps$pop_diff(group = "region", statistic = "jost_d")
  expect_silent(ggplot2::ggplotGrob(plot_diff_manhattan(pd)))                 # max across pairs
  expect_silent(ggplot2::ggplotGrob(plot_diff_manhattan(pd, combine = "mean")))
  pr <- pd$pairs[1, ]
  expect_silent(ggplot2::ggplotGrob(plot_diff_manhattan(pd, pair = c(pr$a, pr$b))))
  expect_error(plot_diff_manhattan(pd, pair = c("nope", "nada")), "no such pair")
})

test_that("pop_diff accepts a genotype override (full unpruned set)", {
  ps <- example_pop_structure("africa", umap = FALSE)
  G  <- ps$genotype()
  # a wider matrix (duplicated columns) stands in for an 'unpruned' set here
  wide <- cbind(G, G)
  colnames(wide) <- c(colnames(G), paste0(colnames(G), "_b"))
  pd_default <- ps$pop_diff(group = "region")
  pd_over    <- ps$pop_diff(group = "region", genotype = wide)
  expect_equal(length(pd_over$snp), ncol(wide))               # used the override matrix
  expect_gt(length(pd_over$snp), length(pd_default$snp))
  # a load_genotypes-style list is accepted too
  lst <- list(genotype = wide, sample.id = rownames(wide))
  expect_s3_class(ps$pop_diff(group = "region", genotype = lst), "pop_diff")
})

# ggplot orders guides by their `order`, and with the default 0 on all of them it falls back
# to a hash of the guide, which depends on the labels. Read the built gtable, not the scales.
diff_legend_titles <- function(p) {
  g <- if (inherits(p, "patchwork")) patchwork::patchworkGrob(p) else ggplot2::ggplotGrob(p)
  out <- list()
  for (k in grep("guide-box", g$layout$name)) {
    gb <- g$grobs[[k]]
    if (is.null(gb$layout)) next
    for (j in seq_along(gb$grobs)) {
      gr <- gb$grobs[[j]]
      if (is.null(gr$layout)) next
      ti <- which(grepl("^title", gr$layout$name))[1]
      if (is.na(ti)) next
      lab <- tryCatch(gr$grobs[[ti]]$children[[1]]$label, error = function(e) NULL)
      if (!is.null(lab) && length(lab) == 1)
        out[[length(out) + 1]] <- list(t = gb$layout$t[j], lab = gsub("\n", " ", lab))
    }
  }
  if (!length(out)) return(character(0))
  vapply(out[order(vapply(out, `[[`, 0, "t"))], `[[`, "", "lab")
}

test_that("the heatmap legends follow the order `annotate` was given", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("patchwork")
  testthat::skip_if_not_installed("ggnewscale")
  ps <- example_pop_structure("africa", umap = FALSE)
  meta <- ps$get_meta()
  meta$batch <- factor(ifelse(seq_len(nrow(meta)) %% 2 == 0, "run_A", "run_B"))
  pd <- jost_d(ps, group = "site")

  # the statistic first, then the annotations as asked for -- reversing the request
  # reverses the stack, which is what tells us the order is honoured and not hashed
  expect_equal(diff_legend_titles(
    plot_diff_heatmap(pd, annotate = c("country", "region", "batch"), meta = meta)),
    c("Jost's D (mean)", "country", "region", "batch"))
  expect_equal(diff_legend_titles(
    plot_diff_heatmap(pd, annotate = c("batch", "country"), meta = meta)),
    c("Jost's D (mean)", "batch", "country"))

  # and it does not move when the data behind the labels changes
  sites <- unique(as.character(meta$site))
  stacks <- lapply(list(sites, sites[1:4], sites[1:3]), function(keep) {
    sub <- ps$subset(site = keep)
    m <- sub$get_meta()
    m$batch <- factor(ifelse(seq_len(nrow(m)) %% 2 == 0, "run_A", "run_B"))
    diff_legend_titles(plot_diff_heatmap(jost_d(sub, group = "site"),
                                         annotate = c("country", "region", "batch"), meta = m))
  })
  expect_equal(length(unique(stacks)), 1L)

  # legends in the margin are collected panel by panel, so they follow the strip layout:
  # the annotations top to bottom, then the heatmap's own scale
  expect_equal(diff_legend_titles(
    plot_diff_heatmap(pd, annotate = c("country", "region"), meta = meta,
                      legend_inside = FALSE)),
    c("country", "region", "Jost's D (mean)"))
})

test_that("annotation strips reuse the object's shared colour map", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("patchwork")
  ps <- example_pop_structure("africa", umap = FALSE)
  regs <- sort(unique(as.character(ps$get_meta()$region)))
  mine <- stats::setNames(c("#E20134", "#003C86")[seq_along(regs)], regs)
  ps$set_colors(list(region = mine))

  fills <- function(p) {
    if (inherits(p, "patchwork")) { q <- p; class(q) <- setdiff(class(q), "patchwork"); p <- q }
    b <- ggplot2::ggplot_build(p)
    unique(unlist(lapply(b$data, function(d) c(d$colour, d$fill))))
  }
  # a level keeps the colour it has in the UMAP and the admixture bars
  expect_true(all(mine %in% fills(ps$plot_diff_heatmap(group = "site", annotate = "region"))))

  # an explicit palette still wins over the shared map
  green <- stats::setNames(rep("#00FF00", length(regs)), regs)
  got <- fills(ps$plot_diff_heatmap(group = "site", annotate = "region",
                                    annotate_colours = list(region = green)))
  expect_true("#00FF00" %in% got)
  expect_false(any(mine %in% got))
})
