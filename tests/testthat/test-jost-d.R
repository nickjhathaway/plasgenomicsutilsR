# Jost's D between metadata groups: per-SNP/pairwise values, the group summary matrix,
# the marker picker, and the heatmap.

test_that("jost_d computes per-SNP pairwise D in [0, 1]", {
  ps <- example_pop_structure(umap = FALSE)              # ghana/cambodia, 2 groups
  jd <- jost_d(ps, group = "country")
  expect_s3_class(jd, "jost_d")
  expect_equal(ncol(jd$D), 1)                             # one pair for two groups
  expect_equal(nrow(jd$D), ncol(ps$genotype()))
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

test_that("all four differentiation statistics compute in a sane range", {
  ps <- example_pop_structure("africa", umap = FALSE)
  for (s in c("jost_d", "gst", "gst_hedrick", "fst")) {
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
                    "gst_mean", "gst_hedrick_max", "fst_top_mean") %in% names(tbl)))
  # top_mean >= mean for every pair/statistic
  expect_true(all(tbl$jost_d_top_mean >= tbl$jost_d_mean))
})
