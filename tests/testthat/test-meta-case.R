# Metadata comes from whoever assembled it, and `Sample` / `sample` / `SAMPLE` are the same
# column to everyone except a string comparison. Every entry point that takes `meta` accepts
# any of them.

.caps_meta <- function(case = "Sample") {
  m <- data.frame(x = paste0("s", 1:8),
                  region = rep(c("north", "south"), each = 4),
                  stringsAsFactors = FALSE)
  names(m)[1] <- case
  m
}

.small_geno <- function() {
  set.seed(4)
  G <- matrix(stats::rbinom(8 * 30, 2, 0.4), 8, 30)
  dimnames(G) <- list(paste0("s", 1:8),
                      paste0("Pf3D7_07_v3:", seq(1000, by = 900, length.out = 30)))
  G
}

test_that("the sample column is found whatever its case", {
  n <- plasgenomicsutilsR:::.normalise_meta
  expect_equal(names(n(data.frame(sample = 1, region = 1))), c("sample", "region"))
  for (case in c("Sample", "SAMPLE", "SaMpLe")) {
    m <- data.frame(a = 1, region = 1); names(m)[1] <- case
    expect_message(out <- n(m), "as `sample`")
    expect_equal(names(out), c("sample", "region"))
  }
  # an exact match wins outright, so this is not ambiguous
  expect_equal(names(n(data.frame(sample = 1, Sample = 1))), c("sample", "Sample"))
  # two non-exact variants are, and guessing which holds the ids is not ours to do
  expect_error(n(data.frame(Sample = 1, SAMPLE = 1)), "differing only in case")
  # a table with no sample column is left for the caller to complain about
  expect_equal(names(n(data.frame(id = 1))), "id")
  expect_null(n(NULL))
})

test_that("PopStructure and IbdResults take a capitalised sample column", {
  G <- .small_geno()
  ps <- suppressMessages(PopStructure$new(G, meta = .caps_meta("Sample")))
  expect_true("sample" %in% names(ps$get_meta()))
  expect_setequal(ps$get_meta()$sample, rownames(G))
  # and it behaves like any other object afterwards
  expect_length(ps$subset(region = "north")$get_samples(), 4)

  pf <- data.frame(sample1 = "s1", sample2 = "s2", ibd_fraction_accessible = 0.3)
  obj <- suppressMessages(ibd_results(pair_fraction = pf, meta = .caps_meta("SAMPLE"),
                                      group_col_in_meta = "region"))
  expect_true("sample" %in% names(obj$get_meta()))

  # add_meta too
  ps2 <- PopStructure$new(G)
  suppressMessages(ps2$add_meta(.caps_meta("Sample")))
  expect_true("sample" %in% names(ps2$get_meta()))
})

test_that("every exported function taking `meta` normalises it", {
  ns <- asNamespace("plasgenomicsutilsR")
  takes_meta <- Filter(function(nm) {
    f <- get(nm, envir = ns)
    is.function(f) && "meta" %in% names(formals(f))
  }, getNamespaceExports(ns))
  # the guard: a new meta-taking function that skips the normaliser fails here
  missing <- Filter(function(nm) {
    txt <- paste(deparse(body(get(nm, envir = ns))), collapse = " ")
    !grepl(".normalise_meta(meta)", txt, fixed = TRUE)
  }, takes_meta)
  expect_equal(missing, character(0))
  expect_gt(length(takes_meta), 15)      # it really is checking a set, not an empty one
})

test_that("results are identical whichever case the metadata used", {
  G <- .small_geno()
  lower <- .caps_meta("sample")
  a <- pop_diversity(G, group = "region", meta = lower)
  b <- suppressMessages(pop_diversity(G, group = "region", meta = .caps_meta("Sample")))
  expect_equal(a, b)

  d1 <- jost_d(G, group = "region", meta = lower)
  d2 <- suppressMessages(jost_d(G, group = "region", meta = .caps_meta("SAMPLE")))
  expect_equal(d1$D, d2$D)
})


# ---- the remaining public surface -------------------------------------------------

test_that("color_palette grows through the colour-blind-friendly sets, then interpolates", {
  expect_length(color_palette(0), 0)
  for (n in c(1, 8, 12, 15, 30)) {
    p <- color_palette(n)
    expect_length(p, n)
    expect_false(anyDuplicated(p) > 0)                 # distinct colours at every size
    expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", p)))
  }
  # the small palettes are prefixes of themselves, so a level keeps its colour as a
  # cohort grows by one group
  expect_equal(color_palette(4), color_palette(8)[1:4])
})

test_that("the deprecated aliases still reach their replacements", {
  ps <- example_pop_structure("africa", umap = FALSE)
  pd <- jost_d(ps, group = "site")
  skip_if_not_installed("ggplot2")
  # plot_jost_d_heatmap() is plot_diff_heatmap() fixed to Jost's D
  a <- plot_jost_d_heatmap(pd, dendrogram = FALSE, cluster = FALSE)
  b <- plot_diff_heatmap(pd, stat = "mean", dendrogram = FALSE, cluster = FALSE)
  expect_equal(class(a), class(b))
  expect_equal(ggplot2::ggplot_build(a)$data, ggplot2::ggplot_build(b)$data)
})

test_that("ld_index reports Ia and rbarD per group", {
  ps <- example_pop_structure("africa", umap = FALSE)
  res <- ld_index(ps, group = "region", max_snps = 60)
  expect_true(all(c("group", "n_samples", "n_loci", "ia", "rbar_d") %in% names(res)))
  expect_setequal(as.character(res$group), unique(as.character(ps$get_meta()$region)))
  # rbar_d is Ia standardised to [-1, 1]; ia itself is unbounded above
  expect_true(all(res$rbar_d >= -1 & res$rbar_d <= 1, na.rm = TRUE))
  expect_true(all(res$n_loci > 0))
  # `min_samples` is the floor a group has to clear to be estimated at all
  expect_true(all(res$n_samples >= 4))
  expect_warning(none <- ld_index(ps, group = "region", max_snps = 60, min_samples = 1e6),
                 "no group had a usable")
  expect_equal(nrow(none), 0)
})
