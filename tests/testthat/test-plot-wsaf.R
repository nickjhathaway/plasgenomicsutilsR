test_that("plot_wsaf draws one panel per sample", {
  skip_if_not_installed("ggplot2")
  set.seed(1)
  sites <- rbind(
    data.frame(sample = "dominant", minor_frac = stats::rbeta(300, 2, 40)),
    data.frame(sample = "even_mix", minor_frac = stats::rbeta(300, 40, 42) / 2))

  p <- plot_wsaf(sites)
  expect_s3_class(p, "ggplot")
  expect_setequal(levels(p$data$.panel), c("dominant", "even_mix"))
  # built without warnings: nothing silently clipped out of the histogram
  expect_silent(invisible(ggplot2::ggplotGrob(p)))
  expect_named(attr(p, "plasgenomics_dims"), c("width", "height"))
})

test_that("a profile orders and labels the panels by class", {
  skip_if_not_installed("ggplot2")
  set.seed(2)
  sites <- rbind(
    data.frame(sample = "mix", minor_frac = stats::rbeta(200, 40, 42) / 2),
    data.frame(sample = "dom", minor_frac = stats::rbeta(200, 2, 40)),
    data.frame(sample = "cln", minor_frac = stats::rbeta(200, 1, 90)))
  prof <- data.frame(
    sample = c("mix", "dom", "cln"),
    class = c("mixed", "dominant_clone", "monoclonal"),
    min_freq_needed = c(NA, 0.20, 0.02))

  p <- plot_wsaf(sites, profile = prof)
  # clonal first, mixtures last, so a cohort reads best-to-worst
  expect_equal(as.character(levels(p$data$.panel)),
               c("cln\nmonoclonal", "dom\ndominant_clone", "mix\nmixed"))
  # each class gets its own fill
  expect_equal(length(unique(p$data$.fill)), 3L)
})

test_that("samples selects and orders the panels as given", {
  skip_if_not_installed("ggplot2")
  sites <- data.frame(sample = rep(c("a", "b", "c"), each = 50),
                      minor_frac = stats::runif(150, 0, 0.5))
  p <- plot_wsaf(sites, samples = c("c", "a"))
  expect_equal(levels(p$data$.panel), c("c", "a"))
  expect_setequal(unique(p$data$sample), c("c", "a"))
})

test_that("colours override the class palette, either spelling", {
  skip_if_not_installed("ggplot2")
  sites <- data.frame(sample = rep("s", 50), minor_frac = stats::runif(50, 0, 0.5))
  prof <- data.frame(sample = "s", class = "dominant_clone")
  expect_equal(unique(plot_wsaf(sites, prof,
                                colours = c(dominant_clone = "#123456"))$data$.fill),
               "#123456")
  expect_equal(unique(plot_wsaf(sites, prof,
                                colors = c(dominant_clone = "#654321"))$data$.fill),
               "#654321")
})

test_that("the wsmaf view spans the whole range and marks both sides", {
  skip_if_not_installed("ggplot2")
  sites <- data.frame(sample = rep("s", 100),
                      minor_frac = stats::runif(100, 0, 0.5),
                      wsmaf = stats::runif(100, 0, 1))
  p <- plot_wsaf(sites, value = "wsmaf", comparable_at = 0.35)
  expect_true(max(p$data$.v) > 0.5)                       # not clipped to the minor half
  vlines <- vapply(p$layers, function(l) inherits(l$geom, "GeomVline"), logical(1))
  expect_equal(sum(vlines), 2L)                            # 0.35 and 0.65
})

test_that("plot_wsaf refuses input it cannot read", {
  skip_if_not_installed("ggplot2")
  sites <- data.frame(sample = rep("s", 20), minor_frac = stats::runif(20, 0, 0.5))
  expect_error(plot_wsaf(data.frame(minor_frac = 1:5)), "needs a `sample` column")
  expect_error(plot_wsaf(sites, value = "wsmaf"), "no `wsmaf` column")
  expect_error(plot_wsaf(sites, samples = "nope"), "none of those samples")
  expect_error(plot_wsaf(sites, profile = data.frame(sample = "s")),
               "`sample` and `class`")
  expect_error(plot_wsaf(data.frame(sample = "s", minor_frac = NA_real_)),
               "no finite")
})

test_that("scales toggles between per-panel and shared y axes", {
  skip_if_not_installed("ggplot2")
  # one sample with many het calls, one with few: only a shared axis shows the difference
  sites <- rbind(
    data.frame(sample = "many", minor_frac = stats::runif(2000, 0, 0.5)),
    data.frame(sample = "few", minor_frac = stats::runif(20, 0, 0.5)))

  expect_true(plot_wsaf(sites)$facet$params$free$y)                  # the default
  expect_false(plot_wsaf(sites, scales = "fixed")$facet$params$free$y)
  # a shared axis puts both panels on one range, so the sparse one reads as nearly empty
  built <- ggplot2::ggplot_build(plot_wsaf(sites, scales = "fixed"))
  rng <- lapply(built$layout$panel_params, function(p) range(p$y.range))
  expect_equal(rng[[1]], rng[[2]])
  expect_error(plot_wsaf(sites, scales = "free_x"), "should be one of")
})
