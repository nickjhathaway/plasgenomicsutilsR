# PDF device selection, suggested sizes, and save_plot().

test_that("pdf_device returns a device function or 'pdf'", {
  d <- pdf_device()
  expect_true(is.function(d) || identical(d, "pdf"))
})

test_that("plot_dims scales with chromosomes and groups", {
  testthat::skip_if_not_installed("ggplot2")
  ex <- example_ibd_results()
  full <- plot_dims(ex, "tugofwar")
  sub  <- plot_dims(ex, "tugofwar", chroms = c("7", "8"))
  expect_named(full, c("width", "height"))
  expect_lt(sub[["width"]], full[["width"]])            # fewer chromosomes -> narrower
  expect_gt(plot_dims(ex, "heatmap")[["height"]], 0)
  expect_gt(plot_dims(ex, "triangles")[["width"]], 0)
})

test_that("plots carry a suggested size and save_plot uses it", {
  testthat::skip_if_not_installed("ggplot2")
  ex <- example_ibd_results()
  p <- plot_ibd_tugofwar(ex, group = "Tanzania")
  d <- attr(p, "plasgenomics_dims")
  expect_named(d, c("width", "height"))
  f <- tempfile(fileext = ".png")
  save_plot(f, p)
  expect_true(file.exists(f) && file.info(f)$size > 0)
})

test_that("save_plot writes a list of plots as a multi-page PDF", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("scales")
  ex <- example_ibd_results()
  pw <- ex$get_pairwise_group()
  ids <- unique(paste0(pw$chr, ":", pw$pos))[1:3]
  plots <- plot_pairwise_ibd_for_genes(ex, snps = ids, individual = TRUE)
  expect_length(plots, 3)

  out <- tempfile(fileext = ".pdf")
  expect_equal(save_plot(out, plots, device = "pdf"), out)   # base pdf -> greppable catalog
  expect_true(file.exists(out))
  raw <- readBin(out, "raw", file.info(out)$size)
  expect_gt(length(grepRaw(charToRaw("/Count 3"), raw, fixed = TRUE)), 0)   # one page per plot
  # a list of plots can only be a multi-page PDF
  expect_error(save_plot(tempfile(fileext = ".png"), plots), "\\.pdf")
  # a single plot still saves as before
  expect_silent(save_plot(tempfile(fileext = ".pdf"), plots[[1]], device = "pdf"))
})

test_that("save_plot fits the canvas to a fixed-aspect panel", {
  testthat::skip_if_not_installed("ggplot2")
  df <- data.frame(x = c(0, 4), y = c(0, 1))
  # coord_fixed with a 4:1 data range -> the panel must be 4x wider than tall
  p <- ggplot2::ggplot(df, ggplot2::aes(x, y)) + ggplot2::geom_point() +
    ggplot2::coord_fixed()
  m <- .plot_metrics(p)
  expect_equal(m$aspect, 0.25, tolerance = 1e-6)
  expect_true(m$fixed_w > 0 && m$fixed_h > 0)   # axes + margins take real inches

  # width given -> height solved so the panel exactly fills it
  got <- .fit_plot_dims(p, width = 8)
  expect_equal(got$width, 8)
  expect_equal(got$height, round(0.25 * (8 - m$fixed_w) + m$fixed_h, 2))
  # height given -> width solved the other way
  got <- .fit_plot_dims(p, height = 4)
  expect_equal(got$height, 4)
  expect_equal(got$width, round((4 - m$fixed_h) / 0.25 + m$fixed_w, 2))
  # round-trip: solving for height then feeding it back gives the same width
  h <- .fit_plot_dims(p, width = 8)$height
  expect_equal(.fit_plot_dims(p, height = h)$width, 8, tolerance = 0.05)

  # both supplied are honoured verbatim, and fit = FALSE opts out entirely
  expect_equal(.fit_plot_dims(p, width = 9, height = 9), list(width = 9, height = 9))
  attr(p, "plasgenomics_dims") <- c(width = 5, height = 12)
  expect_equal(.fit_plot_dims(p, fit = FALSE), list(width = 5, height = 12))
  # fitting overrides an attached height that would leave blank margin
  expect_lt(.fit_plot_dims(p)$height, 12)
  expect_equal(.fit_plot_dims(p)$width, 5)      # attached width is kept as the anchor
})

test_that("a free-coord plot keeps its requested / attached size", {
  testthat::skip_if_not_installed("ggplot2")
  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) + ggplot2::geom_point()
  expect_null(.plot_metrics(p)$aspect)          # nothing locks the shape
  attr(p, "plasgenomics_dims") <- c(width = 8, height = 4)
  expect_equal(.fit_plot_dims(p), list(width = 8, height = 4))
  # one side given: scale the other by the attached aspect rather than inventing a shape
  expect_equal(.fit_plot_dims(p, width = 4), list(width = 4, height = 2))
  expect_equal(.fit_plot_dims(p, height = 8), list(width = 16, height = 8))
})

test_that("multi-page output uses one page size that fits the tallest page", {
  testthat::skip_if_not_installed("ggplot2")
  sq <- function(n) ggplot2::ggplot(data.frame(x = c(0, n), y = c(0, 1)),
                                    ggplot2::aes(x, y)) +
    ggplot2::geom_point() + ggplot2::coord_fixed()
  plots <- list(sq(4), sq(1))                   # a wide page and a square one
  for (i in seq_along(plots)) attr(plots[[i]], "plasgenomics_dims") <- c(width = 6, height = 6)
  f <- tempfile(fileext = ".pdf")
  save_plot(f, plots)
  expect_true(file.exists(f) && file.size(f) > 0)
  # the square page needs the most height at width 6, so that is the page height
  tall <- .fit_plot_dims(plots[[2]], width = 6)$height
  expect_gt(tall, .fit_plot_dims(plots[[1]], width = 6)$height)
})
