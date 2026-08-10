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

test_that("output is still written when cairo claims to exist but cannot open", {
  skip_if_not_installed("ggplot2")
  # The condition that broke macOS R-devel CI: `capabilities("cairo")` reports TRUE because
  # R was built with cairo, but the shared object's X11 dependencies are missing, so opening
  # the device warns ("failed to load cairo DLL") and writes nothing. `pdf_device()` must
  # judge cairo by opening one, not by the build-time capability.
  local_mocked_bindings(.cairo_pdf_works = function() FALSE)
  expect_identical(pdf_device(), "pdf")

  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) + ggplot2::geom_point()
  single <- tempfile(fileext = ".pdf")
  multi <- tempfile(fileext = ".pdf")
  on.exit(unlink(c(single, multi)), add = TRUE)

  suppressMessages(save_plot(single, p))
  expect_true(file.exists(single) && file.size(single) > 0)

  suppressMessages(save_plot(multi, list(p, p)))
  expect_true(file.exists(multi) && file.size(multi) > 0)
})

test_that("sizing a plot leaves no device open and writes no stray Rplots.pdf", {
  skip_if_not_installed("ggplot2")
  # Measuring a gtable asks the device for font metrics; with no device open R opens the
  # default one, which in a script is `pdf` (leaving Rplots.pdf behind) and in a notebook
  # chunk is the chunk's own device (surfacing as an empty figure under the block).
  dir <- tempfile("devcheck"); dir.create(dir)
  old <- setwd(dir); on.exit(setwd(old), add = TRUE)
  before <- length(grDevices::dev.list())

  x <- example_ibd_results()
  p <- plot_selection_manhattan(x)
  m <- .plot_metrics(p)
  expect_true(is.finite(m$fixed_w))
  expect_length(grDevices::dev.list(), before)
  expect_false(file.exists("Rplots.pdf"))

  save_plot("out.pdf", p)
  expect_true(file.exists("out.pdf"))
  expect_length(grDevices::dev.list(), before)
  expect_false(file.exists("Rplots.pdf"))
})

test_that("plot_admixture sizes itself without touching the current device", {
  skip_if_not_installed("ggplot2")
  dir <- tempfile("devcheck2"); dir.create(dir)
  old <- setwd(dir); on.exit(setwd(old), add = TRUE)
  before <- length(grDevices::dev.list())

  q <- matrix(c(0.7, 0.3, 0.2, 0.8, 0.5, 0.5), ncol = 2, byrow = TRUE,
              dimnames = list(c("a", "b", "c"), NULL))
  p <- plot_admixture(q)
  expect_s3_class(p, "ggplot")
  expect_length(grDevices::dev.list(), before)
  expect_false(file.exists("Rplots.pdf"))
})
