# PDF device selection, suggested sizes, and save_plot().

test_that("pdf_device returns a device function or 'pdf'", {
  d <- pdf_device()
  expect_true(is.function(d) || identical(d, "pdf"))
})

test_that("plot_dims scales with chromosomes and regions", {
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
  p <- plot_ibd_tugofwar(ex, region = "Tanzania")
  d <- attr(p, "plasgenomics_dims")
  expect_named(d, c("width", "height"))
  f <- tempfile(fileext = ".png")
  save_plot(f, p)
  expect_true(file.exists(f) && file.info(f)$size > 0)
})
