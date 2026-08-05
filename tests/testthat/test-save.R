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
