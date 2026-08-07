test_that("reference registry returns Pf3D7 facts", {
  ref <- get_reference("pf3d7")
  expect_equal(ref$species, "Plasmodium falciparum")
  expect_length(ref$core_chrom_lengths_bp, 14)
  expect_equal(unname(ref$bp_per_cm), 15000)
  expect_true("pf3d7" %in% available_references())
  expect_error(get_reference("nope"))
})

test_that("normalise_chr strips assorted spellings", {
  expect_equal(normalise_chr("Pf3D7_07_v3"), "7")
  expect_equal(normalise_chr("chr07"), "7")
  expect_equal(normalise_chr("14"), "14")
  expect_equal(normalise_chr(3), "3")
  expect_equal(normalise_chr(c("Pf3D7_01_v3", "chr2")), c("1", "2"))
})

test_that("_pkgdown.yml indexes every vignette and every export", {
  # pkgdown errors the site build on an unindexed vignette, and that failure only surfaces
  # when the tag fires the deploy -- long after the change that caused it
  skip_if_not_installed("yaml")
  yml <- file.path("..", "..", "_pkgdown.yml")
  skip_if_not(file.exists(yml))
  cfg <- yaml::read_yaml(yml)

  listed <- unlist(lapply(cfg$articles, `[[`, "contents"), use.names = FALSE)
  present <- tools::file_path_sans_ext(
    list.files(file.path("..", "..", "vignettes"), pattern = "[.]Rmd$"))
  expect_setequal(present, listed)

  ref <- unlist(lapply(cfg$reference, `[[`, "contents"), use.names = FALSE)
  ref <- sub("^starts_with\\(\"(.*)\"\\)$", "\\1", trimws(ref))
  exports <- getNamespaceExports("plasgenomicsutilsR")
  missing <- exports[!vapply(exports, function(f)
    any(f == ref) || any(startsWith(f, ref[!ref %in% exports])), logical(1))]
  expect_equal(missing, character(0))
})
