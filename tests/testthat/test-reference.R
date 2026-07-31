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
