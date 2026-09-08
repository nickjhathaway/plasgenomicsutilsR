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
  # the organelles collapse to the keys the bundled datasets use, so a FASTA-named table
  # (a tandem-repeat BED, say) joins PF3D7_GENES on every sequence
  expect_equal(normalise_chr(c("Pf3D7_API_v3", "API", "Pf_M76611", "Pf3D7_MIT_v3", "MIT")),
               c("API", "API", "MIT", "MIT", "MIT"))
  expect_true(all(unique(PF3D7_GENES$chrom) %in% c(1:14, "API", "MIT")))
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


test_that("the bundled gene tables record the strand each gene is read on", {
  for (nm in c("PF3D7_GENES", "PF_EXAMPLE_DRUG_GENES", "PF3D7_PARALOG_GENES")) {
    g <- get(nm)
    expect_true("strand" %in% names(g), info = nm)
    expect_setequal(unique(g$strand), c("+", "-"))
    expect_false(anyNA(g$strand), info = nm)
  }
  # coordinates stay low-to-high whatever the strand -- the column is the only orientation
  expect_true(all(PF3D7_GENES$end > PF3D7_GENES$start))
  # a known minus-strand gene, so a sign flip in the build would not pass silently
  expect_equal(PF3D7_GENES$strand[PF3D7_GENES$gene_id == "PF3D7_1343700"], "-")
  expect_equal(PF3D7_GENES$strand[PF3D7_GENES$gene_id == "PF3D7_0709000"], "+")
})

test_that("adding strand did not disturb what already reads these tables", {
  # .gene_track() selects the columns it needs by name and drops the rest, which is what
  # makes appending a column safe
  expect_equal(names(plasgenomicsutilsR:::.gene_track(PF3D7_GENES)),
               c("name", "chr", "start", "end"))
})


test_that("the paralog table's strands agree with PF3D7_GENES where the two overlap", {
  # they are joined from the same GFF, so a disagreement would mean the two datasets were
  # built from different releases -- which nothing downstream would notice
  m <- merge(PF3D7_PARALOG_GENES[, c("gene_id", "strand")],
             PF3D7_GENES[, c("gene_id", "strand")], by = "gene_id",
             suffixes = c("_par", "_genes"))
  expect_gt(nrow(m), 0)
  expect_equal(m$strand_par, m$strand_genes)
  # about a quarter are pseudogenes, so they are absent from PF3D7_GENES and still carry one
  extra <- setdiff(PF3D7_PARALOG_GENES$gene_id, PF3D7_GENES$gene_id)
  expect_gt(length(extra), 0)
  expect_false(anyNA(PF3D7_PARALOG_GENES$strand[
    PF3D7_PARALOG_GENES$gene_id %in% extra]))
})
