# `private$panel_list[[private$pick_panel(panel)]]$x` is a trap, and four accessors fell into
# it. R evaluates `private$panel_list` **before** the index expression, and `pick_panel()` is
# what calls `ensure_panels()` to build that list -- so on an object whose panels have not been
# materialised yet, the first call reads the list as it was (empty) and returns the `%||%`
# default or NULL. The second call, the list now built, returns the truth.
#
# `$encoding()` is the dangerous one: it returns `"dosage"` for an allele-index panel, and it
# is the field every guard keys on. An object that answers "dosage" once can launder an index
# panel into a dosage analysis. Measured on a round-tripped object before the fix:
#
#     call 1: encoding = dosage       | sites = NULL | pruned = NULL
#     call 2: encoding = allele_index | sites = 5    | pruned = FALSE
#
# Objects reach that state through `$save()` / `load_pop_structure()`, which is why
# `.refresh_pop_structure()` exists -- it clears the materialised panels on purpose.

.index_ps <- function() {
  skip_if_not_installed("SeqArray")
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=Pf3D7_01_v3,length=640851>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4",
    "Pf3D7_01_v3\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t1/1\t0/0\t1/1",
    "Pf3D7_01_v3\t400\t.\tC\tT,G\t.\tPASS\t.\tGT\t0/0\t1/1\t2/2\t1/1",
    "Pf3D7_01_v3\t500\t.\tG\tA\t.\tPASS\t.\tGT\t1/1\t1/1\t1/1\t0/0"), vcf)
  g <- load_genotypes(vcf, gds = tempfile(fileext = ".gds"), prune = FALSE,
                      variants = "all", encoding = "allele_index")
  PopStructure$new(g, meta = data.frame(sample = g$sample.id, country = "X"))
}

.round_tripped <- function() {
  ps <- .index_ps()
  f <- tempfile(fileext = ".rds")
  ps$save(f, compress = FALSE)
  load_pop_structure(f)
}

test_that("encoding is right on the FIRST call, not only the second", {
  ps <- .round_tripped()
  first <- ps$encoding()
  second <- ps$encoding()
  expect_equal(first, "allele_index")
  expect_equal(first, second)
})

test_that("sites, allele and pruned answer on the first call too", {
  ps <- .round_tripped()
  expect_equal(nrow(ps$sites()), 3L)
  expect_equal(ps$allele(), "alt")
  expect_false(ps$pruned())
})

test_that("what the first call says matches what the first call hands back", {
  # the consequence, not the mechanism: every dosage guard in the package keys off
  # `$encoding()`, and a matrix of allele indices described as "dosage" is read as
  # alternate copy numbers -- allele 2 becomes "two copies of the alternate"
  ps <- .round_tripped()
  enc <- ps$encoding()
  g <- ps$genotype()
  expect_equal(enc, "allele_index")
  expect_true(max(g, na.rm = TRUE) > 2 || identical(enc, "allele_index"))
  # and the derived dosage panel is reached by asking for it, not by being told a lie
  d <- ps$genotype(needs = "dosage")
  expect_equal(ps$encoding("biallelic_dosage"), "dosage")
  expect_true(all(d %in% c(0L, 2L, NA_integer_)))
})
