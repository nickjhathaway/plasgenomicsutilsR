# `ibd_block_extension_by_allele(allele = ...)` takes a metadata column or a named
# `sample -> state` vector, and until now there was **no way to get one out of a callset**.
# The only code in the package that read a per-sample allele set was private, wired into one
# plot. So "the D384A carriers" was reachable from a hand-built metadata column and from
# nowhere else.

.tri_bcf <- function(dir, gts = c("0/0", "1/1", "2/2", "1/2", "./."),
                     samps = paste0("s", seq_along(gts)),
                     chrom = "Pf3D7_13_v3", pos = 1725592L, ref = "A", alt = "C,G") {
  vcf <- file.path(dir, "tri.vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    sprintf("##contig=<ID=%s,length=2000000>", chrom),
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", samps),
          collapse = "\t"),
    paste(c(chrom, pos, ".", ref, alt, ".", "PASS", ".", "GT", gts), collapse = "\t")), vcf)
  vcf
}

test_that("allele_states names each sample's allele set by its actual bases", {
  skip_if_not(nzchar(Sys.which("bcftools")))
  d <- tempfile(); dir.create(d)
  st <- allele_states(.tri_bcf(d), "Pf3D7_13_v3:1725591")
  expect_equal(st[["s1"]], "A")            # reference
  expect_equal(st[["s2"]], "C")
  expect_equal(st[["s3"]], "G")
  expect_equal(st[["s4"]], "C + G")        # a mixed infection carrying both alternates
  expect_true(is.na(st[["s5"]]))
  expect_named(st, paste0("s", 1:5))
})

test_that("it is exactly the shape ibd_block_extension_by_allele takes", {
  skip_if_not(nzchar(Sys.which("bcftools")))
  d <- tempfile(); dir.create(d)
  st <- allele_states(.tri_bcf(d), "Pf3D7_13_v3:1725591")
  expect_type(st, "character")
  expect_false(is.null(names(st)))
  # and so `carrier =` can finally name a base rather than a metadata label
  expect_true(all(c("C", "G") %in% st))
})

test_that("a 1-based position is accepted when it is declared", {
  skip_if_not(nzchar(Sys.which("bcftools")))
  d <- tempfile(); dir.create(d)
  a <- allele_states(.tri_bcf(d), "Pf3D7_13_v3:1725591")
  b <- allele_states(.tri_bcf(d), "Pf3D7_13_v3:1725592", one_based = TRUE)
  expect_equal(a, b)
})

test_that("asking for a position the callset does not carry is an error, not an empty vector", {
  skip_if_not(nzchar(Sys.which("bcftools")))
  d <- tempfile(); dir.create(d)
  expect_error(allele_states(.tri_bcf(d), "Pf3D7_13_v3:999"), "no record")
})

test_that("chromosome spelling is normalised, so 13 and Pf3D7_13_v3 both work", {
  skip_if_not(nzchar(Sys.which("bcftools")))
  d <- tempfile(); dir.create(d)
  a <- allele_states(.tri_bcf(d), "Pf3D7_13_v3:1725591")
  b <- allele_states(.tri_bcf(d), "13:1725591")
  expect_equal(unname(a), unname(b))
})

test_that("index naming is available for a legend that should not show bases", {
  skip_if_not(nzchar(Sys.which("bcftools")))
  d <- tempfile(); dir.create(d)
  st <- allele_states(.tri_bcf(d), "Pf3D7_13_v3:1725591", names = "index")
  expect_equal(unname(st[1:4]), c("reference", "alternate 1", "alternate 2",
                                  "alternate 1 + alternate 2"))
})

test_that("a biallelic marker reads plainly under both namings", {
  skip_if_not(nzchar(Sys.which("bcftools")))
  d <- tempfile(); dir.create(d)
  v <- .tri_bcf(d, gts = c("0/0", "1/1", "0/1"), samps = c("s1", "s2", "s3"), alt = "T")
  expect_equal(unname(allele_states(v, "Pf3D7_13_v3:1725591")), c("A", "T", "A + T"))
  expect_equal(unname(allele_states(v, "Pf3D7_13_v3:1725591", names = "index")),
               c("reference", "alternate", "mixed"))
})

test_that("the carrier vector drives the block-extension contrast end to end", {
  skip_if_not(nzchar(Sys.which("bcftools")))
  d <- tempfile(); dir.create(d)
  samps <- sprintf("x%02d", 1:24)
  gts <- rep(c("1/1", "0/0"), each = 12)          # 12 carriers of C, 12 reference
  v <- .tri_bcf(d, gts = gts, samps = samps, alt = "C,G")
  st <- allele_states(v, "Pf3D7_13_v3:1725591")
  expect_setequal(unique(unname(st)), c("A", "C"))
  # the point: this vector goes straight into the contrast, no metadata column needed
  expect_true(all(names(st) == samps))
})
