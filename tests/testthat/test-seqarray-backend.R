# The genotype backend, and the reason it had to change.
#
# SNPRelate's SNP-GDS holds a *dosage*: how many copies of one allele a sample has. That is
# two numbers where a multiallelic site needs k, so `biallelic.only` refused to read those
# records at all (measured: 0 of 200 on a real Pf7 slice) and `copy.num.of.ref` read them and
# gave every alternate the same number.
#
# SeqArray's SeqVarGDS holds allele *indices*, so a third allele is representable rather than
# approximated. These tests pin the three encodings that fall out of it, and pin that a
# biallelic panel still comes back exactly as it did.

.mixed_vcf <- function() {
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=Pf3D7_01_v3,length=640851>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4",
    "Pf3D7_01_v3\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t1/1\t0/0\t1/1",   # biallelic SNV
    "Pf3D7_01_v3\t200\t.\tT\t.\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0",   # no ALT
    "Pf3D7_01_v3\t300\t.\tAT\tA\t.\tPASS\t.\tGT\t0/0\t1/1\t0/0\t0/0",  # indel
    "Pf3D7_01_v3\t400\t.\tC\tT,G\t.\tPASS\t.\tGT\t0/0\t1/1\t2/2\t1/1", # multiallelic SNV
    "Pf3D7_01_v3\t500\t.\tG\tA\t.\tPASS\t.\tGT\t1/1\t1/1\t1/1\t1/1"),  # invariant
    vcf)
  vcf
}

skip_backend <- function() {
  testthat::skip_if_not_installed("SeqArray")
}

test_that("the default panel is the biallelic SNVs, coded as alt dosage, as before", {
  skip_backend()
  g <- load_genotypes(.mixed_vcf(), gds = tempfile(fileext = ".gds"), prune = FALSE)
  expect_equal(colnames(g$genotype), c("Pf3D7_01_v3:99", "Pf3D7_01_v3:499"))
  expect_equal(g$encoding, "dosage")
  expect_equal(unname(g$genotype[, "Pf3D7_01_v3:99"]), c(0L, 2L, 0L, 2L))
  expect_true(all(g$genotype[, "Pf3D7_01_v3:499"] == 2L))   # invariant, every sample alt
})

test_that("allele_index reads every record and keeps the alternates apart", {
  skip_backend()
  g <- load_genotypes(.mixed_vcf(), gds = tempfile(fileext = ".gds"), prune = FALSE,
                      variants = "all", encoding = "allele_index")
  expect_equal(g$encoding, "allele_index")
  expect_equal(ncol(g$genotype), 5L)
  # the whole point: 1/1 and 2/2 are different numbers, not the same one
  expect_equal(unname(g$genotype[, "Pf3D7_01_v3:399"]), c(0L, 1L, 2L, 1L))
  expect_equal(g$sites$alt[[which(g$sites$site_key == "Pf3D7_01_v3:399")]], c("T", "G"))
})

test_that("the old dosage reading collapses the alternates, and says so", {
  skip_backend()
  expect_warning(
    g <- load_genotypes(.mixed_vcf(), gds = tempfile(fileext = ".gds"), prune = FALSE,
                        variants = "all", encoding = "dosage"),
    "cannot say which")
  # s2 carries T and s3 carries G, and a dosage gives them the same number
  expect_equal(unname(g$genotype[c(2, 3), "Pf3D7_01_v3:399"]), c(2L, 2L))
})

test_that("a multiallelic panel refuses to be read as dosage without variants = all", {
  skip_backend()
  # the default panel has no multiallelic records in it, so there is nothing to collapse
  expect_no_warning(
    load_genotypes(.mixed_vcf(), gds = tempfile(fileext = ".gds"), prune = FALSE))
})

test_that("the skipped records are counted from the callset itself, not guessed", {
  skip_backend()
  msgs <- paste(testthat::capture_messages(
    load_genotypes(.mixed_vcf(), gds = tempfile(fileext = ".gds"), prune = FALSE)),
    collapse = " ")
  expect_match(msgs, "2 biallelic SNVs")
  expect_match(msgs, "1 multiallelic")
  expect_match(msgs, "1 indel")
  expect_match(msgs, "1 no_alt")
})

test_that("a heterozygous call is a mixed infection and reads as missing under allele_index", {
  skip_backend()
  # the package's standing convention: the parasite is haploid, so 1/2 is two clones rather
  # than a diploid genotype, and a single allele index cannot name both
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=Pf3D7_01_v3,length=640851>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3",
    "Pf3D7_01_v3\t400\t.\tC\tT,G\t.\tPASS\t.\tGT\t1/1\t1/2\t2/2"), vcf)
  g <- load_genotypes(vcf, gds = tempfile(fileext = ".gds"), prune = FALSE,
                      variants = "all", encoding = "allele_index")
  expect_equal(unname(g$genotype[, 1]), c(1L, NA, 2L))
})

test_that("an allele_index panel is refused by the dosage statistics", {
  skip_backend()
  g <- load_genotypes(.mixed_vcf(), gds = tempfile(fileext = ".gds"), prune = FALSE,
                      variants = "all", encoding = "allele_index")
  expect_error(pop_diversity(g), "allele_index")
  expect_error(pop_structure(g), "allele_index")
})

test_that("positions stay 0-based and the sites table stays aligned", {
  skip_backend()
  g <- load_genotypes(.mixed_vcf(), gds = tempfile(fileext = ".gds"), prune = FALSE,
                      variants = "all", encoding = "allele_index")
  expect_equal(g$positions, "0-based")
  expect_equal(g$sites$site_key, colnames(g$genotype))
  expect_equal(g$sites$pos, c(99, 199, 299, 399, 499))
  expect_equal(g$sites$n_alt, c(1L, 0L, 1L, 2L, 1L))
})

test_that("ref dosage still flips, and only for the dosage encoding", {
  skip_backend()
  a <- load_genotypes(.mixed_vcf(), gds = tempfile(fileext = ".gds"), prune = FALSE)
  r <- load_genotypes(.mixed_vcf(), gds = tempfile(fileext = ".gds"), prune = FALSE,
                      allele = "ref")
  expect_equal(unname(r$genotype[, 1]), 2L - unname(a$genotype[, 1]))
  # an allele index is not a count of anything, so there is nothing to flip
  expect_error(
    load_genotypes(.mixed_vcf(), gds = tempfile(fileext = ".gds"), prune = FALSE,
                   variants = "all", encoding = "allele_index", allele = "ref"),
    "allele index")
})

test_that("alt stays positionally aligned to the allele indices, star included", {
  skip_backend()
  # This one is not negotiable: an allele index names a slot in the record's own ALT list, so
  # dropping `*` from `alt` would make index 1 point at the wrong base. `A > *,T` has T at
  # index 2, and `alt[[i]][2]` has to be "T".
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=Pf3D7_01_v3,length=640851>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3",
    "Pf3D7_01_v3\t400\t.\tA\t*,T\t.\tPASS\t.\tGT\t0/0\t1/1\t2/2"), vcf)
  g <- load_genotypes(vcf, gds = tempfile(fileext = ".gds"), prune = FALSE,
                      variants = "all", encoding = "allele_index")
  expect_equal(g$sites$alt[[1]], c("*", "T"))
  expect_equal(g$sites$alt[[1]][g$genotype["s3", 1]], "T")
  # s2's call IS the `*`, and `*` is not a base -- the default blanks it (see
  # test-star-allele.R). What matters here is that blanking a call does not RENUMBER the
  # alleles beside it: T is still index 2, not index 1.
  expect_equal(unname(g$genotype[, 1]), c(0L, NA, 2L))
  g2 <- load_genotypes(vcf, gds = tempfile(fileext = ".gds"), prune = FALSE,
                       variants = "all", encoding = "allele_index", star = "allele")
  expect_equal(unname(g2$genotype[, 1]), c(0L, 1L, 2L))
})

test_that("a record carrying a spanning deletion is skipped by the biallelic panel", {
  skip_backend()
  # `A > *,T` is not a clean SNP site: part of the cohort has no base there to compare. The
  # Python package's `--snps-only` drops it for the same reason, and the two classifications
  # have to agree or a panel and the callset it came from disagree about what a SNP is.
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=Pf3D7_01_v3,length=640851>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3",
    "Pf3D7_01_v3\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t1/1\t1/1",
    "Pf3D7_01_v3\t400\t.\tA\t*,T\t.\tPASS\t.\tGT\t0/0\t2/2\t2/2"), vcf)
  msgs <- paste(testthat::capture_messages(
    g <- load_genotypes(vcf, gds = tempfile(fileext = ".gds"), prune = FALSE)),
    collapse = " ")
  expect_equal(colnames(g$genotype), "Pf3D7_01_v3:99")
  expect_match(msgs, "1 spanning_del")

  # under `variants = "all"` it is there, and the ALT column is still described faithfully
  a <- load_genotypes(vcf, gds = tempfile(fileext = ".gds"), prune = FALSE,
                      variants = "all", encoding = "allele_index")
  r <- a$sites[a$sites$site_key == "Pf3D7_01_v3:399", ]
  expect_equal(r$n_alt, 2L)          # the ALT column really does list two
  expect_equal(r$n_alt_real, 1L)     # only one of them is an allele
  expect_true(r$has_spanning_del)
})

test_that("the record classifier agrees with the companion Python package's", {
  # These are the same cases `tests/test_filter_pipeline.py`'s
  # `test_a_record_is_one_class_and_a_mixed_one_is_not_a_snp` pins, with the same answers.
  # If the two drift, a panel and the callset it was built from disagree about what a SNP is
  # -- and the disagreement shows up as records vanishing between the pipeline and R, which
  # is exactly the kind of loss this work exists to make impossible.
  f <- function(ref, alt) plasgenomicsutilsR:::.classify_alleles(ref, list(alt))
  expect_equal(f("A", "T"), "biallelic_snv")
  expect_equal(f("A", c("T", "G")), "multiallelic")
  expect_equal(f("A", c("T", "ATT")), "mixed")
  expect_equal(f("AT", "A"), "indel")
  expect_equal(f("AT", "GC"), "mnp")
  expect_equal(f("A", character(0)), "no_alt")
  # a padded substitution is one SNV, not a multi-base change -- counting it as an MNP
  # invents a population of them that is not in the data
  expect_equal(f("TTATA", "CTATA"), "biallelic_snv")
  expect_equal(f("ATCG", "GTCA"), "mnp")
  # and a `*` disqualifies the record whatever else is there
  expect_equal(f("T", "*"), "spanning_del")
  expect_equal(f("A", c("*", "T")), "spanning_del")
  expect_equal(f("A", c("T", "*")), "spanning_del")
})
