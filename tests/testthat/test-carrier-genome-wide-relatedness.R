# Background relatedness off the focal chromosome, split by carriage. The focal chromosome is
# 7; carriers share a lot on the other chromosomes, reference samples share little, so the
# fold should come out well above 1 and the rank-sum should separate them.

seg <- function(pairs, chr, start, len) {
  data.frame(sample1 = pairs$sample1, sample2 = pairs$sample2, chr = chr,
             start = start, end = start + len - 1,
             different = 0, Nsnp = 40, stringsAsFactors = FALSE)
}
prs <- function(a, b) data.frame(sample1 = a, sample2 = b, stringsAsFactors = FALSE)

CARR <- sprintf("c%02d", 1:6)
REFS <- sprintf("r%02d", 1:6)

# every within-group carrier pair shares `car_len` off chr 7, every reference pair `ref_len`;
# a couple of on-chr-7 segments are added to prove they are excluded
mk <- function(car_len = 300000, ref_len = 40000) {
  cpairs <- utils::combn(CARR, 2)
  rpairs <- utils::combn(REFS, 2)
  cp <- prs(cpairs[1, ], cpairs[2, ])
  rp <- prs(rpairs[1, ], rpairs[2, ])
  bl <- rbind(
    seg(cp, "3", 100000, car_len),               # carrier background, off focal
    seg(rp, "3", 100000, ref_len),               # reference background, off focal
    seg(cp, "7", 100000, 2000000))               # on the focal chr: must be ignored
  meta <- data.frame(sample = c(CARR, REFS), region = "north",
                     pin = rep(c("Present", "Absent"), each = 6),
                     stringsAsFactors = FALSE)
  ibd_results(blocks = bl, meta = meta, group_col_in_meta = "region")
}

test_that("carriers come out more related off the focal chromosome", {
  r <- carrier_genome_wide_relatedness(mk(), allele = "pin", focal_chr = "Pf3D7_07_v3",
                                       carrier = "Present", reference = "Absent")
  cc <- r[r$class == "carrier/carrier", ]
  rr <- r[r$class == "reference/reference", ]
  expect_equal(cc$n_pairs, choose(6, 2))        # every within-group carrier pair, sharers or not
  expect_equal(rr$n_pairs, choose(6, 2))
  # 300 kb vs 40 kb of background, the focal-chromosome 2 Mb segment excluded from both
  expect_equal(cc$mean_ibd_mb, 0.3, tolerance = 1e-6)
  expect_equal(rr$mean_ibd_mb, 0.04, tolerance = 1e-6)
  expect_equal(cc$fold_vs_ref, 7.5, tolerance = 1e-4)
  expect_true(is.na(rr$p_vs_ref))               # reference row never tests against itself
  expect_lt(cc$p_vs_ref, 0.05)
  expect_equal(attr(r, "focal_chr"), "7")
})

test_that("the fraction denominator is the callable genome minus the focal chromosome", {
  r <- carrier_genome_wide_relatedness(mk(), allele = "pin", focal_chr = "Pf3D7_07_v3",
                                       carrier = "Present", reference = "Absent")
  off <- sum(PF3D7_CORE_CHROM_LENGTHS_BP[names(PF3D7_CORE_CHROM_LENGTHS_BP) != "7"])
  expect_equal(attr(r, "off_focal_bp"), off)
  cc <- r[r$class == "carrier/carrier", ]
  expect_equal(cc$mean_ibd_frac, 0.3 * 1e6 / off, tolerance = 1e-9)
})

test_that("pairs with an end outside the carrier/reference sets are dropped, not pooled", {
  ibd <- mk()
  # a named-vector allele with a third state: those pairs must vanish from every stratum
  st <- stats::setNames(rep(c("Present", "Absent"), each = 6), c(CARR, REFS))
  st[REFS[1]] <- "Other"                         # move one reference sample to a third state
  r <- carrier_genome_wide_relatedness(ibd, allele = st, focal_chr = "7",
                                       carrier = "Present", reference = "Absent")
  rr <- r[r$class == "reference/reference", ]
  expect_equal(rr$n_pairs, choose(5, 2))         # r01 removed, five reference samples left
})

test_that("a set-valued carrier lets several alleles count against one reference", {
  ibd <- mk()
  st <- stats::setNames(c(rep("A675V", 3), rep("R561H", 3), rep("REF", 6)), c(CARR, REFS))
  r <- carrier_genome_wide_relatedness(ibd, allele = st, focal_chr = "7",
                                       carrier = c("A675V", "R561H"), reference = "REF")
  cc <- r[r$class == "carrier/carrier", ]
  # both mutant alleles pooled into carrier: all six carriers, choose(6, 2) pairs
  expect_equal(cc$n_pairs, choose(6, 2))
  expect_equal(unname(attr(r, "states")$carrier), c("A675V", "R561H"))
})

test_that("it refuses to imply a second state at a multiallelic site", {
  ibd <- mk()
  st <- stats::setNames(c(rep("A675V", 3), rep("R561H", 3), rep("REF", 6)), c(CARR, REFS))
  expect_error(
    carrier_genome_wide_relatedness(ibd, allele = st, focal_chr = "7", carrier = "A675V"),
    "cannot be inferred")
  expect_error(
    carrier_genome_wide_relatedness(ibd, allele = st, focal_chr = "7"),
    "name `carrier")
})
