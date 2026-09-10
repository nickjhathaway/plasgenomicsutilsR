blk <- function(pairs, chr, start, len) {
  data.frame(sample1 = pairs$sample1, sample2 = pairs$sample2, chr = chr,
             start = start, end = start + len - 1,
             different = 0, Nsnp = 40, stringsAsFactors = FALSE)
}
pr2 <- function(a, b) data.frame(sample1 = a, sample2 = b, stringsAsFactors = FALSE)
LOCUS <- data.frame(name = "sweep", chr = "7", start = 400000, end = 410000,
                    stringsAsFactors = FALSE)
spread <- function(len, n) len + (seq_len(n) - (n + 1) / 2) * 1000

# `n` disjoint pairs drawn from `samples`, each sharing `locus_len` over the locus and 20 kb
# on chrs 1 and 2
share <- function(samples, locus_len) {
  n <- length(samples) / 2
  p <- pr2(samples[seq(1, length(samples), 2)], samples[seq(2, length(samples), 2)])
  len <- spread(locus_len, n)
  rbind(blk(p, "7", 405000 - len / 2, len),
        blk(p, "1", 100000, 20000), blk(p, "2", 100000, 20000))
}

# 12 carriers (c01..c12) and 12 reference samples (r01..r12) in one region
CARRIERS <- sprintf("c%02d", 1:12)
REFS <- sprintf("r%02d", 1:12)

mk_allele_ibd <- function(carrier_len = 80000, ref_len = 40000, n_car = 6, n_ref = 6,
                          extra = NULL) {
  bl <- rbind(share(CARRIERS[seq_len(2 * n_car)], carrier_len),
              share(REFS[seq_len(2 * n_ref)], ref_len))
  if (!is.null(extra)) bl <- rbind(bl, extra)
  meta <- data.frame(sample = c(CARRIERS, REFS),
                     region = "north",
                     pin = rep(c("Present", "Absent"), each = 12),
                     stringsAsFactors = FALSE)
  ibd_results(blocks = bl, meta = meta, group_col_in_meta = "region")
}

test_that("the strata are split by carriage and each gets its own extension statistic", {
  r <- ibd_block_extension_by_allele(mk_allele_ibd(), LOCUS, allele = "pin",
                                     carrier = "Present", reference = "Absent")

  expect_equal(nrow(r), 1L)
  expect_equal(r$n_carrier, 6L)
  expect_equal(r$n_reference, 6L)
  expect_equal(r$n_discordant, 0L)
  # six pairs, so the median log ratio is the geometric mean of the middle two and lands a
  # hair off the round number; the arithmetic is still exact
  expect_equal(r$ratio_carrier, 4, tolerance = 1e-4)     # 80 kb over a 20 kb baseline
  expect_equal(r$ratio_reference, 2, tolerance = 1e-4)   # 40 kb over the same
  expect_equal(r$ratio_contrast, 2, tolerance = 1e-4)
  expect_lt(r$p_length, 0.05)                   # six against six, cleanly separated
  expect_equal(unname(attr(r, "states")), c("Present", "Absent"))
  # the per-stratum detail is kept, and agrees
  st <- attr(r, "strata")
  expect_setequal(st$stratum, c("carrier", "reference"))
  expect_equal(st$paired_ratio[st$stratum == "carrier"], 4, tolerance = 1e-4)
})

test_that("the fraction statistic uses every possible pair, not only the sharing ones", {
  # 12 carriers give choose(12, 2) = 66 possible pairs, of which 6 share; same for reference
  r <- ibd_block_extension_by_allele(mk_allele_ibd(), LOCUS, allele = "pin",
                                     carrier = "Present", reference = "Absent")
  expect_equal(r$n_carrier_possible, 66)
  expect_equal(r$n_reference_possible, 66)
  expect_equal(r$frac_carrier, 6 / 66)
  expect_equal(r$frac_reference, 6 / 66)
  expect_equal(r$odds_ratio, 1, tolerance = 1e-6)
  expect_equal(r$p_fraction, 1)                 # identical fractions

  # now let the carriers share far more often than the reference samples
  more <- share(c(CARRIERS[1], CARRIERS[3], CARRIERS[5], CARRIERS[7],
                  CARRIERS[9], CARRIERS[11], CARRIERS[2], CARRIERS[4]), 80000)
  r2 <- ibd_block_extension_by_allele(mk_allele_ibd(extra = more), LOCUS, allele = "pin",
                                      carrier = "Present", reference = "Absent")
  expect_gt(r2$n_carrier, r$n_carrier)
  expect_gt(r2$frac_carrier, r2$frac_reference)
  expect_gt(r2$odds_ratio, 1)
  expect_lt(r2$p_fraction, r$p_fraction)
})

test_that("a discordant pair is counted separately and kept out of both strata", {
  # one carrier and one reference sample sharing a long segment over the locus
  dis <- share(c(CARRIERS[1], REFS[1]), 200000)
  r <- ibd_block_extension_by_allele(mk_allele_ibd(extra = dis), LOCUS, allele = "pin",
                                     carrier = "Present", reference = "Absent")
  expect_equal(r$n_discordant, 1L)
  expect_equal(r$n_carrier, 6L)                 # unchanged
  expect_equal(r$n_reference, 6L)
  expect_equal(r$ratio_carrier, 4, tolerance = 1e-4)   # the discordant pair moves neither
  expect_equal(r$ratio_reference, 2, tolerance = 1e-4)
  expect_true("discordant" %in% attr(r, "pair_ratios")$stratum)
})

test_that("equal strata report no allele effect, which is the point of the split", {
  r <- ibd_block_extension_by_allele(mk_allele_ibd(carrier_len = 80000, ref_len = 80000),
                                     LOCUS, allele = "pin",
                                     carrier = "Present", reference = "Absent")
  expect_equal(r$ratio_carrier, 4, tolerance = 1e-4)
  expect_equal(r$ratio_reference, 4, tolerance = 1e-4)
  expect_equal(r$ratio_contrast, 1)
  expect_gt(r$p_length, 0.05)
  # and the locus still looks extended when the pairs are pooled, which is exactly the
  # statement the split refutes
  pooled <- ibd_block_extension_test(mk_allele_ibd(carrier_len = 80000, ref_len = 80000),
                                     LOCUS)
  expect_equal(pooled$paired_ratio, 4, tolerance = 1e-4)
  expect_lt(pooled$p_paired, 0.05)
})

test_that("thin strata get no length test but keep their fraction", {
  r <- ibd_block_extension_by_allele(mk_allele_ibd(n_ref = 2), LOCUS, allele = "pin",
                                     carrier = "Present", reference = "Absent",
                                     min_pairs = 5)
  expect_equal(r$n_reference, 2L)
  expect_true(is.na(r$ratio_reference))         # below min_pairs
  expect_true(is.na(r$p_length))
  expect_false(is.na(r$frac_reference))         # the denominator is there either way
  expect_false(is.na(r$p_fraction))
})

test_that("states default sensibly, samples with no state drop out, and inputs are checked", {
  ibd <- mk_allele_ibd()
  expect_message(d <- ibd_block_extension_by_allele(ibd, LOCUS, allele = "pin"),
                 "as reference")
  expect_setequal(unname(attr(d, "states")), c("Present", "Absent"))

  # a named vector works as well as a column, and NA drops the pair from every stratum
  v <- stats::setNames(rep(c("Present", "Absent"), each = 12), c(CARRIERS, REFS))
  v[CARRIERS[1:2]] <- NA
  n <- ibd_block_extension_by_allele(ibd, LOCUS, allele = v,
                                     carrier = "Present", reference = "Absent")
  expect_equal(n$n_carrier, 5L)
  expect_equal(n$n_carrier_possible, choose(10, 2))

  expect_error(ibd_block_extension_by_allele(data.frame(a = 1), LOCUS, allele = "pin"),
               "must be an IbdResults")
  expect_error(ibd_block_extension_by_allele(ibd, LOCUS, allele = "nope"),
               "must name a column")
  expect_error(ibd_block_extension_by_allele(ibd, LOCUS, allele = "pin",
                                             carrier = "Present", reference = "Present"),
               "two different states")
})

test_that("the object method matches the function", {
  ibd <- mk_allele_ibd()
  expect_equal(ibd$ibd_block_extension_by_allele(LOCUS, allele = "pin", carrier = "Present",
                                                 reference = "Absent"),
               ibd_block_extension_by_allele(ibd, LOCUS, allele = "pin", carrier = "Present",
                                             reference = "Absent"))
})
