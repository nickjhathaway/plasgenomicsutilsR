# A cohort where the sharing at one place on chr 7 is genuinely longer than everywhere else,
# and a flat background of segments spread along chr 7 and chr 8 to tile windows over.
blk <- function(pairs, chr, start, len) {
  data.frame(sample1 = pairs$sample1, sample2 = pairs$sample2, chr = chr,
             start = start, end = start + len - 1,
             different = 0, Nsnp = 40, stringsAsFactors = FALSE)
}
mk_pairs <- function(n, offset = 0) {
  i <- seq_len(n)
  data.frame(sample1 = paste0("s", offset + 2 * i - 1),
             sample2 = paste0("s", offset + 2 * i), stringsAsFactors = FALSE)
}
REGIONS <- data.frame(name = c("r7", "r8"), chr = c("7", "8"),
                      start = 0, end = 1000000, stringsAsFactors = FALSE)

mk_scan_ibd <- function() {
  p <- mk_pairs(6)
  # background: each pair carries a 20 kb segment at a different place on each chromosome,
  # so windows all along both chromosomes have something to see
  bg <- do.call(rbind, lapply(seq_len(nrow(p)), function(i) rbind(
    blk(p[i, ], "7", 20000 + (i - 1) * 30000, 20000),
    blk(p[i, ], "7", 600000 + (i - 1) * 30000, 20000),
    blk(p[i, ], "8", 20000 + (i - 1) * 30000, 20000),
    blk(p[i, ], "8", 500000 + (i - 1) * 30000, 20000))))
  # the hot spot: every pair shares a long segment over 400-410 kb on chr 7
  hot <- blk(p, "7", 405000 - (80000 + (seq_len(6) - 3.5) * 1000) / 2,
             80000 + (seq_len(6) - 3.5) * 1000)
  bl <- rbind(bg, hot)
  meta <- data.frame(sample = unique(c(bl$sample1, bl$sample2)), region = "north",
                     stringsAsFactors = FALSE)
  ibd_results(blocks = bl, meta = meta, group_col_in_meta = "region")
}

test_that("the scan tiles the regions and reports a percentile per window", {
  s <- ibd_block_extension_scan(mk_scan_ibd(), regions = REGIONS, width = 500, step = 50000,
                                min_pairs = 1, min_windows = 5)

  expect_true(all(s$source == "window"))
  expect_true(all(s$span_bp == 500))                  # uniform width is the whole point
  expect_equal(unname(attr(s, "window")), c(500, 50000))
  expect_true(all(s$percentile >= 0 & s$percentile <= 100))
  expect_equal(max(s$percentile), 100)
  expect_equal(s$n_windows[1], sum(s$source == "window"))
  # ordered by descending ratio, and the hot spot is the top of the distribution
  expect_equal(s$paired_ratio, sort(s$paired_ratio, decreasing = TRUE))
  expect_equal(as.character(s$chr[1]), "7")
  expect_true(s$start[1] >= 350000 && s$start[1] <= 450000)
  # a scan does not carry the per-pair table
  expect_null(attr(s, "pair_ratios"))
})

test_that("named loci are placed against the window distribution", {
  ibd <- mk_scan_ibd()
  hot <- data.frame(name = "hot", chr = "7", start = 400000, end = 410000,
                    stringsAsFactors = FALSE)
  cold <- data.frame(name = "cold", chr = "8", start = 30000, end = 40000,
                     stringsAsFactors = FALSE)
  s <- ibd_block_extension_scan(ibd, loci = rbind(hot, cold), regions = REGIONS,
                                width = 500, step = 50000, min_pairs = 1, min_windows = 5)

  loci <- s[s$source == "locus", ]
  expect_setequal(as.character(loci$locus), c("hot", "cold"))
  expect_gt(loci$percentile[loci$locus == "hot"], loci$percentile[loci$locus == "cold"])
  expect_gt(loci$percentile[loci$locus == "hot"], 90)
  # loci are listed first, and the windows are still all there behind them
  expect_equal(as.character(head(s$source, 2)), c("locus", "locus"))
  expect_equal(sum(s$source == "window"), s$n_windows[1])
  # the locus rows agree with running the test on them directly
  d <- ibd_block_extension_test(ibd, hot, min_pairs = 1)
  expect_equal(loci$paired_ratio[loci$locus == "hot"], d$paired_ratio)
})

test_that("percentiles are withheld from a group with too few windows", {
  ibd <- mk_scan_ibd()
  expect_warning(s <- ibd_block_extension_scan(ibd, regions = REGIONS, width = 500,
                                               step = 300000, min_pairs = 1,
                                               min_windows = 50), "fewer than 50")
  expect_true(all(is.na(s$percentile)))
  expect_true(all(s$n_windows < 50))
})

test_that("windows carry the sharing rule through, and the inputs are checked", {
  ibd <- mk_scan_ibd()
  cp <- ibd_block_extension_scan(ibd, regions = REGIONS, width = 500, step = 50000,
                                 sharing = "complete", min_pairs = 1, min_windows = 5)
  expect_equal(attr(cp, "sharing"), "complete")
  # a 500 bp window is spanned by any segment overlapping it well inside, so complete and
  # overlap agree at the top; they need not agree at a window clipped by a segment end
  expect_equal(as.character(cp$chr[1]), "7")

  expect_error(ibd_block_extension_scan(data.frame(a = 1)), "must be an IbdResults")
  expect_error(ibd_block_extension_scan(ibd, regions = REGIONS, width = 0), "positive")
  expect_error(ibd_block_extension_scan(ibd, regions = REGIONS, step = 0), "positive")
  expect_error(ibd_block_extension_scan(ibd, regions = data.frame(
    name = "tiny", chr = "7", start = 0, end = 100), width = 500), "at least")
})

test_that("the object method matches the function", {
  ibd <- mk_scan_ibd()
  expect_equal(ibd$ibd_block_extension_scan(regions = REGIONS, width = 500, step = 50000,
                                            min_pairs = 1, min_windows = 5),
               ibd_block_extension_scan(ibd, regions = REGIONS, width = 500, step = 50000,
                                        min_pairs = 1, min_windows = 5))
})
