# Synthetic segments with a known answer. An IbdResults shifts `end` to make the interval
# half-open, so these helpers take a *length* and the stored segment is exactly that long.
# Lengths at a locus are spread by 1 kb per pair, symmetrically about the nominal length, so
# the per-pair ratios are distinct -- the signed-rank test is then exact rather than a
# tie-corrected approximation, and its p-value is 1 / 2^n_pairs.
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

spread <- function(len, n) len + (seq_len(n) - (n + 1) / 2) * 1000   # median is `len`

LOCUS <- data.frame(name = "sweep", chr = "7", start = 400000, end = 410000,
                    stringsAsFactors = FALSE)

# `n` pairs, each sharing `locus_len` across the locus on chr 7 and `ref_len` on chrs 1 and 2.
# Locus segments are centred on the locus, so every pair spans it whatever its length.
LOCUS_MID <- 405000

mk_group_blocks <- function(n, offset = 0, locus_len = 80000, ref_len = 20000) {
  p <- mk_pairs(n, offset)
  len <- spread(locus_len, n)
  rbind(blk(p, "7", LOCUS_MID - len / 2, len),
        blk(p, "1", 100000, ref_len),
        blk(p, "2", 100000, ref_len))
}

mk_ibd <- function(n = 5, locus_len = 80000, ref_len = 20000, extra_blocks = NULL) {
  bl <- mk_group_blocks(n, locus_len = locus_len, ref_len = ref_len)
  if (!is.null(extra_blocks)) bl <- rbind(bl, extra_blocks)
  meta <- data.frame(sample = unique(c(bl$sample1, bl$sample2)), region = "north",
                     stringsAsFactors = FALSE)
  ibd_results(blocks = bl, meta = meta, group_col_in_meta = "region")
}

test_that("paired_ratio is the locus segment over the pair's own baseline", {
  r <- ibd_block_extension_test(mk_ibd(n = 5, locus_len = 80000, ref_len = 20000), LOCUS)

  expect_equal(nrow(r), 1L)
  expect_equal(as.character(r$locus), "sweep")
  expect_equal(r$group, "north")
  expect_equal(r$n_pairs, 5L)
  expect_equal(r$locus_median, 80000)
  expect_equal(r$ref_median, 20000)
  expect_equal(r$gw_median, 20000)
  expect_equal(r$paired_ratio, 4)              # 80 kb over 20 kb, exactly
  expect_equal(r$naive_ratio, 4)
  expect_equal(r$span_bp, 10000)
  # five distinct positive log ratios: the exact one-sided signed-rank p is 1 / 2^5
  expect_equal(r$p_paired, 1 / 32)
  expect_equal(r$p_two_sided, 2 / 32)
  expect_equal(r$q_paired, r$p_paired)         # one locus, one group: BH cannot move it
  # the per-pair ratios it was built from come back on the result
  pr <- attr(r, "pair_ratios")
  expect_equal(nrow(pr), 5L)
  expect_true(all(pr$log2_ratio > 0))
  expect_equal(stats::median(pr$log2_ratio), 2)
  expect_true(all(c("locus", "group", "sample1", "sample2", "locus_len", "ref_n",
                    "ref_median", "log2_ratio") %in% names(pr)))
})

test_that("the pair's baseline excludes the whole locus chromosome", {
  # the same pairs plus a long extra segment on chr 7 nowhere near the locus. A genome-wide
  # baseline would take it in and flatten the ratio; the sweep's own chromosome has to go
  # entirely, not just the locus interval.
  far <- blk(mk_pairs(5), "7", 1500000, 500000)
  r_plain <- ibd_block_extension_test(mk_ibd(), LOCUS)
  r_extra <- ibd_block_extension_test(mk_ibd(extra_blocks = far), LOCUS)

  expect_equal(r_extra$ref_median, 20000)
  expect_equal(r_extra$paired_ratio, r_plain$paired_ratio)
  expect_equal(r_extra$p_paired, r_plain$p_paired)
  expect_equal(r_extra$locus_median, 80000)    # the off-locus segment is not the locus one
  # it does move the naive background, which is why naive_ratio is not the statistic
  expect_lt(r_extra$naive_ratio, r_plain$naive_ratio)
})

test_that("naive_ratio counts related pairs where paired_ratio does not", {
  # Twenty background pairs sharing only short segments away from the locus, and five
  # unusually related pairs carrying long segments everywhere, the locus included. The
  # locus set is all related pairs, so the naive ratio reads their relatedness as a locus
  # signal while the paired ratio, measuring each against itself, barely moves.
  bg <- mk_pairs(20, offset = 100)
  related <- mk_pairs(5)
  bl <- rbind(
    blk(bg, "1", 100000, 20000), blk(bg, "2", 100000, 20000),
    blk(bg, "7", 900000, 20000),                       # chr 7, but not over the locus
    blk(related, "7", 200000, spread(500000, 5)),      # spans the locus
    blk(related, "1", 100000, 400000), blk(related, "2", 100000, 400000))
  meta <- data.frame(sample = unique(c(bl$sample1, bl$sample2)), region = "north",
                     stringsAsFactors = FALSE)
  r <- ibd_block_extension_test(ibd_results(blocks = bl, meta = meta,
                                            group_col_in_meta = "region"), LOCUS)

  expect_equal(r$n_pairs, 5L)
  expect_equal(r$gw_median, 20000)
  expect_equal(r$locus_median, 500000)
  expect_equal(r$naive_ratio, 25)
  # ref_median far above gw_median is the tell that the sharing pairs are unusually related
  expect_equal(r$ref_median, 400000)
  expect_gt(r$ref_median, r$gw_median)
  expect_equal(r$paired_ratio, 1.25)
  expect_gt(r$naive_ratio, r$paired_ratio)
})

test_that("a locus no longer than the pairs' own background gets no p-value", {
  # every pair's locus segment is exactly as long as its baseline, so no spread here
  p <- mk_pairs(5)
  bl <- rbind(blk(p, "7", 395000, 20000), blk(p, "1", 100000, 20000),
              blk(p, "2", 100000, 20000))
  meta <- data.frame(sample = unique(c(bl$sample1, bl$sample2)), region = "north",
                     stringsAsFactors = FALSE)
  r <- ibd_block_extension_test(ibd_results(blocks = bl, meta = meta,
                                            group_col_in_meta = "region"), LOCUS)
  # every log ratio is exactly 0, so the signed-rank test has nothing to rank and says so
  # rather than erroring or leaning on the R version's handling of dropped zeros
  expect_equal(r$paired_ratio, 1)
  expect_equal(r$p_paired, 1)
  expect_equal(r$p_two_sided, 1)
})

test_that("min_ref_blocks drops pairs whose baseline is too thin", {
  # a sixth pair sharing at the locus with a single reference segment
  tp <- data.frame(sample1 = "t1", sample2 = "t2", stringsAsFactors = FALSE)
  thin <- rbind(blk(tp, "7", 365000, 80000), blk(tp, "1", 100000, 20000))
  ibd <- mk_ibd(extra_blocks = thin)

  expect_equal(ibd_block_extension_test(ibd, LOCUS, min_ref_blocks = 1)$n_pairs, 6L)
  expect_equal(ibd_block_extension_test(ibd, LOCUS, min_ref_blocks = 2)$n_pairs, 5L)
})

test_that("min_pairs drops a stratum rather than reporting it untested", {
  ibd <- mk_ibd(n = 4)
  expect_equal(ibd_block_extension_test(ibd, LOCUS, min_pairs = 4)$n_pairs, 4L)
  expect_equal(nrow(ibd_block_extension_test(ibd, LOCUS, min_pairs = 5)), 0L)
  # an empty result still has the full set of columns, so several loci can be bound together
  empty <- ibd_block_extension_test(ibd, LOCUS, min_pairs = 5)
  expect_true(all(c("locus", "group", "n_pairs", "paired_ratio", "q_paired") %in% names(empty)))
})

test_that("only within-group pairs count, and each group is tested on its own", {
  cross <- blk(data.frame(sample1 = "s1", sample2 = "s101"), "7", 300000, 600000)
  bl <- rbind(mk_group_blocks(5),                                  # north, ratio 4
              mk_group_blocks(5, offset = 100, locus_len = 40000), # south, ratio 2
              cross)
  meta <- data.frame(sample = unique(c(bl$sample1, bl$sample2)), stringsAsFactors = FALSE)
  meta$region <- ifelse(as.integer(sub("^s", "", meta$sample)) > 100, "south", "north")
  r <- ibd_block_extension_test(ibd_results(blocks = bl, meta = meta,
                                            group_col_in_meta = "region"), LOCUS)

  expect_equal(nrow(r), 2L)
  expect_setequal(r$group, c("north", "south"))
  expect_equal(r$paired_ratio[r$group == "north"], 4)
  expect_equal(r$paired_ratio[r$group == "south"], 2)
  expect_equal(sum(r$n_pairs), 10L)
  expect_equal(r$group, c("north", "south"))   # ordered by descending paired_ratio
  # the 600 kb cross-group segment belongs to neither stratum
  expect_false(any(attr(r, "pair_ratios")$locus_len == 600000))
})

test_that("segments come from get_blocks(), so the short-segment floor is already applied", {
  # A 5 kb segment across the locus is under the 15 kb floor and must not become a pair's
  # locus length. Mixing a filtered numerator with an unfiltered denominator is what
  # inflates these ratios several-fold.
  up <- data.frame(sample1 = "u1", sample2 = "u2", stringsAsFactors = FALSE)
  extra <- rbind(blk(up, "7", 402000, 5000),
                 blk(up, "1", 100000, 20000), blk(up, "2", 100000, 20000))
  r <- ibd_block_extension_test(mk_ibd(extra_blocks = extra), LOCUS)

  expect_equal(r$n_pairs, 5L)
  expect_false(any(attr(r, "pair_ratios")$locus_len == 5000))
})

test_that("locus chromosomes are normalised and the overlap rule matches gene_ibd_pairs", {
  ibd <- mk_ibd()
  long <- transform(LOCUS, chr = "Pf3D7_07_v3")
  expect_equal(ibd_block_extension_test(ibd, long)$paired_ratio,
               ibd_block_extension_test(ibd, LOCUS)$paired_ratio)
  expect_equal(ibd_block_extension_test(ibd, long)$chr, "7")

  # the longest segment ends at 446000, so an interval starting there does not overlap it:
  # touching is not overlapping
  abut <- data.frame(name = "abut", chr = "7", start = 446000, end = 450000)
  expect_equal(nrow(ibd_block_extension_test(ibd, abut, min_pairs = 1)), 0L)
  inside <- data.frame(name = "inside", chr = "7", start = 445999, end = 450000)
  expect_equal(ibd_block_extension_test(ibd, inside, min_pairs = 1)$n_pairs, 1L)
  # `within` pads the interval, bringing the abutting segment back
  expect_equal(ibd_block_extension_test(ibd, abut, within = 1, min_pairs = 1)$n_pairs, 1L)

  # exactly the pairs gene_ibd_pairs() finds at the same interval
  gp <- gene_ibd_pairs(ibd, genes = LOCUS)
  pr <- attr(ibd_block_extension_test(ibd, LOCUS), "pair_ratios")
  expect_setequal(paste(gp$sample1, gp$sample2), paste(pr$sample1, pr$sample2))
})

test_that("rows of an interval frame that share a name become one locus", {
  ibd <- mk_ibd()
  # three marker intervals of one haplotype, as an uncollapsed aa_intervals() table looks
  parts <- data.frame(name = "pin", chr = "7",
                      start = c(401000, 404000, 409000),
                      end   = c(401003, 404003, 409003), stringsAsFactors = FALSE)
  r <- ibd_block_extension_test(ibd, parts)

  expect_equal(nrow(r), 1L)
  expect_equal(r$start, 401000)
  expect_equal(r$end, 409003)                  # spans min start to max end
  expect_equal(r$span_bp, 8003)
  expect_equal(r$n_pairs, 5L)

  # a name repeating on two chromosomes is two features, not one interval spanning both --
  # a gene track full of paralogues has to survive being passed as a data frame
  para <- data.frame(name = "dup", chr = c("7", "8"), start = 400000, end = 410000,
                     stringsAsFactors = FALSE)
  expect_warning(r2 <- ibd_block_extension_test(ibd, para, min_pairs = 1), "repeat")
  expect_equal(nrow(r2), 1L)                   # only the chr 7 copy has any sharing
  expect_equal(r2$chr, "7")
  expect_equal(levels(r2$locus), "dup")

  # a gene_id keeps the parts of one feature together and the copies apart
  ided <- data.frame(name = "dup", gene_id = c("g1", "g1", "g2"), chr = c("7", "7", "8"),
                     start = c(401000, 409000, 400000), end = c(401003, 409003, 410000),
                     stringsAsFactors = FALSE)
  expect_warning(r3 <- ibd_block_extension_test(ibd, ided, min_pairs = 1), "repeat")
  expect_equal(nrow(r3), 1L)
  expect_equal(r3$gene_id, "g1")
  expect_equal(c(r3$start, r3$end), c(401000, 409003))   # the two chr 7 rows merged
})

test_that("adjust chooses the family the q-values are computed over", {
  # Three loci shared by 5, 6 and 7 of the same seven pairs, so their exact p-values are
  # 1/32, 1/64 and 1/128 and BH has something to reorder.
  p7 <- mk_pairs(7)
  bl <- rbind(
    blk(p7[1:5, ], "7", 365000, spread(80000, 5)),
    blk(p7[1:6, ], "8", 365000, spread(80000, 6)),
    blk(p7, "9", 365000, spread(80000, 7)),
    blk(p7, "1", 100000, 20000), blk(p7, "2", 100000, 20000))
  meta <- data.frame(sample = unique(c(bl$sample1, bl$sample2)), region = "north",
                     stringsAsFactors = FALSE)
  ibd <- ibd_results(blocks = bl, meta = meta, group_col_in_meta = "region")
  loci <- data.frame(name = c("a", "b", "c"), chr = c("7", "8", "9"),
                     start = 400000, end = 410000, stringsAsFactors = FALSE)

  per_locus <- ibd_block_extension_test(ibd, loci, adjust = "per_locus")
  per_group <- ibd_block_extension_test(ibd, loci, adjust = "per_group")
  all_rows <- ibd_block_extension_test(ibd, loci, adjust = "all")
  none <- ibd_block_extension_test(ibd, loci, adjust = "none")

  expect_equal(levels(per_locus$locus), c("a", "b", "c"))   # input order, not ratio order
  expect_equal(per_locus$n_pairs, c(5L, 6L, 7L))
  expect_equal(per_locus$p_paired, c(1 / 32, 1 / 64, 1 / 128))
  # each locus's baseline excludes only its own chromosome, so the other two loci's long
  # segments sit in it: the ratio is well above 1 but not the naive 4
  expect_true(all(per_locus$paired_ratio > 1))

  # one group per locus, so a per-locus family is a family of one: q is p untouched
  expect_equal(per_locus$q_paired, per_locus$p_paired)
  # all three loci sit in one group, so the scan family has three members
  expect_equal(per_group$q_paired, c(0.03125, 0.0234375, 0.0234375))
  expect_true(any(per_group$q_paired > per_group$p_paired))
  expect_equal(all_rows$q_paired, per_group$q_paired)       # a single group: same family
  expect_true(all(is.na(none$q_paired)))
  expect_equal(attr(none, "adjust"), "none")
})

test_that("a scan-sized call says per-locus is the wrong default family", {
  ibd <- mk_ibd()
  many <- data.frame(name = paste0("g", 1:30), chr = "7",
                     start = 400000, end = 410000, stringsAsFactors = FALSE)
  expect_message(ibd_block_extension_test(ibd, many), 'adjust = "per_group"')
  expect_silent(ibd_block_extension_test(ibd, many, adjust = "per_group"))
  expect_silent(ibd_block_extension_test(ibd, LOCUS))       # a handful of loci says nothing
})

test_that("samples with no group are reported and left out", {
  ibd <- mk_ibd()
  meta <- ibd$get_meta()
  meta$region[meta$sample %in% c("s1", "s2")] <- NA
  ibd$set_meta(meta)
  expect_message(r <- ibd_block_extension_test(ibd, LOCUS, min_pairs = 1),
                 "no `region`")
  expect_equal(r$n_pairs, 4L)
})

test_that("the object method matches the function, and inputs are checked", {
  ibd <- mk_ibd()
  expect_equal(ibd$ibd_block_extension_test(LOCUS), ibd_block_extension_test(ibd, LOCUS))

  expect_error(ibd_block_extension_test(data.frame(a = 1), LOCUS), "must be an IbdResults")
  expect_error(ibd_block_extension_test(ibd, LOCUS, group = "nope"), "no column 'nope'")
  expect_error(ibd_block_extension_test(ibd, data.frame()), "required column")
  expect_error(ibd_block_extension_test(
    ibd, data.frame(name = character(), chr = character(),
                    start = numeric(), end = numeric())), "empty")
})

test_that('sharing = "complete" requires the segment to span the whole locus', {
  # three pairs whose segments cover the locus outright, and two that only clip its edge
  span <- mk_pairs(3)
  clip <- mk_pairs(2, offset = 100)
  bl <- rbind(
    blk(span, "7", LOCUS_MID - spread(80000, 3) / 2, spread(80000, 3)),  # cover it
    blk(clip, "7", 392000, 16000),                                       # ends inside it
    blk(rbind(span, clip), "1", 100000, 20000),
    blk(rbind(span, clip), "2", 100000, 20000))
  meta <- data.frame(sample = unique(c(bl$sample1, bl$sample2)), region = "north",
                     stringsAsFactors = FALSE)
  ibd <- ibd_results(blocks = bl, meta = meta, group_col_in_meta = "region")

  ov <- ibd_block_extension_test(ibd, LOCUS, min_pairs = 1)
  cp <- ibd_block_extension_test(ibd, LOCUS, sharing = "complete", min_pairs = 1)

  expect_equal(ov$n_pairs, 5L)                 # everything that touches the locus
  expect_equal(cp$n_pairs, 3L)                 # only the segments that span it
  expect_equal(cp$locus_median, 80000)
  expect_equal(cp$paired_ratio, 4)
  # the clipping pairs drag the overlap median down, so complete reads higher here
  expect_lt(ov$paired_ratio, cp$paired_ratio)
  expect_equal(attr(cp, "sharing"), "complete")
  expect_equal(attr(ov, "sharing"), "overlap")
})

test_that('under "complete" the locus width is a floor on the segments counted', {
  # this is why a "complete" scan cannot rank loci of different widths against each other:
  # widen the locus past a pair's segment and that pair stops counting at all
  ibd <- mk_ibd()                              # segments 78-82 kb, centred on the locus
  narrow <- data.frame(name = "narrow", chr = "7", start = 404000, end = 406000)
  wide <- data.frame(name = "wide", chr = "7", start = LOCUS_MID - 40000,
                     end = LOCUS_MID + 40000)  # 80 kb, as wide as the segments themselves

  n <- ibd_block_extension_test(ibd, narrow, sharing = "complete", min_pairs = 1)
  w <- ibd_block_extension_test(ibd, wide, sharing = "complete", min_pairs = 1)

  expect_equal(n$n_pairs, 5L)                  # every segment spans a 2 kb locus
  expect_equal(w$n_pairs, 3L)                  # only segments of 80 kb or more span an 80 kb one
  expect_equal(n$locus_median, 80000)          # the whole 78-82 kb spread
  expect_equal(w$locus_median, 81000)          # the 78 and 79 kb segments are gone
  expect_gt(w$locus_median, n$locus_median)    # the width has cut off the short end
  expect_gt(w$paired_ratio, n$paired_ratio)
  # under "overlap" the same two intervals keep every pair, so width costs nothing there
  expect_equal(ibd_block_extension_test(ibd, narrow, min_pairs = 1)$n_pairs, 5L)
  expect_equal(ibd_block_extension_test(ibd, wide, min_pairs = 1)$n_pairs, 5L)
})

test_that('`within` pads the interval under either sharing rule', {
  ibd <- mk_ibd()
  # 80 kb wide, so only the segments of at least 80 kb cover it; padding widens it further
  wide <- data.frame(name = "wide", chr = "7", start = LOCUS_MID - 40000,
                     end = LOCUS_MID + 40000)
  expect_equal(ibd_block_extension_test(ibd, wide, sharing = "complete",
                                        min_pairs = 1)$n_pairs, 3L)
  # padded to 84 kb it is wider than every segment, so nothing spans it any more
  expect_equal(nrow(ibd_block_extension_test(ibd, wide, sharing = "complete",
                                             within = 2000, min_pairs = 1)), 0L)
})
