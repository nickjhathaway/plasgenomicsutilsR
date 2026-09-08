# Tandem repeats: parsing a finder's BED, the period of a unit, the slippage thresholds at
# their boundaries, merging runs that flow into one another, and the bundled Pf3D7 table.
# The thresholds are the ones the HEOME Pf3D7 design used, so their edges are asserted
# exactly: a 10 bp homopolymer is kept and an 11 bp one is not.

.tr_bed <- function(...) {
  rec <- list(...)
  data.frame(chrom = vapply(rec, `[[`, "", 1), start = as.numeric(vapply(rec, `[[`, "", 2)),
             end = as.numeric(vapply(rec, `[[`, "", 3)),
             name = vapply(rec, function(r) paste0(r[1], "-", r[2], "-", r[3], "__", r[4]), ""),
             stringsAsFactors = FALSE)
}

test_that("the name field yields the unit, its period and the finder's copy number", {
  bed <- .tr_bed(c("Pf3D7_07_v3", 100, 120, "A_x20"),
                 c("Pf3D7_07_v3", 130, 160, "ATAT_x7.5"),
                 c("Pf3D7_07_v3", 500, 560, "TATT_x15"),
                 c("Pf3D7_07_v3", 600, 618, "AATAAT_x3"),
                 c("Pf3D7_07_v3", 700, 760, "TAACTCATTTATTACT_x2.17778"))
  tr <- tandem_repeats(bed)
  expect_s3_class(tr, "tbl_df")
  expect_equal(tr$repeat_unit, c("A", "ATAT", "TATT", "AATAAT", "TAACTCATTTATTACT"))
  expect_equal(tr$unit_size, c(1L, 4L, 4L, 6L, 16L))
  # the period is the true repeat length, whatever spelling the finder chose
  expect_equal(tr$period, c(1L, 2L, 4L, 3L, 16L))
  expect_equal(tr$copies, c(20, 7.5, 15, 3, 2.17778))
  expect_equal(tr$width, c(20, 30, 60, 18, 60))
  expect_equal(tr$chr, rep("7", 5))
  expect_equal(tr$chrom, rep("Pf3D7_07_v3", 5))
})

test_that("the location prefix is optional and only UNIT_xCOPIES is required", {
  bed <- data.frame(chrom = "Pf3D7_07_v3", start = c(100, 130), end = c(120, 160),
                    name = c("A_x20", "atat_x7.5"))          # bare, and lower case
  tr <- tandem_repeats(bed)
  expect_equal(tr$repeat_unit, c("A", "ATAT"))
  expect_equal(tr$period, c(1L, 2L))
  expect_equal(tr$copies, c(20, 7.5))
  with_prefix <- bed; with_prefix$name <- paste0("Pf3D7_07_v3-", bed$start, "-", bed$end, "__", bed$name)
  expect_equal(tandem_repeats(with_prefix)[, -6], tr[, -6])   # everything but `name`
})

test_that("a name that is not UNIT_xCOPIES gives NA and is flagged on width alone", {
  bed <- data.frame(chrom = "Pf3D7_07_v3", start = c(100, 200, 300, 400),
                    end = c(160, 210, 330, 440),
                    name = c("Pf3D7_07_v3-100-160", "Pf3D7_07_v3-200-210",
                             "Pf3D7_07_v3-300-330__AT", "AT_x"))
  expect_message(tr <- tandem_repeats(bed), "not `UNIT_xCOPIES`")
  expect_true(all(is.na(tr$repeat_unit)))
  expect_true(all(is.na(tr$period)))
  expect_true(all(is.na(tr$copies)))
  fl <- flag_tandem_repeats(tr)
  expect_equal(fl$avoid, c(TRUE, FALSE, FALSE, FALSE))   # 60 bp >= 50; the rest kept
})

test_that("a BED file is read the same way, gzipped or not, and positional columns work", {
  bed <- .tr_bed(c("Pf3D7_07_v3", 100, 120, "A_x20"), c("Pf3D7_07_v3", 130, 160, "AT_x15"))
  f <- tempfile(fileext = ".bed")
  utils::write.table(cbind(bed, 20, "+"), f, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = FALSE)         # a 6-column BED
  from_file <- tandem_repeats(f)
  expect_equal(from_file$repeat_unit, c("A", "AT"))
  expect_equal(from_file$start, c(100, 130))

  fz <- tempfile(fileext = ".bed.gz")
  con <- gzfile(fz, "w"); writeLines(readLines(f), con); close(con)
  expect_equal(tandem_repeats(fz), from_file)

  # unnamed columns are taken by position
  pos <- setNames(bed, c("V1", "V2", "V3", "V4"))
  expect_equal(tandem_repeats(pos), from_file)

  expect_error(tandem_repeats(bed[, 1:3]), "name")
  expect_error(tandem_repeats(tempfile()), "no such file")
  three <- tempfile(); writeLines("Pf3D7_07_v3\t1\t2", three)
  expect_error(tandem_repeats(three), "needs chrom, start, end and a name")
})

test_that("the slippage thresholds bite at exactly the lengths the design used", {
  bed <- .tr_bed(c("Pf3D7_07_v3", 100, 110, "A_x10"),      # 10: kept
                 c("Pf3D7_07_v3", 200, 211, "A_x11"),      # 11: flagged
                 c("Pf3D7_07_v3", 300, 311, "AT_x5.5"),    # 11: kept
                 c("Pf3D7_07_v3", 400, 412, "AT_x6"),      # 12: flagged
                 c("Pf3D7_07_v3", 500, 520, "AAT_x6.67"),  # 20: kept
                 c("Pf3D7_07_v3", 600, 621, "AAT_x7"),     # 21: flagged
                 c("Pf3D7_07_v3", 700, 749, "TATT_x12"),   # 49, period 4: kept
                 c("Pf3D7_07_v3", 800, 850, "TATT_x12.5")) # 50: flagged
  fl <- flag_tandem_repeats(bed)
  expect_equal(fl$avoid, rep(c(FALSE, TRUE), 4))

  # a dinucleotide spelled as a 4-mer is still judged as a dinucleotide
  quad <- .tr_bed(c("Pf3D7_07_v3", 100, 112, "ATAT_x3"))
  expect_true(flag_tandem_repeats(quad)$avoid)

  # thresholds can be given unnamed (periods 1, 2, 3, ...), extended, or dropped
  expect_equal(flag_tandem_repeats(bed, c(11, 12, 21))$avoid, fl$avoid)
  expect_true(flag_tandem_repeats(bed, c("1" = 11, "2" = 12, "3" = 21, "4" = 40))$avoid[7])
  expect_equal(flag_tandem_repeats(bed, NULL)$avoid, c(rep(FALSE, 7), TRUE))
  expect_equal(flag_tandem_repeats(bed, NULL, min_width = 20)$avoid,
               c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE))
  expect_error(flag_tandem_repeats(bed, c(a = 11)), "periods")
  expect_error(flag_tandem_repeats(bed, c("1" = 11, "1" = 12)), "twice")
  expect_error(flag_tandem_repeats(bed, min_width = "50"), "one number")
})

test_that("merging joins overlapping, abutting and nearby repeats and records what went in", {
  bed <- .tr_bed(c("Pf3D7_07_v3", 100, 116, "A_x16"),      # A run ...
                 c("Pf3D7_07_v3", 115, 130, "AT_x7.5"),    # ... flowing into an AT run
                 c("Pf3D7_07_v3", 130, 145, "AAT_x5"),     # abutting (end == start)
                 c("Pf3D7_07_v3", 150, 170, "TATT_x5"),    # 5 bp away
                 c("Pf3D7_07_v3", 100, 116, "AA_x8"),      # the same span, spelled by the other finder
                 c("Pf3D7_08_v3", 100, 116, "A_x16"))      # another chromosome
  m <- merge_tandem_repeats(bed)
  expect_equal(nrow(m), 3L)
  expect_equal(m$start, c(100, 150, 100))
  expect_equal(m$end, c(145, 170, 116))
  expect_equal(m$chr, c("7", "7", "8"))
  expect_equal(m$n_repeats, c(4L, 1L, 1L))
  expect_equal(m$repeat_units[1], "A,AA,AT,AAT")
  expect_equal(m$periods[1], "1,2,3")
  expect_equal(m$min_period, c(1L, 4L, 1L))
  expect_equal(m$name[1], "Pf3D7_07_v3-100-145")
  expect_equal(m$width, m$end - m$start)

  # a gap reaches the repeat 5 bp away
  expect_equal(nrow(merge_tandem_repeats(bed, gap = 5)), 2L)
  expect_equal(merge_tandem_repeats(bed, gap = 5)$end[1], 170)
  expect_equal(nrow(merge_tandem_repeats(bed, gap = 4)), 3L)
  expect_error(merge_tandem_repeats(bed, gap = -1), "gap")

  # a flag on any member flags the block
  fl <- flag_tandem_repeats(bed)
  mf <- merge_tandem_repeats(fl)
  expect_equal(mf$avoid, c(TRUE, FALSE, TRUE))
})

test_that("a short repeat nested inside a long one does not break the running merge", {
  bed <- .tr_bed(c("Pf3D7_07_v3", 100, 400, "AACCCTA_x43"),
                 c("Pf3D7_07_v3", 110, 130, "CCTAAAC_x3"),
                 c("Pf3D7_07_v3", 380, 420, "AACCCTA_x6"))
  m <- merge_tandem_repeats(bed)
  expect_equal(nrow(m), 1L)
  expect_equal(c(m$start, m$end), c(100, 420))
})

test_that("tandem_repeats_to_avoid keeps flagged runs and long chains of short ones", {
  bed <- .tr_bed(c("Pf3D7_07_v3", 100, 111, "A_x11"),       # flagged on its own
                 c("Pf3D7_07_v3", 111, 117, "AT_x3"),       # harmless, but absorbed into it
                 c("Pf3D7_07_v3", 300, 320, "AAT_x6.67"),   # three harmless repeats ...
                 c("Pf3D7_07_v3", 320, 340, "ATT_x6.67"),   # ... that chain to 60 bp
                 c("Pf3D7_07_v3", 340, 360, "TTA_x6.67"),
                 c("Pf3D7_07_v3", 500, 510, "A_x10"))       # harmless and alone
  av <- tandem_repeats_to_avoid(bed)
  expect_equal(nrow(av), 2L)
  expect_equal(av$start, c(100, 300))
  expect_equal(av$end, c(117, 360))
  expect_equal(av$n_repeats, c(2L, 3L))
  expect_false("avoid" %in% names(av))

  # the chain only counts once merged: without merging none of its members is flagged
  expect_false(any(flag_tandem_repeats(bed)$avoid[3:5]))
  # and a gap can be what makes the chain
  spaced <- bed; spaced$start[4:5] <- spaced$start[4:5] + 3; spaced$end[4:5] <- spaced$end[4:5] + 3
  spaced$name <- paste0(spaced$chrom, "-", spaced$start, "-", spaced$end, "__",
                        sub("^.*__", "", spaced$name))
  expect_equal(nrow(tandem_repeats_to_avoid(spaced)), 1L)
  expect_equal(nrow(tandem_repeats_to_avoid(spaced, gap = 3)), 2L)

  expect_equal(nrow(tandem_repeats_to_avoid(bed[0, ])), 0L)
})

test_that("the bundled Pf3D7 table loads, is cached, and reproduces the design's mask", {
  tr <- pf3d7_tandem_repeats()
  expect_s3_class(tr, "tbl_df")
  expect_equal(names(tr), c("chrom", "chr", "start", "end", "width", "name", "repeat_unit",
                            "unit_size", "period", "copies"))
  expect_gt(nrow(tr), 270000)
  expect_equal(sort(unique(tr$chr)), sort(c(as.character(1:14), "API", "MIT")))
  expect_false(anyNA(tr$repeat_unit))
  expect_true(all(tr$end > tr$start))
  expect_true(all(tr$period <= tr$unit_size))
  expect_identical(pf3d7_tandem_repeats(), tr)           # cached
  expect_equal(attr(tr, "source_version"), "2021-04-06")

  # the first record of the file, and the design's own numbers: 127,399 records flagged,
  # and a mask of about 76k blocks over about 4 Mb
  expect_equal(tr$name[1], "Pf3D7_01_v3-2-356__AACCCTA_x50.5714")
  fl <- flag_tandem_repeats(tr)
  expect_equal(sum(fl$avoid), 127399L)
  av <- tandem_repeats_to_avoid(fl)
  expect_true(nrow(av) > 75000 && nrow(av) < 77000)
  expect_true(sum(av$width) > 4.0e6 && sum(av$width) < 4.3e6)
})

test_that("bed_merge merges plain intervals, with an optional gap", {
  iv <- data.frame(chr = c("7", "7", "7", "7", "Pf3D7_08_v3"),
                   start = c(100, 150, 300, 320, 100), end = c(200, 250, 310, 330, 200))
  m <- bed_merge(iv)
  expect_equal(nrow(m), 4L)
  expect_equal(m$start, c(100, 300, 320, 100))
  expect_equal(m$end, c(250, 310, 330, 200))
  expect_equal(m$n, c(2L, 1L, 1L, 1L))
  expect_equal(m$chr, c("7", "7", "7", "8"))
  expect_equal(m$chrom[4], "Pf3D7_08_v3")                # the input's own spelling
  g <- bed_merge(iv, gap = 10)
  expect_equal(nrow(g), 3L)
  expect_equal(g$end[2], 330)
  expect_equal(bed_merge(iv, chrom = "chr")$width, m$width)
  expect_equal(nrow(bed_merge(iv[0, ])), 0L)
  # abutting intervals merge, as bedtools merge does
  expect_equal(nrow(bed_merge(data.frame(chr = "1", start = c(0, 10), end = c(10, 20)))), 1L)
})

test_that("bed_subtract pads the mask and the mask can be the tandem-repeat blocks", {
  gene <- data.frame(name = "g", chr = "Pf3D7_07_v3", start = 1000, end = 2000)
  mask <- data.frame(chr = "7", start = c(1200, 1500), end = c(1210, 1520))
  plain <- bed_subtract(gene, mask)
  expect_equal(plain$start, c(1000, 1210, 1520))
  expect_equal(plain$end, c(1200, 1500, 2000))
  padded <- bed_subtract(gene, mask, pad = 10)
  expect_equal(padded$start, c(1000, 1220, 1530))
  expect_equal(padded$end, c(1190, 1490, 2000))
  expect_error(bed_subtract(gene, mask, pad = -1), "pad")
  # padding cannot run off the start of the chromosome
  expect_equal(nrow(bed_subtract(data.frame(chr = "7", start = 0, end = 5),
                                 data.frame(chr = "7", start = 2, end = 3), pad = 10)), 0L)

  # a real design step: genes minus the bundled mask, written with the reference's names
  genes <- PF3D7_GENES[PF3D7_GENES$name %in% c("pfcrt", "pfdhfr"), ]
  av <- tandem_repeats_to_avoid(pf3d7_tandem_repeats())
  cut <- bed_subtract(genes, av, pad = 10)
  expect_true(all(c("gene_id", "Pf3D7_chrom", "piece", "width") %in% names(cut)))
  expect_true(all(cut$width > 0))
  expect_lt(sum(cut$width), sum(genes$end - genes$start))
  # nothing left overlaps the padded mask, on any of the pieces
  expect_equal(nrow(bed_intersect(cut, av)$overlap), 0L)
  f <- tempfile(fileext = ".bed")
  write_bed(cut, f)
  first <- strsplit(readLines(f)[1], "\t")[[1]]
  expect_equal(first[1], "Pf3D7_04_v3")                   # not the short "4"
})

test_that("write_bed prefers the reference's own chromosome spelling", {
  # PF3D7_GENES-style: `chrom` is the short key, `Pf3D7_chrom` the reference's name
  x <- data.frame(Pf3D7_chrom = "Pf3D7_07_v3", chrom = "7", start = 1, end = 2)
  f <- tempfile(fileext = ".bed")
  write_bed(x, f)
  expect_equal(strsplit(readLines(f), "\t")[[1]][1], "Pf3D7_07_v3")
  # aa_intervals-style: `chrom` is the reference's name, `chr` the key
  y <- data.frame(chrom = "Pf3D7_07_v3", chr = "7", start = 1, end = 2)
  write_bed(y, f)
  expect_equal(strsplit(readLines(f), "\t")[[1]][1], "Pf3D7_07_v3")
  # an explicit column still wins
  write_bed(x, f, chrom = "chrom")
  expect_equal(strsplit(readLines(f), "\t")[[1]][1], "7")
})

test_that("dedupe drops a repeat inside another with the same unit up to rotation, and nothing else", {
  bed <- .tr_bed(c("Pf3D7_07_v3", 100, 130, "ATA_x10"),     # the run ...
                 c("Pf3D7_07_v3", 101, 128, "TAT_x9"),      # ... reported again, rotated, inside it
                 c("Pf3D7_07_v3", 100, 130, "ATAATA_x5"),   # identical span, doubled spelling
                 c("Pf3D7_07_v3", 110, 122, "A_x12"),       # a different period inside it: kept
                 c("Pf3D7_07_v3", 300, 320, "AAT_x6.67"),
                 c("Pf3D7_08_v3", 101, 128, "TAT_x9"))      # same span on another chromosome
  m <- merge_tandem_repeats(bed)
  d <- merge_tandem_repeats(bed, dedupe = TRUE)
  expect_equal(d[, c("chr", "start", "end", "width")], m[, c("chr", "start", "end", "width")])
  expect_equal(m$n_repeats, c(4L, 1L, 1L))
  expect_equal(d$n_repeats, c(2L, 1L, 1L))
  expect_equal(m$repeat_units[1], "ATA,ATAATA,TAT,A")
  expect_equal(d$repeat_units[1], "ATA,A")
  expect_equal(d$periods[1], "1,3")

  # flags survive: the inner A run is what flags the first block, so it cannot be dropped;
  # the 27 bp TAT run on chr 8 is flagged on its own
  fl <- flag_tandem_repeats(bed)
  expect_equal(merge_tandem_repeats(fl, dedupe = TRUE)$avoid, c(TRUE, FALSE, TRUE))
  expect_equal(tandem_repeats_to_avoid(bed, dedupe = TRUE)[, c("start", "end")],
               tandem_repeats_to_avoid(bed)[, c("start", "end")])

  # the canonical unit is the period core in the smallest rotation of itself or its
  # reverse complement (TAT is the reverse complement of ATA, not a rotation)
  expect_equal(plasgenomicsutilsR:::.canonical_unit(c("ATA", "TAT", "ATAATA", "TA", "ATAT", "A", "T", NA),
                                                    c(3L, 3L, 3L, 2L, 2L, 1L, 1L, NA)),
               c("AAT", "AAT", "AAT", "AT", "AT", "A", "A", NA))
})

test_that("on the bundled table dedupe leaves the mask block for block", {
  tr <- pf3d7_tandem_repeats()
  a <- tandem_repeats_to_avoid(tr)
  b <- tandem_repeats_to_avoid(tr, dedupe = TRUE)
  expect_equal(b[, c("chr", "start", "end")], a[, c("chr", "start", "end")])
  expect_lt(sum(b$n_repeats), sum(a$n_repeats))
})
