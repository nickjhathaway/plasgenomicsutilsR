# Subtracting a callset's SNPs from the codons you want, to find what was never called.
# The behaviour that matters is base resolution: a codon with one of its three bases in the
# panel has to come back as the other two, since a residue cannot be read from one base.

.codons <- function() {
  cds <- read_gff_cds(system.file("extdata", "pf3d7_drug_gene_cds.gff",
                                  package = "plasgenomicsutilsR"))
  aa_intervals(data.frame(transcript_id = c("pfcrt", "pfcrt", "pfkelch13"),
                          aa_position = c(72, 76, 580)), cds)
}

test_that("a partly covered interval comes back as the bases that are not", {
  want <- .codons()
  mid <- want$start[want$aa_position == 76] + 1        # the middle base only
  gaps <- bed_subtract(want, paste0("Pf3D7_07_v3:", mid))

  k76 <- gaps[gaps$aa_position == 76, ]
  expect_equal(nrow(k76), 2L)                          # split either side of the covered base
  expect_equal(sum(k76$width), 2)
  expect_equal(k76$piece, 1:2)
  expect_false(any(k76$start <= mid & k76$end > mid))   # neither piece contains it
  # the untouched codons are unchanged, and row order is kept
  expect_equal(sum(gaps$width[gaps$aa_position != 76]), 6)
  expect_equal(as.character(unique(gaps$aa_position)), c("72", "76", "580"))
})

test_that("a fully covered interval disappears and an untouched one survives whole", {
  want <- .codons()
  s72 <- want$start[want$aa_position == 72]
  gaps <- bed_subtract(want, paste0("Pf3D7_07_v3:", s72 + 0:2))
  expect_false(72 %in% gaps$aa_position)
  expect_equal(nrow(gaps), 2L)
  expect_true(all(gaps$width == 3))

  # nothing to subtract, and a different chromosome, both leave every base
  expect_equal(sum(bed_subtract(want, character(0))$width), 9)
  expect_equal(sum(bed_subtract(want, "Pf3D7_01_v3:1000")$width), 9)
})

test_that("the pieces never overlap what was subtracted", {
  want <- .codons()
  set.seed(2)
  bases <- unlist(Map(function(c_, s, e) paste0(c_, ":", seq(s, e - 1)),
                      want$chrom, want$start, want$end))
  have <- sample(bases, 4)
  gaps <- bed_subtract(want, have)
  # subtracting the same set again is a no-op, which it can only be if nothing overlaps
  expect_equal(nrow(bed_subtract(gaps, have)), nrow(gaps))
  expect_equal(sum(gaps$width), length(bases) - length(have))
})

test_that("locs2 is accepted as a table, a matrix, or a PopStructure", {
  want <- .codons()
  s <- want$start[want$aa_position == 76]
  tbl <- data.frame(chr = "Pf3D7_07_v3", start = s, end = s + 2)
  by_tbl <- bed_subtract(want, tbl)
  by_ids <- bed_subtract(want, paste0("Pf3D7_07_v3:", s + 0:1))
  expect_equal(by_tbl, by_ids)

  ps <- example_pop_structure(umap = FALSE)
  expect_s3_class(bed_subtract(want, ps), "tbl_df")            # PopStructure
  expect_s3_class(bed_subtract(want, ps$genotype()), "tbl_df") # bare matrix

  expect_error(bed_subtract(want, data.frame(chr = "7")), "no column")
  expect_error(bed_subtract(data.frame(start = 1, end = 2), tbl), "no 'chr' column")
})

test_that("min_width drops slivers", {
  want <- .codons()
  s <- want$start[want$aa_position == 76]
  blk <- data.frame(chr = "Pf3D7_07_v3", start = s, end = s + 2)   # leaves 1 base
  expect_equal(nrow(bed_subtract(want, blk)), 3L)
  expect_equal(nrow(bed_subtract(want, blk, min_width = 2)), 2L)
})

test_that("write_bed writes the contig name the reference uses, not the normalised one", {
  want <- .codons()
  f <- tempfile(fileext = ".bed")
  write_bed(bed_subtract(want, character(0)), f)
  lines <- strsplit(readLines(f), "\t", fixed = TRUE)

  # `chr` is "7" for matching; a BED read against a real FASTA needs "Pf3D7_07_v3"
  expect_true(all(grepl("^Pf3D7_", vapply(lines, `[`, "", 1))))
  expect_equal(length(lines), 3L)
  expect_equal(sum(vapply(lines, function(l) as.numeric(l[3]) - as.numeric(l[2]), 0)), 9)
  expect_equal(vapply(lines, length, 0L), rep(4L, 3))          # name column carried
  # sorted, and 0-based half-open as written
  starts <- as.numeric(vapply(lines, `[`, "", 2))
  expect_false(is.unsorted(starts[vapply(lines, `[`, "", 1) == "Pf3D7_07_v3"]))

  # an explicit column wins, and a missing one is named
  f2 <- tempfile(fileext = ".bed")
  write_bed(want, f2, chrom = "chr")
  expect_true(all(vapply(strsplit(readLines(f2), "\t"), `[`, "", 1) %in% c("7", "13")))
  expect_error(write_bed(want, f2, chrom = "nope"), "no `nope` column")
  expect_error(write_bed(data.frame(start = 1, end = 2), f2), "needs a `chrom`")
})
