# genomic_range_aa_positions(): the residues an interval covers, per transcript. Checked on
# the real pfcrt (13 exons, codons straddling introns) and the minus-strand pfkelch13, plus a
# hand-built GFF with two isoforms.

.gff_path <- function() system.file("extdata", "pf3d7_drug_gene_cds.gff", package = "plasgenomicsutilsR")
.fa_path <- function() system.file("extdata", "pf3d7_drug_gene_regions.fasta.gz",
                                   package = "plasgenomicsutilsR")

test_that("a whole gene covers every codon, from 1, whatever is upstream of it", {
  cds <- read_gff_cds(.gff_path())
  crt <- PF3D7_GENES[PF3D7_GENES$name == "pfcrt", ]
  r <- genomic_range_aa_positions(crt, cds)
  expect_equal(nrow(r), 1L)
  expect_equal(r$transcript_id, "PF3D7_0709000.1")
  expect_equal(r$gene_name, "pfcrt")
  expect_equal(r$aa_start, 1L)
  # 424 residues plus the stop codon, which the CDS includes and which is numbered like any
  # other codon (as snp_aa_positions() does)
  expect_equal(r$aa_end, 425L)
  expect_equal(r$n_aa, 425L)
  expect_false(r$partial_start); expect_false(r$partial_end)
  expect_true(r$coding)
  with_seq <- genomic_range_aa_positions(crt, cds, fasta = .fa_path())
  expect_equal(nchar(with_seq$aa_seq), 425L)
  expect_true(startsWith(with_seq$aa_seq, "M"))
  expect_true(endsWith(with_seq$aa_seq, "*"))
  expect_equal(substr(with_seq$aa_seq, 76, 76), "K")
  # starting 1 kb upstream and ending 1 kb downstream changes nothing
  wide <- crt; wide$start <- wide$start - 1000; wide$end <- wide$end + 1000
  expect_equal(genomic_range_aa_positions(wide, cds)$aa_end, 425L)
  expect_equal(genomic_range_aa_positions(wide, cds)$aa_start, 1L)
  # the gene's own columns come along
  expect_true(all(c("gene_id", "Pf3D7_chrom", "name") %in% names(r)))
})

test_that("only codons held whole count unless partial = TRUE, and the flags say which end", {
  cds <- read_gff_cds(.gff_path())
  # pfcrt codon 76 is [403623, 403626); codon 77 follows at [403626, 403629)
  k76 <- data.frame(chr = "Pf3D7_07_v3", start = 403623, end = 403626)
  r <- genomic_range_aa_positions(k76, cds)
  expect_equal(c(r$aa_start, r$aa_end, r$n_aa), c(76L, 76L, 1L))

  mid <- data.frame(chr = "Pf3D7_07_v3", start = 403624, end = 403625)   # middle base only
  expect_false(genomic_range_aa_positions(mid, cds)$coding)
  expect_true(is.na(genomic_range_aa_positions(mid, cds)$aa_start))
  expect_equal(nrow(genomic_range_aa_positions(mid, cds, keep = "hits")), 0L)
  p <- genomic_range_aa_positions(mid, cds, partial = TRUE)
  expect_equal(c(p$aa_start, p$aa_end), c(76L, 76L))
  expect_true(p$partial_start); expect_true(p$partial_end)

  # from the second base of codon 76 through codon 78: 77-78 whole, 76 partial
  span <- data.frame(chr = "Pf3D7_07_v3", start = 403624, end = 403632)
  r <- genomic_range_aa_positions(span, cds)
  expect_equal(c(r$aa_start, r$aa_end, r$n_aa), c(77L, 78L, 2L))
  expect_false(r$partial_start)
  p <- genomic_range_aa_positions(span, cds, partial = TRUE)
  expect_equal(c(p$aa_start, p$aa_end, p$n_aa), c(76L, 78L, 3L))
  expect_true(p$partial_start); expect_false(p$partial_end)
})

test_that("a codon straddling an intron is whole only when both exons' bases are covered", {
  cds <- read_gff_cds(.gff_path())
  # pfcrt codon 178: bases 404109, 404110 (exon 3 end) and 404283 (exon 4 start), 1-based.
  # Exon 4 is 404283-404415, so codon 179 starts at 1-based 404284 = 0-based 404283.
  from_exon4 <- data.frame(chr = "Pf3D7_07_v3", start = 404282, end = 404415)
  r <- genomic_range_aa_positions(from_exon4, cds)
  expect_equal(r$aa_start, 179L)
  p <- genomic_range_aa_positions(from_exon4, cds, partial = TRUE)
  expect_equal(p$aa_start, 178L)
  expect_true(p$partial_start)
  # reaching back over the intron to the end of exon 3 makes 178 whole
  over <- data.frame(chr = "Pf3D7_07_v3", start = 404108, end = 404415)
  expect_equal(genomic_range_aa_positions(over, cds)$aa_start, 178L)
  # the intron alone covers nothing
  intron <- data.frame(chr = "Pf3D7_07_v3", start = 404110, end = 404282)
  expect_false(genomic_range_aa_positions(intron, cds)$coding)
})

test_that("minus-strand genes number from the high end, and the sequence reads in transcript order", {
  cds <- read_gff_cds(.gff_path())
  k13 <- PF3D7_GENES[PF3D7_GENES$name == "pfkelch13", ]
  r <- genomic_range_aa_positions(k13, cds, fasta = .fa_path())
  expect_equal(r$strand, "-")
  expect_equal(r$aa_start, 1L)
  expect_equal(nchar(r$aa_seq), r$n_aa)
  expect_true(startsWith(r$aa_seq, "M"))
  expect_equal(substr(r$aa_seq, 580, 580), "C")            # C580
  # the low-coordinate end of the gene is the C-terminus
  tail3 <- data.frame(chr = "Pf3D7_13_v3", start = k13$start, end = k13$start + 9)
  t <- genomic_range_aa_positions(tail3, cds)
  expect_equal(t$aa_end, r$aa_end)
  expect_gt(t$aa_start, 700)

  # sequence on the plus strand: pfcrt K76 in a 3-codon window
  win <- data.frame(chr = "Pf3D7_07_v3", start = 403620, end = 403629)
  s <- genomic_range_aa_positions(win, cds, fasta = .fa_path())
  expect_equal(c(s$aa_start, s$aa_end), c(75L, 77L))
  expect_equal(nchar(s$aa_seq), 3L)
  expect_equal(substr(s$aa_seq, 2, 2), "K")
  # without sequence there is no aa_seq column
  expect_false("aa_seq" %in% names(genomic_range_aa_positions(win, cds)))
})

test_that("pieces of a masked gene, several transcripts, keep = 'all' vs 'hits', input order", {
  cds <- read_gff_cds(.gff_path())
  crt <- PF3D7_GENES[PF3D7_GENES$name == "pfcrt", ]
  pieces <- bed_subtract(crt, tandem_repeats_to_avoid(pf3d7_tandem_repeats()), pad = 10)
  r <- genomic_range_aa_positions(pieces, cds, keep = "hits")
  expect_true(all(r$aa_start <= r$aa_end))
  expect_true(all(r$n_aa == r$aa_end - r$aa_start + 1))    # a piece covers a contiguous run
  expect_true(all(diff(r$aa_start) > 0))                    # in input (genomic) order
  expect_lt(sum(r$n_aa), 425)                               # the masked repeats cost residues
  expect_true(all(r$piece %in% pieces$piece))

  # two intervals, one non-coding, in the order given
  two <- data.frame(chr = c("Pf3D7_07_v3", "Pf3D7_07_v3", "Pf3D7_01_v3"),
                    start = c(403623, 100, 100), end = c(403626, 200, 200), tag = c("a", "b", "c"))
  all_rows <- genomic_range_aa_positions(two, cds)
  expect_equal(all_rows$tag, c("a", "b", "c"))
  expect_equal(all_rows$coding, c(TRUE, FALSE, FALSE))
  expect_equal(genomic_range_aa_positions(two, cds, keep = "hits")$tag, "a")
  expect_equal(genomic_range_aa_positions(two, cds, genes = NULL)$gene_name, rep(NA_character_, 3))

  # a gene with two isoforms gives one row per transcript
  lines <- c("##gff-version 3",
             "c1\tt\tCDS\t101\t112\t.\t+\t0\tID=a;Parent=T1;gene_id=G",
             "c1\tt\tCDS\t101\t106\t.\t+\t0\tID=b;Parent=T2;gene_id=G",
             "c1\tt\tCDS\t201\t206\t.\t+\t0\tID=c;Parent=T2;gene_id=G")
  f <- tempfile(fileext = ".gff"); writeLines(lines, f)
  iso <- genomic_range_aa_positions(data.frame(chr = "c1", start = 100, end = 300), f, genes = NULL)
  expect_equal(iso$transcript_id, c("T1", "T2"))
  expect_equal(iso$aa_end, c(4L, 4L))
  # a half-open end: [100, 106) is codons 1-2 of both
  expect_equal(genomic_range_aa_positions(data.frame(chr = "c1", start = 100, end = 106), f,
                                          genes = NULL)$aa_end, c(2L, 2L))
  expect_error(genomic_range_aa_positions(data.frame(chr = "c1", start = 10, end = 5), f), "inverted")
})
