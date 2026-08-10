# A hand-built GFF, so the expected coordinates can be worked out by hand: a plus-strand
# transcript whose first exon is 5 bp (so codon 2 straddles the intron) and a minus-strand one.
tiny_gff <- function() {
  lines <- c(
    "##gff-version 3",
    paste("chr1\ttest\tCDS\t101\t105\t.\t+\t0",
          "ID=T1-CDS1;Parent=T1;gene_id=G1", sep = "\t"),
    paste("chr1\ttest\tCDS\t201\t210\t.\t+\t1",
          "ID=T1-CDS2;Parent=T1;gene_id=G1", sep = "\t"),
    paste("chr2\ttest\tCDS\t501\t512\t.\t-\t0",
          "ID=T2-CDS1;Parent=T2;gene_id=G2", sep = "\t"))
  f <- tempfile(fileext = ".gff"); writeLines(lines, f); f
}

test_that("read_gff_cds pulls the coding exons and their transcript", {
  cds <- read_gff_cds(tiny_gff())
  expect_equal(nrow(cds), 3)
  expect_setequal(cds$transcript_id, c("T1", "T1", "T2"))
  expect_setequal(cds$gene_id, c("G1", "G1", "G2"))
  expect_equal(cds$phase, c(0L, 1L, 0L))
  expect_error(read_gff_cds(tempfile()), "path to a GFF")
})

test_that("a codon maps to the bases it is actually made of", {
  cds <- read_gff_cds(tiny_gff())
  aa <- function(tx, k) aa_intervals(data.frame(transcript_id = tx, aa_position = k),
                                     cds, genes = NULL, one_based_output = TRUE)

  # plus strand, first exon 101-105: codon 1 = 101,102,103
  r1 <- aa(  "T1", 1)
  expect_equal(r1$codon_positions, "101,102,103")
  expect_false(r1$spans_intron)
  # codon 2 = 104,105 then the exon runs out, so it continues into the second exon. That exon
  # has phase 1, but phase only matters on the FIRST exon of the transcript -- the codon
  # carries on from where the previous one stopped.
  r2 <- aa("T1", 2)
  expect_equal(r2$codon_positions, "104,105,201")
  expect_true(r2$spans_intron)
  expect_equal(r2$n_exons, 2L)
  # start/end span the intron, so the width is not three bases
  expect_equal(c(r2$start, r2$end), c(104, 201))
  # codon 3 is wholly in the second exon
  expect_equal(aa("T1", 3)$codon_positions, "202,203,204")

  # minus strand: transcript order runs from the highest coordinate down
  m <- aa("T2", 1)
  expect_equal(m$strand, "-")
  expect_equal(m$codon_positions, "510,511,512")
  expect_equal(aa("T2", 2)$codon_positions, "507,508,509")
})

test_that("phase on the first exon shifts the frame", {
  # same first exon, but the transcript starts mid-codon: 2 bases are dropped
  lines <- c("##gff-version 3",
             paste("chr1\tt\tCDS\t101\t112\t.\t+\t2", "ID=P-CDS1;Parent=P;gene_id=GP",
                   sep = "\t"))
  f <- tempfile(fileext = ".gff"); writeLines(lines, f)
  r <- aa_intervals(data.frame(transcript_id = "P", aa_position = 1), read_gff_cds(f),
                    genes = NULL, one_based_output = TRUE)
  expect_equal(r$codon_positions, "103,104,105")
})

test_that("the interval is 0-based half-open by default, like the rest of the package", {
  cds <- read_gff_cds(tiny_gff())
  zero <- aa_intervals(data.frame(transcript_id = "T1", aa_position = 1), cds, genes = NULL)
  one <- aa_intervals(data.frame(transcript_id = "T1", aa_position = 1), cds, genes = NULL,
                      one_based_output = TRUE)
  expect_equal(zero$start, one$start - 1L)
  expect_equal(zero$end, one$end)
  expect_equal(zero$end - zero$start, 3)
  # and it is shaped like an interval table, so the gene-track helper accepts it
  expect_silent(.gene_track(zero))
  expect_equal(zero$name, "T1-AA1")
})

test_that("a transcript can be named by id, gene or symbol", {
  cds <- read_gff_cds(tiny_gff())
  genes <- data.frame(name = "mygene", gene_id = "G1")
  by_tx <- aa_intervals(data.frame(transcript_id = "T1", aa_position = 1), cds, genes = NULL)
  by_gene <- aa_intervals(data.frame(transcript_id = "G1", aa_position = 1), cds, genes = NULL)
  by_sym <- aa_intervals(data.frame(transcript_id = "mygene", aa_position = 1), cds,
                         genes = genes)
  expect_equal(by_gene$codon_positions, by_tx$codon_positions)
  expect_equal(by_sym$codon_positions, by_tx$codon_positions)
  expect_equal(by_sym$name, "T1-AA1")     # named for the transcript actually used
})

test_that("what cannot be placed is reported, not silently dropped", {
  cds <- read_gff_cds(tiny_gff())
  # T1 codes 15 bases = 5 codons once the phase is taken off, so 99 is past the end
  expect_warning(
    expect_error(aa_intervals(data.frame(transcript_id = "T1", aa_position = 99), cds,
                              genes = NULL), "no amino-acid position could be placed"),
    "past the end of the protein")
  expect_warning(aa_intervals(data.frame(transcript_id = c("T1", "nope"),
                                         aa_position = c(1, 1)), cds, genes = NULL),
                 "no transcript found for: nope")
  expect_error(aa_intervals(data.frame(x = 1), cds), "needs `transcript_id`")
  expect_error(aa_intervals(data.frame(transcript_id = "T1", aa_position = -1), cds),
               "whole number 1 or greater")
})

test_that("published resistance codons land on their known SNPs", {
  # the real Pf3D7 GFF is not shipped, so this only runs where it is available
  gff <- "/Users/nhathaway/Documents/tank/data/genomes/plasmodium/genomes/pf/info/gff/Pf3D7.gff"
  skip_if_not(file.exists(gff), "Pf3D7 GFF not available")
  cds <- read_gff_cds(gff)
  aa <- data.frame(
    transcript_id = c("pfcrt", "pfdhps", "pfdhps", "pfkelch13", "pfmdr1", "pfdhfr"),
    aa_position   = c(76,      437,      540,      580,         86,       108))
  r <- aa_intervals(aa, cds, one_based_output = TRUE)
  known <- c(403625, 549685, 549993, 1725259, 958145, 748410)   # published Pf3D7 positions
  for (i in seq_len(nrow(r))) {
    bases <- as.numeric(strsplit(r$codon_positions[i], ",")[[1]])
    expect_true(known[i] %in% bases,
                info = paste(r$name[i], "should cover", known[i]))
  }
  expect_equal(r$strand[r$transcript_id == "PF3D7_1343700.1"], "-")   # kelch13 is minus
})
