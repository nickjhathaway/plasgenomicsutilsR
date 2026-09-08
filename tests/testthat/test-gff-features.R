# read_gff_features(): a gene's exons, introns, CDS or span as 0-based half-open intervals.
# Checked against the real pfcrt annotation (13 exons, plus strand), the minus-strand
# pfkelch13, and a hand-built GFF with two isoforms and Ensembl-style attributes.

.drug_gff <- function() system.file("extdata", "pf3d7_drug_gene_cds.gff",
                                    package = "plasgenomicsutilsR")

test_that("exons come back 0-based half-open, numbered in transcript orientation", {
  ex <- read_gff_features(.drug_gff(), "exon", ids = "pfcrt")
  expect_s3_class(ex, "tbl_df")
  expect_equal(nrow(ex), 13L)
  expect_equal(unique(ex$name), "pfcrt")
  expect_equal(unique(ex$gene_id), "PF3D7_0709000")
  expect_equal(ex$index, 1:13)
  # GFF exon 1 is 402385-403312 (1-based inclusive) -> [402384, 403312)
  expect_equal(c(ex$start[1], ex$end[1]), c(402384, 403312))
  expect_equal(ex$width, ex$end - ex$start)
  expect_equal(unique(ex$chrom), "Pf3D7_07_v3")
  expect_equal(unique(ex$chr), "7")
  expect_false("transcript_id" %in% names(ex))

  # minus strand: exon 1 has the highest coordinates
  k13 <- read_gff_features(.drug_gff(), "exon", ids = "pfkelch13")
  expect_equal(unique(k13$strand), "-")
  expect_equal(k13$index[which.max(k13$start)], 1L)
  expect_equal(k13$index[which.min(k13$start)], nrow(k13))
})

test_that("introns are exactly the gaps between exons, and cds / gene are what they say", {
  ex <- read_gff_features(.drug_gff(), "exon", ids = "pfcrt")
  intr <- read_gff_features(.drug_gff(), "intron", ids = "pfcrt")
  expect_equal(nrow(intr), 12L)
  expect_equal(intr$start, ex$end[-13])
  expect_equal(intr$end, ex$start[-1])
  expect_equal(unique(intr$feature), "intron")
  # exons + introns tile the gene span exactly
  gene <- read_gff_features(.drug_gff(), "gene", ids = "pfcrt")
  expect_equal(nrow(gene), 1L)
  expect_equal(sum(ex$width) + sum(intr$width), gene$width)
  expect_equal(c(gene$start, gene$end), c(min(ex$start), max(ex$end)))
  expect_equal(nrow(bed_intersect(ex, intr)$overlap), 0L)

  # the CDS lies within the exons and is shorter (UTR in exon 1 and exon 13)
  cds <- read_gff_features(.drug_gff(), "cds", ids = "pfcrt")
  expect_equal(nrow(cds), 13L)
  expect_lt(sum(cds$width), sum(ex$width))
  expect_equal(nrow(bed_subtract(cds, ex)), 0L)
  # and agrees with read_gff_cds() once that reader's 1-based starts are shifted
  raw <- read_gff_cds(.drug_gff())
  raw <- raw[raw$gene_id == "PF3D7_0709000", ]
  expect_equal(sort(cds$start), sort(raw$start - 1))
  expect_equal(sort(cds$end), sort(raw$end))

  # a single-exon gene has no introns
  expect_error(read_gff_features(.drug_gff(), "intron", ids = "pfdhfr"), "single exon")
})

test_that("ids resolve gene ids, GFF Names, package names and transcript ids", {
  by_id <- read_gff_features(.drug_gff(), "gene", ids = "PF3D7_0709000")
  by_gff_name <- read_gff_features(.drug_gff(), "gene", ids = "crt")
  by_pkg_name <- read_gff_features(.drug_gff(), "gene", ids = "pfcrt")
  by_tx <- read_gff_features(.drug_gff(), "gene", ids = "PF3D7_0709000.1")
  expect_equal(by_gff_name, by_id)
  expect_equal(by_pkg_name, by_id)
  expect_equal(by_tx, by_id)
  # the GFF's own Name is used when no gene table is given
  expect_equal(read_gff_features(.drug_gff(), "gene", ids = "CRT", genes = NULL)$name, "CRT")
  expect_error(read_gff_features(.drug_gff(), "gene", ids = "pfnope"), "pfnope")
  expect_error(read_gff_features(.drug_gff(), "gene", ids = "pfnope", genes = NULL),
               "not a gene id, gene `Name` or transcript id")

  all6 <- read_gff_features(.drug_gff(), "gene")
  expect_equal(nrow(all6), 6L)
  expect_setequal(all6$name, c("pfdhfr", "pfmdr1", "pfaat1", "pfcrt", "pfdhps", "pfkelch13"))
  expect_equal(all6$gene_id[1], "PF3D7_0417200")               # GFF order
})

test_that("isoforms merge per gene and stay apart per transcript; Ensembl attributes work", {
  # G1 has two isoforms: T1 exons 101-200 + 301-400, T2 exons 101-200 + 251-400 (1-based),
  # so per gene the exons are [100,200) and [250,400) and the only intron is [200,250)
  lines <- c("##gff-version 3",
             "c1\tt\tgene\t101\t400\t.\t+\t.\tID=gene:G1;Name=Alpha",
             "c1\tt\tmRNA\t101\t400\t.\t+\t.\tID=transcript:T1;Parent=gene:G1",
             "c1\tt\tmRNA\t101\t400\t.\t+\t.\tID=transcript:T2;Parent=gene:G1",
             "c1\tt\texon\t101\t200\t.\t+\t.\tParent=transcript:T1",
             "c1\tt\texon\t301\t400\t.\t+\t.\tParent=transcript:T1",
             "c1\tt\texon\t101\t200\t.\t+\t.\tParent=transcript:T2",
             "c1\tt\texon\t251\t400\t.\t+\t.\tParent=transcript:T2",
             "c1\tt\tgene\t901\t950\t.\t-\t.\tID=gene:G2",
             "c1\tt\tmRNA\t901\t950\t.\t-\t.\tID=transcript:T3;Parent=gene:G2",
             "c1\tt\texon\t901\t950\t.\t-\t.\tParent=transcript:T3")
  f <- tempfile(fileext = ".gff"); writeLines(lines, f)

  ex <- read_gff_features(f, "exon", genes = NULL)
  expect_equal(ex$gene_id, c("G1", "G1", "G2"))
  expect_equal(ex$name, c("Alpha", "Alpha", NA))
  expect_equal(ex$start, c(100, 250, 900))
  expect_equal(ex$end, c(200, 400, 950))
  intr <- read_gff_features(f, "intron", ids = "G1", genes = NULL)
  expect_equal(c(intr$start, intr$end), c(200, 250))

  tx <- read_gff_features(f, "exon", per = "transcript", genes = NULL)
  expect_equal(tx$transcript_id, c("T1", "T1", "T2", "T2", "T3"))
  expect_equal(tx$start, c(100, 300, 100, 250, 900))
  ti <- read_gff_features(f, "intron", per = "transcript", ids = "G1", genes = NULL)
  expect_equal(ti$transcript_id, c("T1", "T2"))
  expect_equal(ti$width, c(100, 50))
  # naming a transcript restricts the gene to it; naming the gene keeps both
  expect_equal(read_gff_features(f, "exon", per = "transcript", ids = "T2", genes = NULL)$transcript_id,
               c("T2", "T2"))
  expect_equal(nrow(read_gff_features(f, "exon", per = "transcript", ids = "G1", genes = NULL)), 4L)
  expect_equal(read_gff_features(f, "gene", genes = NULL)$name, c("Alpha", NA))
  expect_error(read_gff_features(f, "cds", genes = NULL), "no CDS features")

  # a GFF with CDS but no exon rows uses the CDS as exons, and says so
  cds_only <- tempfile(fileext = ".gff")
  writeLines(sub("\texon\t", "\tCDS\t", lines), cds_only)
  expect_message(cx <- read_gff_features(cds_only, "exon", genes = NULL), "using its CDS")
  expect_equal(cx$start, ex$start)
})

test_that("the design flow: gene -> exons -> minus the tandem repeats", {
  exons <- read_gff_features(.drug_gff(), "exon", ids = c("pfcrt", "pfkelch13"))
  mask <- tandem_repeats_to_avoid(pf3d7_tandem_repeats())
  targets <- bed_subtract(exons, mask, pad = 10)
  expect_true(all(c("gene_id", "name", "index", "piece", "width") %in% names(targets)))
  expect_lt(sum(targets$width), sum(exons$width))
  expect_equal(nrow(bed_intersect(targets, mask)$overlap), 0L)
  f <- tempfile(fileext = ".bed")
  write_bed(targets, f)
  expect_true(all(grepl("^Pf3D7_(07|13)_v3\t", readLines(f))))
})
