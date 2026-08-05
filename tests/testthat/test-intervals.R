# bed_intersect() + the bundled genome-region datasets + duplicate-Name handling.

test_that("bed_intersect partitions overlap / only1 / only2", {
  a <- data.frame(id = c("x", "y", "z"), chr = c("7", "7", "8"),
                  start = c(100, 5000, 10), end = c(200, 5200, 20))
  b <- data.frame(chr = c("7", "8"), start = c(50, 900), end = c(1000, 1000))
  hit <- bed_intersect(a, b)
  expect_setequal(hit$overlap$id, "x")                 # x (100-200) inside 7:50-1000
  expect_setequal(hit$only1$id, c("y", "z"))           # y off the region, z on 8:10-20 (< 900)
  expect_equal(nrow(hit$only2), 1)                     # the 8:900-1000 region hit nothing
  expect_true(all(c("chr.2", "overlap_start", "overlap_end", "overlap_bp") %in% names(hit$overlap)))
  expect_equal(hit$overlap$overlap_start, 100)         # intersection = max(starts)
  expect_equal(hit$overlap$overlap_end, 200)
})

test_that("bed_intersect matches chromosome spellings and accepts a chrom alias", {
  genes <- data.frame(name = "g", chrom = "7", start = 100, end = 200)   # 'chrom' alias
  region <- data.frame(chr = "Pf3D7_07_v3", start = 50, end = 1000)      # long spelling
  hit <- bed_intersect(genes, region)
  expect_equal(nrow(hit$overlap), 1)                   # 7 == Pf3D7_07_v3 via normalise_chr
})

test_that("empty overlap still returns well-formed pieces", {
  a <- data.frame(chr = "7", start = 1, end = 2)
  b <- data.frame(chr = "7", start = 100, end = 200)
  hit <- bed_intersect(a, b)
  expect_equal(nrow(hit$overlap), 0)
  expect_equal(nrow(hit$only1), 1)
  expect_equal(nrow(hit$only2), 1)
})

test_that("core-region / paralog datasets classify genes as expected", {
  expect_true(all(c("Pf3D7_chrom", "start", "end", "chrom") %in% names(PF3D7_CORE_REGIONS)))
  expect_true(all(c("gene_id", "description") %in% names(PF3D7_PARALOG_GENES)))
  core <- bed_intersect(PF3D7_GENES, PF3D7_CORE_REGIONS)
  expect_true("PF3D7_0709000" %in% core$overlap$gene_id)     # pfcrt is core
  expect_true("PF3D7_0100100" %in% core$only1$gene_id)       # chr1 telomere var gene is not
  expect_gt(nrow(core$overlap), nrow(core$only1))            # most genes are core
})

test_that("repeated gene Names are disambiguated by gene_id (not collapsed)", {
  blocks <- data.frame(sample1 = c("s1", "s2"), sample2 = c("s2", "s3"),
                       chr = "Pf3D7_01_v3", start = c(1, 10), end = c(1e6, 20),
                       different = c(0, 1))
  meta <- data.frame(sample = c("s1", "s2", "s3"), region = c("A", "A", "B"))
  gtab <- PF3D7_GENES[PF3D7_GENES$name == "pfvar", ][1:3, ]   # 3 gene_ids, one Name
  ibd <- ibd_results(blocks = blocks, meta = meta, reference = "pf3d7")
  ov <- NULL
  expect_warning(ov <- gene_ibd_overlap(ibd, genes = gtab, group = "region"), "repeat")
  expect_length(levels(ov$gene), 3)                          # 3 distinct features, not 1
  expect_setequal(unique(ov$name), "pfvar")                  # raw name kept
  expect_setequal(unique(ov$gene_id), gtab$gene_id)          # unique ids kept
})
