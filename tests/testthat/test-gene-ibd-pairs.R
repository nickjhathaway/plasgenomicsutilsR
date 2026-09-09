
test_that("gene_ibd_pairs puts metadata on both ends of each pair", {
  blocks <- data.frame(
    sample1 = c("s1", "s1", "s2"), sample2 = c("s2", "s3", "s3"),
    chr = "7", start = c(395000, 395000, 395000), end = c(415000, 415000, 415000),
    different = 0, Nsnp = 40)
  meta <- data.frame(sample = paste0("s", 1:3),
                     region = factor(c("north", "north", "south"),
                                     levels = c("south", "north")),
                     country = c("UG", "TZ", "KE"), stringsAsFactors = FALSE)
  ibd <- ibd_results(blocks = blocks, meta = meta, group_col_in_meta = "region",
                     genes = PF_EXAMPLE_DRUG_GENES)

  p <- gene_ibd_pairs(ibd, genes = "pfcrt", add_meta_cols = c("region", "country"))
  expect_true(all(c("sample1_region", "sample2_region",
                    "sample1_country", "sample2_country") %in% names(p)))
  expect_equal(as.character(p$sample1_region), as.character(meta$region[match(p$sample1, meta$sample)]))
  expect_equal(p$sample2_country, meta$country[match(p$sample2, meta$sample)])
  expect_equal(levels(p$sample1_region), c("south", "north"))
  # the pairs themselves are unchanged by asking for the labels
  expect_equal(p[, seq_len(ncol(gene_ibd_pairs(ibd, genes = "pfcrt")))],
               gene_ibd_pairs(ibd, genes = "pfcrt"))
  # meta is taken from the object, and can be overridden
  expect_equal(ibd$gene_ibd_pairs(genes = "pfcrt", add_meta_cols = "region"), p[, -(ncol(p) - 1:0)])
  expect_error(gene_ibd_pairs(ibd, genes = "pfcrt", add_meta_cols = "nope"), "no column 'nope'")
})

test_that("an empty gene_ibd_pairs result still carries the requested columns", {
  # a gene nothing is IBD over: the columns must line up so an rbind of several genes works
  blocks <- data.frame(sample1 = "s1", sample2 = "s2", chr = "7",
                       start = 100000, end = 130000, different = 0, Nsnp = 40)
  meta <- data.frame(sample = c("s1", "s2"), region = c("north", "south"))
  ibd <- ibd_results(blocks = blocks, meta = meta, group_col_in_meta = "region",
                     genes = PF_EXAMPLE_DRUG_GENES)
  e <- gene_ibd_pairs(ibd, genes = "pfcrt", add_meta_cols = "region")
  expect_equal(nrow(e), 0)
  expect_true(all(c("sample1_region", "sample2_region") %in% names(e)))
})
