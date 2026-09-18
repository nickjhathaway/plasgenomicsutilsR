
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

test_that("gene_ibd_pairs clusters over the blocks `sharing` selects, not the whole gene's", {
  # pfcrt is Pf3D7_07_v3:403221-406317. Two groups whose members share the whole gene, plus
  # one partial block bridging them -- the edge plot_ibd_network(sharing = "complete") never
  # draws, and the one that decides whether the two groups are one cluster or two.
  # blocks have to clear ibd_results()'s 15 kb / 15 SNP floor to survive at all
  spanning <- function(a, b) data.frame(sample1 = a, sample2 = b, chr = "7",
                                        start = 395000, end = 420000)
  blocks <- rbind(
    spanning("s1", "s2"), spanning("s1", "s3"), spanning("s2", "s3"),   # group A
    spanning("s4", "s5"),                                              # group B
    data.frame(sample1 = "s3", sample2 = "s4", chr = "7",              # the partial bridge
               start = 390000, end = 405000))                           # stops inside pfcrt
  blocks$different <- 0
  blocks$Nsnp <- 40
  ibd <- ibd_results(blocks = blocks, genes = PF_EXAMPLE_DRUG_GENES)

  p <- gene_ibd_pairs(ibd, genes = "pfcrt")
  expect_equal(sum(p$coverage == "complete"), 4L)
  expect_equal(sum(p$coverage == "partial"), 1L)
  # the stored ids are the overlap clustering: the partial bridge merges both groups
  expect_equal(unique(p$gene_cluster_id), 1L)
  expect_equal(unique(p$gene_cluster_size), 5L)
  # and filtering the table afterwards does NOT recompute them -- still one merged cluster
  expect_equal(unique(subset(p, coverage == "complete")$gene_cluster_id), 1L)
  expect_equal(unique(subset(p, coverage == "complete")$gene_cluster_size), 5L)

  # asking for complete sharing up front re-clusters over just those blocks: two groups
  cp <- gene_ibd_pairs(ibd, genes = "pfcrt", sharing = "complete")
  expect_true(all(cp$coverage == "complete"))
  expect_equal(nrow(cp), 4L)
  # same pairs as the filtered table, only the cluster columns differ
  drop_cl <- function(d) d[, setdiff(names(d), c("gene_cluster_id", "gene_cluster_size"))]
  expect_equal(drop_cl(cp), drop_cl(subset(p, coverage == "complete")))
  expect_equal(sort(unique(cp$gene_cluster_id)), c(1L, 2L))
  gp <- unique(data.frame(id = cp$gene_cluster_id, size = cp$gene_cluster_size))
  expect_equal(gp$size[gp$id == 1], 3L)   # s1, s2, s3
  expect_equal(gp$size[gp$id == 2], 2L)   # s4, s5
  # the bridged pair is gone, so s3 and s4 are no longer in the same cluster
  expect_equal(nrow(subset(cp, sample1 == "s3" & sample2 == "s4")), 0L)
  expect_false(cp$gene_cluster_id[cp$sample1 == "s1"][1] ==
                 cp$gene_cluster_id[cp$sample1 == "s4"][1])
})

test_that("sharing = 'complete' asks a block to span the padded interval", {
  # a block covering pfcrt (403221-406317) but not pfcrt +- 20 kb
  blocks <- data.frame(sample1 = "s1", sample2 = "s2", chr = "7",
                       start = 400000, end = 420000, different = 0, Nsnp = 40)
  ibd <- ibd_results(blocks = blocks, genes = PF_EXAMPLE_DRUG_GENES)

  expect_equal(nrow(gene_ibd_pairs(ibd, genes = "pfcrt", sharing = "complete")), 1L)
  # `coverage` is measured against the gene alone, so padding does not change it ...
  w <- gene_ibd_pairs(ibd, genes = "pfcrt", within = 20000)
  expect_equal(w$coverage, "complete")
  # ... but the edge rule, like plot_ibd_network()'s, is against the padded interval
  n <- gene_ibd_pairs(ibd, genes = "pfcrt", within = 20000, sharing = "complete")
  expect_equal(nrow(n), 0L)
  expect_true(all(c("gene_cluster_id", "gene_cluster_size") %in% names(n)))
})
