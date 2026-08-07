
# --------------------------------------------------------------------------- #
#  subset_groups() / restrict_groups()                                         #
# --------------------------------------------------------------------------- #

.group_fixture <- function() {
  blocks <- data.frame(
    sample1 = c("a", "a", "b", "c", "d", "e"),
    sample2 = c("b", "c", "c", "d", "e", "f"),
    chr = "Pf3D7_07_v3", start = 400000, end = 460000,
    Nsnp = 40L, different = 0L, stringsAsFactors = FALSE)
  meta <- data.frame(sample = letters[1:6], region = c("N", "N", "N", "S", "S", "S"),
                     stringsAsFactors = FALSE)
  ibd_results(blocks = blocks, meta = meta, group_col_in_meta = "region",
              genes = PF_EXAMPLE_DRUG_GENES, min_block_snp = 0, min_block_kb = 0)
}

test_that("keep and drop each narrow the aggregated tables, and compose", {
  ibd <- example_ibd_results()
  all_g <- levels(factor(ibd$get_selection()$group))
  skip_if(length(all_g) < 4)

  d <- ibd$subset_groups(drop = all_g[1])
  expect_setequal(levels(factor(d$get_selection()$group)), all_g[-1])
  expect_setequal(levels(factor(d$get_per_snp_group()$group)), all_g[-1])

  k <- ibd$subset_groups(keep = all_g[1:2])
  expect_setequal(levels(factor(k$get_selection()$group)), all_g[1:2])

  # drop is applied after keep, so the two compose rather than conflict
  kd <- ibd$subset_groups(keep = all_g[1:3], drop = all_g[1])
  expect_setequal(levels(factor(kd$get_selection()$group)), all_g[2:3])

  # subset_groups() copies; the receiver is untouched
  expect_setequal(levels(factor(ibd$get_selection()$group)), all_g)
  # ...and restrict_groups() is the in-place form
  r <- ibd$clone(deep = TRUE)
  r$restrict_groups(drop = all_g[1])
  expect_setequal(levels(factor(r$get_selection()$group)), all_g[-1])
})

test_that("a group pair survives only when both of its groups do", {
  ibd <- example_ibd_results()
  all_g <- levels(factor(ibd$get_selection()$group))
  skip_if(is.null(ibd$get_pairwise_group()) || length(all_g) < 3)

  keep <- all_g[1:2]
  pw <- ibd$subset_groups(keep = keep)$get_pairwise_group()
  expect_true(all(as.character(pw$group_a) %in% keep))
  expect_true(all(as.character(pw$group_b) %in% keep))
  # k groups give k(k+1)/2 unordered pairs, so the row count scales by that ratio
  n_all <- length(all_g)
  expect_equal(nrow(pw) / nrow(ibd$get_pairwise_group()),
               (2 * 3 / 2) / (n_all * (n_all + 1) / 2), tolerance = 1e-6)
})

test_that("dropping a group narrows its samples too, and trims the group order", {
  ibd <- .group_fixture()
  expect_equal(ibd$get_group_order(), c("N", "S"))

  n <- ibd$subset_groups(drop = "S")
  expect_setequal(n$get_analyzed_samples(), c("a", "b", "c"))
  expect_equal(nrow(n$get_blocks()), 3)            # a-b, a-c, b-c
  expect_equal(nrow(n$get_meta()), 3)
  # the dropped group leaves the order rather than lingering as an unused level
  expect_equal(n$get_group_order(), "N")
  # block-derived output follows
  expect_lt(nrow(gene_ibd_pairs(n)), nrow(gene_ibd_pairs(ibd)))
  expect_equal(nrow(ibd$get_blocks()), 6)          # original untouched
})

test_that("an unknown group warns and an empty result errors", {
  ibd <- .group_fixture()
  expect_warning(ibd$subset_groups(drop = "Nowhere"), "not in these results")
  expect_error(suppressWarnings(ibd$subset_groups(keep = "Nowhere")), "no groups")
  expect_error(ibd$subset_groups(drop = c("N", "S")), "no groups")
  # neither argument is a no-op rather than an error
  expect_equal(nrow(ibd$subset_groups()$get_blocks()), 6)
})

test_that("a group subset still plots", {
  skip_if_not_installed("ggplot2")
  ibd <- example_ibd_results()
  sub <- ibd$subset_groups(keep = levels(factor(ibd$get_selection()$group))[1:2])
  expect_s3_class(plot_selection_manhattan(sub), "ggplot")
  expect_s3_class(plot_ibd_sharing_manhattan(sub), "ggplot")
  expect_s3_class(plot_ibd_pairwise_group_heatmap(sub), "ggplot")
})

test_that("a group only `meta` knows about does not have to be ordered", {
  # meta is routinely a superset of the groups the result tables carry: a group with a
  # single sample contributes no within-group pair, so no per-SNP or selection row
  sel <- data.frame(group = rep(c("N", "S"), each = 3), chr = "Pf3D7_07_v3",
                    pos = rep(c(1000, 2000, 3000), 2), neg_log10_p = seq(1, 8, length.out = 6),
                    stringsAsFactors = FALSE)
  meta <- data.frame(sample = letters[1:5], region = c("N", "N", "S", "S", "W"),
                     stringsAsFactors = FALSE)
  ibd <- ibd_results(selection = sel, meta = meta, group_col_in_meta = "region")

  # ordering only the groups that appear in the results is accepted...
  expect_warning(ibd$set_group_order(c("S", "N")), "does not name meta group")
  expect_equal(levels(ibd$get_selection()$group)[1:2], c("S", "N"))
  # ...and the meta-only group keeps its level, so its sample is not silently NA'd out of
  # the block-derived output
  expect_true("W" %in% ibd$get_group_order())
  expect_false(anyNA(ibd$get_meta()$region))

  # a group that IS in the results still may not be omitted
  expect_error(ibd$set_group_order("N"), "missing from")
})

test_that("subsetting sees groups that only `meta` carries", {
  ibd <- .group_fixture()                      # blocks + meta only, no aggregated tables
  expect_silent(sub <- ibd$subset_groups(drop = "S"))
  expect_equal(sub$get_group_order(), "N")
  expect_setequal(sub$get_analyzed_samples(), c("a", "b", "c"))
  expect_warning(ibd$subset_groups(drop = "Nowhere"), "not in these results")
})
