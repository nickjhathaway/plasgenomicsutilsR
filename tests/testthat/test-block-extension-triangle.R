blk <- function(pairs, chr, start, len) {
  data.frame(sample1 = pairs$sample1, sample2 = pairs$sample2, chr = chr,
             start = start, end = start + len - 1,
             different = 0, Nsnp = 40, stringsAsFactors = FALSE)
}
spread <- function(len, n) len + (seq_len(n) - (n + 1) / 2) * 1000
LOCUS <- data.frame(name = "sweep", chr = "7", start = 400000, end = 410000,
                    stringsAsFactors = FALSE)

# north pairs share 80 kb over the locus, south pairs 40 kb, and the north-south pairs that
# do share get 80 kb -- so the triangle has three distinct cells with known values
mk_two_group <- function(n_cross = 5) {
  n <- data.frame(sample1 = c("n1", "n3", "n5", "n7", "n9"),
                  sample2 = c("n2", "n4", "n6", "n8", "n10"), stringsAsFactors = FALSE)
  s <- data.frame(sample1 = c("s1", "s3", "s5", "s7", "s9"),
                  sample2 = c("s2", "s4", "s6", "s8", "s10"), stringsAsFactors = FALSE)
  x <- data.frame(sample1 = c("n1", "n3", "n5", "n7", "n9")[seq_len(n_cross)],
                  sample2 = c("s1", "s3", "s5", "s7", "s9")[seq_len(n_cross)],
                  stringsAsFactors = FALSE)
  mk <- function(p, len) {
    k <- nrow(p)
    rbind(blk(p, "7", 405000 - spread(len, k) / 2, spread(len, k)),
          blk(p, "1", 100000, 20000), blk(p, "2", 100000, 20000))
  }
  bl <- rbind(mk(n, 80000), mk(s, 40000), mk(x, 80000))
  meta <- data.frame(sample = unique(c(bl$sample1, bl$sample2)), stringsAsFactors = FALSE)
  meta$region <- ifelse(grepl("^n", meta$sample), "north", "south")
  ibd_results(blocks = bl, meta = meta, group_col_in_meta = "region")
}

test_that('pairs = "all" adds the cross-group cells and keeps the diagonal identical', {
  ibd <- mk_two_group()
  w <- ibd_block_extension_test(ibd, LOCUS)
  a <- ibd_block_extension_test(ibd, LOCUS, pairs = "all")
  b <- ibd_block_extension_test(ibd, LOCUS, pairs = "between")

  expect_true(all(c("group_a", "group_b") %in% names(a)))
  expect_false("group" %in% names(a))
  expect_equal(nrow(a), 3L)                      # north-north, south-south, north-south
  expect_equal(nrow(b), 1L)
  expect_true(all(b$group_a != b$group_b))

  # the diagonal of "all" is exactly the "within" result
  diag <- a[a$group_a == a$group_b, ]
  expect_equal(sort(diag$paired_ratio), sort(w$paired_ratio))
  expect_equal(sum(diag$n_pairs), sum(w$n_pairs))
  expect_equal(a$paired_ratio[a$group_a == "north" & a$group_b == "north"], 4)
  expect_equal(a$paired_ratio[a$group_a == "south" & a$group_b == "south"], 2)
  expect_equal(b$paired_ratio, 4)                # cross pairs share the long haplotype
  expect_equal(b$n_pairs, 5L)
})

test_that("a cross-group pair is measured against its own baseline, like any other", {
  ibd <- mk_two_group()
  a <- ibd_block_extension_test(ibd, LOCUS, pairs = "all")
  pr <- attr(a, "pair_ratios")
  expect_true(all(c("group_a", "group_b") %in% names(pr)))
  cross <- pr[pr$group_a != pr$group_b, ]
  expect_equal(nrow(cross), 5L)
  expect_equal(unique(cross$ref_median), 20000)  # its own segments off chromosome 7
  expect_equal(stats::median(cross$log2_ratio), 2)
})

test_that("the triangle draws every group pair, greying the ones with too few pairs", {
  skip_if_not_installed("ggplot2")
  ibd <- mk_two_group()
  p <- plot_pairwise_block_extension(ibd, LOCUS, min_pairs = 1)
  expect_s3_class(p, "ggplot")
  d <- p$data
  expect_equal(nrow(d), 3L)                      # both diagonals plus the one off-diagonal
  expect_setequal(as.character(d$group_a), c("north", "south"))
  expect_true(all(is.finite(d$paired_ratio)))

  # a cohort where only 3 pairs span the two groups: the diagonal is testable, the
  # off-diagonal is not, and it must be drawn grey rather than dropped or shown as ~1
  thin <- mk_two_group(n_cross = 3)
  d2 <- plot_pairwise_block_extension(thin, LOCUS, min_pairs = 5)$data
  expect_equal(nrow(d2), 3L)                     # the cell is still drawn
  cross <- d2[d2$group_a != d2$group_b, ]
  expect_equal(nrow(cross), 1L)
  expect_true(is.na(cross$paired_ratio))         # blank, not a number
  expect_true(is.na(cross$n_pairs))
  expect_equal(sum(!is.na(d2$paired_ratio)), 2L) # both diagonal cells survive
})

test_that("the fill limits are symmetric about 1 so the midpoint means no extension", {
  skip_if_not_installed("ggplot2")
  ibd <- mk_two_group()
  p <- plot_pairwise_block_extension(ibd, LOCUS, min_pairs = 1)
  # ggplot2 stores the limits already through the transform, so they are on the log2 scale
  lim <- p$scales$get_scales("fill")$limits
  expect_equal(lim[1], -lim[2])
  expect_gt(lim[2], 0)
  # a ratio of 4 against a ratio of 2 puts the wider one at the edge
  expect_equal(2^lim[2], 4)
})

test_that("individual = TRUE gives one plot per locus, and the method matches", {
  skip_if_not_installed("ggplot2")
  ibd <- mk_two_group()
  two <- rbind(LOCUS, data.frame(name = "other", chr = "1", start = 100000, end = 110000))
  l <- plot_pairwise_block_extension(ibd, two, min_pairs = 1, individual = TRUE)
  expect_named(l, c("sweep", "other"))
  expect_s3_class(l[[1]], "ggplot")
  expect_s3_class(ibd$plot_pairwise_block_extension(LOCUS, min_pairs = 1), "ggplot")
})
