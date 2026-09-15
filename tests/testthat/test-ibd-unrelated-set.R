
# a pair table in the shape .pair_edges() expects, from an edge list of "related" pairs
.pairs_from_edges <- function(samples, from, to, w = 0.5, background = 0) {
  all_pairs <- t(utils::combn(samples, 2))
  df <- data.frame(sample1 = all_pairs[, 1], sample2 = all_pairs[, 2],
                   ibd_fraction_accessible = background, stringsAsFactors = FALSE)
  for (i in seq_along(from)) {
    hit <- (df$sample1 == from[i] & df$sample2 == to[i]) |
           (df$sample1 == to[i] & df$sample2 == from[i])
    df$ibd_fraction_accessible[hit] <- w
  }
  df
}

test_that("the selected set is valid and maximal", {
  set.seed(3)
  s <- paste0("s", 1:40)
  all_pairs <- t(utils::combn(s, 2))
  df <- data.frame(sample1 = all_pairs[, 1], sample2 = all_pairs[, 2],
                   ibd_fraction_accessible = stats::runif(nrow(all_pairs), 0, 0.06),
                   stringsAsFactors = FALSE)

  u <- ibd_unrelated_set(df, max_ibd = 0.02, restarts = 50)
  sel <- u$sample[u$selected]

  # valid: no pair inside the set is above the cutoff
  inside <- df[df$sample1 %in% sel & df$sample2 %in% sel, ]
  expect_true(all(inside$ibd_fraction_accessible <= 0.02))
  expect_lte(attr(u, "max_ibd_within_set"), 0.02)
  # and it reports the real sharing inside the set, not just that the cutoff held
  expect_equal(attr(u, "max_ibd_within_set"), max(inside$ibd_fraction_accessible))
  expect_equal(attr(u, "set_size"), length(sel))

  # maximal: every unselected sample is blocked by something that was selected
  expect_true(all(!is.na(u$blocked_by[!u$selected])))
  expect_true(all(u$blocked_by[!u$selected] %in% sel))
  expect_true(all(is.na(u$blocked_by[u$selected])))
})

test_that("a path of five finds the three-node independent set", {
  s <- paste0("s", 1:5)
  df <- .pairs_from_edges(s, from = c("s1", "s2", "s3", "s4"), to = c("s2", "s3", "s4", "s5"))
  u <- ibd_unrelated_set(df, max_ibd = 0.01, restarts = 30)
  expect_equal(sum(u$selected), 3L)
  expect_setequal(u$sample[u$selected], c("s1", "s3", "s5"))
})

test_that("a star beats keeping one sample per single-linkage cluster", {
  # the hub joins every leaf, so this is ONE component: one-per-cluster would keep 1 sample.
  # the unrelated set is the 10 leaves, which is the whole reason this function exists.
  s <- c("hub", paste0("leaf", 1:10))
  df <- .pairs_from_edges(s, from = rep("hub", 10), to = paste0("leaf", 1:10))
  u <- ibd_unrelated_set(df, max_ibd = 0.01, restarts = 30)

  expect_equal(sum(u$selected), 10L)
  expect_false("hub" %in% u$sample[u$selected])
  # and the clustering really does call it a single component
  cl <- ibd_pair_clusters(df, min_ibd = 0.01)
  expect_equal(length(unique(stats::na.omit(cl$cluster_id))), 1L)
})

test_that("samples related to nobody are always selected", {
  s <- c("a", "b", "loner1", "loner2")
  df <- .pairs_from_edges(s, from = "a", to = "b")
  u <- ibd_unrelated_set(df, max_ibd = 0.01, restarts = 10)
  expect_true(all(u$selected[u$sample %in% c("loner1", "loner2")]))
  expect_equal(u$n_links[u$sample == "loner1"], 0L)
  expect_true(is.na(u$max_ibd[u$sample == "loner1"]))
  # exactly one of the linked pair survives
  expect_equal(sum(u$selected[u$sample %in% c("a", "b")]), 1L)
})

test_that("a looser cutoff never gives a smaller set", {
  set.seed(11)
  s <- paste0("s", 1:30)
  all_pairs <- t(utils::combn(s, 2))
  df <- data.frame(sample1 = all_pairs[, 1], sample2 = all_pairs[, 2],
                   ibd_fraction_accessible = stats::runif(nrow(all_pairs), 0, 0.1),
                   stringsAsFactors = FALSE)
  sizes <- vapply(c(0.01, 0.02, 0.05, 0.08),
                  function(t) sum(ibd_unrelated_set(df, max_ibd = t, restarts = 40)$selected),
                  integer(1))
  expect_false(is.unsorted(sizes))
})

test_that("the result reproduces and leaves the caller's random state alone", {
  set.seed(5)
  s <- paste0("s", 1:25)
  all_pairs <- t(utils::combn(s, 2))
  df <- data.frame(sample1 = all_pairs[, 1], sample2 = all_pairs[, 2],
                   ibd_fraction_accessible = stats::runif(nrow(all_pairs), 0, 0.05),
                   stringsAsFactors = FALSE)

  a <- ibd_unrelated_set(df, max_ibd = 0.02, restarts = 25, seed = 7)
  b <- ibd_unrelated_set(df, max_ibd = 0.02, restarts = 25, seed = 7)
  expect_equal(a, b)

  set.seed(99)
  before <- stats::runif(1)
  set.seed(99)
  invisible(ibd_unrelated_set(df, max_ibd = 0.02, restarts = 25, seed = 7))
  expect_equal(stats::runif(1), before)
})

test_that("metadata is carried by sample and bad columns error", {
  s <- paste0("s", 1:5)
  df <- .pairs_from_edges(s, from = c("s1", "s3"), to = c("s2", "s4"))
  meta <- data.frame(sample = s,
                     region = factor(c("north", "north", "south", "south", "north"),
                                     levels = c("south", "north")),
                     stringsAsFactors = FALSE)
  u <- ibd_unrelated_set(df, max_ibd = 0.01, meta = meta, add_meta_cols = "region")
  expect_true("region" %in% names(u))
  expect_equal(levels(u$region), c("south", "north"))
  expect_equal(as.character(u$region), as.character(meta$region[match(u$sample, meta$sample)]))
  expect_error(ibd_unrelated_set(df, max_ibd = 0.01, meta = meta, add_meta_cols = "nope"),
               "no column 'nope'")
})

test_that("an IbdResults carrying a pair table works, and meta comes from it", {
  s <- paste0("s", 1:5)
  df <- .pairs_from_edges(s, from = c("s1", "s3"), to = c("s2", "s4"))
  meta <- data.frame(sample = s, region = c(rep("north", 3), rep("south", 2)),
                     stringsAsFactors = FALSE)
  ibd <- ibd_results(pair_fraction = df, meta = meta, group_col_in_meta = "region")

  expect_equal(ibd$ibd_unrelated_set(max_ibd = 0.01),
               ibd_unrelated_set(df, max_ibd = 0.01))
  u <- ibd$ibd_unrelated_set(max_ibd = 0.01, add_meta_cols = "region")
  expect_equal(as.character(u$region), meta$region[match(u$sample, meta$sample)])
})
