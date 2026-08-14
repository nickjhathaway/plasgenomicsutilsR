# Group-level summaries of the genome-wide per-pair IBD fraction. The numbers are checked
# against the pairs by hand, since the whole point of the function is the denominator: which
# sample pairs a group-pair row is averaging over.

.tiny_pairs <- function() {
  data.frame(sample1 = c("s1", "s1", "s1", "s2", "s2", "s3"),
             sample2 = c("s2", "s3", "s4", "s3", "s4", "s4"),
             ibd_fraction_accessible = c(0.02, 0.31, 0.04, 0.28, 0.03, 0.05),
             stringsAsFactors = FALSE)
}
.tiny_meta <- function() {
  data.frame(sample = paste0("s", 1:4),
             region = c("north", "north", "south", "south"), stringsAsFactors = FALSE)
}

test_that("each row summarises exactly the pairs spanning its two groups", {
  res <- pair_fraction_summary(.tiny_pairs(), group = "region", meta = .tiny_meta())
  expect_equal(nrow(res), 3L)                       # north/north, north/south, south/south

  between <- res[res$group_a == "north" & res$group_b == "south", ]
  cross <- c(0.31, 0.04, 0.28, 0.03)                # s1-s3, s1-s4, s2-s3, s2-s4
  expect_equal(between$n_pairs, 4L)
  expect_equal(between$mean, mean(cross))
  expect_equal(between$median, stats::median(cross))
  expect_equal(between$sd, stats::sd(cross))
  expect_equal(between$min, min(cross))
  expect_equal(between$max, max(cross))
  expect_equal(c(between$q25, between$q75),
               unname(stats::quantile(cross, c(0.25, 0.75))))

  # a within-group row is that group's own pairs, and nothing else
  expect_equal(res$n_pairs[res$group_a == "north" & res$group_b == "north"], 1L)
  expect_equal(res$mean[res$group_a == "north" & res$group_b == "north"], 0.02)
  expect_true(is.na(res$sd[res$group_a == "north" & res$group_b == "north"]))  # n = 1
})

test_that("a pair counts once whichever way round its samples are stored", {
  flipped <- .tiny_pairs()
  flipped[] <- lapply(flipped, identity)
  s1 <- flipped$sample1; flipped$sample1 <- flipped$sample2; flipped$sample2 <- s1
  a <- pair_fraction_summary(.tiny_pairs(), group = "region", meta = .tiny_meta())
  b <- pair_fraction_summary(flipped, group = "region", meta = .tiny_meta())
  expect_equal(a, b)
})

test_that("the summary recovers sharing that was built into the data", {
  set.seed(7)
  s <- paste0("s", 1:12)
  g <- t(utils::combn(s, 2))
  meta <- data.frame(sample = s,
                     region = factor(rep(c("west", "east", "north"), each = 4),
                                     levels = c("north", "west", "east")),
                     stringsAsFactors = FALSE)
  gm <- stats::setNames(as.character(meta$region), meta$sample)
  # within-group pairs share an order of magnitude more than between-group ones
  base <- ifelse(gm[g[, 1]] == gm[g[, 2]], 0.30, 0.03)
  pf <- data.frame(sample1 = g[, 1], sample2 = g[, 2],
                   ibd_fraction_accessible = base + stats::rnorm(nrow(g), 0, 0.005),
                   stringsAsFactors = FALSE)
  obj <- ibd_results(pair_fraction = pf, meta = meta, group_col_in_meta = "region")
  res <- obj$pair_fraction_summary()

  expect_equal(nrow(res), 6L)                       # 3 within + 3 between
  w <- res$group_a == res$group_b
  expect_true(all(res$mean[w] > 0.25))
  expect_true(all(res$mean[!w] < 0.05))
  # the object's group order drives the rows, not an alphabetical one
  expect_equal(levels(res$group_a), c("north", "west", "east"))
  expect_equal(as.character(res$group_a), c("north", "north", "north", "west", "west", "east"))

  # n_pairs is the full pairwise count, which is what the mean divides by
  expected <- ifelse(w, choose(res$n_samples_a, 2), res$n_samples_a * res$n_samples_b)
  expect_equal(res$n_pairs, as.integer(expected))
})

test_that("pairs that cannot contribute are dropped, and said so", {
  p <- .tiny_pairs()
  m <- .tiny_meta()[1:3, ]                           # s4 has no metadata
  expect_message(pair_fraction_summary(p, group = "region", meta = m),
                 "touch a sample with no `region`")
  res <- suppressMessages(pair_fraction_summary(p, group = "region", meta = m))
  expect_equal(sum(res$n_pairs), 3L)                 # s1-s2, s1-s3, s2-s3

  p2 <- .tiny_pairs(); p2$ibd_fraction_accessible[2] <- NA
  expect_message(pair_fraction_summary(p2, group = "region", meta = .tiny_meta()),
                 "no `ibd_fraction_accessible` value")
  res2 <- suppressMessages(pair_fraction_summary(p2, group = "region", meta = .tiny_meta()))
  expect_equal(res2$n_pairs[res2$group_a == "north" & res2$group_b == "south"], 3L)
})

test_that("quantile columns are named for the percent asked for", {
  res <- pair_fraction_summary(.tiny_pairs(), group = "region", meta = .tiny_meta(),
                               probs = c(0.05, 0.5, 0.95))
  expect_true(all(c("q5", "q50", "q95") %in% names(res)))
  # 50% must not be trimmed to "5" and collide with a genuine 5% column
  expect_false(anyDuplicated(names(res)) > 0)
  expect_equal(res$q50, res$median)
  expect_named(pair_fraction_summary(.tiny_pairs(), group = "region", meta = .tiny_meta(),
                                     probs = c(0.025, 0.975))[, 11:12], c("q2.5", "q97.5"))
})

test_that("the edges of the group set behave", {
  p <- .tiny_pairs()
  one <- data.frame(sample = paste0("s", 1:4), region = "only", stringsAsFactors = FALSE)
  # one group still has a within-group row ...
  expect_equal(nrow(pair_fraction_summary(p, group = "region", meta = one)), 1L)
  # ... and nothing to compare between
  expect_error(pair_fraction_summary(p, group = "region", meta = one, within = FALSE),
               "no between-group comparison")

  expect_equal(nrow(pair_fraction_summary(p, group = "region", meta = .tiny_meta(),
                                          within = FALSE)), 1L)
  expect_message(pair_fraction_summary(p, group = "region", meta = .tiny_meta(),
                                       min_pairs = 4), "fewer than 4")
  # the default drops empty combinations, and says that rather than "fewer than 1"
  expect_message(pair_fraction_summary(p, group = "region", meta = .tiny_meta()[1:3, ]),
                 "have no pairs in the table")

  expect_error(pair_fraction_summary(p, group = "region", meta = .tiny_meta(),
                                     value = "nope"), "no column 'nope'")
  expect_error(pair_fraction_summary(p, group = "nope", meta = .tiny_meta()),
               "no metadata column")
  expect_error(pair_fraction_summary(p, group = "region", meta = .tiny_meta(),
                                     probs = c(0.5, 2)), "proportions")
  expect_error(pair_fraction_summary(p[0, ], group = "region", meta = .tiny_meta()),
               "empty")
})

test_that("the class method needs a pair table and honours `value`", {
  obj <- ibd_results(meta = .tiny_meta(), group_col_in_meta = "region")
  expect_error(obj$pair_fraction_summary(), "no pair table")

  p <- .tiny_pairs()
  p$ibd_fraction_full_genome <- p$ibd_fraction_accessible * 0.5
  obj$set_pair_fraction(p)
  res <- obj$pair_fraction_summary(value = "ibd_fraction_full_genome")
  expect_equal(attr(res, "value"), "ibd_fraction_full_genome")
  expect_equal(res$mean,
               suppressMessages(obj$pair_fraction_summary())$mean * 0.5)
})

test_that("summarising by a column other than the declared group one works", {
  set.seed(1)
  s <- paste0("s", 1:30)
  g <- t(utils::combn(s, 2))
  site <- paste0("site", rep(1:6, each = 5))
  meta <- data.frame(sample = s, uganda_site = site,
                     region = ifelse(site %in% paste0("site", 1:3), "north", "south"),
                     stringsAsFactors = FALSE)
  pf <- data.frame(sample1 = g[, 1], sample2 = g[, 2],
                   ibd_fraction_accessible = stats::runif(nrow(g), 0.01, 0.5),
                   stringsAsFactors = FALSE)
  obj <- ibd_results(pair_fraction = pf, meta = meta, group_col_in_meta = "region")

  # the stored group order is the region one; asking for a different column must use that
  # column's own levels rather than intersecting site names with region labels
  res <- obj$pair_fraction_summary(group = "uganda_site")
  expect_equal(nrow(res), choose(6, 2) + 6)
  expect_setequal(levels(res$group_a), unique(site))
  expect_equal(unique(res$n_pairs[res$group_a == res$group_b]), choose(5, 2))
  expect_equal(unique(res$n_pairs[res$group_a != res$group_b]), 5L * 5L)

  expect_equal(nrow(obj$pair_fraction_summary(group = "region")), 3L)
})

test_that("a group order that names only some levels still summarises them all", {
  p <- .tiny_pairs()
  meta <- .tiny_meta()
  obj <- ibd_results(pair_fraction = p, meta = meta, group_col_in_meta = "region")
  suppressWarnings(obj$set_group_order("south"))       # north not named
  res <- obj$pair_fraction_summary()
  expect_equal(nrow(res), 3L)                          # nothing dropped
  expect_equal(levels(res$group_a), c("south", "north"))   # named first, rest appended
})
