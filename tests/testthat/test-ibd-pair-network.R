# Genome-wide IBD network from the per-pair fraction table.

skip_if_no_graph <- function() {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("igraph")
  testthat::skip_if_not_installed("ggraph")
}

# s1-s2-s3 form a chain of related samples, s4/s5 share a little, s6 shares with nobody
make_pairs <- function() {
  s <- paste0("s", 1:6)
  grid <- t(utils::combn(s, 2))
  df <- data.frame(sample1 = grid[, 1], sample2 = grid[, 2],
                   ibd_fraction_accessible = 0.001, stringsAsFactors = FALSE)
  set <- function(df, a, b, v) { df$ibd_fraction_accessible[df$sample1 == a & df$sample2 == b] <- v; df }
  df <- set(df, "s1", "s2", 0.80)
  df <- set(df, "s2", "s3", 0.40)
  df <- set(df, "s4", "s5", 0.05)
  df$pair <- paste0(df$sample1, "__", df$sample2)
  df
}

make_meta <- function() data.frame(
  sample = paste0("s", 1:6),
  region = factor(c("A", "A", "B", "B", "C", "C"), levels = c("C", "B", "A")),
  marker = c("x", "y", "x", NA, "y", NA), stringsAsFactors = FALSE)

edge_layer <- function(p) {
  b <- ggplot2::ggplot_build(p)
  i <- which(vapply(b$plot$layers, function(l) inherits(l$geom, "GeomSegment"), logical(1)))
  if (!length(i)) return(NULL)
  b$data[[i[1]]]
}

test_that("edges are the pairs above min_ibd, and the rest are drawn as unconnected", {
  skip_if_no_graph()
  p <- plot_ibd_pair_network(make_pairs(), meta = make_meta(), min_ibd = 0.01)
  e <- edge_layer(p)
  expect_equal(nrow(e), 3L)                       # s1-s2, s2-s3, s4-s5
  expect_true(grepl("6 samples, 3 pairs", p$labels$subtitle))
  expect_true(grepl("1 unconnected", p$labels$subtitle))   # s6

  # raising the cutoff drops the weak pair and isolates two more samples
  p2 <- plot_ibd_pair_network(make_pairs(), meta = make_meta(), min_ibd = 0.10)
  expect_equal(nrow(edge_layer(p2)), 2L)
  expect_true(grepl("3 unconnected", p2$labels$subtitle))
})

test_that("edge width maps to the IBD fraction", {
  skip_if_no_graph()
  p <- plot_ibd_pair_network(make_pairs(), meta = make_meta(), min_ibd = 0.01)
  e <- edge_layer(p)
  # ggplot renames the mapped aesthetic, so compare against the weights that were kept,
  # in the order the edges were built
  kept <- sort(c(0.80, 0.40, 0.05), decreasing = TRUE)
  expect_equal(length(unique(round(e$linewidth, 6))), 3L)
  expect_equal(order(e$linewidth), order(kept))          # width tracks the fraction
  expect_equal(which.max(e$linewidth), which.max(kept))

  s <- ggplot2::ggplot_build(p)$plot$scales$get_scales("linewidth")
  expect_equal(s$name, "IBD")
  # breaks come back in log2 space; in data space they are powers of two
  br <- 2^s$get_breaks()
  br <- br[is.finite(br)]
  expect_true(all(br > 0))
  expect_equal(br, 2^round(log2(br)))
})

test_that("colour and shape groups, the title and the subtitle all toggle", {
  skip_if_no_graph()
  pairs <- make_pairs(); meta <- make_meta()
  p <- plot_ibd_pair_network(pairs, meta = meta, color_group = "region",
                             shape_group = "marker", title = "custom", subtitle = FALSE)
  expect_equal(p$labels$title, "custom")
  expect_null(p$labels$subtitle)
  # the metadata factor's own level order is kept
  expect_equal(ggplot2::ggplot_build(p)$plot$scales$get_scales("colour")$get_limits(),
               c("C", "B", "A"))
  # NA in the shape column does not borrow a real level's shape
  s <- ggplot2::ggplot_build(p)$plot$scales$get_scales("shape")
  lv <- s$get_limits(); lv <- lv[!is.na(lv)]
  expect_false(s$na.value %in% unname(s$map(lv)))

  expect_null(plot_ibd_pair_network(pairs, meta = meta, title = FALSE)$labels$title)
  expect_equal(plot_ibd_pair_network(pairs, meta = meta, subtitle = "mine")$labels$subtitle,
               "mine")
})

test_that("a missing endpoint column or weight column says what to do", {
  skip_if_no_graph()
  bad <- make_pairs()[, c("pair", "ibd_fraction_accessible")]
  expect_error(plot_ibd_pair_network(bad), "no sample1/sample2 column")
  expect_error(plot_ibd_pair_network(make_pairs(), weight = "nope"),
               "no column 'nope'")
})

test_that("an IbdResults carries the pair table and narrows it with its groups", {
  skip_if_no_graph()
  ibd <- ibd_results(pair_fraction = make_pairs(), meta = make_meta(),
                     group_col_in_meta = "region", reference = "pf3d7")
  expect_equal(nrow(ibd$get_pair_fraction()), 15L)
  expect_s3_class(plot_ibd_pair_network(ibd, color_group = "region", min_ibd = 0.01), "ggplot")

  # dropping a group drops its samples from the pair table too
  sub <- ibd$subset_groups(drop = "C")
  pf <- sub$get_pair_fraction()
  expect_false(any(c("s5", "s6") %in% c(pf$sample1, pf$sample2)))
  expect_equal(nrow(pf), 6L)                       # C(4,2) among s1-s4
})

test_that("samples = restricts the network without touching the file", {
  skip_if_no_graph()
  p <- plot_ibd_pair_network(make_pairs(), meta = make_meta(),
                             samples = paste0("s", 1:3), min_ibd = 0.01)
  expect_true(grepl("^3 samples, 2 pairs", p$labels$subtitle))
})

# ggplot orders guides by their `order`, and with the default 0 it falls back to a hash of the
# guide -- so the stacking changed with the labels and two plots of the same cohort could put
# their legends in different orders. Read the built gtable, not the scales.
legend_titles <- function(p) {
  g <- ggplot2::ggplotGrob(p)
  i <- which(g$layout$name == "guide-box-right")
  if (!length(i)) return(character(0))
  gb <- g$grobs[[i]]
  titles <- vapply(gb$grobs, function(x) {
    v <- tryCatch(x$grobs[[which(grepl("title", x$layout$name))[1]]]$children[[1]]$label,
                  error = function(e) NA_character_)
    if (length(v) == 1) as.character(v) else NA_character_
  }, character(1))
  keep <- !is.na(titles)
  titles[keep][order(gb$layout$t[keep])]
}

test_that("the legends stack in the same order whatever the data", {
  skip_if_no_graph()
  mk <- function(n, regions, markers) {
    s <- paste0("s", seq_len(n)); grid <- t(utils::combn(s, 2))
    set.seed(1)
    list(df = data.frame(sample1 = grid[, 1], sample2 = grid[, 2],
                         ibd_fraction_accessible = stats::runif(nrow(grid), 0.02, 0.9)),
         meta = data.frame(sample = s, region = rep(regions, length.out = n),
                           marker = rep(markers, length.out = n)))
  }
  seen <- lapply(list(list(6, c("A", "B"), c("x", "y")),
                      list(8, c("A", "B", "C"), c("x", "y")),
                      list(10, c("N", "S"), c("wt", "mut", "na"))),
                 function(cfg) {
                   d <- mk(cfg[[1]], cfg[[2]], cfg[[3]])
                   legend_titles(plot_ibd_pair_network(
                     d$df, meta = d$meta, color_group = "region", shape_group = "marker",
                     min_ibd = 0.01))
                 })
  expect_equal(seen[[1]], c("region", "marker", "IBD"))
  expect_equal(seen[[2]], seen[[1]])          # a third region must not reshuffle them
  expect_equal(seen[[3]], seen[[1]])          # nor a third marker
})

test_that("the per-gene network stacks colour above shape, whatever the data", {
  skip_if_no_graph()
  mk <- function(n, regions, markers) {
    genes <- data.frame(name = "pfcrt", chr = "7", start = 403000, end = 406000)
    blocks <- data.frame(sample1 = paste0("s", 1:(n - 1)), sample2 = paste0("s", 2:n),
                         chr = "Pf3D7_07_v3", start = 403500, end = 405500, different = 0)
    meta <- data.frame(sample = paste0("s", seq_len(n)),
                       region = rep(regions, length.out = n),
                       marker = rep(markers, length.out = n))
    ibd_results(genes = genes, blocks = blocks, meta = meta,
                min_block_snp = 0, min_block_kb = 0, reference = "pf3d7")
  }
  for (cfg in list(list(8, c("A", "B"), c("x", "y")),
                   list(10, c("A", "B", "C"), c("x", "y")),
                   list(12, c("A", "B"), c("x", "y", "z")))) {
    p <- plot_ibd_network(mk(cfg[[1]], cfg[[2]], cfg[[3]]), gene = "pfcrt",
                          color_group = "region", shape_group = "marker")
    expect_equal(legend_titles(p), c("region", "marker"),
                 info = sprintf("%d samples", cfg[[1]]))
  }
})

test_that("legend breaks are the powers of two inside the data, not rounded past it", {
  f <- plasgenomicsutilsR:::.ibd_weight_breaks

  # the case from a real run: the callable-map denominator puts a fully-shared pair a hair
  # over 1, and rounding the top outwards used to add a censored `2` break whose octave
  # pushed the smallest real break off the bottom of the legend
  expect_equal(f(c(0.03001, 1.000000615142868)),
               c(0.03125, 0.0625, 0.125, 0.25, 0.5, 1))

  # every break lies within the data
  for (w in list(c(0.031, 0.9), c(0.002, 0.5), c(0.03001, 1.0000006))) {
    br <- f(w)
    expect_true(all(br >= min(w) & br <= max(w)),
                info = paste(range(w), collapse = ".."))
  }

  # a wide range is thinned but still spans it, keeping the largest break
  wide <- f(c(0.0005, 1))
  expect_lte(length(wide), 6L)
  expect_equal(max(wide), 1)
  expect_true(min(wide) < 0.01)

  # too narrow to label in octaves -> let ggplot choose rather than invent one key
  expect_s3_class(f(c(0.4, 0.6)), "waiver")
  expect_s3_class(f(numeric(0)), "waiver")
})

test_that("the class exposes the pair network the same way it exposes the per-gene one", {
  skip_if_no_graph()
  ibd <- ibd_results(pair_fraction = make_pairs(), meta = make_meta(),
                     group_col_in_meta = "region", reference = "pf3d7")
  expect_true("plot_ibd_pair_network" %in% names(IbdResults$public_methods))
  p <- ibd$plot_ibd_pair_network(min_ibd = 0.01, color_group = "region")
  expect_s3_class(p, "ggplot")
  # ... and it errors helpfully when no pair table was attached
  bare <- ibd_results(meta = make_meta(), reference = "pf3d7")
  expect_error(bare$plot_ibd_pair_network(), "no pair table")
})
