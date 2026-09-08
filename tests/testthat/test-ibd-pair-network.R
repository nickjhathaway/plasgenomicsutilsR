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

test_that("border outlines the nodes, and refuses to fight the shape encoding", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("igraph")
  skip_if_not_installed("ggraph")
  pairs <- make_pairs(); meta <- make_meta()
  grp <- names(meta)[!names(meta) %in% "sample"][1]

  plain <- plot_ibd_pair_network(pairs, meta = meta, min_ibd = 0.01, color_group = grp)
  edged <- plot_ibd_pair_network(pairs, meta = meta, min_ibd = 0.01, color_group = grp, border = "black")
  node_layer <- function(p) {
    b <- ggplot2::ggplot_build(p)
    b$data[[which(vapply(b$data, function(d) "shape" %in% names(d), logical(1)))[1]]]
  }
  # off by default: the group is the mark's colour and nothing outlines it
  expect_true(is.na(formals(plot_ibd_pair_network)$border))
  expect_false("fill" %in% names(node_layer(plain)) &&
                 length(unique(node_layer(plain)$fill)) > 1)

  # on: shape 21, one outline colour, and the group moved to the fill so it can have one
  nl <- node_layer(edged)
  expect_equal(unique(nl$shape), 21)
  expect_equal(unique(nl$colour), "black")
  expect_gt(length(unique(nl$fill)), 1)
  # the same categories, still one legend, just keyed on fill now
  expect_equal(sort(unique(nl$fill)), sort(unique(node_layer(plain)$colour)))

  # only shapes 21-25 have an outline separate from a fill, so honouring both would mean
  # throwing the caller's shapes away -- refuse instead of choosing for them
  expect_error(plot_ibd_pair_network(pairs, meta = meta, min_ibd = 0.01, color_group = grp,
                                       shape_group = grp, border = "black"),
               "cannot be combined with `shape_group`")
  # and with the border off, shapes behave exactly as before
  expect_s3_class(plot_ibd_pair_network(pairs, meta = meta, min_ibd = 0.01, shape_group = grp), "ggplot")
})


# --- ibd_pair_clusters() / ibd_pair_links() -----------------------------------------------
# make_pairs(): s1-s2 0.80, s2-s3 0.40, s4-s5 0.05, everything else 0.001. At min_ibd = 0.03
# that is one chain {s1,s2,s3}, one pair {s4,s5}, and s6 alone.

test_that("the clusters are the components the network draws", {
  cl <- ibd_pair_clusters(make_pairs(), min_ibd = 0.03)
  expect_equal(names(cl), c("sample", "connected", "cluster_id", "cluster_size",
                            "n_links", "max_ibd"))
  expect_equal(nrow(cl), 6)

  expect_setequal(cl$sample[cl$connected], c("s1", "s2", "s3", "s4", "s5"))
  expect_equal(cl$sample[!cl$connected], "s6")

  # ids run largest cluster first, like gene_cluster_id
  expect_equal(sort(cl$sample[cl$cluster_id %in% 1]), c("s1", "s2", "s3"))
  expect_equal(sort(cl$sample[cl$cluster_id %in% 2]), c("s4", "s5"))
  expect_equal(cl$cluster_size[match(c("s1", "s4"), cl$sample)], c(3L, 2L))
  # a sample connected to nobody gets NA, not a cluster of its own
  expect_true(is.na(cl$cluster_id[cl$sample == "s6"]))
  expect_true(is.na(cl$cluster_size[cl$sample == "s6"]))
})

test_that("single linkage joins a chain whose ends never share directly", {
  # s1 and s3 share 0.001, well under the cutoff, yet both sit in cluster 1 through s2
  cl <- ibd_pair_clusters(make_pairs(), min_ibd = 0.03)
  expect_equal(cl$cluster_id[cl$sample == "s1"], cl$cluster_id[cl$sample == "s3"])
  expect_equal(cl$n_links[match(c("s1", "s2", "s3", "s6"), cl$sample)], c(1L, 2L, 1L, 0L))
  expect_equal(cl$max_ibd[match(c("s1", "s2", "s3"), cl$sample)], c(0.80, 0.80, 0.40))
  expect_true(is.na(cl$max_ibd[cl$sample == "s6"]))
})

test_that("connected samples come first, unconnected last", {
  cl <- ibd_pair_clusters(make_pairs(), min_ibd = 0.03)
  expect_equal(cl$sample, c("s1", "s2", "s3", "s4", "s5", "s6"))
  expect_false(any(diff(cl$connected) > 0))   # never unconnected then connected again
})

test_that("min_ibd moves the boundary the same way it does in the plot", {
  # 0.06 drops the s4-s5 edge (they share 0.05); 0.5 leaves only s1-s2
  expect_setequal(with(ibd_pair_clusters(make_pairs(), min_ibd = 0.06),
                       sample[!connected]), c("s4", "s5", "s6"))
  expect_setequal(with(ibd_pair_clusters(make_pairs(), min_ibd = 0.5),
                       sample[connected]), c("s1", "s2"))
  # nobody linked at all: every sample still comes back, all unconnected
  none <- ibd_pair_clusters(make_pairs(), min_ibd = 0.99)
  expect_equal(nrow(none), 6)
  expect_false(any(none$connected))
  expect_true(all(is.na(none$cluster_id)))
})

test_that("the counts match what the plot puts in its subtitle", {
  skip_if_no_graph()
  for (mi in c(0.03, 0.06, 0.5)) {
    cl <- ibd_pair_clusters(make_pairs(), min_ibd = mi)
    p <- plot_ibd_pair_network(make_pairs(), min_ibd = mi)
    expect_match(p$labels$subtitle,
                 sprintf("^%d samples, %d pairs sharing", nrow(cl),
                         nrow(ibd_pair_links(make_pairs(), min_ibd = mi))), info = mi)
    n_iso <- sum(!cl$connected)
    expect_match(p$labels$subtitle, sprintf("\\(%d unconnected\\)", n_iso), info = mi)
  }
})

test_that("ibd_pair_links is the edge list, highest sharing first", {
  e <- ibd_pair_links(make_pairs(), min_ibd = 0.03)
  expect_equal(names(e), c("sample1", "sample2", "ibd_fraction"))
  expect_equal(e$ibd_fraction, c(0.80, 0.40, 0.05))
  expect_equal(paste(e$sample1, e$sample2), c("s1 s2", "s2 s3", "s4 s5"))
  # the samples those edges touch are exactly the connected ones
  cl <- ibd_pair_clusters(make_pairs(), min_ibd = 0.03)
  expect_setequal(unique(c(e$sample1, e$sample2)), cl$sample[cl$connected])
})

test_that("`weight` and `samples` are honoured, and a bad column is named", {
  df <- make_pairs()
  df$ibd_fraction_full_genome <- df$ibd_fraction_accessible / 2
  expect_equal(ibd_pair_links(df, weight = "ibd_fraction_full_genome",
                              min_ibd = 0.03)$ibd_fraction, c(0.40, 0.20))
  sub <- ibd_pair_clusters(df, min_ibd = 0.03, samples = c("s1", "s2", "s6"))
  expect_equal(sub$sample, c("s1", "s2", "s6"))
  expect_equal(sub$sample[!sub$connected], "s6")
  expect_error(ibd_pair_clusters(df, weight = "nope"), "no column 'nope'")
})

test_that("an IbdResults carrying a pair table can be asked directly", {
  ibd <- ibd_results(pair_fraction = make_pairs(), meta = make_meta())
  expect_equal(ibd$ibd_pair_clusters(min_ibd = 0.03),
               ibd_pair_clusters(make_pairs(), min_ibd = 0.03))
  expect_equal(nrow(ibd$ibd_pair_links(min_ibd = 0.03)), 3)
  expect_error(ibd_pair_clusters(ibd_results(meta = make_meta())), "no pair table")
})

test_that("add_meta_cols puts each column on both ends, in the order asked for", {
  meta <- make_meta()
  meta$country <- c("UG", "UG", "TZ", "TZ", "KE", "KE")
  e <- ibd_pair_links(make_pairs(), min_ibd = 0.03,
                      add_meta_cols = c("region", "country"), meta = meta)
  expect_equal(names(e), c("sample1", "sample2", "ibd_fraction",
                           "sample1_region", "sample2_region",
                           "sample1_country", "sample2_country"))
  # s1-s2 are both region A; s2-s3 spans A and B
  expect_equal(as.character(e$sample1_region), c("A", "A", "B"))
  expect_equal(as.character(e$sample2_region), c("A", "B", "C"))
  expect_equal(e$sample1_country, c("UG", "UG", "TZ"))
  # a single column is fine, and the edges themselves are untouched
  one <- ibd_pair_links(make_pairs(), min_ibd = 0.03, add_meta_cols = "region", meta = meta)
  expect_equal(one[, 1:3], ibd_pair_links(make_pairs(), min_ibd = 0.03))
})

test_that("a factor metadata column keeps its level order on both ends", {
  # make_meta() orders region C < B < A, which a merge would turn alphabetical
  meta <- make_meta()
  e <- ibd_pair_links(make_pairs(), min_ibd = 0.03, add_meta_cols = "region", meta = meta)
  expect_s3_class(e$sample1_region, "factor")
  expect_equal(levels(e$sample1_region), c("C", "B", "A"))
  expect_equal(levels(e$sample2_region), c("C", "B", "A"))
})

test_that("meta comes from the IbdResults when there is one", {
  ibd <- ibd_results(pair_fraction = make_pairs(), meta = make_meta())
  e <- ibd$ibd_pair_links(min_ibd = 0.03, add_meta_cols = "region")
  expect_equal(as.character(e$sample1_region), c("A", "A", "B"))
  # an explicit meta still wins over the object's
  other <- make_meta(); other$region <- "Z"
  expect_true(all(ibd$ibd_pair_links(min_ibd = 0.03, add_meta_cols = "region",
                                     meta = other)$sample1_region == "Z"))
})

test_that("a sample missing from meta gets NA rather than losing its edge", {
  meta <- make_meta()[1:2, ]                       # s3..s6 unknown
  e <- ibd_pair_links(make_pairs(), min_ibd = 0.03, add_meta_cols = "region", meta = meta)
  expect_equal(nrow(e), 3)
  expect_equal(is.na(e$sample2_region), c(FALSE, TRUE, TRUE))
})

test_that("add_meta_cols says what is wrong when it cannot be honoured", {
  expect_error(ibd_pair_links(make_pairs(), add_meta_cols = "region"), "needs meta")
  expect_error(ibd_pair_links(make_pairs(), add_meta_cols = "nope", meta = make_meta()),
               "no column 'nope'")
  # the error names what it could have used instead
  expect_error(ibd_pair_links(make_pairs(), add_meta_cols = "nope", meta = make_meta()),
               "region, marker")
})

test_that("add_meta_cols reads a capitalised sample column like everything else does", {
  meta <- make_meta()
  names(meta)[names(meta) == "sample"] <- "Sample"
  e <- suppressMessages(
    ibd_pair_links(make_pairs(), min_ibd = 0.03, add_meta_cols = "region", meta = meta))
  expect_equal(as.character(e$sample1_region), c("A", "A", "B"))
})

test_that("ibd_pair_clusters carries metadata onto each sample", {
  meta <- make_meta()
  meta$country <- c("UG", "UG", "TZ", "TZ", "KE", "KE")
  cl <- ibd_pair_clusters(make_pairs(), min_ibd = 0.03,
                          add_meta_cols = c("region", "country"), meta = meta)
  expect_equal(names(cl), c("sample", "connected", "cluster_id", "cluster_size",
                            "n_links", "max_ibd", "region", "country"))
  expect_equal(as.character(cl$region), c("A", "A", "B", "B", "C", "C"))
  expect_equal(cl$country, c("UG", "UG", "TZ", "TZ", "KE", "KE"))
  # a factor keeps its level order here too, and the rest of the table is untouched
  expect_equal(levels(cl$region), c("C", "B", "A"))
  expect_equal(cl[, 1:6], ibd_pair_clusters(make_pairs(), min_ibd = 0.03))
  # and it reaches an unconnected sample as readily as a clustered one
  expect_equal(as.character(cl$region[cl$sample == "s6"]), "C")
})

test_that("ibd_pair_clusters takes meta from the IbdResults, and refuses a name clash", {
  ibd <- ibd_results(pair_fraction = make_pairs(), meta = make_meta())
  expect_equal(as.character(ibd$ibd_pair_clusters(min_ibd = 0.03,
                                                 add_meta_cols = "region")$region),
               c("A", "A", "B", "B", "C", "C"))
  # a metadata column named like one the table already reports would be silently overwritten
  clash <- make_meta(); clash$connected <- "yes"
  expect_error(ibd_pair_clusters(make_pairs(), add_meta_cols = "connected", meta = clash),
               "would overwrite 'connected'")
  expect_error(ibd_pair_clusters(make_pairs(), add_meta_cols = "nope", meta = make_meta()),
               "no column 'nope'")
})
