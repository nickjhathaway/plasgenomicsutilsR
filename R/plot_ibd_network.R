# Sample-level IBD network at a gene or locus: nodes are samples, an edge joins two
# samples that share an IBD block over the interval. Built from the IBD blocks loaded
# into an IbdResults (ibd_results(blocks=, meta=)).

# Distinguishable point shapes, solid ones first so a small node still reads clearly.
.SHAPE_PALETTE <- c(16, 17, 15, 18, 8, 3, 4, 7, 9, 10, 12, 13, 14)

.shape_palette <- function(n) {
  if (n > length(.SHAPE_PALETTE)) {
    stop("shape_group has ", n, " levels but only ", length(.SHAPE_PALETTE),
         " shapes are distinguishable; pass `shapes` explicitly or colour by this ",
         "column instead", call. = FALSE)
  }
  .SHAPE_PALETTE[seq_len(n)]
}

# Values for a manual scale. A *named* vector maps level -> value in any order and may
# cover only some levels; an unnamed one is taken positionally against `levels`. Anything
# missing falls back to the automatic palette, so a partial mapping is allowed.
.match_scale_values <- function(values, levels, what, palette) {
  auto <- stats::setNames(palette(length(levels)), levels)
  if (is.null(values)) return(auto)
  if (!is.null(names(values)) && any(nzchar(names(values)))) {
    unknown <- setdiff(names(values), levels)
    if (length(unknown)) {
      warning("`", what, "` names not among the levels are ignored: ",
              paste(unknown, collapse = ", "), call. = FALSE)
    }
    keep <- intersect(names(values), levels)
    auto[keep] <- values[keep]
    return(auto)
  }
  if (length(values) < length(levels)) {
    stop("`", what, "` has ", length(values), " value(s) for ", length(levels),
         " levels; give one per level or use a named vector", call. = FALSE)
  }
  stats::setNames(values[seq_along(levels)], levels)
}

# Resolve a gene name / locus string / interval to list(chr, start, end, label).
.resolve_locus <- function(x, gene, locus) {
  if (!is.null(gene)) {
    g <- x$get_genes()
    if (is.null(g)) stop("gene= given but this IbdResults has no gene track", call. = FALSE)
    .check_gene_request(g, gene[1])
    row <- g[tolower(g$name) == tolower(gene[1]), , drop = FALSE][1, ]
    return(list(chr = normalise_chr(row$chr), start = as.numeric(row$start),
                end = as.numeric(row$end), label = as.character(row$name)))
  }
  if (is.null(locus)) {
    stop("give a gene= (name in the track) or locus= (\"chr:pos\", \"chr:start-end\", ",
         "or a data frame with chr/start/end)", call. = FALSE)
  }
  # a bare position is the single base [pos, pos + 1), keeping every interval half-open
  if (is.data.frame(locus)) {
    st <- as.numeric(locus$start[1])
    en <- if ("end" %in% names(locus)) as.numeric(locus$end[1]) else st + 1
    return(list(chr = normalise_chr(locus$chr[1]), start = st, end = en,
                label = sprintf("%s:%s-%s", normalise_chr(locus$chr[1]), st, en)))
  }
  s <- as.character(locus)[1]
  parts <- strsplit(s, ":", fixed = TRUE)[[1]]
  chr <- normalise_chr(paste(parts[-length(parts)], collapse = ":"))
  rng <- parts[length(parts)]
  if (grepl("-", rng, fixed = TRUE)) {
    se <- as.numeric(strsplit(rng, "-", fixed = TRUE)[[1]]); st <- se[1]; en <- se[2]
  } else {
    st <- as.numeric(rng); en <- st + 1
  }
  list(chr = chr, start = st, end = en, label = s)
}

# Force-directed layouts pull every adjacent pair together with the same strength, so a
# densely inter-connected group collapses to a blob in which individual nodes cannot be
# resolved. Weighting each edge by `(1 - J)^spread` -- J being the Jaccard overlap of the
# two endpoints' neighbourhoods -- leaves edges that bridge otherwise separate
# neighbourhoods at full strength (J = 0 gives weight 1 at any `spread`) while edges buried
# inside a near-clique (whose endpoints share nearly all their neighbours) pull only
# weakly, so repulsion opens such a group out into a readable disc. The exponent form is
# unbounded above, so `spread` past 1 keeps loosening the densest groups. Jaccard is
# computed per edge from adjacency lists rather than from a full node x node matrix.
.clique_relax_weights <- function(g, spread) {
  if (is.null(spread) || !is.finite(spread) || spread <= 0) return(NULL)
  el <- igraph::as_edgelist(g, names = FALSE)
  if (!nrow(el)) return(NULL)
  nb <- lapply(igraph::adjacent_vertices(g, igraph::V(g)), as.integer)
  deg <- igraph::degree(g)
  shared <- mapply(function(u, v) length(intersect(nb[[u]], nb[[v]])), el[, 1], el[, 2])
  jac <- shared / pmax(deg[el[, 1]] + deg[el[, 2]] - shared, 1)
  pmax((1 - jac)^spread, 1e-4)
}

#' IBD network at a gene or locus
#'
#' A sample-level network for one gene or locus: each node is a sample and an edge joins
#' two samples whose pair shares an **IBD block overlapping** the interval. Optionally
#' colours and/or shapes nodes by metadata columns, and optionally drops isolated nodes
#' for a cleaner picture -- or keeps them so the total sample count is visible.
#'
#' Needs an [IbdResults] built with `blocks =` (and `meta =` for grouping / the full
#' sample set). See [gene_ibd_overlap()] for the pairwise-fraction summary of the same data.
#'
#' @param x An [IbdResults] with IBD `blocks`.
#' @param gene A single gene name from the object's track (its interval is used).
#' @param locus Alternatively, a locus: `"chr:pos"`, `"chr:start-end"`, or a one-row data
#'   frame with `chr`, `start`, `end`. Give exactly one of `gene` / `locus`. Coordinates
#'   are 0-based half-open (see [plasgenomicsutilsR-coordinates]), so `"chr:1000-2000"` is
#'   the 1000 bases from 0-based 1000 up to but not including 2000, and a bare `"chr:1000"`
#'   is the single base at 0-based 1000.
#' @param color_group Optional metadata column to colour nodes by (needs `meta`).
#' @param colors Colours for `color_group`. A **named** `level -> colour` vector maps by
#'   name and may cover only some levels (the rest keep their automatic colour); an
#'   unnamed vector is taken positionally, in the column's level order. `NULL` (default)
#'   uses [meta_colors()], so the same mapping can be shared with other plots.
#' @param shape_group Optional metadata column to set node **shape** by, so a second
#'   variable can be read off the same plot (needs `meta`). Independent of `color_group`:
#'   use either, both, or neither.
#' @param shapes Shapes for `shape_group`, named or positional exactly like `colors`
#'   (values are \pkg{ggplot2} shape codes). `NULL` picks distinguishable defaults and
#'   errors if the column has more levels than there are distinct shapes.
#' @param within Pad the interval by this many bp on both sides (default `0`).
#' @param sharing What an edge requires of a pair's IBD segment:
#'   \describe{
#'     \item{`"overlap"`}{(default) the segment touches the interval anywhere -- the pair
#'       shares *some* of the gene/locus.}
#'     \item{`"complete"`}{the segment spans the whole interval -- the pair shares the
#'       *entire* gene/locus. A stricter, usually much sparser graph.}
#'   }
#'   `within` applies either way, so with padding `"complete"` asks the segment to cover the
#'   padded interval. Use [gene_ibd_pairs()] to see per-pair which of the two each segment
#'   satisfies (its `coverage` column) and how much of the gene is covered.
#' @param include_isolated Keep samples with no IBD edge (default `FALSE`). `TRUE` also
#'   shows every analyzed sample; the isolated samples are drawn in a grid below the
#'   connected component (sorted by `group` when colouring by one). `FALSE` shows only the
#'   connected nodes.
#' @param layout A \pkg{ggraph}/\pkg{igraph} layout name (default `"fr"`, Fruchterman-
#'   Reingold). Other options include `"kk"` (Kamada-Kawai), `"stress"`, `"drl"`,
#'   `"circle"`, `"nicely"`, `"graphopt"`, and `"lgl"`. The connected component is laid out
#'   on its own, with its aspect preserved.
#' @param spread How strongly to open up densely inter-connected groups (default `1.5`).
#'   Edge attraction is weighted by `(1 - J)^spread`, where `J` is the Jaccard overlap of
#'   the two samples' IBD neighbourhoods: edges inside a near-clique pull only weakly, so
#'   the group expands enough to resolve individual samples, while edges linking otherwise
#'   separate groups keep their full pull and cluster separation is preserved. `0` leaves
#'   edges unweighted; raise it above the default to loosen the densest groups further, at
#'   the cost of them taking up more of the canvas. Applies to weight-aware layouts
#'   (`"fr"`, `"kk"`, `"drl"`, `"stress"`, ...); others are unaffected.
#' @param node_size,node_alpha Node point aesthetics.
#' @param edge_colour,edge_width Edge aesthetics (default width `1`).
#' @param edge_alpha Edge opacity (default `0.5`). `NULL` instead scales opacity down with
#'   edge count (`~120 / n_edges`, clamped to \[0.06, 0.6]).
#' @param title Plot title: `NULL` (default) uses `"IBD network: <label>"`, a string sets a
#'   custom title, and `NA`/`FALSE` draws no title.
#' @param subtitle Show the sample / IBD-pair count line under the title (default `TRUE`).
#' @param seed Random seed for stochastic layouts (e.g. `"fr"`).
#' @return A ggplot (\pkg{ggraph}) object.
#' @examples
#' \dontrun{
#' ibd <- ibd_results(blocks = "hmm.txt", meta = meta, genes = PF_EXAMPLE_DRUG_GENES)
#' plot_ibd_network(ibd, gene = "pfcrt", color_group = "region")
#' # colour by region, shape by year: two variables on one plot
#' plot_ibd_network(ibd, gene = "pfcrt", color_group = "region",
#'                  shape_group = "collection_year")
#' # a named vector pins chosen levels; the rest keep their automatic colour
#' plot_ibd_network(ibd, gene = "pfcrt", color_group = "region",
#'                  colors = c(North = "#1f78b4", East = "#33a02c"))
#' plot_ibd_network(ibd, locus = "Pf3D7_07_v3:403500", include_isolated = TRUE)
#' }
#' @export
plot_ibd_network <- function(x, gene = NULL, locus = NULL,
                             color_group = NULL, colors = NULL,
                             shape_group = NULL, shapes = NULL, within = 0,
                             sharing = c("overlap", "complete"),
                             include_isolated = FALSE, layout = "fr", spread = 1.5,
                             node_size = 3, node_alpha = 0.9,
                             edge_colour = "grey65", edge_alpha = 0.5, edge_width = 1,
                             title = NULL, subtitle = TRUE, seed = 42) {
  .need_package("ggplot2", "plot_ibd_network()")
  .need_package("igraph", "plot_ibd_network()")
  .need_package("ggraph", "plot_ibd_network()")
  sharing <- match.arg(sharing)
  blocks <- x$get_blocks()
  if (is.null(blocks)) {
    stop("this IbdResults has no IBD blocks; build it with ibd_results(blocks = , meta = )",
         call. = FALSE)
  }
  iv <- .resolve_locus(x, gene, locus)

  bl <- blocks
  bl$chr <- normalise_chr(bl$chr)
  lo <- iv$start - within
  hi <- iv$end + within
  m <- bl$chr == iv$chr & if (sharing == "complete") {
    bl$start <= lo & bl$end >= hi          # the segment has to span the whole interval
  } else {
    bl$start < hi & bl$end > lo            # any half-open [start, end) overlap
  }
  sub <- bl[m, , drop = FALSE]
  # distinct undirected pairs sharing IBD over the interval
  a <- pmin(sub$sample1, sub$sample2); b <- pmax(sub$sample1, sub$sample2)
  edges <- unique(data.frame(from = a, to = b, stringsAsFactors = FALSE))

  meta <- x$get_meta()
  analyzed <- x$get_analyzed_samples()
  if (is.null(analyzed)) analyzed <- unique(c(bl$sample1, bl$sample2))
  connected <- unique(c(edges$from, edges$to))
  isolated <- if (include_isolated) setdiff(analyzed, connected) else character(0)
  if (!length(connected) && !length(isolated)) {
    stop("no samples share an IBD block over ", iv$label,
         " (set include_isolated = TRUE to still show the samples)", call. = FALSE)
  }
  # sample -> value lookups for each aesthetic, keeping the column's own level order
  lookup <- function(col, what) {
    if (is.null(col)) return(NULL)
    if (is.null(meta)) stop(what, " needs meta; build with ibd_results(meta = )", call. = FALSE)
    if (!col %in% names(meta)) stop("meta has no column '", col, "'", call. = FALSE)
    v <- .as_group_factor(meta[[col]])
    list(map = stats::setNames(as.character(v), as.character(meta$sample)),
         levels = .levels_of(v))
  }
  col_of <- lookup(color_group, "color_group")
  shp_of <- lookup(shape_group, "shape_group")

  # lay out the connected component on its own, normalised with a single scale (aspect kept)
  set.seed(seed)
  main <- NULL
  if (length(connected)) {
    g <- igraph::graph_from_data_frame(edges, directed = FALSE,
                                       vertices = data.frame(name = connected))
    w <- .clique_relax_weights(g, spread)
    lay <- if (is.null(w)) ggraph::create_layout(g, layout = layout) else
      tryCatch({ set.seed(seed); ggraph::create_layout(g, layout = layout, weights = w) },
               error = function(e) { set.seed(seed); ggraph::create_layout(g, layout = layout) })
    scl <- max(diff(range(lay$x)), diff(range(lay$y)), 1e-9)   # single scale keeps aspect
    main <- data.frame(name = lay$name,
                       x = (lay$x - min(lay$x)) / scl,
                       y = (lay$y - min(lay$y)) / scl,
                       stringsAsFactors = FALSE)
  }

  # isolated samples: a grid below the connected component, sorted by the grouping(s) so
  # like samples sit together -- colour first, then shape
  iso <- NULL
  if (length(isolated)) {
    keys <- list()
    if (!is.null(col_of)) keys <- c(keys, list(match(col_of$map[isolated], col_of$levels)))
    if (!is.null(shp_of)) keys <- c(keys, list(match(shp_of$map[isolated], shp_of$levels)))
    if (length(keys)) isolated <- isolated[do.call(order, c(keys, list(isolated)))]
    xr <- if (!is.null(main)) range(main$x) else c(0, 1)
    span <- diff(xr)
    ncw  <- max(1L, ceiling(sqrt(length(isolated)) * 2.4))
    nrw  <- ceiling(length(isolated) / ncw)
    step <- span / max(1, ncw - 1)                                   # square grid spacing
    y_top <- (if (!is.null(main)) min(main$y) else 0) - 3 * step     # gap below main
    gx <- xr[1] + ((seq_along(isolated) - 1) %% ncw) * step
    gy <- y_top - ((seq_along(isolated) - 1) %/% ncw) * step
    iso <- data.frame(name = isolated, x = gx, y = gy, stringsAsFactors = FALSE)
  }

  nodes <- rbind(main, iso)
  if (!is.null(col_of)) nodes$.colour <- factor(col_of$map[nodes$name], levels = col_of$levels)
  if (!is.null(shp_of)) nodes$.shape <- factor(shp_of$map[nodes$name], levels = shp_of$levels)

  # edges in the connected-component coordinates (isolated nodes have none)
  edf <- NULL
  if (!is.null(main) && nrow(edges)) {
    edf <- data.frame(
      x  = main$x[match(edges$from, main$name)], y  = main$y[match(edges$from, main$name)],
      xe = main$x[match(edges$to,   main$name)], ye = main$y[match(edges$to,   main$name)])
    edf <- edf[stats::complete.cases(edf), , drop = FALSE]
  }
  # scale edge opacity down with edge count when not set explicitly
  if (is.null(edge_alpha)) edge_alpha <- max(0.06, min(0.6, 120 / max(1, nrow(edges))))

  p <- ggplot2::ggplot()
  if (!is.null(edf) && nrow(edf)) {
    p <- p + ggplot2::geom_segment(
      data = edf, ggplot2::aes(.data$x, .data$y, xend = .data$xe, yend = .data$ye),
      colour = edge_colour, alpha = edge_alpha, linewidth = edge_width)
  }
  # nodes: colour and shape are independent, so either, both or neither can be mapped
  aes_args <- list(x = quote(.data$x), y = quote(.data$y))
  if (!is.null(col_of)) aes_args$colour <- quote(.data$.colour)
  if (!is.null(shp_of)) aes_args$shape <- quote(.data$.shape)
  fixed <- list(data = nodes, mapping = do.call(ggplot2::aes, aes_args),
                size = node_size, alpha = node_alpha)
  if (is.null(col_of)) fixed$colour <- "#2166ac"
  p <- p + do.call(ggplot2::geom_point, fixed)

  if (!is.null(col_of)) {
    node_cols <- .match_scale_values(colors, col_of$levels, "colors",
                                     function(n) meta_colors(
                                       data.frame(g = factor(col_of$levels,
                                                             levels = col_of$levels)))[["g"]])
    p <- p + ggplot2::scale_colour_manual(values = node_cols, name = color_group,
                                          na.value = "grey70", drop = FALSE)
  }
  if (!is.null(shp_of)) {
    node_shapes <- .match_scale_values(shapes, shp_of$levels, "shapes",
                                       function(n) .shape_palette(n))
    p <- p + ggplot2::scale_shape_manual(values = node_shapes, name = shape_group,
                                         na.value = 4, drop = FALSE)
  }
  # dashed separator + label between the connected component (y >= 0) and the isolated grid
  if (!is.null(iso)) {
    sep_y <- (min(main$y) + max(iso$y)) / 2
    p <- p +
      ggplot2::annotate("segment", x = min(nodes$x), xend = max(nodes$x),
                        y = sep_y, yend = sep_y, colour = "grey80", linewidth = 0.3,
                        linetype = "dashed") +
      ggplot2::annotate("text", x = min(nodes$x), y = max(iso$y) + (min(main$y) - max(iso$y)) * 0.28,
                        hjust = 0, vjust = 0.5, label = sprintf("%d unconnected", length(isolated)),
                        size = 3, colour = "grey45")
  }

  ttl <- if (is.null(title)) paste0("IBD network: ", iv$label)
         else if (isFALSE(title) || (length(title) == 1 && is.na(title))) NULL
         else as.character(title)
  sub_lab <- if (isTRUE(subtitle)) {
    # name the criterion: the same gene gives very different graphs under the two
    sprintf("%d samples, %d IBD pairs sharing %s%s", nrow(nodes), nrow(edges),
            if (sharing == "complete") "the whole interval" else "part of the interval",
            if (length(isolated)) sprintf(" (%d unconnected)", length(isolated)) else "")
  } else NULL

  p <- p +
    ggplot2::labs(title = ttl, subtitle = sub_lab) +
    ggplot2::coord_fixed() +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   legend.position = "right")
  # suggested output size, scaled to the sample count (used by save_plot())
  sz <- max(4, min(8, sqrt(nrow(nodes)) * 0.35))
  attr(p, "plasgenomics_dims") <- c(width = round(sz, 1),
    height = round(sz * (if (length(isolated)) 1.25 else 1), 1))
  p
}
