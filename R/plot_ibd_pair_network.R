# Genome-wide IBD network: nodes are samples, an edge is the fraction of the genome a pair
# shares IBD. The per-gene network (plot_ibd_network()) answers "who shares this locus";
# this one answers "who is related, and how closely", off the pair table that
# `plasgenomicsutils ibd_fraction_and_snp_density` writes.

#: default edge weight: the fraction of the *callable* genome a pair shares
.PAIR_FRACTION_COL <- "ibd_fraction_accessible"

# The endpoints are columns, not something to recover by splitting the joined `pair` key --
# that is guesswork the moment a sample name contains the separator.
.pair_endpoints <- function(df) {
  miss <- setdiff(c("sample1", "sample2"), names(df))
  if (length(miss))
    stop("the pair table has no ", paste(miss, collapse = "/"), " column. Re-run ",
         "`plasgenomicsutils ibd_fraction_and_snp_density` -- it writes the endpoints ",
         "alongside `pair`.", call. = FALSE)
  df$sample1 <- as.character(df$sample1)
  df$sample2 <- as.character(df$sample2)
  df
}

#' Genome-wide IBD relatedness network
#'
#' Nodes are samples; an edge joins a pair sharing more than `min_ibd` of the genome IBD, and
#' its width is how much. [plot_ibd_network()] asks who shares one locus, where every edge means
#' the same thing; this asks who is related overall, so the amount is the point.
#'
#' Samples with no edge above the cutoff are drawn as a grid underneath, separated by a dashed
#' rule and counted -- an unrelated sample is a result, and dropping it silently would overstate
#' how connected the cohort is.
#'
#' @param pairs The per-pair table from `plasgenomicsutils ibd_fraction_and_snp_density`
#'   (`*.pair_ibd_fraction.tsv.gz`): a path, a data frame, or an [IbdResults] carrying one.
#'   Needs `sample1`/`sample2`, which that command writes.
#' @param meta Sample metadata for `color_group` / `shape_group`; taken from the [IbdResults]
#'   when `pairs` is one.
#' @param weight Column holding the fraction. Defaults to `ibd_fraction_accessible`, the
#'   callable-genome denominator; `ibd_fraction_full_genome` divides by the whole genome
#'   instead, so it reads lower for the same pair.
#' @param min_ibd Draw an edge only above this fraction (default `0.01`). Every pair is in the
#'   table, most of them sharing essentially nothing, so without a cutoff the graph is complete.
#' @param samples Optional subset of samples to keep.
#' @param color_group,colour_group,shape_group Metadata columns for node colour and shape.
#' @param colors,colours,shapes Values for those scales, named or positional (see [plot_ibd_network()]).
#' @param na_shape,na_colour,na_color What a sample with no value in those columns gets.
#'
#'   Legends stack in a fixed order -- colour, shape, then the IBD width -- so plots stay
#'   comparable; override with `+ ggplot2::guides(linewidth = ggplot2::guide_legend(order = 1))`.
#' @param include_isolated Show samples with no edge (default `TRUE`).
#' @param layout,spread,seed Layout algorithm, clique spreading, and the seed that makes it
#'   reproducible.
#' @param node_size,node_alpha,edge_colour,edge_color,edge_alpha Node and edge aesthetics.
#' @param weight_range Narrowest and widest edge, in `linewidth` units.
#' @param weight_breaks Legend breaks; defaults to powers of two spanning the data, since
#'   sharing runs over orders of magnitude.
#' @param title Plot title; `NULL` (default) writes one, `FALSE` or `NA` drops it.
#' @param subtitle `TRUE` (default) counts samples, edges and unconnected samples; a string
#'   replaces it, `FALSE` drops it.
#' @return A ggplot object.
#' @seealso [plot_ibd_network()] for one gene or locus.
#' @examples
#' \dontrun{
#' plot_ibd_pair_network("ibd_fraction.pair_ibd_fraction.tsv.gz", meta = meta,
#'                       color_group = "region", min_ibd = 0.03)
#' }
#' @export
plot_ibd_pair_network <- function(pairs, meta = NULL, weight = NULL, min_ibd = 0.01,
                                  samples = NULL,
                                  color_group = NULL, colors = NULL,
                                  shape_group = NULL, shapes = NULL,
                                  na_shape = .NA_SHAPE, na_colour = "grey70",
                                  include_isolated = TRUE, layout = "fr", spread = 1.5,
                                  node_size = 3, node_alpha = 0.9,
                                  edge_colour = "grey65", edge_alpha = 0.6,
                                  weight_range = c(0.15, 2.6), weight_breaks = NULL,
                                  title = NULL, subtitle = TRUE, seed = 42,
                                  colour_group = NULL, colours = NULL, na_color = NULL, edge_color = NULL) {
  color_group <- .alias_arg("color_group", "colour_group")
  colors <- .alias_arg("colors", "colours")
  na_colour <- .alias_arg("na_colour", "na_color")
  edge_colour <- .alias_arg("edge_colour", "edge_color")
  .need_package("ggplot2", "plot_ibd_pair_network()")
  .need_package("igraph", "plot_ibd_pair_network()")
  .need_package("ggraph", "plot_ibd_pair_network()")

  if (inherits(pairs, "IbdResults")) {
    if (is.null(meta)) meta <- pairs$get_meta()
    pf <- pairs$get_pair_fraction()
    if (is.null(pf))
      stop("this IbdResults has no pair table; build it with ",
           "ibd_results(pair_fraction = ) or call $set_pair_fraction()", call. = FALSE)
    pairs <- pf
  } else if (is.character(pairs) && length(pairs) == 1) {
    pairs <- .read_maybe(pairs, "pair table")
  }
  df <- as.data.frame(pairs, stringsAsFactors = FALSE)
  if (!nrow(df)) stop("the pair table is empty", call. = FALSE)
  df <- .pair_endpoints(df)

  if (is.null(weight)) weight <- .PAIR_FRACTION_COL
  if (!weight %in% names(df))
    stop("the pair table has no column '", weight, "'. Available: ",
         paste(names(df), collapse = ", "), call. = FALSE)
  df$.w <- suppressWarnings(as.numeric(df[[weight]]))

  all_samples <- unique(c(df$sample1, df$sample2))
  if (!is.null(samples)) {
    all_samples <- intersect(all_samples, as.character(samples))
    df <- df[df$sample1 %in% all_samples & df$sample2 %in% all_samples, , drop = FALSE]
  }
  keep <- !is.na(df$.w) & df$.w > min_ibd
  edges <- data.frame(from = df$sample1[keep], to = df$sample2[keep],
                      weight = df$.w[keep], stringsAsFactors = FALSE)
  if (!nrow(edges) && !include_isolated) {
    stop("no pair shares more than min_ibd = ", min_ibd,
         " (lower it, or set include_isolated = TRUE to still show the samples)", call. = FALSE)
  }

  .draw_ibd_network(
    edges = edges, analyzed = all_samples, meta = meta,
    color_group = color_group, colors = colors, shape_group = shape_group, shapes = shapes,
    na_shape = na_shape, na_colour = na_colour, include_isolated = include_isolated,
    layout = layout, spread = spread, node_size = node_size, node_alpha = node_alpha,
    edge_colour = edge_colour, edge_alpha = edge_alpha,
    weight_name = "IBD", weight_range = weight_range, weight_breaks = weight_breaks,
    title = if (is.null(title)) "Genome-wide IBD network" else title,
    subtitle = subtitle,
    subtitle_text = function(n_nodes, n_edges, n_iso)
      sprintf("%d samples, %d pairs sharing > %g of the genome%s", n_nodes, n_edges, min_ibd,
              if (n_iso) sprintf(" (%d unconnected)", n_iso) else ""),
    seed = seed)
}
