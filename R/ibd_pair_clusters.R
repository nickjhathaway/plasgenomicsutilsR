# The genome-wide companion to add_ibd_clusters() / gene_ibd_pairs(), which cluster over one
# gene or locus. Here the edges are whole-genome relatedness, so there is one clustering for
# the cohort rather than one per interval -- what plot_ibd_pair_network() draws.

#' Genome-wide IBD clusters, and which samples are connected
#'
#' Who is related to whom across the genome as a table: one row per sample, saying whether it
#' shares more than `min_ibd` with anyone, and which single-linkage cluster it falls in. These
#' are the components [plot_ibd_pair_network()] draws -- the connected blobs and the grid of
#' unconnected samples beneath them -- so the same `weight` and `min_ibd` give the same answer
#' the picture shows.
#'
#' Single linkage means a sample joins a cluster if it shares with **any** member, so a chain
#' of pairs is one cluster even where its ends never share directly. Ids run largest cluster
#' first, so `cluster_id == 1` is the biggest group; a sample sharing with nobody gets `NA`.
#' [add_ibd_clusters()] and [gene_ibd_pairs()] number their per-gene clusters the same way.
#'
#' @param pairs The per-pair table from `plasgenomicsutils ibd_fraction_and_snp_density`
#'   (`*.pair_ibd_fraction.tsv.gz`): a path, a data frame, or an [IbdResults] carrying one.
#' @param weight Column holding the fraction. Defaults to `ibd_fraction_accessible`, the
#'   callable-genome denominator; `ibd_fraction_full_genome` divides by the whole genome
#'   instead, so it reads lower for the same pair.
#' @param min_ibd Count a pair as linked only above this fraction (default `0.01`). Pass what
#'   you plotted with, or the clusters will not match the components.
#' @param samples Optional subset of samples to keep.
#' @param add_meta_cols Metadata columns to carry onto each row, by sample. A sample missing
#'   from `meta` gets `NA`; factor columns keep their level order, so a region stays in map
#'   order rather than turning alphabetical.
#' @param meta Sample metadata for `add_meta_cols`; taken from the [IbdResults] when `pairs`
#'   is one, so it only needs giving alongside a plain table.
#' @return A tibble with one row per sample, connected samples first (largest cluster first,
#'   then by name) and the unconnected after them:
#'   \describe{
#'     \item{`sample`}{the sample.}
#'     \item{`connected`}{whether it shares more than `min_ibd` with any other sample.}
#'     \item{`cluster_id`, `cluster_size`}{its single-linkage cluster and how many samples are
#'       in it; `NA` when it is connected to nobody.}
#'     \item{`n_links`}{how many samples it is linked to.}
#'     \item{`max_ibd`}{the largest fraction it shares with any of them, `NA` when none.}
#'   }
#'   plus one column for each of `add_meta_cols`.
#' @seealso [plot_ibd_pair_network()] for the picture, [ibd_pair_links()] for the edges
#'   themselves, [add_ibd_clusters()] and [gene_ibd_pairs()] for the per-gene equivalents.
#' @examples
#' \dontrun{
#' cl <- ibd_pair_clusters(ibd_all, min_ibd = 0.03)
#'
#' # the connected samples, and the unconnected ones
#' cl$sample[cl$connected]
#' cl$sample[!cl$connected]
#'
#' # the biggest cluster
#' subset(cl, cluster_id == 1)
#'
#' # with where each sample came from, to see whether a cluster spans sites
#' cl <- ibd_all$ibd_pair_clusters(min_ibd = 0.03, add_meta_cols = c("region", "country"))
#' with(subset(cl, connected), table(cluster_id, region))
#'
#' # colour the network by cluster, the way add_ibd_clusters() does per gene
#' ibd_all$set_meta(merge(ibd_all$get_meta(), cl[, c("sample", "cluster_id")], by = "sample"))
#' ibd_all$plot_ibd_pair_network(min_ibd = 0.03, color_group = "cluster_id")
#' }
#' @export
ibd_pair_clusters <- function(pairs, weight = NULL, min_ibd = 0.01, samples = NULL,
                              add_meta_cols = NULL, meta = NULL) {
  if (!is.null(add_meta_cols) && is.null(meta) && inherits(pairs, "IbdResults"))
    meta <- pairs$get_meta()
  meta <- .normalise_meta(meta)
  pe <- .pair_edges(pairs, weight, min_ibd, samples)
  edges <- pe$edges
  analyzed <- sort(pe$analyzed)

  out <- data.frame(sample = analyzed, connected = FALSE,
                    cluster_id = NA_integer_, cluster_size = NA_integer_,
                    n_links = 0L, max_ibd = NA_real_, stringsAsFactors = FALSE)
  if (nrow(edges)) {
    cl <- .single_linkage(edges$from, edges$to)
    i <- match(analyzed, names(cl$id))
    out$cluster_id <- unname(cl$id[i])
    out$cluster_size <- unname(cl$size[i])
    out$connected <- !is.na(out$cluster_id)
    # both endpoints of every edge, so a sample's links and its closest relative are one
    # tabulation over the stacked ends rather than two passes over the pair table
    ends <- c(edges$from, edges$to)
    w <- c(edges$weight, edges$weight)
    out$n_links <- unname(table(factor(ends, levels = analyzed))[analyzed])
    best <- vapply(split(w, factor(ends, levels = analyzed)), function(z)
      if (length(z)) max(z) else NA_real_, numeric(1))
    out$max_ibd <- unname(best[analyzed])
  }
  out$n_links <- as.integer(out$n_links)
  # connected first, largest cluster first, then by name; unconnected after them
  out <- out[order(!out$connected, out$cluster_id, out$sample), , drop = FALSE]
  out <- .add_sample_meta(out, add_meta_cols, meta)
  rownames(out) <- NULL
  tibble::as_tibble(out)
}

#' The IBD pairs a network draws
#'
#' The edge list behind [plot_ibd_pair_network()]: one row per sample pair sharing more than
#' `min_ibd` of the genome. [ibd_pair_clusters()] answers which samples those edges connect;
#' this is the sharing itself.
#'
#' `add_meta_cols` carries metadata onto both ends of each edge, which is what turns the list
#' into something to summarise over -- whether a pair is within a site or between two, say:
#'
#' ```
#' e <- ibd_all$ibd_pair_links(min_ibd = 0.03, add_meta_cols = "region")
#' table(within_region = e$sample1_region == e$sample2_region)
#' ```
#'
#' @inheritParams ibd_pair_clusters
#' @param add_meta_cols Metadata columns to attach to both endpoints. Each `col` becomes
#'   `sample1_col` and `sample2_col`, in the order asked for. A sample missing from `meta`
#'   gets `NA`; factor columns keep their level order, so a region stays in map order rather
#'   than turning alphabetical.
#' @param meta Sample metadata for `add_meta_cols`; taken from the [IbdResults] when `pairs`
#'   is one, so it only needs giving alongside a plain table.
#' @return A tibble of `sample1`, `sample2` and `ibd_fraction` (the `weight` column's value),
#'   highest sharing first, plus a pair of columns for each of `add_meta_cols`.
#' @seealso [ibd_pair_clusters()], [plot_ibd_pair_network()],
#'   [pair_fraction_summary()] for the sharing already summarised by group.
#' @examples
#' \dontrun{
#' ibd_pair_links(ibd_all, min_ibd = 0.03)
#'
#' # the same edges, labelled with where each end came from
#' ibd_all$ibd_pair_links(min_ibd = 0.03, add_meta_cols = c("region", "country"))
#' }
#' @export
ibd_pair_links <- function(pairs, weight = NULL, min_ibd = 0.01, samples = NULL,
                           add_meta_cols = NULL, meta = NULL) {
  if (!is.null(add_meta_cols) && is.null(meta) && inherits(pairs, "IbdResults"))
    meta <- pairs$get_meta()
  meta <- .normalise_meta(meta)
  e <- .pair_edges(pairs, weight, min_ibd, samples)$edges
  out <- data.frame(sample1 = e$from, sample2 = e$to, ibd_fraction = e$weight,
                    stringsAsFactors = FALSE)
  out <- out[order(-out$ibd_fraction, out$sample1, out$sample2), , drop = FALSE]
  out <- .add_endpoint_meta(out, add_meta_cols, meta)
  rownames(out) <- NULL
  tibble::as_tibble(out)
}
