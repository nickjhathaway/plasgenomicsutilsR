# Single-linkage IBD clusters, written into the stored metadata so any plot that colours or
# shapes by a metadata column can use them.

# Sample -> cluster id / size over one interval, using the same edge rule as
# plot_ibd_network(): `within` pads the interval, `sharing` decides whether a block has to
# span the whole thing or merely touch it. Shares `.single_linkage()` with gene_ibd_pairs(),
# so the ids are the same numbers that table reports.
.cluster_over_interval <- function(blocks, iv, within, sharing) {
  bl <- blocks
  bl$chr <- normalise_chr(bl$chr)
  lo <- iv$start - within
  hi <- iv$end + within
  m <- bl$chr == iv$chr & if (sharing == "complete") {
    bl$start <= lo & bl$end >= hi
  } else {
    bl$start < hi & bl$end > lo
  }
  d <- bl[m, , drop = FALSE]
  if (!nrow(d)) return(NULL)
  a <- pmin(d$sample1, d$sample2); b <- pmax(d$sample1, d$sample2)
  .single_linkage(a, b)
}

#' Add IBD cluster ids to the stored metadata
#'
#' For each gene or locus, works out the single-linkage clusters of samples sharing IBD over
#' it and writes them into the object's metadata as `<label>_cluster_id` (and
#' `<label>_cluster_size` when `size = TRUE`). Any plot that reads a metadata column can then
#' use them -- notably [plot_ibd_network()], where the clusters *are* the connected components
#' being drawn:
#'
#' ```
#' ibd$add_ibd_clusters(genes = "pfcrt")
#' plot_ibd_network(ibd, gene = "pfcrt", color_group = "pfcrt_cluster_id")
#' ```
#'
#' Single linkage means a sample joins a cluster if it shares with **any** member, so a chain
#' of pairs is one cluster even where its ends never share directly. Ids run largest cluster
#' first and are numbered per interval, so cluster 1 at one gene is unrelated to cluster 1 at
#' another. A sample that shares with nobody over the interval gets `NA`.
#'
#' `within` and `sharing` mean exactly what they do in [plot_ibd_network()] and
#' [gene_ibd_pairs()], and they decide which edges exist -- so pass the same values you plot
#' with, or the colours will not match the components.
#'
#' @param x An [IbdResults] with `blocks` and `meta` loaded.
#' @param genes Gene names in the object's track, a gene-interval data frame (`name`,
#'   `chr`/`chrom`, `start`, `end`) to use instead of that track, or `NULL` for every gene
#'   in it. Ignored
#'   when `locus` is given.
#' @param locus A locus instead of a gene: `"chr:pos"`, `"chr:start-end"`, or a data frame
#'   with `chr`/`start`/`end`. Its `name` (or the coordinate string) labels the column.
#' @param within Pad the interval by this many bp on both sides before deciding whether a
#'   block overlaps it (default `0`).
#' @param sharing `"overlap"` (default) counts a pair when any IBD block touches the
#'   interval; `"complete"` only when a block spans the whole of it.
#' @param size Also add `<label>_cluster_size`, the number of samples in the cluster
#'   (default `FALSE`).
#' @param prefix Optional string put before each column name, for keeping several settings
#'   side by side (e.g. `prefix = "complete_"`).
#' Re-running replaces any column of the same name, so changing `within` or `sharing` and
#' calling again re-clusters rather than accumulating stale columns. Use `prefix` when you
#' want two settings side by side.
#'
#' @return Invisibly `x`, with the metadata extended. `x$get_meta()` shows the new columns.
#' @examples
#' \dontrun{
#' ibd <- ibd_results(blocks = "blocks.hmm.txt", meta = meta,
#'                    genes = PF_EXAMPLE_DRUG_GENES)
#' add_ibd_clusters(ibd)                       # one column per gene in the track
#' plot_ibd_network(ibd, gene = "pfcrt", color_group = "pfcrt_cluster_id")
#' }
#' @export
add_ibd_clusters <- function(x, genes = NULL, locus = NULL, within = 0,
                             sharing = c("overlap", "complete"), size = FALSE,
                             prefix = "") {
  sharing <- match.arg(sharing)
  blocks <- x$get_blocks()
  if (is.null(blocks))
    stop("this IbdResults has no IBD blocks; build it with ibd_results(blocks = )",
         call. = FALSE)
  meta <- x$get_meta()
  if (is.null(meta))
    stop("add_ibd_clusters() writes into the metadata; build with ibd_results(meta = )",
         call. = FALSE)

  ivs <- if (!is.null(locus)) {
    list(.resolve_locus(x, NULL, locus))
  } else {
    # names select from the object's track; a table replaces it, so a gene the object was
    # not built with can be clustered without rebuilding it
    g <- .gene_track_for(x, genes)
    supplied <- if (is.null(genes) || is.character(genes)) NULL else g
    lapply(as.character(g$name), function(nm) .resolve_locus(x, nm, NULL, supplied))
  }

  samples <- as.character(meta$sample)
  had <- names(meta)
  added <- character(0)
  for (iv in ivs) {
    cl <- .cluster_over_interval(blocks, iv, within, sharing)
    col <- paste0(prefix, iv$label, "_cluster_id")
    # a factor ordered 1, 2, 3 ... so a plot's legend and palette follow cluster size rather
    # than the alphabet ("10" before "2")
    ids <- if (is.null(cl)) rep(NA_integer_, length(samples)) else unname(cl$id[samples])
    meta[[col]] <- factor(ids, levels = sort(unique(stats::na.omit(ids))))
    added <- c(added, col)
    if (size) {
      scol <- paste0(prefix, iv$label, "_cluster_size")
      meta[[scol]] <- if (is.null(cl)) NA_integer_ else unname(cl$size[samples])
      added <- c(added, scol)
    }
  }
  x$set_meta(meta)
  # a deterministic column name means a re-run replaces rather than accumulates
  message(if (any(added %in% had)) "updated " else "added ", length(added),
          " column(s): ", paste(added, collapse = ", "))
  invisible(x)
}
