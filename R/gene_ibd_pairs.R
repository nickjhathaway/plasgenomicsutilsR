# The pair-level companion to gene_ibd_overlap(): instead of a fraction per group-pair,
# one row per sample pair whose IBD segment touches a gene, with how much of the gene that
# segment covers.

# Connected components over an edge list: single-linkage clustering, where a sample joins a
# cluster if it shares IBD with ANY member of it. Union-find in base R rather than igraph,
# which is only Suggested and should not become required to build a table.
# Ids are assigned largest cluster first, so `gene_cluster_id == 1` is the biggest group at
# that gene and the numbering means the same thing across genes.
.single_linkage <- function(a, b) {
  nodes <- unique(c(a, b))
  parent <- seq_along(nodes)
  find <- function(i) {
    while (parent[i] != i) i <- parent[i]
    i
  }
  ia <- match(a, nodes); ib <- match(b, nodes)
  for (k in seq_along(ia)) {
    ra <- find(ia[k]); rb <- find(ib[k])
    if (ra != rb) parent[ra] <- rb
  }
  root <- vapply(seq_along(nodes), find, integer(1))
  size <- tabulate(root, nbins = length(nodes))
  # rank by size, then by the first member's name, so the id is stable run to run
  first <- vapply(split(nodes, root)[as.character(sort(unique(root)))],
                  function(z) sort(z)[1], character(1))
  ord <- sort(unique(root))[order(-size[sort(unique(root))], first)]
  id <- match(root, ord)
  list(id = stats::setNames(id, nodes),
       size = stats::setNames(size[root], nodes))
}


#' Sample pairs sharing IBD over each gene
#'
#' An adjacency list of the pairs that are IBD across one or more genes: one row per
#' sample pair x IBD block x gene, saying how much of the gene the block covers. Pairs
#' with no IBD over a gene are simply absent -- use [gene_ibd_overlap()] when you need the
#' denominator (all pairs compared) rather than just the sharing ones.
#'
#' A pair appears more than once for a gene only if it has several separate IBD segments
#' spanning it.
#'
#' @param x An [IbdResults] built with `blocks =`.
#' @param genes Gene names from the object's track (default all), or a gene-interval data
#'   frame (`name`, `chr`/`chrom`, `start`, `end`). Several genes come back in one table.
#' @param within Pad each gene interval by this many bp on both sides when deciding
#'   whether a block overlaps it (default `0`). Coverage is always measured against the
#'   gene's own span, so a block that reaches only into the padding covers `0`.
#' @return A tibble with one row per pair x block x gene:
#'   \describe{
#'     \item{`sample1`, `sample2`}{the IBD pair, ordered so `sample1 < sample2`.}
#'     \item{`chr`}{chromosome of the block and gene.}
#'     \item{`block_start`, `block_end`}{the IBD segment, 0-based half-open.}
#'     \item{`gene`, `name`, `gene_id`}{gene labels (`gene` is unique, see
#'       [gene_ibd_overlap()]).}
#'     \item{`gene_start`, `gene_end`}{the gene interval, 0-based half-open.}
#'     \item{`coverage`}{`"complete"` when the block spans the whole gene, else
#'       `"partial"`.}
#'     \item{`covered_start`, `covered_end`}{the covered portion of the gene -- the gene's
#'       own bounds when `coverage` is `"complete"`.}
#'     \item{`gene_cluster_id`, `gene_cluster_size`}{single-linkage cluster of samples
#'       sharing at this gene, and how many samples are in it. A sample joins a cluster if
#'       it shares with **any** member, so a chain of pairs is one cluster even where its
#'       ends never share directly -- which is what [plot_ibd_network()] draws as a
#'       connected component. Ids run largest first, so `1` is the biggest group at that
#'       gene; they are per gene, so cluster 1 at `pfcrt` and cluster 1 at `pfdhps` are
#'       unrelated.}
#'     \item{`covered_bp`, `percent_covered`}{width of that portion, and it as a percentage
#'       of the gene's length.}
#'   }
#' @seealso [gene_ibd_overlap()] for the per-group-pair fractions,
#'   [plasgenomicsutilsR-coordinates] for the interval convention.
#' @examples
#' \dontrun{
#' ibd <- ibd_results(blocks = "hmm.txt", genes = PF_EXAMPLE_DRUG_GENES)
#' pairs <- gene_ibd_pairs(ibd, genes = c("pfcrt", "pfdhps"))
#' subset(pairs, coverage == "complete")
#' }
#' @export
gene_ibd_pairs <- function(x, genes = NULL, within = 0) {
  blocks <- x$get_blocks()
  if (is.null(blocks)) {
    stop("this IbdResults has no IBD blocks; build it with ibd_results(blocks = )",
         call. = FALSE)
  }
  gtrack <- .gene_track_for(x, genes)
  gid <- if ("gene_id" %in% names(gtrack)) as.character(gtrack$gene_id) else NULL

  blocks$chr <- normalise_chr(blocks$chr)
  by_chr <- split(seq_len(nrow(blocks)), blocks$chr)

  out <- vector("list", nrow(gtrack))
  for (i in seq_len(nrow(gtrack))) {
    f <- gtrack[i, ]
    gs <- as.numeric(f$start); ge <- as.numeric(f$end)
    idx <- by_chr[[f$chr]]
    if (is.null(idx) || !length(idx)) next
    d <- blocks[idx, , drop = FALSE]
    m <- d$start < (ge + within) & d$end > (gs - within)
    if (!any(m)) next
    d <- d[m, , drop = FALSE]

    cs <- pmax(d$start, gs)                       # intersection with the gene itself
    ce <- pmin(d$end, ge)
    covered <- pmax(0, ce - cs)
    complete <- covered >= (ge - gs)
    s1 <- pmin(d$sample1, d$sample2); s2 <- pmax(d$sample1, d$sample2)
    out[[i]] <- data.frame(
      sample1 = s1, sample2 = s2, chr = f$chr,
      block_start = d$start, block_end = d$end,
      gene = gtrack$.label[i], name = as.character(f$name),
      gene_id = if (is.null(gid)) NA_character_ else gid[i],
      gene_start = gs, gene_end = ge,
      coverage = ifelse(complete, "complete", "partial"),
      covered_start = ifelse(covered > 0, cs, NA_real_),
      covered_end = ifelse(covered > 0, ce, NA_real_),
      covered_bp = covered,
      percent_covered = if (ge > gs) 100 * covered / (ge - gs) else NA_real_,
      stringsAsFactors = FALSE)
    cl <- .single_linkage(s1, s2)
    out[[i]]$gene_cluster_id <- unname(cl$id[s1])
    out[[i]]$gene_cluster_size <- unname(cl$size[s1])
  }
  out <- out[!vapply(out, is.null, logical(1))]
  if (!length(out)) {
    return(tibble::tibble(
      sample1 = character(0), sample2 = character(0), chr = character(0),
      block_start = numeric(0), block_end = numeric(0), gene = character(0),
      name = character(0), gene_id = character(0), gene_start = numeric(0),
      gene_end = numeric(0), coverage = character(0), covered_start = numeric(0),
      covered_end = numeric(0), covered_bp = numeric(0), percent_covered = numeric(0),
      gene_cluster_id = integer(0), gene_cluster_size = integer(0)))
  }
  res <- do.call(rbind, out)
  res$gene <- factor(res$gene, levels = gtrack$.label[gtrack$.label %in% res$gene])
  res <- res[order(res$gene, res$sample1, res$sample2), , drop = FALSE]
  rownames(res) <- NULL
  tibble::as_tibble(res)
}
