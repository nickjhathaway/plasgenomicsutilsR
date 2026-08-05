# Per-gene IBD-block overlap between sample groups, computed in R from the blocks loaded
# into an IbdResults (the ad-hoc counterpart of `plasgenomicsutils ibd_gene_overlap`).

# GFF `Name` repeats across gene families (var, rifin, ...) while `gene_id` is unique.
# Build a unique display label per feature: the plain name where it is unique, and
# `name (gene_id)` (or make.unique) where it collides -- so tables and plot facets stay
# 1:1 with genes.
.disambiguate_gene_labels <- function(name, gene_id = NULL) {
  name <- as.character(name)
  dup <- name %in% name[duplicated(name)]
  if (!any(dup)) return(name)
  if (!is.null(gene_id) && any(!is.na(gene_id))) {
    lab <- ifelse(dup & !is.na(gene_id), paste0(name, " (", gene_id, ")"), name)
    if (anyDuplicated(lab)) lab <- make.unique(lab)
    lab
  } else {
    make.unique(name)
  }
}

#' Per-gene IBD-block overlap between groups
#'
#' For each gene and each pair of groups, the fraction of sample pairs that share an IBD
#' **block overlapping the gene** -- a pair counts when any of its IBD segments overlaps the
#' (optionally padded) gene interval, so a segment spanning the gene counts even with no
#' genotyped SNP inside it. The denominator is *all* pairs compared in that group-pair (from
#' the analyzed-sample set), so pairs that are never IBD still count against the fraction.
#'
#' Needs an [IbdResults] built with `blocks =` and `meta =` (or use the Python
#' `ibd_gene_overlap` tool and pass its table as `gene_overlap =`).
#'
#' @param x An [IbdResults] with IBD `blocks` and `meta`.
#' @param genes Gene names from the object's track (default all), or a gene-interval data
#'   frame (`name`, `chr`/`chrom`, `start`, `end`) to use instead of the track.
#' @param group Metadata column defining the groups (default: the first non-`sample`
#'   column of `meta`).
#' @param within Pad each gene interval by this many bp on both sides (default `0`).
#' @return A tibble, one row per gene x group-pair: `gene` (a unique display label -- the
#'   `name`, disambiguated as `name (gene_id)` where GFF Names repeat across gene
#'   families), `name`, `gene_id`, `chr`, `start`, `end`, `group_a`, `group_b`,
#'   `n_pairs_ibd`, `n_pairs_total`, `frac_pairs_ibd`.
#' @export
gene_ibd_overlap <- function(x, genes = NULL, group = NULL, within = 0) {
  blocks <- x$get_blocks()
  if (is.null(blocks)) {
    stop("this IbdResults has no IBD blocks; build it with ibd_results(blocks = , meta = )",
         call. = FALSE)
  }
  meta <- x$get_meta()
  if (is.null(meta) || !"sample" %in% names(meta)) {
    stop("block grouping needs meta with a 'sample' column; pass ibd_results(meta = )",
         call. = FALSE)
  }
  if (is.null(group)) group <- setdiff(names(meta), "sample")[1]
  if (!group %in% names(meta)) stop("meta has no column '", group, "'", call. = FALSE)
  s2g <- stats::setNames(as.character(meta[[group]]), as.character(meta$sample))

  # group counts among the analyzed samples (denominator), known-group only
  ag <- s2g[x$get_analyzed_samples()]
  ag <- ag[!is.na(ag)]
  cnt <- table(ag)
  grp_levels <- sort(names(cnt))
  if (length(grp_levels) < 2) stop("need at least two groups among analyzed samples", call. = FALSE)
  total_pairs <- function(a, b) {
    if (a == b) cnt[[a]] * (cnt[[a]] - 1) / 2 else cnt[[a]] * cnt[[b]]
  }

  # gene intervals: names into the track, or a supplied interval data frame
  if (is.null(genes) || is.character(genes)) {
    gtrack <- x$get_genes()
    if (is.null(gtrack)) stop("no gene track; pass genes= or ibd_results(genes = )", call. = FALSE)
    if (is.character(genes)) {
      .check_gene_request(gtrack, genes)
      gtrack <- gtrack[tolower(gtrack$name) %in% tolower(genes), , drop = FALSE]
    }
  } else {
    gtrack <- .as_gene_track(genes)
  }
  gtrack$chr <- normalise_chr(gtrack$chr)
  # GFF Names repeat across gene families (var/rifin/...); gene_id is unique. Key each
  # feature by a unique display label so duplicate names don't collapse in the table/facets.
  gid <- if ("gene_id" %in% names(gtrack)) as.character(gtrack$gene_id) else NULL
  gtrack$.label <- .disambiguate_gene_labels(gtrack$name, gid)
  if (!identical(gtrack$.label, as.character(gtrack$name))) {
    warning("some gene Names repeat (var/rifin/...); disambiguated by gene_id in the ",
            "`gene` column (`name`/`gene_id` kept as columns)", call. = FALSE)
  }

  blocks$chr <- normalise_chr(blocks$chr)
  by_chr <- split(seq_len(nrow(blocks)), blocks$chr)

  gp <- expand.grid(a = grp_levels, b = grp_levels, stringsAsFactors = FALSE)
  gp <- gp[gp$a <= gp$b, , drop = FALSE]
  gp_key <- paste(gp$a, gp$b, sep = "\r")

  out <- vector("list", nrow(gtrack))
  for (i in seq_len(nrow(gtrack))) {
    f <- gtrack[i, ]
    lo <- as.numeric(f$start) - within
    hi <- as.numeric(f$end) + within
    counts <- stats::setNames(integer(nrow(gp)), gp_key)
    idx <- by_chr[[f$chr]]
    if (!is.null(idx) && length(idx)) {
      d <- blocks[idx, , drop = FALSE]
      m <- d$start < hi & d$end > lo                      # half-open [start, end) overlap
      if (any(m)) {
        s1 <- d$sample1[m]; s2 <- d$sample2[m]
        key <- ifelse(s1 < s2, paste(s1, s2), paste(s2, s1))   # distinct pairs
        keep <- !duplicated(key)
        ga <- s2g[s1[keep]]; gb <- s2g[s2[keep]]
        ok <- !is.na(ga) & !is.na(gb) & ga %in% grp_levels & gb %in% grp_levels
        ga <- ga[ok]; gb <- gb[ok]
        a2 <- ifelse(ga <= gb, ga, gb); b2 <- ifelse(ga <= gb, gb, ga)
        tab <- table(paste(a2, b2, sep = "\r"))
        counts[names(tab)] <- as.integer(tab)
      }
    }
    tp <- mapply(total_pairs, gp$a, gp$b)
    num <- as.integer(counts[gp_key])
    out[[i]] <- data.frame(
      gene = gtrack$.label[i], name = as.character(f$name),
      gene_id = if (is.null(gid)) NA_character_ else gid[i],
      chr = f$chr, start = as.numeric(f$start), end = as.numeric(f$end),
      group_a = gp$a, group_b = gp$b, n_pairs_ibd = num, n_pairs_total = as.numeric(tp),
      frac_pairs_ibd = ifelse(tp > 0, num / tp, NA_real_), stringsAsFactors = FALSE)
  }
  res <- do.call(rbind, out)
  res$gene <- factor(res$gene, levels = gtrack$.label)
  tibble::as_tibble(res)
}
