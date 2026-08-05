# Map above-threshold IBD-selection SNPs onto a gene track: which genes sit under (or
# near) a positive-selection signal. Reads the selection table and the per-group
# significance threshold from an IbdResults object; the gene track is the object's, or
# one passed in (so you can scan the whole genome without attaching 5,000 genes to the
# object and disturbing the plot defaults).

# Coerce a gene table (path or data frame) to chr/start/end(/name/gene_id), accepting
# `chrom` as a `chr` alias (as in PF3D7_GENES) and normalising chromosome spellings.
.as_gene_track <- function(genes) {
  genes <- .read_maybe(genes, "genes")
  if (is.null(genes)) return(NULL)
  if (!"chr" %in% names(genes) && "chrom" %in% names(genes)) genes$chr <- genes$chrom
  .require_cols(genes, c("chr", "start", "end"), "genes")
  genes$chr <- normalise_chr(genes$chr)
  if (!"name" %in% names(genes)) {
    genes$name <- paste0(genes$chr, ":", genes$start, "-", genes$end)
  }
  genes
}

#' Genes under positive selection
#'
#' Intersects the significant SNPs of the IBD selection statistic (those at or above the
#' per-group Bonferroni threshold on `neg_log10_p`) with a gene track, returning the genes
#' hit by a selection signal. Because VCF filtering can leave the peak SNP just outside a
#' gene, a SNP counts for a gene when it falls within `within` bp of the gene interval
#' (default 2 kb), not only strictly inside it.
#'
#' Thresholding is always on `neg_log10_p` (the only metric the Bonferroni threshold is
#' defined for); the returned row reports **every** metric column of the selection table
#' at that gene's peak SNP (`peak_*`), so you still see `maf` / `z_score` / `chi2_stat`,
#' etc.
#'
#' @param x An [IbdResults] object (needs a `selection` table; per-group thresholds come
#'   from the object, or pass `threshold`).
#' @param within Maximum distance in bp between a significant SNP and the gene interval
#'   for the SNP to count (default `2000`). `0` requires the SNP to be strictly inside.
#' @param genes Optional gene track to scan (path or data frame with `chr`/`chrom`,
#'   `start`, `end`, and optionally `name`, `gene_id`), overriding the object's track.
#'   Pass [PF3D7_GENES] here to scan every gene without attaching it to the object (which
#'   would make the plot functions draw a line per gene). Default: the object's `genes`.
#' @param groups Optional subset of groups to scan (needs a `group` column); default all.
#' @param threshold Optional threshold override on the `neg_log10_p` scale: a single number
#'   for every group, or a named vector (`group -> threshold`). Default uses the object's
#'   thresholds.
#' @return A tibble, one row per (group, gene) hit, sorted by group then descending peak
#'   `neg_log10_p`: `group` (dropped for a group-less selection table), `gene_id`, `name`,
#'   `chr`, `gene_start`, `gene_end`, `n_snps` (significant SNPs in the window),
#'   `min_distance` (0 when a SNP is inside the gene), `peak_pos`, and a `peak_<metric>`
#'   column for every metric in the selection table (e.g. `peak_neg_log10_p`, `peak_maf`).
#' @examples
#' ibd <- example_ibd_results()
#' pos_selection_genes(ibd)                          # object's drug-gene track, 2 kb window
#' pos_selection_genes(ibd, within = 0)              # strictly-inside only
#' pos_selection_genes(ibd, genes = PF3D7_GENES)     # scan every gene, track left untouched
#' @export
pos_selection_genes <- function(x, within = 2000, genes = NULL,
                                groups = NULL, threshold = NULL) {
  sel <- x$get_selection()
  if (is.null(sel)) stop("this IbdResults has no selection table", call. = FALSE)
  if (!"neg_log10_p" %in% names(sel)) {
    stop("selection table has no 'neg_log10_p' column (pos_selection_genes thresholds ",
         "on that metric)", call. = FALSE)
  }
  gtrack <- if (!is.null(genes)) .as_gene_track(genes) else x$get_genes()
  if (is.null(gtrack) || !nrow(gtrack)) {
    stop("no gene track: pass genes= (e.g. genes = PF3D7_GENES) or build the object with ",
         "ibd_results(..., genes = )", call. = FALSE)
  }
  has_gid <- "gene_id" %in% names(gtrack)
  # every non-structural selection column is a "metric" reported at the peak SNP
  metric_cols <- setdiff(names(sel), c("chr", "pos", "group", "cum_pos"))

  thr_tbl <- x$get_thresholds()
  grp_levels <- if ("group" %in% names(sel)) {
    if (is.null(groups)) sort(unique(stats::na.omit(sel$group))) else groups
  } else {
    NA_character_                      # one global scan
  }

  thr_for <- function(g) {
    if (!is.null(threshold)) {
      if (!is.null(names(threshold)) && !is.na(g)) return(as.numeric(threshold[[as.character(g)]]))
      return(as.numeric(threshold[1]))
    }
    if (is.null(thr_tbl)) return(NA_real_)
    if (is.na(g) || !"group" %in% names(thr_tbl) || all(is.na(thr_tbl$group))) {
      return(thr_tbl$threshold[1])
    }
    v <- thr_tbl$threshold[thr_tbl$group == g]
    if (length(v)) v[1] else NA_real_
  }

  out <- list()
  for (g in grp_levels) {
    tv <- thr_for(g)
    if (!is.finite(tv)) {
      warning(sprintf("no finite threshold for group '%s'; pass threshold=",
                      if (is.na(g)) "global" else g), call. = FALSE)
      next
    }
    sub  <- if (is.na(g)) sel else sel[!is.na(sel$group) & sel$group == g, , drop = FALSE]
    hits <- sub[is.finite(sub$neg_log10_p) & sub$neg_log10_p >= tv, , drop = FALSE]
    if (!nrow(hits)) next
    for (ch in unique(hits$chr)) {
      gc <- gtrack[gtrack$chr == ch, , drop = FALSE]
      if (!nrow(gc)) next
      hc <- hits[hits$chr == ch, , drop = FALSE]
      for (i in seq_len(nrow(gc))) {
        # `pos` and the gene interval are both 0-based, the interval half-open
        m <- hc$pos >= gc$start[i] - within & hc$pos < gc$end[i] + within
        if (!any(m)) next
        win  <- hc[m, , drop = FALSE]
        dist <- pmax(0, gc$start[i] - win$pos, win$pos - (gc$end[i] - 1))   # 0 when inside
        prow <- win[which.max(win$neg_log10_p), , drop = FALSE]       # peak by neg_log10_p
        base <- data.frame(
          group = if (is.na(g)) NA_character_ else as.character(g),
          gene_id = if (has_gid) gc$gene_id[i] else NA_character_,
          name = gc$name[i], chr = ch,
          gene_start = gc$start[i], gene_end = gc$end[i],
          n_snps = nrow(win), min_distance = min(dist), peak_pos = prow$pos,
          stringsAsFactors = FALSE)
        peak <- stats::setNames(prow[, metric_cols, drop = FALSE], paste0("peak_", metric_cols))
        out[[length(out) + 1L]] <- cbind(base, peak)
      }
    }
  }
  if (!length(out)) {
    empty <- data.frame(group = character(), gene_id = character(), name = character(),
                        chr = character(), gene_start = integer(), gene_end = integer(),
                        n_snps = integer(), min_distance = numeric(), peak_pos = integer(),
                        stringsAsFactors = FALSE)
    for (mc in paste0("peak_", metric_cols)) empty[[mc]] <- numeric()
    return(tibble::as_tibble(empty))
  }
  res <- do.call(rbind, out)
  res <- res[order(res$group, -res$peak_neg_log10_p), , drop = FALSE]
  rownames(res) <- NULL
  if (all(is.na(res$group))) res$group <- NULL   # group-less selection table
  tibble::as_tibble(res)
}
