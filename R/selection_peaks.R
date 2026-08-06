# Collapse a per-SNP selection scan into the loci it actually implicates.
#
# Every scan in this package reports one row per SNP, and a sweep spans many SNPs, so a
# count of significant SNPs is not a count of findings -- it mostly reflects how densely
# the region was genotyped. Merging adjacent hits into peaks gives the number that belongs
# in a sentence, and the interval to look up.
#
# This is also the honest way to read a scan whose p-values are not well calibrated (see
# `lambda_gc` from `ibd_selection_statistic`): the *ranking* survives a wrong null, the
# per-SNP tail probability does not.

#: SNPs closer together than this are treated as one peak.
PEAK_GAP_BP <- 20000L

.peak_value_col <- function(df, metric) {
  if (!is.null(metric)) {
    if (!metric %in% names(df))
      stop(sprintf("no `%s` column in the scan", metric), call. = FALSE)
    return(metric)
  }
  hit <- intersect(c("neg_log10_p", "beta", "value", "ihs"), names(df))
  if (!length(hit)) stop("cannot tell which column to rank by; set `metric`", call. = FALSE)
  hit[1]
}

# Which rows count as hits, given the criterion.
.peak_hits <- function(df, criterion, metric, cutoff, top, thresholds) {
  v <- df[[metric]]
  if (criterion == "value") {
    if (is.null(cutoff)) stop("`criterion = \"value\"` needs a `cutoff`", call. = FALSE)
    return(is.finite(v) & v >= cutoff)
  }
  if (criterion == "top") {
    keep <- rep(FALSE, nrow(df))
    for (g in unique(df$.grp)) {
      i <- which(df$.grp == g & is.finite(v))
      if (!length(i)) next
      keep[i] <- v[i] >= stats::quantile(v[i], 1 - top, names = FALSE)
    }
    return(keep)
  }
  # bonferroni / fdr: prefer the flag the scan already carries, else the stored line
  flag <- if (criterion == "fdr") "significant_fdr" else "significant"
  if (flag %in% names(df)) return(.flag_true(df[[flag]]))
  if (is.null(thresholds))
    stop(sprintf("the scan has no `%s` column and no thresholds were supplied", flag),
         call. = FALSE)
  col <- if (criterion == "fdr") "neg_log10_p_fdr_threshold" else "threshold"
  if (!col %in% names(thresholds))
    stop(sprintf("the thresholds have no `%s`; regenerate with a current %s", col,
                 "`ibd_selection_statistic`"), call. = FALSE)
  line <- thresholds[[col]][match(as.character(df$.grp), as.character(thresholds$group))]
  if (all(is.na(line)) && nrow(thresholds) == 1L) line <- thresholds[[col]][1]
  is.finite(v) & is.finite(line) & v >= line
}

# a flag column may be logical, character "True"/"False" (from a Python TSV), or NA
.flag_true <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  !is.na(x) & tolower(as.character(x)) %in% c("true", "1", "yes")
}

#' Merge a selection scan's significant SNPs into peaks
#'
#' Turns a per-SNP scan into one row per implicated locus. A sweep leaves a run of
#' neighbouring significant SNPs, so the SNP count answers "how dense is my panel here"
#' while the peak count answers "how many things did I find" -- and only the second belongs
#' in a result.
#'
#' Works on any of the per-SNP scans in the package: an [IbdResults] (its selection table),
#' or the tibble from [run_ihs()], [run_rsb()], [run_xpehh()] or [beta_score()].
#'
#' @param x An [IbdResults], or a scan data frame with `chr`, `pos` and a value column.
#' @param criterion How a SNP qualifies:
#'   * `"bonferroni"` (default) / `"fdr"` — the scan's own `significant` /
#'     `significant_fdr` flag, or the matching stored threshold line.
#'   * `"top"` — the top `top` fraction within each group, the convention for scans whose
#'     null is not trustworthy.
#'   * `"value"` — at or above `cutoff`.
#' @param metric Column to rank by; guessed from the scan when `NULL`
#'   (`neg_log10_p`, then `beta`, `value`, `ihs`).
#' @param cutoff Value cutoff for `criterion = "value"`.
#' @param top Tail fraction for `criterion = "top"` (default 0.01).
#' @param gap Merge hits separated by at most this many bp (default `r PEAK_GAP_BP`).
#' @param min_snps Drop peaks supported by fewer than this many significant SNPs. Raising
#'   it to 2 or 3 is the cheapest way to suppress isolated single-SNP hits.
#' @param pad Widen each peak by this many bp on both sides before reporting and before
#'   matching genes.
#' @param genes Optional gene table (`name`, `chr` or `chrom`, `start`, `end`) to annotate
#'   each peak with the genes it overlaps; [PF3D7_GENES] is a sensible argument.
#' @param thresholds Threshold table, when `x` is a bare data frame carrying no flag column.
#' @return A tibble with one row per peak: the grouping column, `chr`, `start`, `end`,
#'   `width`, `n_snps` (significant SNPs in the peak), `peak_pos` (the single
#'   highest-scoring SNP), `peak_value`, `mean_value`, and -- when a gene table was given
#'   -- four annotation columns. Sorted by `peak_value`.
#'
#'   Everything is anchored on `peak_pos`, not the interval's midpoint: the midpoint is an
#'   artefact of where merging started and stopped, and can land in a gap between genes.
#'
#'   * `peak_genes` — the gene(s) whose span **contains** the peak SNP, comma-separated.
#'     Usually one, since gene spans rarely overlap, and **empty when the peak SNP is
#'     intergenic** — about a third of peaks on a real cohort.
#'   * `nearest_gene`, `distance_to_gene` — of the genes the peak **covers**, the closest
#'     one to the peak SNP and the gap to it in bp, `0` when the SNP is inside it. This is
#'     what to read when `peak_genes` is empty; those intergenic peaks are typically within
#'     a kb or two of a gene in the same peak. Candidates are restricted to the interval,
#'     so a peak that covers no gene reports none rather than pointing at something far
#'     outside it — which matters when `genes` is a short list, where most peaks cover
#'     nothing from it.
#'   * `n_genes` — how many genes the interval spans. A peak of a few hundred kb, as real
#'     IBD sharing regions are, can cross dozens; naming them all helps nobody, so this
#'     counts them. Narrow `gap` if the intervals are wider than you want.
#'
#'   The three agree by construction: `n_genes == 0` implies both name columns are empty,
#'   and a non-empty `peak_genes` implies `distance_to_gene == 0`.
#' @seealso [ihs_genes()] and [beta_genes()] for the per-gene view; this is the per-locus
#'   one, which does not need a gene to exist where the signal is.
#' @examples
#' ps <- example_pop_structure(umap = FALSE)
#' b <- beta_score(ps, group = "country", window = 300000, min_window_snps = 1)
#' selection_peaks(b, criterion = "top", top = 0.2, genes = PF_EXAMPLE_DRUG_GENES)
#' @export
selection_peaks <- function(x, criterion = c("bonferroni", "fdr", "top", "value"),
                            metric = NULL, cutoff = NULL, top = 0.01,
                            gap = PEAK_GAP_BP, min_snps = 1L, pad = 0,
                            genes = NULL, thresholds = NULL) {
  criterion <- match.arg(criterion)
  if (inherits(x, "IbdResults")) {
    if (is.null(thresholds)) thresholds <- x$get_thresholds()
    df <- x$get_selection()
    if (is.null(df)) stop("this IbdResults has no selection table", call. = FALSE)
  } else {
    df <- as.data.frame(x)
  }
  if (!nrow(df)) return(.empty_peaks())
  need <- c("chr", "pos")
  if (!all(need %in% names(df)))
    stop("the scan needs `chr` and `pos` columns", call. = FALSE)

  by <- intersect(c("group", "pair"), names(df))[1]
  df$.grp <- if (is.na(by)) factor("all") else .as_group_factor(df[[by]])
  metric <- .peak_value_col(df, metric)
  df$.v <- df[[metric]]

  hits <- .peak_hits(df, criterion, metric, cutoff, top, thresholds)
  df <- df[hits, , drop = FALSE]
  if (!nrow(df)) {
    warning("no SNP met the criterion", call. = FALSE)
    return(.empty_peaks())
  }

  rows <- list()
  for (g in levels(df$.grp)) {
    sub_g <- df[df$.grp == g, , drop = FALSE]
    if (!nrow(sub_g)) next
    for (chr in unique(sub_g$chr)) {
      s <- sub_g[sub_g$chr == chr, , drop = FALSE]
      s <- s[order(s$pos), , drop = FALSE]
      id <- cumsum(c(TRUE, diff(s$pos) > gap))
      for (k in unique(id)) {
        p <- s[id == k, , drop = FALSE]
        if (nrow(p) < min_snps) next
        top_i <- which.max(p$.v)
        rows[[length(rows) + 1L]] <- data.frame(
          grp = g, chr = as.character(chr),
          start = min(p$pos) - pad, end = max(p$pos) + 1 + pad,
          n_snps = nrow(p), peak_pos = p$pos[top_i],
          peak_value = p$.v[top_i], mean_value = mean(p$.v, na.rm = TRUE),
          stringsAsFactors = FALSE)
      }
    }
  }
  if (!length(rows)) {
    warning("no peak had at least `min_snps` SNPs", call. = FALSE)
    return(.empty_peaks())
  }
  out <- do.call(rbind, rows)
  out$start <- pmax(0, out$start)
  out$width <- out$end - out$start
  out <- out[c("grp", "chr", "start", "end", "width", "n_snps", "peak_pos",
               "peak_value", "mean_value")]
  names(out)[1] <- if (is.na(by)) "group" else by
  if (!is.null(genes)) {
    ann <- .annotate_peaks(out, genes)
    out$peak_genes <- ann$peak_genes
    out$nearest_gene <- ann$nearest_gene
    out$distance_to_gene <- ann$distance_to_gene
    out$n_genes <- ann$n_genes
  }
  out <- out[order(-out$peak_value), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "peak_metric") <- metric
  attr(out, "peak_criterion") <- criterion
  tibble::as_tibble(out)
}

.empty_peaks <- function() {
  tibble::as_tibble(data.frame(
    group = character(), chr = character(), start = numeric(), end = numeric(),
    width = numeric(), n_snps = integer(), peak_pos = numeric(),
    peak_value = numeric(), mean_value = numeric(), stringsAsFactors = FALSE))
}

# What sits at the peak, and how wide a net the interval casts.
#
# Everything is anchored on `peak_pos` -- the single highest-scoring SNP -- rather than the
# interval's midpoint, which is an artefact of where merging happened to start and stop and
# can easily land in a gap between genes. Listing every gene in the interval is useless
# once a peak runs to a few hundred kb, which real IBD sharing regions do, so the interval
# gets a count instead of a roll-call.
#
# The peak SNP is often intergenic (about a third of the time on a real cohort, usually
# within a kb or two of something), so the nearest gene *within the peak* is reported
# alongside. `peak_genes` and `nearest_gene` are kept apart on purpose: "the SNP is in this
# gene" and "the SNP is near this gene" are different claims and collapsing them would
# quietly overstate one. Both stay inside the peak, so all three columns describe the same
# interval and a peak covering no gene says so in every one of them.
.annotate_peaks <- function(peaks, genes) {
  g <- .gene_track(genes)
  key <- normalise_chr(peaks$chr)
  at_peak <- character(nrow(peaks))
  nearest <- character(nrow(peaks))
  dist <- rep(NA_real_, nrow(peaks))
  n_genes <- integer(nrow(peaks))

  for (i in seq_len(nrow(peaks))) {
    same <- which(g$chr == key[i])
    at_peak[i] <- ""
    # Only genes the peak actually covers are candidates. Searching the whole chromosome
    # would name a gene the peak has nothing to do with -- with a short gene list that is
    # the common case, and it contradicts `n_genes`: 0 genes in the interval must not come
    # back with one named beside it.
    within <- same[g$start[same] < peaks$end[i] & g$end[same] > peaks$start[i]]
    n_genes[i] <- length(within)
    if (!length(within)) next

    p <- peaks$peak_pos[i]
    inside <- within[g$start[within] <= p & g$end[within] > p]
    at_peak[i] <- paste(sort(unique(g$name[inside])), collapse = ",")
    # 0 bp when the peak SNP is inside the gene; otherwise the gap to the closer edge
    d <- pmax(g$start[within] - p, p - (g$end[within] - 1), 0)
    best <- which(d == min(d))
    nearest[i] <- paste(sort(unique(g$name[within][best])), collapse = ",")
    dist[i] <- min(d)
  }
  list(peak_genes = at_peak, nearest_gene = nearest, distance_to_gene = dist,
       n_genes = n_genes)
}
