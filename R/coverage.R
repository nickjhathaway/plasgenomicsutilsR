# Reading and plotting the sequencing-depth tables written by
# `plasgenomicsutils coverage_depth_stats` and `coverage_dropout_regions`.
#
# The compute side lives in the Python package because it reads BAMs; everything here
# takes its TSV output and turns it into the QC calls and figures.

#: Depth a base must reach to count towards the headline breadth column.
COVERAGE_QC_THRESHOLD <- 10L
#: Default floors for a sample to pass QC: mean depth, and percent of bases at threshold.
COVERAGE_MIN_MEAN <- 5
COVERAGE_MIN_BREADTH <- 80

#' Read a coverage table
#'
#' Identifier columns are read as text. Sample names are very often bare digits -- a
#' sequencing id, or a BAM named after one -- and guessing a type would turn those into
#' doubles, printing `4089106922` as `4.09e+09` and silently failing to join against the
#' metadata.
#'
#' @param path TSV(.gz) written by `plasgenomicsutils coverage_depth_stats`
#'   (per-sample/per-chromosome), or by `coverage_dropout_regions`.
#' @return A tibble.
#' @examples
#' f <- tempfile(fileext = ".tsv")
#' write.table(data.frame(sample = c("s1", "s2"), chrom = "genome",
#'                        mean = c(48, 4), pct_ge_10x = c(95, 12)),
#'             f, sep = "\t", quote = FALSE, row.names = FALSE)
#' read_coverage(f)
#' @export
read_coverage <- function(path) {
  .need_package("readr", "read_coverage()")
  head <- readr::read_tsv(path, n_max = 0, show_col_types = FALSE, progress = FALSE)
  as_text <- intersect(c("sample", "chrom", "engine", "genes", "top_partner", "flag"),
                       names(head))
  spec <- stats::setNames(rep(list(readr::col_character()), length(as_text)), as_text)
  readr::read_tsv(path, show_col_types = FALSE, progress = FALSE,
                  col_types = do.call(readr::cols, spec))
}

#' Per-sample coverage QC verdict
#'
#' Reduces a coverage table to one row per sample -- the genome-wide row -- and applies
#' the two floors a sample has to clear: enough average depth, and enough of the genome
#' actually reaching a usable depth. The second matters more: selective whole-genome
#' amplification can give a respectable mean while leaving much of the genome at zero,
#' and only the breadth column shows that.
#'
#' @param cov A coverage table from [read_coverage()].
#' @param threshold Depth whose breadth column is used (default `r COVERAGE_QC_THRESHOLD`;
#'   the table must have been produced with that threshold).
#' @param min_mean,min_breadth Floors for mean depth and percent of bases at `threshold`.
#' @return A tibble with one row per sample: `sample`, `mean`, `median`, `sd`,
#'   `pct_ge_<threshold>x`, `pct_zero`, `pass`, and `fail_reason`.
#' @examples
#' cov <- data.frame(sample = c("a", "b"), chrom = "genome", mean = c(45, 4),
#'                   median = c(44, 0), sd = c(9, 7), pct_zero = c(1, 62),
#'                   pct_ge_10x = c(96, 30))
#' coverage_qc(cov)
#' @export
coverage_qc <- function(cov, threshold = COVERAGE_QC_THRESHOLD,
                        min_mean = COVERAGE_MIN_MEAN,
                        min_breadth = COVERAGE_MIN_BREADTH) {
  col <- paste0("pct_ge_", threshold, "x")
  if (!col %in% names(cov))
    stop(sprintf("no `%s` column -- rerun coverage_depth_stats with --thresholds %d",
                 col, threshold), call. = FALSE)
  g <- as.data.frame(cov[cov$chrom == "genome", , drop = FALSE])
  if (!nrow(g)) stop("no `genome` rows in the coverage table", call. = FALSE)

  low_mean <- g$mean < min_mean
  low_breadth <- g[[col]] < min_breadth
  reason <- ifelse(low_mean & low_breadth, "low depth and breadth",
            ifelse(low_mean, "low depth",
            ifelse(low_breadth, "low breadth", NA_character_)))
  keep <- intersect(c("sample", "mean", "median", "sd", "pct_zero", col), names(g))
  out <- g[keep]
  out$pass <- !(low_mean | low_breadth)
  out$fail_reason <- reason
  tibble::as_tibble(out[order(out$mean), , drop = FALSE])
}

#' Per-sample coverage overview
#'
#' Mean depth against the fraction of the genome reaching a usable depth, one point per
#' sample, with the QC floors drawn. The failures separate along whichever axis they
#' failed on, which is the quickest way to tell a shallow run from an uneven one.
#'
#' @inheritParams coverage_qc
#' @param label_failures Name the samples that fail.
#' @return A ggplot object.
#' @examples
#' cov <- data.frame(sample = paste0("s", 1:4), chrom = "genome",
#'                   mean = c(60, 45, 30, 3), pct_ge_10x = c(96, 91, 41, 7))
#' plot_coverage_summary(cov)
#' @export
plot_coverage_summary <- function(cov, threshold = COVERAGE_QC_THRESHOLD,
                                  min_mean = COVERAGE_MIN_MEAN,
                                  min_breadth = COVERAGE_MIN_BREADTH,
                                  label_failures = TRUE) {
  .need_package("ggplot2", "plot_coverage_summary()")
  col <- paste0("pct_ge_", threshold, "x")
  qc <- as.data.frame(coverage_qc(cov, threshold, min_mean, min_breadth))
  qc$.breadth <- qc[[col]]

  p <- ggplot2::ggplot(qc, ggplot2::aes(.data$mean, .data$.breadth)) +
    ggplot2::geom_hline(yintercept = min_breadth, colour = "firebrick",
                        linetype = "dashed", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = min_mean, colour = "firebrick",
                        linetype = "dashed", linewidth = 0.4) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$pass), size = 1.8, alpha = 0.85) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "grey30", `FALSE` = "firebrick"),
                                 name = "passes QC") +
    ggplot2::scale_x_continuous(trans = "log10") +
    ggplot2::labs(x = "mean depth (log scale)",
                  y = sprintf("%% of bases at >= %dx", threshold)) +
    ggplot2::theme_bw(base_size = 10)
  if (label_failures && any(!qc$pass)) {
    .need_package("ggplot2", "plot_coverage_summary()")
    p <- p + ggplot2::geom_text(data = qc[!qc$pass, , drop = FALSE],
                                ggplot2::aes(label = .data$sample),
                                size = 2.4, vjust = -0.9, colour = "firebrick")
  }
  attr(p, "plasgenomics_dims") <- c(width = 6, height = 4.5)
  p
}

#' Coverage per chromosome, per sample
#'
#' A sample-by-chromosome tile of mean depth (or any other column in the table). A
#' chromosome that is systematically low across samples is a reference or amplification
#' problem, not a sample problem -- the two read very differently here.
#'
#' @param cov A coverage table from [read_coverage()].
#' @param metric Column to fill by (default `"mean"`).
#' @param relative Divide each sample's value by its own genome-wide value, so the tile
#'   shows relative rather than absolute depth and deep samples do not dominate.
#' @param sample_order Optional sample order; defaults to increasing genome-wide mean.
#' @return A ggplot object.
#' @examples
#' cov <- expand.grid(sample = paste0("s", 1:4),
#'                    chrom = sprintf("Pf3D7_%02d_v3", 1:3), stringsAsFactors = FALSE)
#' cov$mean <- c(40, 38, 12, 44, 41, 39, 15, 43, 42, 40, 13, 45)
#' cov <- rbind(cov, data.frame(sample = paste0("s", 1:4), chrom = "genome",
#'                              mean = c(41, 40, 13, 44)))
#' plot_coverage_by_chrom(cov)
#' @export
plot_coverage_by_chrom <- function(cov, metric = "mean", relative = TRUE,
                                   sample_order = NULL) {
  .need_package("ggplot2", "plot_coverage_by_chrom()")
  df <- as.data.frame(cov)
  if (!metric %in% names(df))
    stop(sprintf("no `%s` column in the coverage table", metric), call. = FALSE)
  genome <- df[df$chrom == "genome", , drop = FALSE]
  df <- df[df$chrom != "genome", , drop = FALSE]
  if (!nrow(df)) stop("the table has only genome-wide rows", call. = FALSE)

  df$.v <- df[[metric]]
  lab <- metric
  if (relative) {
    ref <- stats::setNames(genome[[metric]], genome$sample)
    df$.v <- df$.v / ref[df$sample]
    lab <- paste(metric, "/ genome-wide")
  }
  ord <- if (!is.null(sample_order)) sample_order else
    genome$sample[order(genome[[metric]])]
  df$sample <- factor(df$sample, levels = ord)
  df$chrom <- factor(normalise_chr(df$chrom),
                     levels = .levels_of(normalise_chr(df$chrom)))

  p <- ggplot2::ggplot(df, ggplot2::aes(.data$chrom, .data$sample, fill = .data$.v)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(name = lab, option = "magma") +
    ggplot2::labs(x = "chromosome", y = NULL) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 5),
                   panel.grid = ggplot2::element_blank())
  attr(p, "plasgenomics_dims") <- c(
    width = 6, height = round(max(3, 0.09 * length(unique(df$sample)) + 1.2), 1))
  p
}

#' Coverage dropouts along the genome
#'
#' The fraction of samples below depth in each window, across the genome, with the merged
#' dropout regions marked. Regions where the line sits near 1 are amplified in almost
#' nobody -- they will read as invariant rather than as missing unless they are excluded.
#'
#' @param windows Per-window table from `coverage_depth_stats --windows-output`, or the
#'   already-merged regions from `coverage_dropout_regions`.
#' @param min_depth A sample counts as uncovered in a window below this mean depth
#'   (ignored when `windows` is already merged).
#' @param min_frac_samples Draw the flag line at this fraction.
#' @param genes,highlight_genes,label_genes Optional gene track, as in [plot_ihs()].
#' @param chroms,skip_chr,reference As in [plot_ihs()].
#' @return A ggplot object.
#' @examples
#' win <- expand.grid(sample = paste0("s", 1:5), chrom = "Pf3D7_01_v3",
#'                    start = seq(0, 3e5, by = 1e5), stringsAsFactors = FALSE)
#' win$end <- win$start + 1e5
#' win$mean_depth <- c(rep(30, 15), rep(1, 5))    # the last window is empty in everyone
#' plot_coverage_dropout(win)
#' @export
plot_coverage_dropout <- function(windows, min_depth = 5, min_frac_samples = 0.9,
                                  genes = NULL, highlight_genes = NULL,
                                  label_genes = NULL, chroms = NULL, skip_chr = NULL,
                                  reference = DEFAULT_REFERENCE) {
  .need_package("ggplot2", "plot_coverage_dropout()")
  df <- as.data.frame(windows)
  if ("frac_samples_uncovered" %in% names(df)) {
    agg <- df
  } else {
    need <- c("sample", "chrom", "start", "end", "mean_depth")
    if (!all(need %in% names(df)))
      stop("`windows` needs sample, chrom, start, end and mean_depth columns",
           call. = FALSE)
    df$.un <- df$mean_depth < min_depth
    agg <- stats::aggregate(list(frac_samples_uncovered = df$.un),
                            by = list(chrom = df$chrom, start = df$start, end = df$end),
                            FUN = mean)
  }
  agg$chr <- normalise_chr(agg$chrom)
  agg$pos <- (agg$start + agg$end) / 2

  p <- .scan_manhattan(agg, "frac_samples_uncovered", ylab = "fraction of samples below depth",
                       reference = reference, chroms = chroms, skip_chr = skip_chr,
                       genes = genes, highlight_genes = highlight_genes,
                       label_genes = label_genes, threshold = min_frac_samples,
                       point_size = 0.4, point_alpha = 0.5)
  p + ggplot2::coord_cartesian(ylim = c(0, 1))
}
