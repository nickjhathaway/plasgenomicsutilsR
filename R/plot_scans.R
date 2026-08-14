# Genome-wide plots for the per-SNP and per-window scans (iHS, Rsb, XP-EHH, beta,
# windowed diversity), all sharing the chromosome layout, banding, gene track and
# auto-sizing already used by the IBD Manhattan plots.

# Any per-SNP or per-window table -> a Manhattan over the cumulative genome axis.
.scan_manhattan <- function(df, value, facet = NULL, ylab = value, title = NULL,
                            reference = DEFAULT_REFERENCE, chroms = NULL, skip_chr = NULL,
                            zoom = NULL, zoom_pad = 0.05,
                            genes = NULL, genes_for_track = NULL, gene_label_angle = 0,
                            highlight_genes = NULL, label_genes = NULL,
                            threshold = NULL, hline = NULL, point_size = 0.6,
                            point_alpha = 0.7, colours = NULL) {
  .need_package("ggplot2", "the scan plots")
  if (!nrow(df)) stop("nothing to plot", call. = FALSE)
  if (!value %in% names(df))
    stop(sprintf("no `%s` column in the scan", value), call. = FALSE)

  full <- .chrom_layout(reference)
  layout <- .select_layout(full, chroms, skip_chr)
  gtrack <- if (is.null(genes)) NULL else .gene_track(genes)
  .check_gene_request(gtrack, highlight_genes)
  z <- .zoom_setup(zoom, gtrack, layout, label_genes, zoom_pad, reference,
                   genes_for_track = genes_for_track)
  # zoomed, the gene names go in the track stacked underneath, not above the panel
  label_genes <- if (!is.null(z)) FALSE
    else if (is.null(label_genes)) !is.null(highlight_genes) else isTRUE(label_genes)

  df$chr <- normalise_chr(df$chr)
  df <- .attach_band(.recum(df, layout), layout)
  df$.y <- df[[value]]
  df <- df[is.finite(df$.y), , drop = FALSE]
  if (!nrow(df)) stop(sprintf("every `%s` value is missing", value), call. = FALSE)
  df <- .crop_to_window(df, z, "SNPs")

  gsel <- if (is.null(z)) .genes_for_layout(gtrack, layout, highlight_genes) else z$genes
  faceted <- !is.null(facet) && facet %in% names(df)
  if (faceted) df$.facet <- .as_group_factor(df[[facet]])
  top_level <- if (faceted) .first_level(df$.facet) else NULL
  genome_frac <- sum(layout$len) / sum(full$len)

  p <- ggplot2::ggplot(df, ggplot2::aes(.data$cum_pos, .data$.y)) +
    .chr_band_layer(layout) +
    .gene_line_layer(gsel) +
    (if (label_genes) .gene_label_layer(gsel, if (faceted) ".facet" else NULL, top_level,
       footprint = .label_footprint_bp(nchar(gsel$name), .fp_layout(z, layout),
                                       .fp_frac(z, genome_frac)))) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$.band), size = point_size,
                        alpha = point_alpha, stroke = 0) +
    ggplot2::scale_colour_manual(values = .band_colours(colours)) +
    .chr_axis(layout) +
    ggplot2::labs(x = "chromosome", y = ylab, title = title) +
    .manhattan_theme() +
    .gene_label_space(label_genes)

  if (!is.null(hline))
    p <- p + ggplot2::geom_hline(yintercept = hline, colour = "grey40",
                                 linetype = "dotted", linewidth = 0.3)
  if (!is.null(threshold))
    p <- p + ggplot2::geom_hline(yintercept = threshold, colour = "firebrick",
                                 linetype = "dashed", linewidth = 0.4)
  if (faceted)
    p <- p + ggplot2::facet_wrap(~ .data$.facet, ncol = 1, strip.position = "right")

  n_panels <- if (faceted) nlevels(df$.facet) else 1L
  p <- .apply_zoom(p, z, layout, n_panels, gene_label_angle)
  attr(p, "plasgenomics_dims") <- if (is.null(z))
    .dims_genome(n_panels, genome_frac, label_genes = label_genes)
    else .dims_zoom(n_panels, attr(p, "plasgenomics_track_in"))
  p
}

# Accept any gene table and hand back the columns the gene track layers expect.
.gene_track <- function(genes) {
  g <- as.data.frame(genes)
  if (!"chr" %in% names(g)) {
    for (alt in c("chrom", "Pf3D7_chrom")) if (alt %in% names(g)) { g$chr <- g[[alt]]; break }
  }
  need <- c("name", "chr", "start", "end")
  if (!all(need %in% names(g)))
    stop("`genes` needs name, chr (or chrom), start and end columns", call. = FALSE)
  g$chr <- normalise_chr(g$chr)
  g[need]
}

#' Manhattan plot of a haplotype-homozygosity scan
#'
#' Draws `neg_log10_p` (or the raw statistic) along the genome, one panel per group or
#' population pair, over the chromosome bands and optional gene markers used by the other
#' genome-wide plots in the package.
#'
#' @param scan The tibble from [run_ihs()], [run_rsb()] or [run_xpehh()].
#' @param metric `"neg_log10_p"` (default) or the statistic itself (`"ihs"` / `"value"`).
#' @param threshold Draw a dashed significance line at this height. `NULL` uses the
#'   `-log10(p)` of a 1% two-sided tail when plotting `neg_log10_p`, the convention these
#'   scans are usually read at; `NA` draws none.
#' @param genes Gene table to mark (e.g. [PF_EXAMPLE_DRUG_GENES]); `NULL` for none.
#' @param highlight_genes,label_genes Which genes to mark and whether to name them.
#' @param chroms,skip_chr Chromosomes to keep or drop.
#' @param zoom Optional single interval to crop to, keeping the same data and the same
#'   coordinates as the genome-wide plot: a chromosome (`"7"`), a range
#'   (`"7:728,081-988,719"`), a gene name from `genes`, or a one-row data frame with
#'   chr/start/end. Every gene in the window is drawn and named unless
#'   `label_genes = FALSE`.
#' @param zoom_pad Context to add around `zoom`, clamped to the chromosome. One value pads
#'   both sides (default 5%); two pad the left and the right, either in that order or named
#'   -- `c(left = 5000, right = 40000)`, and naming only one side pads only that side. Each
#'   side is read on its own: below 1 it is a fraction of the interval's span, at or above 1
#'   it is base pairs, so `c(0.1, 20000)` is legal.
#' @param genes_for_track Optional gene table for the track drawn beneath a zoomed plot
#'   (e.g. [PF3D7_GENES]), so every gene in the window is shown and named while the plot's
#'   own short track still supplies the marked positions inside the panel. Without it the
#'   track and the marks come from the same genes, which means marking a whole annotation
#'   just to see the neighbours.
#' @param gene_label_angle Rotation for the gene names in that track, in degrees. `0`
#'   (default) centres each name under its gene; `45` or `90` runs it down to the left,
#'   which is what keeps long systematic ids from colliding over a dense annotation.
#' @param reference Reference id for the chromosome layout.
#' @param point_size,point_alpha,colours,colors Point and band aesthetics.
#' @return A ggplot object.
#' @examples
#' ps <- example_pop_structure(umap = FALSE)
#' hap <- parasite_haplotypes(ps, maf = 0.05)
#' plot_ihs(run_ihs(hap, group = "country"), genes = PF_EXAMPLE_DRUG_GENES)
#' @export
plot_ihs <- function(scan, metric = c("neg_log10_p", "ihs", "value"), threshold = NULL,
                     genes = NULL, highlight_genes = NULL, label_genes = NULL,
                     chroms = NULL, skip_chr = NULL, zoom = NULL, zoom_pad = 0.05,
                     genes_for_track = NULL, gene_label_angle = 0,
                     reference = DEFAULT_REFERENCE,
                     point_size = 0.6, point_alpha = 0.7, colours = NULL,
                     colors = NULL) {
  colours <- .alias_arg("colours", "colors")
  metric <- match.arg(metric)
  if (metric == "ihs" && !"ihs" %in% names(scan)) metric <- "value"
  facet <- if ("group" %in% names(scan)) "group" else if ("pair" %in% names(scan)) "pair"
  lab <- switch(metric, neg_log10_p = expression(-log[10](italic(p))), metric)
  thr <- if (metric == "neg_log10_p") {
    if (is.null(threshold)) -log10(0.01) else threshold
  } else threshold
  .scan_manhattan(scan, metric, facet = facet, ylab = lab,
                  reference = reference, chroms = chroms, skip_chr = skip_chr,
                  zoom = zoom, zoom_pad = zoom_pad,
                  genes = genes, genes_for_track = genes_for_track,
                  gene_label_angle = gene_label_angle,
                  highlight_genes = highlight_genes, label_genes = label_genes,
                  threshold = if (is.null(thr) || is.na(thr)) NULL else thr,
                  hline = if (metric != "neg_log10_p") 0 else NULL,
                  point_size = point_size, point_alpha = point_alpha, colours = colours)
}

#' Manhattan plot of beta scores
#'
#' @param b The tibble from [beta_score()].
#' @inheritParams plot_ihs
#' @param threshold Draw a dashed line at this beta; `NULL` uses the 99th percentile of
#'   the scores actually plotted, the empirical tail these scans are read at.
#' @return A ggplot object.
#' @examples
#' ps <- example_pop_structure(umap = FALSE)
#' b <- beta_score(ps, group = "country", window = 300000, min_window_snps = 1)
#' plot_beta(b)
#' @export
plot_beta <- function(b, threshold = NULL, genes = NULL, highlight_genes = NULL,
                      label_genes = NULL, chroms = NULL, skip_chr = NULL,
                      zoom = NULL, zoom_pad = 0.05,
                      genes_for_track = NULL, gene_label_angle = 0,
                      reference = DEFAULT_REFERENCE, point_size = 0.6,
                      point_alpha = 0.7, colours = NULL,
                      colors = NULL) {
  colours <- .alias_arg("colours", "colors")
  thr <- if (is.null(threshold)) {
    v <- b$beta[is.finite(b$beta)]
    if (length(v) > 20) stats::quantile(v, 0.99, names = FALSE) else NULL
  } else if (is.na(threshold)) NULL else threshold
  .scan_manhattan(b, "beta", facet = if ("group" %in% names(b)) "group" else NULL,
                  ylab = expression(beta), reference = reference, chroms = chroms,
                  skip_chr = skip_chr, zoom = zoom, zoom_pad = zoom_pad,
                  genes = genes, genes_for_track = genes_for_track,
                  gene_label_angle = gene_label_angle,
                  highlight_genes = highlight_genes,
                  label_genes = label_genes, threshold = thr, hline = 0,
                  point_size = point_size, point_alpha = point_alpha, colours = colours)
}

#' Windowed diversity along the genome
#'
#' One panel per group, with the statistic plotted at the window mid-points.
#'
#' @param div A [pop_diversity()] result computed with `by = "window"`.
#' @param metric Which column to draw: `"pi"` (default), `"he"`, `"theta_w"`,
#'   `"tajima_d"`, `"hap_div"`, `"seg_sites"`, ...
#' @inheritParams plot_ihs
#' @return A ggplot object.
#' @examples
#' ps <- example_pop_structure(umap = FALSE)
#' d <- pop_diversity(ps, group = "country", by = "window", window = 500000)
#' plot_diversity(d, metric = "he")
#' @export
plot_diversity <- function(div, metric = "pi", genes = NULL, highlight_genes = NULL,
                           label_genes = NULL, chroms = NULL, skip_chr = NULL,
                           zoom = NULL, zoom_pad = 0.05,
                           genes_for_track = NULL, gene_label_angle = 0,
                           reference = DEFAULT_REFERENCE, point_size = 0.6,
                           point_alpha = 0.7, colours = NULL,
                           colors = NULL) {
  colours <- .alias_arg("colours", "colors")
  if (!"start" %in% names(div) || all(is.na(div$start)))
    stop("plot_diversity() needs a windowed result: pop_diversity(by = \"window\")",
         call. = FALSE)
  df <- as.data.frame(div)
  df$pos <- (df$start + df$end) / 2
  .scan_manhattan(df, metric, facet = if ("group" %in% names(df)) "group" else NULL,
                  ylab = metric, reference = reference, chroms = chroms,
                  skip_chr = skip_chr, zoom = zoom, zoom_pad = zoom_pad,
                  genes = genes, genes_for_track = genes_for_track,
                  gene_label_angle = gene_label_angle,
                  highlight_genes = highlight_genes, label_genes = label_genes,
                  hline = if (metric == "tajima_d") 0 else NULL,
                  point_size = point_size, point_alpha = point_alpha, colours = colours)
}

#' Linkage-disequilibrium decay curve
#'
#' Mean r-squared against the distance between SNP pairs, one line per group, with the
#' half-decay distance marked.
#'
#' @param ld The tibble from [read_ld_decay()] (or `plasgenomicsutils ld_decay`).
#' @param show_half_decay Mark each group's half-decay distance with a vertical segment.
#' @param colours,colors Optional named colour vector for the groups.
#' @return A ggplot object.
#' @examples
#' # the curve itself comes from `plasgenomicsutils ld_decay`; read_ld_decay() loads it
#' ld <- data.frame(group = rep(c("a", "b"), each = 4),
#'                  bin_mid = rep(c(2500, 7500, 12500, 17500), 2),
#'                  n_pairs = 100,
#'                  mean_r2 = c(0.40, 0.22, 0.14, 0.11, 0.25, 0.24, 0.23, 0.22))
#' ld$group <- factor(ld$group)
#' attr(ld, "ld_half_decay") <- data.frame(group = c("a", "b"),
#'                                         half_decay_bp = c(7100, NA))
#' plot_ld_decay(ld)
#' @export
plot_ld_decay <- function(ld, show_half_decay = TRUE, colours = NULL,
                          colors = NULL) {
  colours <- .alias_arg("colours", "colors")
  .need_package("ggplot2", "plot_ld_decay()")
  if (!nrow(ld)) stop("nothing to plot", call. = FALSE)
  df <- as.data.frame(ld)
  df <- df[is.finite(df$mean_r2), , drop = FALSE]
  p <- ggplot2::ggplot(df, ggplot2::aes(.data$bin_mid, .data$mean_r2,
                                        colour = .data$group)) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::labs(x = "distance between SNPs (bp)", y = expression(mean~r^2),
                  colour = NULL) +
    ggplot2::theme_bw(base_size = 10)
  if (!is.null(colours)) p <- p + ggplot2::scale_colour_manual(values = colours)

  hd <- attr(ld, "ld_half_decay")
  if (show_half_decay && !is.null(hd)) {
    hd <- hd[is.finite(hd$half_decay_bp), , drop = FALSE]
    if (nrow(hd)) {
      hd$group <- factor(hd$group, levels = levels(df$group))
      p <- p + ggplot2::geom_vline(data = hd,
        ggplot2::aes(xintercept = .data$half_decay_bp, colour = .data$group),
        linetype = "dashed", linewidth = 0.4, show.legend = FALSE)
    }
  }
  attr(p, "plasgenomics_dims") <- c(width = 6, height = 4)
  p
}
