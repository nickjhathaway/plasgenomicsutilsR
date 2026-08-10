# A single locus in detail: IBD sharing and a selection scan on one pair of axes.

# Candidate columns for the optional point aesthetics, tried in order. Nothing is mapped
# unless one of them is actually in the scan, so a table without them plots plain points.
.MUTATION_COLS <- c("mutation_type", "mutation_class", "effect", "consequence")

# Resolve an optional aesthetic: NULL auto-detects from `candidates`, NA/FALSE disables it,
# a name is taken as given (and must exist).
.optional_aes <- function(arg, df, candidates, what) {
  if (isFALSE(arg) || (length(arg) == 1 && is.na(arg))) return(NULL)
  if (is.null(arg)) {
    hit <- intersect(candidates, names(df))
    return(if (length(hit)) hit[1] else NULL)
  }
  if (!is.character(arg) || length(arg) != 1)
    stop("`", what, "` must be a column name, NULL to detect one, or NA for none",
         call. = FALSE)
  if (!arg %in% names(df))
    stop("`", what, " = \"", arg, "\"` is not a column of the scan", call. = FALSE)
  arg
}

#' One locus in detail: IBD sharing against a selection scan
#'
#' A zoomed panel over a single interval with two vertical axes: the fraction of pairs
#' sharing each SNP by IBD as a step curve on the left, and a per-SNP selection scan as
#' points on the right, with a gene track underneath. Reading both signals against one
#' another over a few hundred kilobases is what separates a shared haplotype that is also
#' under selection from one that is merely common.
#'
#' The two tracks keep their own units. The scan is drawn on the secondary axis, so a point
#' and a curve at the same height are not the same number: compare shapes and positions,
#' not heights.
#'
#' @param x An [IbdResults] object (needs a per_snp_group table).
#' @param locus The interval to draw: a gene name from the object's track, a range
#'   (`"7:380,000-430,000"`), a whole chromosome (`"7"`), or a one-row data frame with
#'   chr/start/end.
#' @param scan Optional per-SNP scan for the right axis: a [run_ihs()] result, a
#'   [beta_score()] table, or any table with chr, pos and `metric`. `NULL` (default) uses
#'   the object's own IBD selection statistic.
#' @param metric Column of `scan` to draw (default `"neg_log10_p"`).
#' @param groups Optional groups to keep; `NULL` draws every group in the object.
#' @param pad Context to add around `locus`, clamped to the chromosome. One value pads both
#'   sides (default 10%); two pad the left and the right, either in that order or named --
#'   `c(left = 5000, right = 40000)`, and naming only one side pads only that side. Each side
#'   is read on its own: below 1 it is a fraction of the interval's span, at or above 1 it is
#'   base pairs, so `c(0.1, 20000)` is legal.
#' @param min_span Widen the window to at least this many base pairs (default 50 kb). A
#'   single gene padded by a fraction of its own length is usually too narrow to contain
#'   more than a handful of SNPs.
#' @param threshold Line(s) to draw on the scan axis: a height, or the name of a threshold
#'   the selection run wrote -- `"bonferroni"`, `"fdr"`, `"permutation"`, `"empirical"`,
#'   `"both"` or `"all"`, the same vocabulary and colours [plot_selection_manhattan()] uses,
#'   and per group where the run wrote one per group. `NULL` (default) draws the object's
#'   Bonferroni line, or the 1% tail when the right axis is an external `scan`; `NA` draws
#'   none. Naming a kind the run did not write is an error, and naming one at all needs the
#'   object's own statistic, since the threshold table says nothing about an external scan.
#' @param scan_abs Plot the magnitude of a signed statistic. `NULL` (default) takes the
#'   absolute value when the metric actually has negative values -- iHS, Rsb, a z-score and
#'   beta are all read by distance from zero, and both tails mean selection -- and labels the
#'   axis `|metric|`. `FALSE` keeps the sign and extends the panel below zero instead, which
#'   is readable but means a point level with the IBD curve's zero is a value of zero rather
#'   than the floor of the axis.
#' @param size_by Optional scan column to map to point size (e.g. a per-SNP allele
#'   frequency difference). `NULL` (default) draws every point the same size.
#' @param shape_by Optional scan column to map to point shape. `NULL` (default) uses a
#'   mutation-class column (`mutation_type`, `mutation_class`, `effect`, `consequence`) if
#'   the scan has one, and plain points otherwise; `NA` never maps shape.
#' @param point_size Size of the scan points when `size_by` is not used.
#' @param ibd_colour,scan_colour Colours for the IBD curve and the scan points.
#' @param gene_track Draw the gene track under the panel (default `TRUE`).
#' @param genes_for_track Optional gene table for the track drawn under the panel
#'   (e.g. [PF3D7_GENES]), so every gene in the window is shown and named while the object's
#'   own track still supplies the marked positions inside the panel. Without it the
#'   track and the marks come from the same genes, which means marking a whole annotation
#'   just to see the neighbours.
#' @param gene_label_angle Rotation for the gene names in that track, in degrees. `0`
#'   (default) centres each name under its gene; `45` or `90` runs it down to the left,
#'   which is what keeps long systematic ids from colliding over a dense annotation.
#' @return A ggplot object, or a patchwork of the panel over the gene track when
#'   `gene_track = TRUE` and the window contains genes.
#' @examples
#' x <- example_ibd_results()
#' plot_ibd_locus(x, "pfcrt", pad = 50000)
#' @export
plot_ibd_locus <- function(x, locus, scan = NULL, metric = "neg_log10_p", groups = NULL,
                           pad = 0.1, min_span = 50000, threshold = NULL, scan_abs = NULL,
                           size_by = NULL, shape_by = NULL,
                           point_size = 1.4, ibd_colour = "#2166ac",
                           scan_colour = "#d95f02", gene_track = TRUE,
                           genes_for_track = NULL, gene_label_angle = 0) {
  .need_package("ggplot2", "plot_ibd_locus()")
  .need_package("scales", "plot_ibd_locus()")
  ibd <- x$get_per_snp_group()
  if (is.null(ibd)) stop("this IbdResults has no per_snp_group table", call. = FALSE)
  tt <- .top_track(x, scan, metric)
  sel <- tt$df
  if (!metric %in% names(sel))
    stop(sprintf("the scan has no '%s' column", metric), call. = FALSE)

  layout <- x$chrom_layout()
  z <- .zoom_setup(locus, x$get_genes(), layout, label_genes = gene_track, pad = pad,
                   reference = x$reference_id(), min_span = min_span,
                   genes_for_track = genes_for_track)
  if (is.null(z)) stop("`locus` is required", call. = FALSE)

  ibd <- .filter_group(ibd, groups)
  if ("group" %in% names(sel) && !is.null(groups))
    sel <- sel[as.character(sel$group) %in% as.character(groups), , drop = FALSE]
  ibd <- .in_window(.recum(ibd, layout), z)
  sel <- .in_window(.recum(sel, layout), z)
  if (!nrow(ibd))
    stop(sprintf("no IBD SNPs in %s; widen the window with `pad` / `min_span`",
                 .window_label(z)), call. = FALSE)
  sel <- sel[is.finite(sel[[metric]]), , drop = FALSE]

  # A signed statistic (iHS, Rsb, a z-score, beta) is read by how far it is from zero, and
  # both tails mean selection. The IBD track starts at 0, so anchoring the scan there too
  # would put every negative value off the bottom of the panel: take the magnitude instead,
  # and say so on the axis. `scan_abs = FALSE` keeps the sign and lets the axis go negative --
  # readable, but then a point level with the IBD curve's zero is a value of zero, not a floor.
  signed <- nrow(sel) > 0 && any(sel[[metric]] < 0, na.rm = TRUE)
  if (is.null(scan_abs)) scan_abs <- signed
  scan_lab <- metric
  if (isTRUE(scan_abs) && signed) {
    sel[[metric]] <- abs(sel[[metric]])
    scan_lab <- paste0("|", metric, "|")
  }

  # The scan rides a secondary axis, so it is drawn in IBD units and read back through the
  # inverse. Both tracks share zero so the axes agree on where "nothing" is.
  ibd_max <- max(ibd$frac_pairs_ibd, na.rm = TRUE)
  if (!is.finite(ibd_max) || ibd_max <= 0) ibd_max <- 1
  sel_max <- if (nrow(sel)) max(abs(sel[[metric]]), na.rm = TRUE) else NA_real_
  if (!is.finite(sel_max) || sel_max <= 0) sel_max <- 1
  ratio <- ibd_max / sel_max
  if (nrow(sel)) sel$.y <- sel[[metric]] * ratio

  faceted <- "group" %in% names(ibd) && length(unique(ibd$group)) > 1
  if (faceted) ibd$group <- .as_group_factor(ibd$group)
  # a scan without groups is the same evidence for every panel, so repeat it in each
  if (faceted && nrow(sel)) {
    sel <- if ("group" %in% names(sel)) {
      sel$group <- factor(as.character(sel$group), levels = levels(ibd$group))
      sel[!is.na(sel$group), , drop = FALSE]
    } else {
      do.call(rbind, lapply(levels(ibd$group), function(g) transform(sel, group = g)))
    }
  }

  size_col <- .optional_aes(size_by, sel, character(0), "size_by")
  shape_col <- .optional_aes(shape_by, sel, .MUTATION_COLS, "shape_by")
  pt_aes <- ggplot2::aes(x = .data$cum_pos, y = .data$.y)
  if (!is.null(size_col)) pt_aes$size <- ggplot2::aes(size = .data[[size_col]])$size
  if (!is.null(shape_col)) pt_aes$shape <- ggplot2::aes(shape = .data[[shape_col]])$shape

  thr <- .locus_thresholds(x, tt, threshold, metric, ratio, ibd, faceted, scan_colour)
  # A threshold above everything in the window is the informative answer ("nothing here is
  # significant"), so the scale is extended to hold it. Left out, ggplot silently drops the
  # line for being outside the limits and the plot looks as if none was asked for.
  y_top <- max(ibd_max, if (nrow(sel)) max(sel$.y, na.rm = TRUE) else 0,
               if (length(thr$values)) max(thr$values, na.rm = TRUE) * ratio else 0)
  # keeping the sign means the panel has to reach below zero, or the negative half is clipped
  y_bottom <- min(0, if (nrow(sel)) min(sel$.y, na.rm = TRUE) else 0)

  # `size` has to be left out of the call, not passed as NULL: a NULL parameter still
  # counts as set and ggplot then drops the mapped aesthetic ("Ignoring empty aesthetic")
  pt_args <- list(data = sel, mapping = pt_aes, colour = scan_colour, alpha = 0.8)
  if (is.null(size_col)) pt_args$size <- point_size

  p <- ggplot2::ggplot() +
    ggplot2::geom_step(data = ibd, ggplot2::aes(.data$cum_pos, .data$frac_pairs_ibd),
                       colour = ibd_colour, linewidth = 0.4, direction = "mid") +
    (if (nrow(sel)) do.call(ggplot2::geom_point, pt_args)) +
    thr$layers +
    .gene_line_layer(z$genes) +
    ggplot2::scale_y_continuous(
      name = "pairs IBD", labels = scales::percent_format(accuracy = 1),
      limits = c(y_bottom, y_top),
      sec.axis = ggplot2::sec_axis(~ . / ratio,
                                   name = paste0(tt$label, " (", scan_lab, ")"))) +
    .region_axis(z$interval, z$offset) +
    ggplot2::coord_cartesian(xlim = z$xlim) +
    .manhattan_theme() +
    ggplot2::theme(axis.title.y.left = ggplot2::element_text(colour = ibd_colour),
                   axis.text.y.left = ggplot2::element_text(colour = ibd_colour),
                   axis.title.y.right = ggplot2::element_text(colour = scan_colour),
                   axis.text.y.right = ggplot2::element_text(colour = scan_colour))
  if (faceted)
    p <- p + ggplot2::facet_wrap(~ .data$group, ncol = 1, strip.position = "right")

  n_panels <- if (faceted) nlevels(ibd$group) else 1L
  if (gene_track) p <- .stack_gene_track(p, z, layout, n_panels, gene_label_angle)
  attr(p, "plasgenomics_dims") <- .dims_zoom(n_panels, attr(p, "plasgenomics_track_in"))
  p
}

# Threshold line(s) for the scan axis. `threshold` takes a height, or the name of a kind the
# selection run wrote -- the same vocabulary and the same colours [plot_selection_manhattan()]
# uses, resolved by the same helpers so the two plots cannot drift apart -- or NA for none.
# Every line is mapped through the same transform as the points, since the scan is drawn in
# IBD units on the secondary axis.
.locus_thresholds <- function(x, tt, threshold, metric, ratio, ibd, faceted,
                              default_colour) {
  vals <- numeric(0)
  flat <- function(v, colour, linetype = "dashed") {
    if (!length(v) || !is.finite(v[1])) return(NULL)
    vals <<- c(vals, v[1])
    ggplot2::geom_hline(yintercept = v[1] * ratio, colour = colour, linetype = linetype,
                        linewidth = 0.4)
  }
  done <- function(layers) list(layers = layers, values = vals)
  if (length(threshold) == 1 && is.na(threshold)) return(done(NULL))

  if (is.numeric(threshold)) return(done(flat(threshold, default_colour)))

  if (is.character(threshold)) {
    if (!tt$own)
      stop("`threshold = \"", threshold[1], "\"` names a kind in the object's own threshold ",
           "table, which says nothing about an external `scan`; pass a height instead",
           call. = FALSE)
    kinds <- .threshold_kinds(.resolve_threshold(threshold))
    out <- lapply(kinds, function(kind) {
      line <- .threshold_line(x$get_thresholds(), kind, strict = length(kinds) == 1)
      if (is.null(line)) return(NULL)
      sty <- .THRESHOLD_STYLE[[kind]]
      # a per-group threshold gets a line per panel; keep only the panels being drawn
      has_g <- faceted && "group" %in% names(line) && any(!is.na(line$group))
      if (!has_g) return(flat(line$threshold, sty$colour, sty$linetype))
      line$group <- factor(as.character(line$group), levels = levels(ibd$group))
      line <- line[!is.na(line$group), , drop = FALSE]
      if (!nrow(line)) return(NULL)
      line$.y <- line$threshold * ratio
      vals <<- c(vals, line$threshold)
      ggplot2::geom_hline(data = line, ggplot2::aes(yintercept = .data$.y),
                          inherit.aes = FALSE, colour = sty$colour,
                          linetype = sty$linetype, linewidth = 0.4)
    })
    out <- out[!vapply(out, is.null, logical(1))]
    return(done(if (length(out)) out else NULL))
  }

  # the default: the object's own Bonferroni line, or the tail an external scan is read at
  if (metric != "neg_log10_p") return(done(NULL))
  if (!tt$own) return(done(flat(-log10(0.01), default_colour)))
  line <- .threshold_line(x$get_thresholds(), "bonferroni", strict = FALSE)
  if (is.null(line)) return(done(NULL))
  done(flat(max(line$threshold), .THRESHOLD_STYLE$bonferroni$colour,
            .THRESHOLD_STYLE$bonferroni$linetype))
}
