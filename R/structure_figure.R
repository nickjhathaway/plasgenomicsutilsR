# Composite population-structure figure: a UMAP scatter combined in one plot with a
# region-faceted sNMF admixture, sharing one theme (so fonts match), one region colour
# map (so UMAP points and admixture strips match), and collected legends.

# a single region's admixture panel: stacked bars + a colour strip on top (no text)
.admix_panel <- function(q_g, samples_g, fill_vals, header_col, base_size,
                         border = TRUE, border_colour = "black", border_linewidth = 0.15) {
  long <- data.frame(
    sample  = factor(rep(samples_g, times = ncol(q_g)), levels = samples_g),
    cluster = factor(rep(colnames(q_g), each = length(samples_g)), levels = colnames(q_g)),
    q       = as.vector(q_g),
    stringsAsFactors = FALSE)
  ggplot2::ggplot(long, ggplot2::aes(.data$sample, .data$q, fill = .data$cluster)) +
    ggplot2::geom_col(width = 1, colour = if (border) border_colour else NA,
                      linewidth = border_linewidth) +
    ggplot2::annotate("rect", xmin = 0.5, xmax = length(samples_g) + 0.5,
                      ymin = 1.02, ymax = 1.13, fill = header_col) +
    ggplot2::scale_fill_manual(values = fill_vals, name = "cluster", drop = FALSE) +
    ggplot2::scale_y_continuous(limits = c(0, 1.13), expand = c(0, 0)) +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(1, 1, 1, 1))
}

#' Combined UMAP + admixture figure
#'
#' Draws a UMAP scatter and an sNMF admixture bar plot as **one** figure that shares a
#' single theme (matching font sizes), a single region colour map (so UMAP points and the
#' admixture colour strips match), and collected legends. The admixture is faceted by a
#' metadata column with a colour strip -- and no text label -- over each region, samples
#' ordered once and reused, laid out over one or more rows.
#'
#' @param x A [PopStructure] with a UMAP (`run_umap()`) and an sNMF fit (`run_snmf()`).
#' @param group Metadata column to facet/colour the admixture by (default the first
#'   non-`sample` metadata column).
#' @param colour,color Metadata column supplying the shared region colours (default `group`).
#' @param K Number of ancestral populations: an integer for a specific K, or `NULL` /
#'   `"best_k"` to use the cross-entropy best K ([PopStructure]'s `best_k()`).
#' @param rows A list of character vectors giving which `group` levels sit on each
#'   admixture row, e.g. `list(c("DRC", "Kenya"), c("Tanzania", "Uganda"))`. Default: all
#'   levels on one row. Levels not listed are dropped.
#' @param orientation `"vertical"` puts the UMAP above the admixture; `"horizontal"` puts
#'   it to the left.
#' @param sample_order Optional explicit sample order (see [admixture_order()]); by
#'   default computed once (within group) and reused across the figure.
#' @param umap_colour Metadata column colouring the UMAP points (default `colour`).
#' @param region_colours,cluster_colours Optional named colour vectors overriding the
#'   region strip / K-cluster fills.
#' @param region_label Legend title for the region colours (default `colour`).
#' @param base_size Base font size shared by every panel.
#' @param border Outline each sample's admixture bar (default `TRUE`) so neighbours with
#'   nearly identical ancestry stay distinct.
#' @param border_colour,border_color,border_linewidth Colour and width of the per-sample outline.
#' @param legend Where to collect the shared legends (`"right"`, `"left"`, `"bottom"`,
#'   `"top"`).
#' @param legend_point_size Size of the region dots in the UMAP legend (default `3.5`, so
#'   the key reads clearly next to the small plotted points).
#' @param point_size,point_alpha Size and opacity of the UMAP scatter points.
#' @param umap_ratio Relative size of the UMAP vs the admixture block (a single number;
#'   default 1 means roughly equal).
#' @param file Optional path to save to (via [save_plot()]) at the auto-computed size.
#' @param width,height Optional output size in inches (default auto from sample counts,
#'   row count, and orientation); also attached to the result as attributes.
#' @return A \pkg{patchwork} object (invisibly if `file` is given), with `width`/`height`
#'   attributes carrying the suggested output size.
#' @examples
#' \dontrun{
#' ps <- example_pop_structure("africa")
#' ps$run_snmf(K = 1:9)
#' plot_structure_figure(ps, group = "site",
#'   rows = list(c("DRC", "Ethiopia", "Sudan"),
#'               c("Kenya_East", "Kenya_West", "Tanzania_East", "Tanzania_West")))
#' }
#' @export
plot_structure_figure <- function(x, group = NULL, colour = group, K = NULL, rows = NULL,
                                  orientation = c("vertical", "horizontal"),
                                  sample_order = NULL, umap_colour = colour,
                                  region_colours = NULL, cluster_colours = NULL,
                                  region_label = NULL, base_size = 11,
                                  border = TRUE, border_colour = "black",
                                  border_linewidth = 0.15, legend = "right",
                                  legend_point_size = 3.5, point_size = 1.6,
                                  point_alpha = 0.8, umap_ratio = 1, file = NULL,
                                  width = NULL, height = NULL,
                                  color = NULL, border_color = NULL) {
  colour <- .alias_arg("colour", "color")
  border_colour <- .alias_arg("border_colour", "border_color")
  .need_package("ggplot2", "plot_structure_figure()")
  .need_package("patchwork", "plot_structure_figure()")
  if (!inherits(x, "PopStructure")) stop("`x` must be a PopStructure", call. = FALSE)
  orientation <- match.arg(orientation)
  or <- function(a, b) if (is.null(a)) b else a

  meta <- x$get_meta()
  if (is.null(meta)) stop("this PopStructure has no metadata", call. = FALSE)
  group <- or(group, setdiff(names(meta), "sample")[1])
  colour <- or(colour, group)
  umap_colour <- or(umap_colour, colour)
  region_label <- or(region_label, colour)
  if (is.null(x$umap_df())) stop("run_umap() first", call. = FALSE)

  if (is.null(K) || identical(K, "best_k") || identical(K, "best")) K <- x$best_k()
  q <- x$q(K)
  samples <- rownames(q)

  region_cols <- or(region_colours,
                    or(x$get_colors()[[colour]], meta_colors(meta, cols = colour)[[colour]]))
  cluster_cols <- or(cluster_colours,
                     stats::setNames(.pick_palette(ncol(q)), colnames(q)))

  if (is.null(sample_order)) sample_order <- admixture_order(q, meta = meta, group = group)
  grp_of <- stats::setNames(as.character(meta[[group]]), meta$sample)

  levs <- .levels_of(meta[[group]])
  if (is.null(rows)) rows <- list(levs)

  # one admixture panel per group level, samples in the shared order
  panel_for <- function(g) {
    ss <- sample_order[grp_of[sample_order] == g]
    ss <- ss[!is.na(ss)]
    .admix_panel(q[ss, , drop = FALSE], ss, cluster_cols, unname(region_cols[g]), base_size,
                 border = border, border_colour = border_colour,
                 border_linewidth = border_linewidth)
  }
  counts <- vapply(levs, function(g) sum(grp_of[sample_order] == g, na.rm = TRUE),
                   numeric(1))

  # assemble each row (widths proportional to sample counts), then stack the rows
  row_blocks <- lapply(rows, function(gr) {
    gr <- gr[gr %in% levs]
    patchwork::wrap_plots(lapply(gr, panel_for), nrow = 1, widths = counts[gr])
  })
  admix <- if (length(row_blocks) == 1) row_blocks[[1]]
           else patchwork::wrap_plots(row_blocks, ncol = 1)

  umap <- plot_umap(x, colour = umap_colour, colors = region_cols,
                    point_size = point_size, point_alpha = point_alpha,
                    legend_point_size = legend_point_size) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::labs(colour = region_label)

  fig <- if (orientation == "vertical") {
    patchwork::wrap_plots(umap, admix, ncol = 1, heights = c(umap_ratio, 1),
                          guides = "collect")
  } else {
    patchwork::wrap_plots(umap, admix, nrow = 1, widths = c(umap_ratio, 1.4),
                          guides = "collect")
  }
  fig <- fig & ggplot2::theme(legend.position = legend)

  # auto output size (inches) from sample counts, rows, orientation
  max_row <- max(vapply(rows, function(gr) sum(counts[gr[gr %in% levs]]), numeric(1)))
  n_rows <- length(rows)
  admix_w <- max_row * 0.045 + 1.6
  admix_h <- n_rows * 1.5 + 0.5
  umap_side <- 5
  if (orientation == "vertical") {
    W <- max(admix_w, umap_side + 1.5); H <- umap_side + admix_h
  } else {
    W <- umap_side + admix_w; H <- max(umap_side, admix_h)
  }
  W <- or(width, round(W, 1)); H <- or(height, round(H, 1))
  attr(fig, "width") <- W
  attr(fig, "height") <- H

  if (!is.null(file)) {
    save_plot(file, fig, width = W, height = H)
    return(invisible(fig))
  }
  fig
}
