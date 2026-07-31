# Saving plots to PDF with good font handling, and size suggestions that scale with
# how much is being drawn (chromosomes across, regions / features down).

.GENOME_FULL_WIDTH <- 16   # inches for a full 14-chromosome genome-wide track
.MIN_PLOT_WIDTH <- 6

# genome-wide tracks (manhattan / selection / tug-of-war): width scales with the
# genome fraction shown, height with the number of region panels.
.dims_genome <- function(n_panels, genome_frac, per_panel = 1.7, label_genes = FALSE) {
  w <- max(.MIN_PLOT_WIDTH, .GENOME_FULL_WIDTH * genome_frac)
  h <- 0.9 + per_panel * max(1, n_panels) + if (label_genes) 0.7 else 0
  c(width = round(w, 1), height = round(h, 1))
}

# region x region heatmap: as above for width; height also grows with the region
# count (each of the n_panels anchor facets has n_regions rows).
.dims_heatmap <- function(n_panels, n_regions, genome_frac, label_genes = FALSE) {
  w <- max(.MIN_PLOT_WIDTH, .GENOME_FULL_WIDTH * genome_frac) + 1.2
  h <- 0.6 + n_panels * (n_regions * 0.28 + 0.35) + if (label_genes) 0.7 else 0
  c(width = round(w, 1), height = round(h, 1))
}

# triangles: a grid of n_features square panels, each n_regions x n_regions.
.dims_triangles <- function(n_features, n_regions, ncol = NULL) {
  if (is.null(ncol) || ncol < 1) ncol <- max(1L, min(n_features, ceiling(sqrt(n_features))))
  nrows <- ceiling(n_features / ncol)
  w <- ncol * (n_regions * 0.5 + 0.6) + 1.3
  h <- nrows * (n_regions * 0.5 + 0.7)
  c(width = round(w, 1), height = round(h, 1))
}

#' Suggested output dimensions for an IBD plot
#'
#' Returns a `c(width, height)` (inches) that scales with how much a plot will draw:
#' the genome fraction shown across (chromosomes kept), and the number of region
#' panels / triangle features down. The `plot_*()` functions already attach this to
#' the plot they return, so [save_plot()] uses it automatically; call `plot_dims()`
#' yourself only to inspect or override the numbers.
#'
#' @param x An [IbdResults] object.
#' @param type One of `"manhattan"`, `"selection"`, `"tugofwar"`, `"heatmap"`,
#'   `"triangles"`.
#' @param regions,chroms,skip_chr,genes,snps,ncol,label_genes The same selection
#'   arguments you pass to the plot, so the counts match what will be drawn.
#' @return A named numeric vector `c(width = , height = )` in inches.
#' @examples
#' plot_dims(example_ibd_results(), "tugofwar")
#' @export
plot_dims <- function(x, type = c("manhattan", "selection", "tugofwar",
                                  "heatmap", "triangles"),
                      regions = NULL, chroms = NULL, skip_chr = NULL,
                      genes = NULL, snps = NULL, ncol = NULL, label_genes = FALSE) {
  type <- match.arg(type)
  full <- x$chrom_layout()
  lay <- .select_layout(full, chroms, skip_chr)
  gfrac <- sum(lay$len) / sum(full$len)

  if (type == "triangles") {
    pw <- x$get_pairwise_region()
    n_regions <- if (is.null(pw)) 2 else length(unique(c(pw$region_a, pw$region_b)))
    gt <- x$get_genes()
    n_gene <- if (!is.null(gt) && (is.null(snps) || !is.null(genes))) {
      if (!is.null(genes)) sum(tolower(gt$name) %in% tolower(genes)) else nrow(gt)
    } else 0
    n_snp <- if (is.null(snps)) 0 else if (is.data.frame(snps)) nrow(snps) else length(snps)
    return(.dims_triangles(max(1, n_gene + n_snp), max(2, n_regions), ncol))
  }
  if (type == "heatmap") {
    pw <- x$get_pairwise_region()
    regs <- if (is.null(pw)) 1 else length(unique(c(pw$region_a, pw$region_b)))
    return(.dims_heatmap(regs, regs, gfrac, label_genes))
  }
  tbl <- if (type == "selection") x$get_selection() else x$get_per_snp_region()
  regs <- if (!is.null(tbl) && "region" %in% names(tbl)) {
    r <- unique(tbl$region)
    if (!is.null(regions)) r <- intersect(r, regions)
    length(r)
  } else 1
  .dims_genome(regs, gfrac, label_genes = label_genes)
}


#' The preferred PDF graphics device
#'
#' Returns [grDevices::cairo_pdf()] when cairo is available and the platform is not
#' Windows (cairo embeds fonts more reliably, but its PDF output can be unreliable on
#' Windows), otherwise the string `"pdf"`. Use it with [ggplot2::ggsave()]:
#' `ggsave(file, plot, device = pdf_device())`.
#'
#' @return A device function ([grDevices::cairo_pdf]) or the string `"pdf"`.
#' @examples
#' pdf_device()
#' @export
pdf_device <- function() {
  cairo_ok <- isTRUE(capabilities("cairo")) && .Platform$OS.type != "windows"
  if (cairo_ok) grDevices::cairo_pdf else "pdf"
}

#' Save a plot, preferring the cairo PDF device
#'
#' A thin wrapper around [ggplot2::ggsave()] that, for `.pdf` output, defaults to the
#' cairo PDF device (better font embedding) via [pdf_device()], falling back to the
#' standard `pdf` device where cairo is unavailable or unreliable (e.g. Windows).
#' Non-PDF outputs, and any explicit `device`, pass straight through to `ggsave()`.
#'
#' @param filename Output path.
#' @param plot Plot to save (defaults to the last plot displayed).
#' @param device Graphics device. `NULL` (default) auto-selects for `.pdf` output
#'   (see [pdf_device()]); pass `grDevices::cairo_pdf`, `"pdf"`, or any device
#'   [ggplot2::ggsave()] accepts to force a choice.
#' @param width,height Output size in inches. `NULL` (default) uses the size the
#'   `plot_*()` function attached to the plot (see [plot_dims()]), falling back to
#'   `ggsave()`'s default if none is present.
#' @param ... Passed to [ggplot2::ggsave()] (e.g. `dpi`, `units`).
#' @return `filename`, invisibly.
#' @examples
#' \dontrun{
#' # size is chosen automatically from the number of regions / chromosomes drawn
#' save_plot("ibd_manhattan.pdf", plot_ibd_manhattan(example_ibd_results()))
#' }
#' @export
save_plot <- function(filename, plot = ggplot2::last_plot(), device = NULL,
                      width = NULL, height = NULL, ...) {
  .need_package("ggplot2", "save_plot()")
  d <- attr(plot, "plasgenomics_dims")
  if (is.null(width) && !is.null(d)) width <- unname(d["width"])
  if (is.null(height) && !is.null(d)) height <- unname(d["height"])
  if (is.null(device) && grepl("\\.pdf$", filename, ignore.case = TRUE)) device <- pdf_device()
  args <- list(filename = filename, plot = plot, device = device, ...)
  if (!is.null(width)) args$width <- width
  if (!is.null(height)) args$height <- height
  do.call(ggplot2::ggsave, args)
  invisible(filename)
}
