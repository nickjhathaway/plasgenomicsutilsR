# Saving plots to PDF with good font handling, and size suggestions that scale with
# how much is being drawn (chromosomes across, groups / features down).

.GENOME_FULL_WIDTH <- 16   # inches for a full 14-chromosome genome-wide track
.MIN_PLOT_WIDTH <- 6

# genome-wide tracks (manhattan / selection / tug-of-war): width scales with the
# genome fraction shown, height with the number of group panels.
.dims_genome <- function(n_panels, genome_frac, per_panel = 1.7, label_genes = FALSE) {
  w <- max(.MIN_PLOT_WIDTH, .GENOME_FULL_WIDTH * genome_frac)
  h <- 0.9 + per_panel * max(1, n_panels) + if (label_genes) 0.7 else 0
  c(width = round(w, 1), height = round(h, 1))
}

# group x group heatmap: as above for width; height also grows with the group
# count (each of the n_panels anchor facets has n_groups rows).
.dims_heatmap <- function(n_panels, n_groups, genome_frac, label_genes = FALSE) {
  w <- max(.MIN_PLOT_WIDTH, .GENOME_FULL_WIDTH * genome_frac) + 1.2
  h <- 0.6 + n_panels * (n_groups * 0.28 + 0.35) + if (label_genes) 0.7 else 0
  c(width = round(w, 1), height = round(h, 1))
}

# triangles: a grid of n_features square panels, each n_groups x n_groups.
.dims_triangles <- function(n_features, n_groups, ncol = NULL) {
  if (is.null(ncol) || ncol < 1) ncol <- max(1L, min(n_features, ceiling(sqrt(n_features))))
  nrows <- ceiling(n_features / ncol)
  w <- ncol * (n_groups * 0.5 + 0.6) + 1.3
  h <- nrows * (n_groups * 0.5 + 0.7)
  c(width = round(w, 1), height = round(h, 1))
}

#' Suggested output dimensions for an IBD plot
#'
#' Returns a `c(width, height)` (inches) that scales with how much a plot will draw:
#' the genome fraction shown across (chromosomes kept), and the number of group
#' panels / triangle features down. The `plot_*()` functions already attach this to
#' the plot they return, so [save_plot()] uses it automatically; call `plot_dims()`
#' yourself only to inspect or override the numbers.
#'
#' @param x An [IbdResults] object.
#' @param type One of `"manhattan"`, `"selection"`, `"tugofwar"`, `"heatmap"`,
#'   `"triangles"`.
#' @param groups,chroms,skip_chr,genes,snps,ncol,label_genes The same selection
#'   arguments you pass to the plot, so the counts match what will be drawn.
#' @return A named numeric vector `c(width = , height = )` in inches.
#' @examples
#' plot_dims(example_ibd_results(), "tugofwar")
#' @export
plot_dims <- function(x, type = c("manhattan", "selection", "tugofwar",
                                  "heatmap", "triangles"),
                      groups = NULL, chroms = NULL, skip_chr = NULL,
                      genes = NULL, snps = NULL, ncol = NULL, label_genes = FALSE) {
  type <- match.arg(type)
  full <- x$chrom_layout()
  lay <- .select_layout(full, chroms, skip_chr)
  gfrac <- sum(lay$len) / sum(full$len)

  if (type == "triangles") {
    pw <- x$get_pairwise_group()
    n_groups <- if (is.null(pw)) 2 else length(unique(c(pw$group_a, pw$group_b)))
    gt <- x$get_genes()
    n_gene <- if (!is.null(gt) && (is.null(snps) || !is.null(genes))) {
      if (!is.null(genes)) sum(tolower(gt$name) %in% tolower(genes)) else nrow(gt)
    } else 0
    n_snp <- if (is.null(snps)) 0 else if (is.data.frame(snps)) nrow(snps) else length(snps)
    return(.dims_triangles(max(1, n_gene + n_snp), max(2, n_groups), ncol))
  }
  if (type == "heatmap") {
    pw <- x$get_pairwise_group()
    regs <- if (is.null(pw)) 1 else length(unique(c(pw$group_a, pw$group_b)))
    return(.dims_heatmap(regs, regs, gfrac, label_genes))
  }
  tbl <- if (type == "selection") x$get_selection() else x$get_per_snp_group()
  regs <- if (!is.null(tbl) && "group" %in% names(tbl)) {
    r <- unique(tbl$group)
    if (!is.null(groups)) r <- intersect(r, groups)
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
#' @param plot A single plot (defaults to the last plot displayed), **or a list/vector of
#'   plots** -- e.g. the list returned by `plot_pairwise_ibd_for_genes(individual = TRUE)`
#'   -- which is written as a **multi-page PDF**, one plot per page (`filename` must end
#'   in `.pdf`).
#' @param device Graphics device. `NULL` (default) auto-selects for `.pdf` output
#'   (see [pdf_device()]); pass `grDevices::cairo_pdf`, `"pdf"`, or any device
#'   [ggplot2::ggsave()] accepts to force a choice.
#' @param width,height Output size in inches. `NULL` (default) uses the size the
#'   `plot_*()` function attached to the plot (see [plot_dims()]), falling back to
#'   `ggsave()`'s default if none is present. For a multi-page list, the page size is the
#'   largest attached size across the plots (one size for every page).
#' @param ... Passed to [ggplot2::ggsave()] (e.g. `dpi`, `units`) for a single plot;
#'   ignored for a multi-page list (which honours only `width` / `height` / `device`).
#' @return `filename`, invisibly.
#' @examples
#' \dontrun{
#' # size is chosen automatically from the number of groups / chromosomes drawn
#' save_plot("ibd_manhattan.pdf", plot_ibd_sharing_manhattan(example_ibd_results()))
#' # a list of plots -> one multi-page PDF
#' save_plot("triangles.pdf",
#'           plot_pairwise_ibd_for_genes(example_ibd_results(), individual = TRUE))
#' }
#' @export
save_plot <- function(filename, plot = ggplot2::last_plot(), device = NULL,
                      width = NULL, height = NULL, ...) {
  .need_package("ggplot2", "save_plot()")
  # a list/vector of plots -> multi-page PDF (one page each)
  if (is.list(plot) && !.is_plot(plot) && length(plot) && all(vapply(plot, .is_plot, logical(1)))) {
    return(.save_plots_multipage(filename, plot, device, width, height))
  }
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

# a single plot object (ggplot, or a patchwork assembly, both inherit "ggplot")
.is_plot <- function(p) ggplot2::is.ggplot(p) || inherits(p, "patchwork")

# largest attached `plasgenomics_dims[key]` across a list of plots, or `default`
.max_dim <- function(plots, key, default) {
  vals <- vapply(plots, function(p) {
    d <- attr(p, "plasgenomics_dims"); if (is.null(d)) NA_real_ else unname(d[key])
  }, numeric(1))
  if (all(is.na(vals))) default else max(vals, na.rm = TRUE)
}

# write a list of plots as one multi-page PDF (one plot per page, uniform page size)
.save_plots_multipage <- function(filename, plots, device, width, height) {
  if (!grepl("\\.pdf$", filename, ignore.case = TRUE)) {
    stop("a list of plots is written as a multi-page PDF, so `filename` must end in '.pdf'; ",
         "got '", filename, "'. Save the plots individually for other formats.", call. = FALSE)
  }
  if (is.null(width))  width  <- .max_dim(plots, "width", default = 7)
  if (is.null(height)) height <- .max_dim(plots, "height", default = 7)
  dev <- if (is.null(device)) pdf_device() else device
  if (is.function(dev)) {
    dev(filename, width = width, height = height, onefile = TRUE)      # e.g. cairo_pdf
  } else if (identical(dev, "pdf") || identical(dev, "cairo_pdf")) {
    grDevices::pdf(filename, width = width, height = height, onefile = TRUE)
  } else {
    stop("multi-page output needs a PDF device; `device` was not a PDF device", call. = FALSE)
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  for (p in plots) print(p)
  invisible(filename)
}
