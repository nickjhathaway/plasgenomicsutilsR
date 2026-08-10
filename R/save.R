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


# `capabilities("cairo")` answers what R was BUILT with, not whether the shared object can
# be loaded now: an R build can ship a cairo.so whose X11 dependencies are absent, so the
# capability reports TRUE and opening the device then fails with "failed to load cairo DLL"
# and writes nothing. Open one for real and see. Probed once per session, since the answer
# cannot change and the probe costs a file.
.cairo_pdf_works <- local({
  known <- NULL
  function() {
    if (!is.null(known)) return(known)
    known <<- isTRUE(capabilities("cairo")) && .Platform$OS.type != "windows" &&
      tryCatch({
        f <- tempfile(fileext = ".pdf")
        on.exit(unlink(f), add = TRUE)
        before <- grDevices::dev.cur()
        # a broken cairo warns rather than errors, so judge it by whether a device opened
        suppressWarnings(grDevices::cairo_pdf(f))
        opened <- !identical(grDevices::dev.cur(), before)
        if (opened) grDevices::dev.off()
        opened
      }, error = function(e) FALSE)
    known
  }
})

#' The preferred PDF graphics device
#'
#' Returns [grDevices::cairo_pdf()] when a cairo device can actually be opened and the
#' platform is not Windows (cairo embeds fonts more reliably, but its PDF output can be
#' unreliable on Windows), otherwise the string `"pdf"`. Availability is settled by opening
#' a throwaway device rather than by `capabilities("cairo")`, which reports what R was built
#' with and so can claim a cairo that fails to load. Use it with [ggplot2::ggsave()]:
#' `ggsave(file, plot, device = pdf_device())`.
#'
#' @return A device function ([grDevices::cairo_pdf]) or the string `"pdf"`.
#' @examples
#' pdf_device()
#' @export
pdf_device <- function() {
  if (.cairo_pdf_works()) grDevices::cairo_pdf else "pdf"
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
#' @param width,height Output size in inches. **Give one and the other is computed** from
#'   the plot's contents; give neither and both are worked out; give both and they are used
#'   as-is. `NULL` for both falls back to the size the `plot_*()` function attached (see
#'   [plot_dims()]), then to `ggsave()`'s default. For a multi-page list, one page size
#'   serves every page: the width every page needs, and the tallest height any page needs at
#'   that width.
#' @param fit Size the canvas to the drawing (default `TRUE`). A plot whose panel has a
#'   locked aspect ratio -- anything using [ggplot2::coord_fixed()], such as
#'   [plot_ibd_network()] or the gene triangles -- only fills a canvas of one particular
#'   shape; on any other shape the leftover appears as blank margin above and below (or
#'   beside) the drawing. Fitting measures the built plot's fixed furniture (titles,
#'   legends, axes, margins) and its panel ratio, then solves for the dimension you did not
#'   supply so there is no leftover. Set `FALSE` to use the requested / attached numbers
#'   verbatim.
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
                      width = NULL, height = NULL, fit = TRUE, ...) {
  .need_package("ggplot2", "save_plot()")
  # a list/vector of plots -> multi-page PDF (one page each)
  if (is.list(plot) && !.is_plot(plot) && length(plot) && all(vapply(plot, .is_plot, logical(1)))) {
    return(.save_plots_multipage(filename, plot, device, width, height, fit = fit))
  }
  dims <- .fit_plot_dims(plot, width, height, fit = fit)
  width <- dims$width; height <- dims$height
  if (is.null(device) && grepl("\\.pdf$", filename, ignore.case = TRUE)) device <- pdf_device()
  args <- list(filename = filename, plot = plot, device = device, ...)
  if (!is.null(width)) args$width <- width
  if (!is.null(height)) args$height <- height
  do.call(ggplot2::ggsave, args)
  invisible(filename)
}

# Height a legend stacked down the side needs, in inches. A side legend sits in its own
# column spanning the canvas, so nothing in the layout forces the canvas to be tall enough
# for it -- if it is not, grid simply clips the keys off the ends. Measuring the guide-box
# grob is the only reliable way to know: key size, text size, titles and how many columns
# the guides wrapped into all feed in.
.guide_box_height <- function(gt, names = c("guide-box-right", "guide-box-left")) {
  i <- which(gt$layout$name %in% names)
  if (!length(i)) return(0)
  hs <- vapply(i, function(k) {
    g <- gt$grobs[[k]]
    if (is.null(g$heights)) return(0)
    h <- try(grid::convertHeight(sum(g$heights), "in", valueOnly = TRUE), silent = TRUE)
    if (inherits(h, "try-error") || !is.finite(h)) 0 else h
  }, numeric(1))
  max(0, max(hs))
}

.side_legend_height <- function(gt) .guide_box_height(gt)

# Measuring a gtable resolves grob and string widths, which asks the graphics device for
# font metrics -- and if no device is open, R opens the DEFAULT one to answer. In a script
# that is `pdf`, which leaves a stray Rplots.pdf in the working directory; inside a knitr /
# Quarto chunk it is the chunk's own device, where the measurement surfaces as an extra
# empty figure under the block. Measuring on a throwaway null device avoids both: it
# produces no file and is not the device anything is being captured from.
.with_null_device <- function(expr) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
}

# Measure a built plot: the space its fixed furniture needs (titles, legends, axes,
# margins), the height any side legend needs, and, when the panel has a locked aspect
# ratio, that ratio. `null` units convert to zero in absolute terms, so summing the
# gtable's widths/heights leaves exactly the non-panel inches.
.plot_metrics <- function(p) .with_null_device(.plot_metrics_here(p))

.plot_metrics_here <- function(p) {
  out <- list(fixed_w = NA_real_, fixed_h = NA_real_, aspect = NULL, n_row = 1L, n_col = 1L,
              legend_h = 0)
  gt <- try(ggplot2::ggplotGrob(p), silent = TRUE)
  if (inherits(gt, "try-error")) return(out)
  out$fixed_w <- grid::convertWidth(sum(gt$widths), "in", valueOnly = TRUE)
  out$fixed_h <- grid::convertHeight(sum(gt$heights), "in", valueOnly = TRUE)
  out$legend_h <- .side_legend_height(gt)
  b <- try(ggplot2::ggplot_build(p), silent = TRUE)
  if (inherits(b, "try-error")) return(out)
  lay <- b$layout$layout
  if (!is.null(lay) && all(c("ROW", "COL") %in% names(lay))) {
    out$n_row <- max(1L, max(lay$ROW)); out$n_col <- max(1L, max(lay$COL))
  }
  # coord_fixed()/coord_equal() report a panel height:width; free coords return NULL
  asp <- try(b$layout$coord$aspect(b$layout$panel_params[[1]]), silent = TRUE)
  if (!inherits(asp, "try-error") && length(asp) == 1 && is.finite(asp) && asp > 0) {
    out$aspect <- as.numeric(asp)
  }
  out
}

# Resolve the output size. A plot whose panel has a locked aspect ratio only fills a canvas
# of one particular shape -- any other shape shows as blank margin -- so given one dimension
# the other is computed rather than guessed, and given neither the attached width is kept
# and the height solved for. Falls back to the attached size (then ggsave's default) for
# plots with no fixed aspect, where any canvas shape is legitimate.
.fit_plot_dims <- function(p, width = NULL, height = NULL, fit = TRUE) {
  d <- attr(p, "plasgenomics_dims")
  att_w <- if (!is.null(d) && "width" %in% names(d)) unname(d[["width"]]) else NULL
  att_h <- if (!is.null(d) && "height" %in% names(d)) unname(d[["height"]]) else NULL
  if (!isTRUE(fit) || (!is.null(width) && !is.null(height))) {
    return(list(width = width %||% att_w, height = height %||% att_h))
  }
  m <- .plot_metrics(p)
  # a side legend is clipped rather than accommodated, so treat its height as a floor
  min_h <- if (is.finite(m$legend_h) && m$legend_h > 0) m$legend_h + 0.25 else 0
  if (is.null(m$aspect) || !is.finite(m$fixed_w) || !is.finite(m$fixed_h)) {
    # no locked aspect: keep whatever was asked for, and take the aspect of the attached
    # size for the missing side so an explicit width still yields a sensible height
    if (!is.null(width) && is.null(height) && !is.null(att_w) && !is.null(att_h)) {
      return(list(width = width, height = round(max(width * att_h / att_w, min_h), 2)))
    }
    if (!is.null(height) && is.null(width) && !is.null(att_w) && !is.null(att_h)) {
      return(list(width = round(height * att_w / att_h, 2), height = height))
    }
    h <- height %||% att_h
    if (!is.null(h)) h <- round(max(h, min_h), 2)
    return(list(width = width %||% att_w, height = h))
  }
  # panel_h_total = aspect * panel_w_each * n_row, panel_w_total = panel_w_each * n_col
  k <- m$aspect * m$n_row / m$n_col
  if (is.null(width) && is.null(height)) width <- att_w %||% max(.MIN_PLOT_WIDTH, m$fixed_w + 4)
  if (!is.null(width)) {
    panel_w <- width - m$fixed_w
    if (!is.finite(panel_w) || panel_w <= 0.5) return(list(width = width, height = height %||% att_h))
    return(list(width = round(width, 2),
                height = round(max(k * panel_w + m$fixed_h, min_h), 2)))
  }
  panel_h <- height - m$fixed_h
  if (!is.finite(panel_h) || panel_h <= 0.5) return(list(width = att_w, height = height))
  list(width = round(panel_h / k + m$fixed_w, 2), height = round(height, 2))
}

# a single plot object (ggplot, or a patchwork assembly, both inherit "ggplot").
# ggplot2 renamed is.ggplot() to is_ggplot() in 3.5.2; accept either installed version.
.is_plot <- function(p) {
  fn <- getNamespace("ggplot2")[["is_ggplot"]]
  if (is.null(fn)) fn <- ggplot2::is.ggplot
  fn(p) || inherits(p, "patchwork")
}

# largest attached `plasgenomics_dims[key]` across a list of plots, or `default`
.max_dim <- function(plots, key, default) {
  vals <- vapply(plots, function(p) {
    d <- attr(p, "plasgenomics_dims"); if (is.null(d)) NA_real_ else unname(d[key])
  }, numeric(1))
  if (all(is.na(vals))) default else max(vals, na.rm = TRUE)
}

# write a list of plots as one multi-page PDF (one plot per page, uniform page size)
.save_plots_multipage <- function(filename, plots, device, width, height, fit = TRUE) {
  if (!grepl("\\.pdf$", filename, ignore.case = TRUE)) {
    stop("a list of plots is written as a multi-page PDF, so `filename` must end in '.pdf'; ",
         "got '", filename, "'. Save the plots individually for other formats.", call. = FALSE)
  }
  # One page size has to serve every page, so take the width every page needs and the
  # tallest height any page needs at that width -- fitting each page separately would give
  # a different size per page, which a single PDF cannot do.
  if (is.null(width)) width <- .max_dim(plots, "width", default = 7)
  if (is.null(height)) {
    if (isTRUE(fit)) {
      hs <- vapply(plots, function(p) {
        h <- .fit_plot_dims(p, width = width, height = NULL, fit = TRUE)$height
        if (is.null(h) || !is.finite(h)) NA_real_ else h
      }, numeric(1))
      height <- if (all(is.na(hs))) .max_dim(plots, "height", default = 7) else max(hs, na.rm = TRUE)
    } else {
      height <- .max_dim(plots, "height", default = 7)
    }
  }
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
