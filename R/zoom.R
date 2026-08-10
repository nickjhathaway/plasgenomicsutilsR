# Zooming a genome-wide plot to one interval.
#
# Zoom is deliberately separate from `chroms` / `skip_chr`: those choose which chromosomes
# make up the cumulative axis, while zoom keeps that axis and only crops it, so a zoomed
# plot and the genome-wide one it came from put the same locus at the same coordinate and
# carry the same data behind it.

# Parse a zoom request into chr/start/end. Accepts a chromosome ("7", "chr7",
# "Pf3D7_07_v3"), a range ("7:728,081-988,719"), a gene name when `genes` is supplied, or
# a one-row data frame with chr/start/end.
.resolve_region <- function(region, genes = NULL, reference = DEFAULT_REFERENCE,
                            fallback = NULL) {
  if (is.null(region)) return(NULL)

  if (is.data.frame(region)) {
    if (!all(c("chr", "start", "end") %in% names(region)))
      stop("a data frame `zoom` needs chr, start and end columns", call. = FALSE)
    return(list(chr = normalise_chr(region$chr[1]), start = as.numeric(region$start[1]),
                end = as.numeric(region$end[1])))
  }
  if (!is.character(region) || length(region) != 1)
    stop("`zoom` must be a single string, or a one-row data frame with chr/start/end",
         call. = FALSE)

  m <- regmatches(region, regexec("^(.*):([0-9,]+)\\s*-\\s*([0-9,]+)$", region))[[1]]
  if (length(m) == 4) {
    num <- function(z) as.numeric(gsub(",", "", z))
    return(list(chr = normalise_chr(m[2]), start = num(m[3]), end = num(m[4])))
  }

  for (src in list(genes, fallback)) {
    if (is.null(src)) next
    g <- .gene_track(src)
    hit <- which(tolower(g$name) == tolower(region))
    # a full annotation carries both a symbol and an id, and either should resolve
    if (!length(hit) && "gene_id" %in% names(as.data.frame(src)))
      hit <- which(toupper(as.data.frame(src)$gene_id) == toupper(region))
    if (length(hit))
      return(list(chr = g$chr[hit[1]], start = as.numeric(g$start[hit[1]]),
                  end = as.numeric(g$end[hit[1]])))
  }

  ch <- normalise_chr(region)
  lens <- .chrom_layout(reference)
  if (ch %in% lens$chr)
    return(list(chr = ch, start = 0, end = lens$len[match(ch, lens$chr)]))

  stop("could not read `zoom = \"", region, "\"`: expected a chromosome, a ",
       "\"chr:start-end\" range",
       if (!is.null(genes) || !is.null(fallback)) ", a gene name in the track" else "",
       ", or a data frame with chr/start/end", call. = FALSE)
}

# How much to add on each side. One number pads both; two pad the left and the right, either
# in that order or named `left` / `right` (and naming just one pads only that side). Each side
# is judged on its own, so a fraction on one and base pairs on the other is legal -- the
# fraction is always of the interval's own span, taken before anything is added, so the two
# sides do not depend on the order they are applied in.
.pad_sides <- function(pad, span) {
  none <- c(0, 0)
  if (is.null(pad)) return(none)
  pad <- pad[!is.na(pad)]
  if (!length(pad)) return(none)
  if (!is.numeric(pad))
    stop("padding must be numeric: one value for both sides, or two for the left and the ",
         "right", call. = FALSE)
  if (length(pad) > 2)
    stop("padding takes at most two values (the left and the right), not ", length(pad),
         call. = FALSE)
  nm <- names(pad)
  if (!is.null(nm) && any(nzchar(nm))) {
    if (!all(nzchar(nm)) || !all(nm %in% c("left", "right")) || anyDuplicated(nm))
      stop("padding names must be \"left\" and/or \"right\", each at most once; name every ",
           "value or none of them", call. = FALSE)
    sides <- c(left = 0, right = 0)
    sides[nm] <- pad
    pad <- unname(sides)
  } else if (length(pad) == 1) {
    pad <- c(pad, pad)
  }
  out <- ifelse(pad < 1, span * pad, pad)
  out[!is.finite(out) | out < 0] <- 0
  out
}

# Pad an interval, clamped to the chromosome, so a gene-sized zoom still shows its
# surroundings. `min_span` widens the window further when the padded interval is still too
# narrow to hold anything -- a 3 kb gene padded by 10% is a 3.6 kb window, which on typical
# SNP density contains no SNPs at all.
.pad_region <- function(iv, pad, layout, min_span = 0) {
  if (is.null(iv)) return(NULL)
  len <- layout$len[match(iv$chr, layout$chr)]
  grow <- function(iv, left, right) {
    iv$start <- max(0, iv$start - left)
    iv$end <- if (is.na(len)) iv$end + right else min(len, iv$end + right)
    iv
  }
  by <- .pad_sides(pad, iv$end - iv$start)
  if (any(by > 0)) iv <- grow(iv, by[1], by[2])
  short <- min_span - (iv$end - iv$start)
  if (is.finite(short) && short > 0) iv <- grow(iv, short / 2, short / 2)
  iv
}

# Inside one chromosome the cumulative axis is just position, so label it in bp/kb/Mb
# rather than leaving the genome-wide chromosome ticks on a 200 kb window.
.region_axis <- function(iv, offset) {
  span <- iv$end - iv$start
  unit <- if (span >= 2e6) 1e6 else if (span >= 1e4) 1e3 else 1
  suffix <- c("bp", "kb", "Mb")[match(unit, c(1, 1e3, 1e6))]
  ggplot2::scale_x_continuous(
    name = paste0("chromosome ", iv$chr, " position (", suffix, ")"),
    # break on round positions within the chromosome, not round values of the cumulative
    # axis, which is offset by every chromosome before this one
    breaks = function(lims) offset + pretty(lims - offset, n = 6),
    labels = function(v) format(round((v - offset) / unit, if (unit == 1) 0 else 1),
                                big.mark = ",", trim = TRUE, scientific = FALSE),
    expand = ggplot2::expansion(mult = 0.01))
}

# Everything a genome-wide plot needs to become a zoomed one. Returns NULL when `zoom` is
# NULL so callers can keep one code path. `genes` is the plot's full gene track: inside a
# window every gene in view is worth drawing and naming, not just the highlighted ones.
.zoom_setup <- function(zoom, genes, layout, label_genes = NULL, pad = 0.05,
                        reference = DEFAULT_REFERENCE, min_span = 0,
                        genes_for_track = NULL) {
  if (is.null(zoom)) return(NULL)
  # a gene name resolves against either source, so `zoom = "pfcrt"` works whether pfcrt is
  # in the plot's own short track or only in the full annotation
  iv <- .resolve_region(zoom, genes, reference, fallback = genes_for_track)
  iv <- .pad_region(iv, pad, layout, min_span)
  offset <- layout$offset[match(iv$chr, layout$chr)]
  if (is.na(offset))
    stop("`zoom` is on chromosome ", iv$chr, ", which this plot is not showing (see ",
         "`chroms` / `skip_chr`)", call. = FALSE)
  xlim <- c(offset + iv$start, offset + iv$end)

  in_window <- function(src) {
    g <- .genes_for_layout(src, layout, NULL)
    if (is.null(g)) return(NULL)
    off <- layout$offset[match(normalise_chr(g$chr), layout$chr)]
    keep <- (off + as.numeric(g$end)) >= xlim[1] & (off + as.numeric(g$start)) <= xlim[2]
    if (any(keep)) g[keep, , drop = FALSE] else NULL
  }
  # `genes` marks positions inside the panel, `genes_for_track` fills the track underneath.
  # Keeping them apart is what lets a plot mark a few genes of interest while still showing
  # every gene in the window, instead of drawing a reference line per gene.
  g <- in_window(genes)
  track <- if (is.null(genes_for_track)) g else in_window(genes_for_track)
  list(interval = iv, offset = offset, xlim = xlim, genes = g, track = track,
       label = if (is.null(label_genes)) !is.null(track) else isTRUE(label_genes),
       # label footprints must be measured against the visible span, not the whole genome
       layout = data.frame(xmin = xlim[1], xmax = xlim[2]))
}

# Rows whose cumulative position falls inside the window. Zoom filters the data as well as
# cropping the axis: several of these plots derive their y scaling from the data (the
# tug-of-war normalises each track to its maximum), so leaving the genome-wide rows in place
# would scale a 50 kb window against a genome-wide peak and flatten it. Filtering does not
# move anything, because the x coordinate is absolute either way.
.in_window <- function(df, z) {
  if (is.null(z) || is.null(df)) return(df)
  df[is.finite(df$cum_pos) & df$cum_pos >= z$xlim[1] & df$cum_pos <= z$xlim[2], ,
     drop = FALSE]
}

.window_label <- function(z) {
  sprintf("%s:%.0f-%.0f", z$interval$chr, z$interval$start, z$interval$end)
}

# Crop to the window, refusing an empty one rather than drawing a blank panel.
.crop_to_window <- function(df, z, what) {
  out <- .in_window(df, z)
  if (!is.null(z) && !nrow(out))
    stop("no ", what, " in ", .window_label(z), "; widen the window with `zoom_pad`",
         call. = FALSE)
  out
}

# Crop and relabel the finished plot, then stack the gene track under it when the window
# has genes to name. The scale and coord are added last so they replace the genome-wide
# ones; ggplot notes both replacements, which is expected here.
.apply_zoom <- function(p, z, layout = NULL, n_panels = 1L, angle = 0) {
  if (is.null(z)) return(p)
  p <- suppressMessages(
    p +
      .region_axis(z$interval, z$offset) +
      # cropping is the point: points outside the window must not spill past the panel
      ggplot2::coord_cartesian(xlim = z$xlim, clip = "on"))
  if (z$label && !is.null(layout)) p <- .stack_gene_track(p, z, layout, n_panels, angle)
  p
}

# Zoomed panels are one window wide, so they do not scale with how much genome is shown.
.ZOOM_WIDTH_IN <- 9

.dims_zoom <- function(n_panels, track_in = 0) {
  if (is.null(track_in) || !is.finite(track_in)) track_in <- 0
  c(width = .ZOOM_WIDTH_IN,
    height = max(2.6, 1.1 + 1.6 * max(1L, n_panels)) + track_in)
}

# Gene-label footprints are estimated from the visible span and the plot's width, both of
# which zoom changes: the window replaces the genome, and the width is fixed rather than
# proportional to how much genome is on screen.
.fp_layout <- function(z, layout) if (is.null(z)) layout else z$layout
.fp_frac <- function(z, genome_frac) if (is.null(z)) genome_frac else .ZOOM_WIDTH_IN / 16

# ---- gene track panel -------------------------------------------------------

# A standalone track of gene boxes with their names, on the same cumulative x axis as the
# plot it sits under. Genome-wide the gene track is reference lines through the panel with
# names above it; inside a window there is room to draw each gene's real extent instead, so
# a zoomed plot gets this panel stacked underneath and keeps its own y axis clean.
#
# Boxes are packed into rows: a gene that would collide with the one before it drops to the
# next row, so overlapping genes stay readable at their true coordinates.
#
# The panel is measured in inches rather than relative units. One row is .ROW_IN, and the
# space under the boxes is whatever the labels need at this angle and length, so text can
# never overflow the strip it was given. The height is returned as an attribute for the
# caller to lay out and to size the page with.
.ROW_IN <- 0.17           # inches per box row
.CHAR_IN <- 2.845 / 72 * 0.5    # per point of text size: pt -> in, ~half-em per character

.gene_track_panel <- function(genes, xlim, height = 0.9, fill = "#4d4d4d", angle = 0,
                              width_in = 9, size = .GENE_LABEL_SIZE, min_width = TRUE) {
  .need_package("ggplot2", "the zoomed gene track")
  if (is.null(genes) || !nrow(genes)) return(NULL)
  g <- genes[order(genes$.gene_xmin), , drop = FALSE]
  span <- diff(xlim)
  # A small gene in a wide window would draw a hairline, so give every box a floor width --
  # but never widen it past a neighbouring box, and never at all when the boxes were already
  # snapped to the SNP columns they own (`min_width = FALSE`), since growing one there would put
  # it back over a SNP it does not contain, which is the thing the snapping avoids.
  short <- min_width & (g$.gene_xmax - g$.gene_xmin) < span * 0.004
  if (any(short)) {
    ord <- order(g$.gene_xmin)
    lo_lim <- rep(-Inf, nrow(g)); hi_lim <- rep(Inf, nrow(g))
    lo_lim[ord][-1] <- g$.gene_xmax[ord][-nrow(g)]
    hi_lim[ord][-nrow(g)] <- g$.gene_xmin[ord][-1]
    mid <- (g$.gene_xmin + g$.gene_xmax) / 2
    g$.gene_xmin[short] <- pmax(mid[short] - span * 0.002, lo_lim[short])
    g$.gene_xmax[short] <- pmin(mid[short] + span * 0.002, hi_lim[short])
  }
  # Anchor on the middle of the VISIBLE part of the gene: a gene straddling the edge of the
  # window has its true mid outside the panel, so a label placed there is cut off entirely.
  g$.mid <- (pmax(g$.gene_xmin, xlim[1]) + pmin(g$.gene_xmax, xlim[2])) / 2

  rad <- angle * pi / 180
  char_in <- size * .CHAR_IN
  bp_per_in <- span / max(2, width_in - 1.2)
  # horizontal room each name takes, in base pairs
  fp <- nchar(g$name) * char_in * cos(rad) * bp_per_in

  # Over a full annotation the boxes are far enough apart while a name like "PF3D7_0720300"
  # is several times wider than the gene under it, so the packing has to know how wide the
  # text is; and hanging rotated names from one baseline is what stops a row's label from
  # crossing the row beneath it.
  # A name is anchored at its gene, so one near the edge of the window would run outside the
  # panel and be cut off. Those flip from centred to starting (or ending) at the gene, which
  # keeps the whole name on the page without moving where it points.
  if (angle == 0) {
    g$.hjust <- 0.5
    g$.hjust[g$.mid - fp / 2 < xlim[1]] <- 0
    g$.hjust[g$.mid + fp / 2 > xlim[2]] <- 1
  } else {
    g$.hjust <- 1                                   # runs down to the left of the gene
    # unless that would leave the window -- and then only if going right actually fits
    flip <- g$.mid - fp < xlim[1] & g$.mid + fp <= xlim[2]
    g$.hjust[flip] <- 0
  }
  lab_lo <- g$.mid - fp * g$.hjust
  lab_hi <- lab_lo + fp

  # Horizontal labels sit under their own box, so a row is only usable if the NAME fits too.
  # Rotated labels all hang from one baseline below the whole track (see below), so there
  # only the boxes decide the packing.
  if (angle == 0) {
    lo <- pmin(g$.gene_xmin, lab_lo)
    hi <- pmax(g$.gene_xmax, lab_hi)
  } else {
    lo <- g$.gene_xmin
    hi <- g$.gene_xmax
  }
  g$.row <- 0L
  last_hi <- numeric(0)
  for (k in seq_len(nrow(g))) {
    r <- 0L
    while (r < length(last_hi) && lo[k] < last_hi[r + 1] + span * 0.01) r <- r + 1L
    g$.row[k] <- r
    last_hi[r + 1] <- hi[k]
  }
  n_rows <- max(g$.row) + 1L
  g$.y <- -g$.row

  if (angle == 0) {
    g$.label_y <- g$.y - height * 0.8
    label_in <- char_in / .CHAR_IN * 0.02 + 0.13          # one line of text
  } else {
    # one baseline for every name, tiered when two rotated footprints still overlap
    ord <- order(lab_lo)
    g <- g[ord, , drop = FALSE]
    lo_o <- lab_lo[ord]; hi_o <- lab_hi[ord]
    # A real gap between names on a tier, not merely "not overlapping": two labels ending a
    # couple of hundred base pairs apart read as one run of text. The clearance cannot be
    # scaled by cos(angle) alone -- at 90 degrees that is zero, so labels would never tier
    # however close together they are (four codons of one gene, say). A rotated name still
    # needs about a line's height of horizontal room from its neighbour.
    line_in <- size * 2.845 / 72 * 1.2
    gap <- max(1.5 * char_in * cos(rad), line_in * sin(rad)) * bp_per_in
    g$.tier <- 0L
    last <- numeric(0)
    for (k in seq_len(nrow(g))) {
      ti <- 0L
      while (ti < length(last) && lo_o[k] < last[ti + 1] + gap) ti <- ti + 1L
      g$.tier[k] <- ti
      last[ti + 1] <- hi_o[k]
    }
    # Each tier drops by the reach of the names on the tier above it, not by the reach of the
    # longest name anywhere, so a tier of short names does not push the next one far down.
    reach <- vapply(split(nchar(g$name), g$.tier), max, numeric(1)) * char_in * sin(rad)
    drop_in <- c(0, cumsum(reach))[seq_along(reach)]
    names(drop_in) <- names(reach)
    g$.label_y <- -(n_rows - 1) - height - drop_in[as.character(g$.tier)] / .ROW_IN
    label_in <- sum(reach) + 0.1
  }

  span_units <- (n_rows - 1) + height + label_in / .ROW_IN + 0.4
  p <- ggplot2::ggplot(g) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = .data$.gene_xmin, xmax = .data$.gene_xmax,
                   ymin = .data$.y - height / 2, ymax = .data$.y + height / 2),
      fill = fill, colour = NA) +
    ggplot2::geom_text(
      ggplot2::aes(x = .data$.mid, y = .data$.label_y, label = .data$name),
      size = size, colour = "grey20", angle = angle, fontface = "italic",
      hjust = g$.hjust, vjust = 1) +
    ggplot2::scale_y_continuous(
      limits = c(height / 2 + 0.2 - span_units, height / 2 + 0.2), expand = ggplot2::expansion(0)) +
    ggplot2::coord_cartesian(xlim = xlim, clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(t = 0, r = 12, b = 0, l = 6))
  attr(p, "track_in") <- span_units * .ROW_IN
  p
}

# Genes with their box extents on the cumulative axis, for .gene_track_panel().
.zoom_gene_boxes <- function(g, layout) {
  if (is.null(g) || !nrow(g)) return(NULL)
  off <- layout$offset[match(normalise_chr(g$chr), layout$chr)]
  g$.gene_xmin <- off + as.numeric(g$start)
  g$.gene_xmax <- off + as.numeric(g$end)
  g
}

# Stack the gene track under a zoomed plot, sharing the x axis: the main panel keeps the
# tick labels, the track just carries the boxes and names.
.stack_gene_track <- function(p, z, layout, n_panels = 1L, angle = 0, min_width = TRUE) {
  track <- .gene_track_panel(.zoom_gene_boxes(z$track, layout), z$xlim, angle = angle,
                             width_in = .ZOOM_WIDTH_IN, min_width = min_width)
  if (is.null(track)) return(p)
  .need_package("patchwork", "the zoomed gene track")
  track_in <- attr(track, "track_in")
  # an absolute height for the strip: the panel above takes whatever is left, so the labels
  # get exactly the inches they were measured for however tall the page ends up
  out <- patchwork::wrap_plots(p, track, ncol = 1) +
    patchwork::plot_layout(heights = grid::unit.c(grid::unit(1, "null"),
                                                  grid::unit(track_in, "in")))
  attr(out, "plasgenomics_track_in") <- track_in
  out
}
