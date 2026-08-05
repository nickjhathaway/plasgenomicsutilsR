# Genome-wide IBD / selection plots over an IbdResults object. Each returns a
# ggplot object (composable / saveable by the caller); ggplot2 and scales are
# optional (Suggests), guarded at call time.

# ---- shared plotting pieces ------------------------------------------------

# alternating grey chromosome bands behind a genome-wide track
.chr_band_layer <- function(layout) {
  bands <- layout[layout$band == "a", , drop = FALSE]
  ggplot2::geom_rect(
    data = bands,
    ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = -Inf, ymax = Inf),
    fill = "#ebebeb", inherit.aes = FALSE
  )
}

# chromosome axis: ticks at each chromosome mid-point, labelled by number
.chr_axis <- function(layout) {
  ggplot2::scale_x_continuous(
    breaks = layout$mid, labels = layout$chr,
    expand = ggplot2::expansion(mult = 0.005)
  )
}

# Validate a highlight_genes / genes request against the track: error (rather than
# silently drawing nothing) when a requested name is not present, so a typo or a gene
# that was never added to the track surfaces instead of vanishing.
.check_gene_request <- function(gtrack, which) {
  if (is.null(which)) return(invisible())
  if (is.null(gtrack) || !nrow(gtrack)) {
    stop("gene(s) requested but this IbdResults has no gene track; ",
         "pass genes= to ibd_results()", call. = FALSE)
  }
  miss <- which[!tolower(which) %in% tolower(gtrack$name)]
  if (length(miss)) {
    stop("gene(s) not in the track: ", paste(miss, collapse = ", "),
         ".\n  Available: ", paste(gtrack$name, collapse = ", "), call. = FALSE)
  }
  invisible()
}

# Genes restricted to a (possibly re-laid-out) layout, with `cum_mid` recomputed so
# reference lines land correctly after skipping chromosomes. `which` optionally keeps
# only named genes to highlight.
.genes_for_layout <- function(genes, layout, which = NULL) {
  if (is.null(genes)) return(NULL)
  if (!is.null(which)) genes <- genes[tolower(genes$name) %in% tolower(which), , drop = FALSE]
  g <- normalise_chr(genes$chr)
  genes <- genes[g %in% layout$chr, , drop = FALSE]
  if (!nrow(genes)) return(NULL)
  g <- normalise_chr(genes$chr)
  genes$cum_mid <- layout$offset[match(g, layout$chr)] +
    (as.numeric(genes$start) + as.numeric(genes$end)) / 2
  genes
}

# vertical reference lines at annotated genes (optional)
.gene_line_layer <- function(genes) {
  if (is.null(genes) || !nrow(genes)) return(NULL)
  ggplot2::geom_vline(
    data = genes, ggplot2::aes(xintercept = .data$cum_mid),
    colour = "grey70", alpha = 0.7, linewidth = 0.3, inherit.aes = FALSE
  )
}

.GENE_LABEL_SIZE <- 2.7   # geom_text size for gene labels (also drives overlap sensing)

# Estimate a gene label's horizontal footprint in base pairs, so overlap detection
# reflects the *rendered* label width rather than a flat fraction of the genome. The
# label is drawn at 45 deg (hjust = 0), so its horizontal extent is the rotated text
# width; we map that from the intended plot width (mirrors .dims_genome) and font size.
.label_footprint_bp <- function(nchar_vec, layout, genome_frac, size = .GENE_LABEL_SIZE) {
  span_bp   <- diff(range(c(layout$xmin, layout$xmax)))
  plot_w_in <- max(6, 16 * genome_frac)          # .GENOME_FULL_WIDTH * frac, floored at 6
  usable_in <- max(2, plot_w_in - 1.2)            # drop axis + margins
  bp_per_in <- span_bp / usable_in
  char_in   <- size * 2.845 / 72 * 0.5            # size -> pt -> in, ~half-em per char
  nchar_vec * char_in * cos(pi / 4) * bp_per_in   # 45 deg -> horizontal component
}

# gene-name labels just above the panel (drawn in the top facet only when
# `facet_col`/`top_level` are given, so names aren't repeated down every panel).
# Every label sits at its TRUE genomic x (the gene mid-point, matching its reference
# line); a label that would overlap the one before it (its estimated text footprint
# reaching past the next label's anchor) is lifted onto a higher vertical tier so they
# clear each other without ever being moved sideways (which would misrepresent the
# position). `footprint` is the per-label width in bp from .label_footprint_bp().
# Needs `coord_cartesian(clip = "off")` + a top margin to sit outside.
.gene_label_layer <- function(genes, facet_col = NULL, top_level = NULL, footprint = NULL) {
  if (is.null(genes) || !nrow(genes) || !"name" %in% names(genes)) return(NULL)
  ord <- order(genes$cum_mid)
  g <- genes[ord, , drop = FALSE]
  fp <- if (is.null(footprint)) rep(0, nrow(g)) else footprint[ord]
  g$.tier <- 0L
  if (nrow(g) >= 2) {
    last_x <- numeric(0); last_fp <- numeric(0)   # last x + its footprint, per tier
    for (i in seq_len(nrow(g))) {
      t <- 0L
      # collide when the previous label on this tier extends right past label i's anchor
      while (t < length(last_x) && g$cum_mid[i] - last_x[t + 1] < last_fp[t + 1]) t <- t + 1L
      g$.tier[i] <- t
      last_x[t + 1] <- g$cum_mid[i]
      last_fp[t + 1] <- fp[i]
    }
  }
  if (!is.null(facet_col) && !is.null(top_level)) g[[facet_col]] <- top_level
  lapply(sort(unique(g$.tier)), function(t) {
    ggplot2::geom_text(
      data = g[g$.tier == t, , drop = FALSE],
      ggplot2::aes(x = .data$cum_mid, y = Inf, label = .data$name),
      inherit.aes = FALSE, vjust = -(0.4 + t * 1.3), hjust = 0, angle = 45,
      size = .GENE_LABEL_SIZE, colour = "grey20")
  })
}

# room above the panel for outside gene labels (only when labelling)
.gene_label_space <- function(label_genes) {
  if (!label_genes) return(NULL)
  list(ggplot2::coord_cartesian(clip = "off"),
       ggplot2::theme(plot.margin = ggplot2::margin(t = 42, r = 12, b = 6, l = 6)))
}

.manhattan_theme <- function() {
  # No grid lines: they read oddly against the alternating grey chromosome bands
  # (visible on the white bands, hidden on the grey). The bands carry the x reference.
  ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = "none",
      strip.background = ggplot2::element_rect(fill = "grey95", colour = NA)
    )
}

# two-tone colour by chromosome parity (the classic Manhattan look)
.band_colours <- function(colours) {
  if (is.null(colours)) c(a = "#3B4CC0", b = "#7E9BD8") else colours
}

# Restrict the chromosome layout to a kept set and re-lay it out contiguously (so
# skipped chromosomes leave no gap). `chroms` keeps only those; `skip_chr` drops
# those; both accept any chromosome spelling (normalised).
.select_layout <- function(layout, chroms = NULL, skip_chr = NULL) {
  keep <- layout$chr
  if (!is.null(chroms)) keep <- keep[keep %in% normalise_chr(chroms)]
  if (!is.null(skip_chr)) keep <- setdiff(keep, normalise_chr(skip_chr))
  if (!length(keep)) stop("no chromosomes left after chroms/skip_chr selection", call. = FALSE)
  sub <- layout[layout$chr %in% keep, , drop = FALSE]
  sub <- sub[order(match(sub$chr, layout$chr)), , drop = FALSE]
  sub$offset <- cumsum(c(0, sub$len))[seq_len(nrow(sub))]
  sub$xmin <- sub$offset
  sub$xmax <- sub$offset + sub$len
  sub$mid <- sub$offset + sub$len / 2
  sub$band <- rep(c("a", "b"), length.out = nrow(sub))
  sub
}

# Recompute `cum_pos` for a table against a (possibly re-laid-out) layout, keeping
# only rows on the layout's chromosomes.
.recum <- function(df, layout) {
  df <- df[df$chr %in% layout$chr, , drop = FALSE]
  df$cum_pos <- layout$offset[match(df$chr, layout$chr)] + as.numeric(df$pos)
  df
}

.filter_group <- function(df, groups) {
  if (is.null(groups) || is.null(df) || !"group" %in% names(df)) return(df)
  df[df$group %in% groups, , drop = FALSE]
}

.attach_band <- function(df, layout) {
  df$.band <- layout$band[match(df$chr, layout$chr)]
  df
}

# white -> pink -> magenta sequential ramp (RdPu-like), the default for the pairwise
# IBD-sharing tile plots (group x group heatmap, drug-gene triangles)
.IBD_FILL_DEFAULT <- c("white", "#fde0dd", "#fa9fb5", "#c51b8a", "#7a0177")

# a log-family fill transform can't take 0 (log2(0) = -Inf, "introduced infinite values");
# map non-positive fills to NA so they land on na.value (the lightest colour) cleanly.
.log_trans <- function(trans) is.character(trans) && length(trans) == 1L && grepl("^log", trans)
.na_nonpositive <- function(v, trans) {
  if (.log_trans(trans)) v[!is.na(v) & v <= 0] <- NA
  v
}

# thin lines at internal chromosome boundaries, to sharpen where one chromosome ends
.chr_boundary_layer <- function(layout, colour = "grey55", linewidth = 0.25) {
  bnd <- layout$offset[-1]
  if (!length(bnd)) return(NULL)
  ggplot2::geom_vline(xintercept = bnd, colour = colour, linewidth = linewidth)
}

# Build the fill scale shared by the tile plots. `fill_scale` fully overrides;
# otherwise a `scale_fill_gradientn` over `colors` with an optional `trans`
# (e.g. "log2", "sqrt") and `limits` (values outside are squished into range, so a
# few extremes near 1 stop crushing the rest of the scale).
.ibd_fill_scale <- function(name = "pairs IBD", trans = "identity", colors = NULL,
                            limits = NULL, fill_scale = NULL) {
  if (!is.null(fill_scale)) return(fill_scale)
  .need_package("scales", "the IBD fill scale")
  if (is.null(colors)) colors <- .IBD_FILL_DEFAULT
  ggplot2::scale_fill_gradientn(
    colours = colors, trans = trans, limits = limits, name = name,
    na.value = colors[1],
    oob = if (is.null(limits)) scales::censor else scales::squish
  )
}

# ---- Manhattan: per-SNP IBD fraction ---------------------------------------

#' IBD Manhattan plot
#'
#' Per-SNP fraction of pairs IBD along the genome. If the table carries a
#' `group` column the plot is faceted one row per group.
#'
#' @param x An [IbdResults] object.
#' @param groups Optional character vector to keep (needs a `group` column).
#' @param chroms Optional chromosomes to keep (any spelling); others are dropped
#'   and the remaining ones re-laid-out contiguously.
#' @param skip_chr Optional chromosomes to drop (complement of `chroms`).
#' @param highlight_genes Optional gene names (from the object's `genes` track) to
#'   draw as reference lines; default all genes in the track. Requesting a name not in
#'   the track is an error (rather than silently drawing nothing).
#' @param label_genes Label the genes with their names. `NULL` (default) labels them
#'   only when `highlight_genes` is given; `TRUE`/`FALSE` forces it.
#' @param point_size,point_alpha Point aesthetics.
#' @param colours Optional length-2 colour vector for the alternating chromosome bands.
#' @return A ggplot object.
#' @export
plot_ibd_sharing_manhattan <- function(x, groups = NULL, chroms = NULL, skip_chr = NULL,
                               highlight_genes = NULL, label_genes = NULL,
                               point_size = 0.5, point_alpha = 0.6, colours = NULL) {
  .need_package("ggplot2", "plot_ibd_sharing_manhattan()")
  .need_package("scales", "plot_ibd_sharing_manhattan()")
  df <- x$get_per_snp_group()
  if (is.null(df)) stop("this IbdResults has no per_snp_group table", call. = FALSE)
  .check_gene_request(x$get_genes(), highlight_genes)
  if (is.null(label_genes)) label_genes <- !is.null(highlight_genes)
  layout <- .select_layout(x$chrom_layout(), chroms, skip_chr)
  df <- .attach_band(.recum(.filter_group(df, groups), layout), layout)
  df <- df[is.finite(df$frac_pairs_ibd), , drop = FALSE]
  genes <- .genes_for_layout(x$get_genes(), layout, highlight_genes)
  faceted <- "group" %in% names(df)
  top_level <- if (faceted) sort(unique(df$group))[1] else NULL
  genome_frac <- sum(layout$len) / sum(x$chrom_layout()$len)

  p <- ggplot2::ggplot(df, ggplot2::aes(.data$cum_pos, .data$frac_pairs_ibd)) +
    .chr_band_layer(layout) +
    .gene_line_layer(genes) +
    (if (label_genes) .gene_label_layer(genes, if (faceted) "group" else NULL, top_level,
       footprint = .label_footprint_bp(nchar(genes$name), layout, genome_frac))) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$.band), size = point_size,
                        alpha = point_alpha, stroke = 0) +
    ggplot2::scale_colour_manual(values = .band_colours(colours)) +
    .chr_axis(layout) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(x = "chromosome", y = "pairs IBD") +
    .manhattan_theme() +
    .gene_label_space(label_genes)
  if (faceted) {
    p <- p + ggplot2::facet_wrap(~ .data$group, ncol = 1, strip.position = "right")
  }
  attr(p, "plasgenomics_dims") <- .dims_genome(
    if (faceted) length(unique(df$group)) else 1L, genome_frac, label_genes = label_genes)
  p
}

# ---- Manhattan: IBD selection statistic ------------------------------------

#' IBD selection-statistic Manhattan plot
#'
#' The selection statistic along the genome, with the Bonferroni threshold drawn
#' as a dashed line when plotting `neg_log10_p` and thresholds are available.
#'
#' @param x An [IbdResults] object.
#' @param metric Which column to plot: `"neg_log10_p"` (default), `"chi2_stat"`,
#'   or `"z_score"`.
#' @param groups Optional character vector to keep (needs a `group` column).
#' @param chroms Optional chromosomes to keep (any spelling); others are dropped and
#'   the remaining ones re-laid-out contiguously.
#' @param skip_chr Optional chromosomes to drop (complement of `chroms`).
#' @param highlight_genes Optional gene names (from the object's `genes` track) to draw
#'   as reference lines; default all genes in the track. Requesting a name not in the
#'   track is an error.
#' @param label_genes Label the genes. `NULL` (default) labels them only when
#'   `highlight_genes` is given; `TRUE`/`FALSE` forces it.
#' @param draw_threshold Draw the significance threshold line (only for `neg_log10_p`).
#' @param point_size,point_alpha Point aesthetics.
#' @param colours Optional length-2 colour vector for the alternating chromosome bands.
#' @return A ggplot object.
#' @export
plot_selection_manhattan <- function(x, metric = c("neg_log10_p", "chi2_stat", "z_score"),
                                     groups = NULL, chroms = NULL, skip_chr = NULL,
                                     highlight_genes = NULL, label_genes = NULL,
                                     draw_threshold = TRUE,
                                     point_size = 0.5, point_alpha = 0.6, colours = NULL) {
  .need_package("ggplot2", "plot_selection_manhattan()")
  metric <- match.arg(metric)
  df <- x$get_selection()
  if (is.null(df)) stop("this IbdResults has no selection table", call. = FALSE)
  if (!metric %in% names(df)) stop(sprintf("selection table has no '%s' column", metric),
                                    call. = FALSE)
  .check_gene_request(x$get_genes(), highlight_genes)
  if (is.null(label_genes)) label_genes <- !is.null(highlight_genes)
  layout <- .select_layout(x$chrom_layout(), chroms, skip_chr)
  df <- .attach_band(.recum(.filter_group(df, groups), layout), layout)
  df$.y <- df[[metric]]
  df <- df[is.finite(df$.y), , drop = FALSE]
  genes <- .genes_for_layout(x$get_genes(), layout, highlight_genes)
  faceted <- "group" %in% names(df)
  top_level <- if (faceted) sort(unique(df$group))[1] else NULL
  genome_frac <- sum(layout$len) / sum(x$chrom_layout()$len)

  p <- ggplot2::ggplot(df, ggplot2::aes(.data$cum_pos, .data$.y)) +
    .chr_band_layer(layout) +
    .gene_line_layer(genes) +
    (if (label_genes) .gene_label_layer(genes, if (faceted) "group" else NULL, top_level,
       footprint = .label_footprint_bp(nchar(genes$name), layout, genome_frac))) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$.band), size = point_size,
                        alpha = point_alpha, stroke = 0) +
    ggplot2::scale_colour_manual(values = .band_colours(colours)) +
    .chr_axis(layout) +
    ggplot2::labs(x = "chromosome", y = metric) +
    .manhattan_theme() +
    .gene_label_space(label_genes)

  thr <- x$get_thresholds()
  if (draw_threshold && metric == "neg_log10_p" && !is.null(thr)) {
    # Drop non-finite thresholds up front: a group with no valid SNPs gets a
    # NaN threshold (written as an empty cell by `ibd_selection_statistic`), and
    # feeding NA to geom_hline draws nothing but warns ("Removed N rows ...").
    thr <- thr[is.finite(thr$threshold), , drop = FALSE]
    has_group <- "group" %in% names(thr) && any(!is.na(thr$group))
    if (has_group && "group" %in% names(df)) {
      thr <- thr[thr$group %in% df$group, , drop = FALSE]
      if (nrow(thr)) {
        p <- p + ggplot2::geom_hline(data = thr,
          ggplot2::aes(yintercept = .data$threshold), colour = "firebrick",
          linetype = "dashed", linewidth = 0.4)
      }
    } else if (nrow(thr)) {
      p <- p + ggplot2::geom_hline(yintercept = thr$threshold[1], colour = "firebrick",
        linetype = "dashed", linewidth = 0.4)
    }
  }
  if ("group" %in% names(df)) {
    p <- p + ggplot2::facet_wrap(~ .data$group, ncol = 1, strip.position = "right")
  }
  attr(p, "plasgenomics_dims") <- .dims_genome(
    if (faceted) length(unique(df$group)) else 1L, genome_frac, label_genes = label_genes)
  p
}

# ---- group x group IBD heatmap along the genome --------------------------

#' Group-by-group IBD heatmap along the genome
#'
#' Per-SNP IBD sharing between group pairs as tiles along the genome, one facet
#' per anchor group. The stored upper-triangle (`group_a <= group_b`) is
#' mirrored so every anchor shows all partner groups. Alternating grey chromosome
#' bands and thin boundary lines mark where each chromosome starts and ends, and the
#' x-axis spans the **full chromosome lengths** so un-genotyped (empty) regions are
#' visible rather than collapsed out.
#'
#' @param x An [IbdResults] object.
#' @param anchor Optional single group to show (one panel) instead of all.
#' @param chroms Optional chromosomes to keep (any spelling); others are dropped and
#'   the remaining ones re-laid-out contiguously.
#' @param skip_chr Optional chromosomes to drop (complement of `chroms`).
#' @param trans Fill-scale transform, e.g. `"identity"` (default), `"log2"`, `"sqrt"`.
#' @param colors Optional colour ramp for the fill (defaults to a single-hue
#'   light-to-dark sequential scale that stays readable when most values are near 0).
#' @param limits Optional `c(lo, hi)` fill limits; values outside are squished into
#'   range, so a few extremes near 1 don't crush the rest of the scale.
#' @param fill_scale Optional ggplot2 fill scale that fully overrides the above.
#' @param highlight_genes Optional gene names from the `genes` track to mark with lines;
#'   requesting a name not in the track is an error.
#' @param label_genes Label the genes above the top panel. `NULL` (default) labels them
#'   only when `highlight_genes` is given; `TRUE`/`FALSE` forces it.
#' @return A ggplot object.
#' @export
plot_ibd_pairwise_group_heatmap <- function(x, anchor = NULL, chroms = NULL, skip_chr = NULL,
                                    trans = "identity", colors = NULL, limits = NULL,
                                    fill_scale = NULL, highlight_genes = NULL,
                                    label_genes = NULL) {
  .need_package("ggplot2", "plot_ibd_pairwise_group_heatmap()")
  df <- x$get_pairwise_group()
  if (is.null(df)) stop("this IbdResults has no pairwise_group table", call. = FALSE)
  .check_gene_request(x$get_genes(), highlight_genes)
  if (is.null(label_genes)) label_genes <- !is.null(highlight_genes)
  layout <- .select_layout(x$chrom_layout(), chroms, skip_chr)

  # mirror the upper triangle so each anchor group sees every partner
  swapped <- df
  swapped$group_a <- df$group_b
  swapped$group_b <- df$group_a
  full <- unique(rbind(df, swapped[swapped$group_a != swapped$group_b, , drop = FALSE]))
  if (!is.null(anchor)) full <- full[full$group_a %in% anchor, , drop = FALSE]
  full <- .recum(full, layout)
  full$frac_pairs_ibd <- .na_nonpositive(full$frac_pairs_ibd, trans)   # 0 -> NA for log fills
  genes <- .genes_for_layout(x$get_genes(), layout, highlight_genes)
  top_a <- sort(unique(full$group_a))[1]
  genome_frac <- sum(layout$len) / sum(x$chrom_layout()$len)

  spacing <- stats::median(diff(sort(unique(full$cum_pos))))
  tile_w <- if (is.finite(spacing) && spacing > 0) spacing * 1.5 else 1

  p <- ggplot2::ggplot(full, ggplot2::aes(.data$cum_pos, .data$group_b,
                                          fill = .data$frac_pairs_ibd)) +
    .chr_band_layer(layout) +                          # alternating grey chromosome bands
    ggplot2::geom_tile(width = tile_w, height = 0.9) +
    .chr_boundary_layer(layout) +                      # + a thin line at each chromosome edge
    .gene_line_layer(genes) +
    (if (label_genes) .gene_label_layer(genes, "group_a", top_a,
       footprint = .label_footprint_bp(nchar(genes$name), layout, genome_frac))) +
    .ibd_fill_scale("pairs IBD", trans = trans, colors = colors, limits = limits,
                    fill_scale = fill_scale) +
    .chr_axis(layout) +
    # span the full chromosome lengths so empty (un-genotyped) regions are shown too
    ggplot2::expand_limits(x = c(0, max(layout$xmax))) +
    ggplot2::labs(x = "chromosome", y = "group") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   strip.background = ggplot2::element_rect(fill = "grey95", colour = NA)) +
    .gene_label_space(label_genes)
  if (is.null(anchor)) {
    p <- p + ggplot2::facet_wrap(~ .data$group_a, ncol = 1, strip.position = "right")
  }
  attr(p, "plasgenomics_dims") <- .dims_heatmap(
    length(unique(full$group_a)), length(unique(c(full$group_a, full$group_b))),
    genome_frac, label_genes = label_genes)
  p
}

# ---- "tug of war" mirror plot ----------------------------------------------

#' IBD / selection "tug-of-war" mirror plot
#'
#' The selection statistic hangs from the top (bars descending) while the per-SNP
#' IBD fraction rises from the bottom (bars ascending), sharing one centred axis so
#' peaks of each line up. A single left axis carries both halves, its tick labels
#' tinted to their track (selection on top, IBD on the bottom).
#'
#' @param x An [IbdResults] object (needs both selection and per_snp_group tables).
#' @param group Group(s) to plot. `NULL` (default) plots every group, faceted one
#'   per row; a single group gives one panel; a vector facets those groups.
#' @param metric Selection metric for the top track (default `"neg_log10_p"`).
#' @param scale `"common"` (default) scales every panel to a shared maximum so they
#'   are directly comparable; `"free"` scales each group to its own maximum for more
#'   per-group detail (the axis then reads as within-group percentages).
#' @param chroms Optional chromosomes to keep (any spelling); others are dropped and
#'   the rest re-laid-out contiguously.
#' @param skip_chr Optional chromosomes to drop (complement of `chroms`).
#' @param highlight_genes Optional gene names from the `genes` track to mark with lines;
#'   requesting a name not in the track is an error.
#' @param label_genes Label the genes. `NULL` (default) labels them only when
#'   `highlight_genes` is given; `TRUE`/`FALSE` forces it.
#' @param draw_threshold Draw the per-group significance threshold (only for
#'   `neg_log10_p` with `scale = "common"`, when thresholds are available).
#' @param selection_colour,ibd_colour Bar and axis colours for the two tracks.
#' @return A ggplot object.
#' @export
plot_ibd_tugofwar <- function(x, group = NULL, metric = "neg_log10_p",
                              scale = c("common", "free"), chroms = NULL, skip_chr = NULL,
                              highlight_genes = NULL, label_genes = NULL,
                              draw_threshold = TRUE,
                              selection_colour = "#fd8d3c", ibd_colour = "#2166ac") {
  .need_package("ggplot2", "plot_ibd_tugofwar()")
  .need_package("scales", "plot_ibd_tugofwar()")
  scale <- match.arg(scale)
  sel <- x$get_selection()
  ibd <- x$get_per_snp_group()
  if (is.null(sel) || is.null(ibd)) {
    stop("plot_ibd_tugofwar() needs both selection and per_snp_group tables", call. = FALSE)
  }
  if (!metric %in% names(sel)) stop(sprintf("selection table has no '%s' column", metric),
                                    call. = FALSE)
  .check_gene_request(x$get_genes(), highlight_genes)
  if (is.null(label_genes)) label_genes <- !is.null(highlight_genes)
  regs <- .resolve_groups(group, sel, ibd)
  if (!is.null(regs)) {
    if ("group" %in% names(sel)) sel <- sel[sel$group %in% regs, , drop = FALSE]
    if ("group" %in% names(ibd)) ibd <- ibd[ibd$group %in% regs, , drop = FALSE]
  }
  layout <- .select_layout(x$chrom_layout(), chroms, skip_chr)
  sel <- .recum(sel[is.finite(sel[[metric]]), , drop = FALSE], layout)
  ibd <- .recum(ibd[is.finite(ibd$frac_pairs_ibd), , drop = FALSE], layout)
  genes <- .genes_for_layout(x$get_genes(), layout, highlight_genes)
  genome_frac <- sum(layout$len) / sum(x$chrom_layout()$len)

  present <- sort(unique(c(sel$group, ibd$group)))
  faceted <- ("group" %in% names(sel)) && length(present) > 1
  top_level <- if (faceted) present[1] else NULL

  # normalise each track so its tallest bar reaches the centre (y = 0). common =
  # one shared max (panels comparable); free = each group to its own max.
  normalized <- scale == "free" && "group" %in% names(sel)
  if (normalized) {
    sm <- tapply(sel[[metric]], sel$group, max, na.rm = TRUE)
    im <- tapply(ibd$frac_pairs_ibd, ibd$group, max, na.rm = TRUE)
    sm[!is.finite(sm) | sm <= 0] <- 1
    im[!is.finite(im) | im <= 0] <- 1
    sel$.tip <- 1 - sel[[metric]] / sm[as.character(sel$group)]
    ibd$.tip <- -1 + ibd$frac_pairs_ibd / im[as.character(ibd$group)]
    sel_max <- 1; ibd_max <- 1
  } else {
    sel_max <- max(sel[[metric]], na.rm = TRUE)
    ibd_max <- max(ibd$frac_pairs_ibd, na.rm = TRUE)
    if (!is.finite(sel_max) || sel_max <= 0) sel_max <- 1
    if (!is.finite(ibd_max) || ibd_max <= 0) ibd_max <- 1
    sel$.tip <- 1 - .safe_scale(sel[[metric]], sel_max)
    ibd$.tip <- -1 + .safe_scale(ibd$frac_pairs_ibd, ibd_max)
  }

  # single left axis: selection breaks in the top half, IBD in the bottom half,
  # each tick label tinted to its track.
  sel_vals <- pretty(c(0, sel_max), 4); sel_vals <- sel_vals[sel_vals >= 0 & sel_vals <= sel_max]
  ibd_vals <- pretty(c(0, ibd_max), 4); ibd_vals <- ibd_vals[ibd_vals >= 0 & ibd_vals <= ibd_max]
  sel_lab <- if (normalized) scales::percent(sel_vals, accuracy = 1) else as.character(sel_vals)
  yax <- data.frame(
    y   = c(1 - sel_vals / sel_max, -1 + ibd_vals / ibd_max),
    lab = c(sel_lab, scales::percent(ibd_vals, accuracy = 1)),
    col = c(rep(selection_colour, length(sel_vals)), rep(ibd_colour, length(ibd_vals))),
    stringsAsFactors = FALSE)
  yax <- yax[order(yax$y), ]

  use_ggtext <- requireNamespace("ggtext", quietly = TRUE)
  if (use_ggtext) {
    # colour each tick label via markdown so we avoid vectorised element_text()
    yax$lab <- paste0("<span style='color:", yax$col, ";'>", yax$lab, "</span>")
    ytitle <- paste0(
      "<span style='color:", selection_colour, ";'>selection (", metric, ", top)</span>",
      " / <span style='color:", ibd_colour, ";'>IBD fraction (bottom)</span>")
    ytext_elem  <- ggtext::element_markdown()
    ytitle_elem <- ggtext::element_markdown(angle = 90)
  } else {
    ytitle <- paste0("selection ", metric, " (top)  /  IBD fraction (bottom)")
    ytext_elem  <- ggplot2::element_text()
    ytitle_elem <- ggplot2::element_text()
  }

  thr_layer <- NULL
  if (draw_threshold && metric == "neg_log10_p" && !normalized) {
    thr <- x$get_thresholds()
    if (!is.null(thr)) {
      thr <- thr[is.finite(thr$threshold), , drop = FALSE]   # drop NaN (no-valid-SNP) groups
      if ("group" %in% names(thr) && any(!is.na(thr$group)) && "group" %in% names(sel)) {
        thr <- thr[thr$group %in% present, , drop = FALSE]
      }
      if (nrow(thr)) {
        thr$.y <- 1 - thr$threshold / sel_max
        thr_layer <- ggplot2::geom_hline(data = thr, ggplot2::aes(yintercept = .data$.y),
          colour = "firebrick", linetype = "dashed", linewidth = 0.4, inherit.aes = FALSE)
      }
    }
  }

  p <- ggplot2::ggplot() +
    .chr_band_layer(layout) +
    .gene_line_layer(genes) +
    (if (label_genes) .gene_label_layer(genes, if (faceted) "group" else NULL, top_level,
       footprint = .label_footprint_bp(nchar(genes$name), layout, genome_frac))) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
    ggplot2::geom_segment(data = sel,
      ggplot2::aes(x = .data$cum_pos, xend = .data$cum_pos, y = 1, yend = .data$.tip),
      colour = selection_colour, linewidth = 0.15) +
    ggplot2::geom_segment(data = ibd,
      ggplot2::aes(x = .data$cum_pos, xend = .data$cum_pos, y = -1, yend = .data$.tip),
      colour = ibd_colour, linewidth = 0.15) +
    thr_layer +
    .chr_axis(layout) +
    ggplot2::scale_y_continuous(name = ytitle, limits = c(-1, 1),
      breaks = yax$y, labels = yax$lab,
      expand = ggplot2::expansion(mult = c(0.02, 0.08))) +
    ggplot2::labs(x = "chromosome") +
    .manhattan_theme() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ytext_elem,
      axis.title.y = ytitle_elem) +
    .gene_label_space(label_genes)
  if (faceted) {
    p <- p + ggplot2::facet_wrap(~ group, ncol = 1, strip.position = "right")
  } else if (!is.null(regs)) {
    p <- p + ggplot2::labs(title = paste0("tug-of-war: ", paste(regs, collapse = ", ")))
  }
  attr(p, "plasgenomics_dims") <- .dims_genome(
    if (faceted) length(present) else 1L, genome_frac, label_genes = label_genes)
  p
}

# ---- drug-gene group x group triangles -----------------------------------

# SNPs from a pairwise-group table overlapping one feature interval
.snps_in_gene <- function(pw, chr, start, end) {
  pw[pw$chr == chr & pw$pos >= start & pw$pos <= end, , drop = FALSE]
}

# Turn user SNP identifiers into single-position features: either "chr:pos" ids or a
# data frame with chr + pos (+ optional name).
.parse_snp_features <- function(snps) {
  if (is.data.frame(snps)) {
    if (!all(c("chr", "pos") %in% names(snps)))
      stop("`snps` data frame needs 'chr' and 'pos' columns", call. = FALSE)
    nm <- if ("name" %in% names(snps)) as.character(snps$name) else paste0(snps$chr, ":", snps$pos)
    return(data.frame(name = nm, chr = normalise_chr(snps$chr),
                      start = as.numeric(snps$pos), end = as.numeric(snps$pos),
                      stringsAsFactors = FALSE))
  }
  ids <- as.character(snps)
  parts <- strsplit(ids, ":", fixed = TRUE)
  chr <- vapply(parts, function(z) paste(z[-length(z)], collapse = ":"), "")
  pos <- suppressWarnings(as.numeric(vapply(parts, function(z) z[length(z)], "")))
  if (any(is.na(pos)))
    stop("`snps` must be 'chr:pos' ids (e.g. 'Pf3D7_07_v3:403222') or a chr/pos data frame",
         call. = FALSE)
  data.frame(name = ids, chr = normalise_chr(chr), start = pos, end = pos,
             stringsAsFactors = FALSE)
}

# per-gene / per-SNP triangle values from the pairwise per-SNP table (the SNP-in-feature
# path: a gene's value aggregates the pairwise IBD of SNPs strictly inside it). Used for
# `snps=` and as the fallback when no IBD blocks / overlap table are available.
.snp_gene_triangle_df <- function(x, genes, snps, agg_fn) {
  pw <- x$get_pairwise_group()
  if (is.null(pw)) stop("this IbdResults has no pairwise_group table", call. = FALSE)
  if (!is.null(genes)) .check_gene_request(x$get_genes(), genes)
  feats <- list()
  gtrack <- x$get_genes()
  if (!is.null(gtrack) && (is.null(snps) || !is.null(genes))) {
    if (!is.null(genes)) gtrack <- gtrack[tolower(gtrack$name) %in% tolower(genes), , drop = FALSE]
    if (nrow(gtrack)) feats[["genes"]] <- data.frame(
      name = gtrack$name, chr = normalise_chr(gtrack$chr),
      start = as.numeric(gtrack$start), end = as.numeric(gtrack$end),
      stringsAsFactors = FALSE)
  }
  if (!is.null(snps)) feats[["snps"]] <- .parse_snp_features(snps)
  feats <- do.call(rbind, feats)
  if (is.null(feats) || !nrow(feats)) {
    stop("nothing to plot: give a genes track (ibd_results(genes=)) and/or snps=", call. = FALSE)
  }
  pw$chr <- normalise_chr(pw$chr)
  parts <- list()
  for (i in seq_len(nrow(feats))) {
    f <- feats[i, ]
    sub <- .snps_in_gene(pw, f$chr, f$start, f$end)
    if (!nrow(sub)) {
      warning(sprintf("feature '%s' has no overlapping SNPs; skipped", f$name), call. = FALSE)
      next
    }
    ag <- stats::aggregate(frac_pairs_ibd ~ group_a + group_b, data = sub, FUN = agg_fn)
    ag$gene <- f$name
    parts[[length(parts) + 1]] <- ag
  }
  if (!length(parts)) stop("none of the requested features had overlapping SNPs", call. = FALSE)
  df <- do.call(rbind, parts)
  df$gene <- factor(df$gene, levels = feats$name[feats$name %in% df$gene])
  df
}

# per-gene triangle values from a precomputed block-overlap table (Python ibd_gene_overlap)
.precomputed_gene_df <- function(tab, gene_names) {
  .require_cols(tab, c("gene", "group_a", "group_b", "frac_pairs_ibd"), "gene_overlap")
  if (!is.null(gene_names)) {
    have <- unique(as.character(tab$gene))
    miss <- gene_names[!tolower(gene_names) %in% tolower(have)]
    if (length(miss)) {
      stop("gene(s) not in the gene_overlap table: ", paste(miss, collapse = ", "), call. = FALSE)
    }
    tab <- tab[tolower(tab$gene) %in% tolower(gene_names), , drop = FALSE]
  }
  # disambiguate repeated gene Names (via gene_id when the table carries it) so distinct
  # genes that share a Name don't collapse into one facet
  gid <- if ("gene_id" %in% names(tab)) as.character(tab$gene_id) else NULL
  key <- if (is.null(gid)) as.character(tab$gene) else paste(tab$gene, gid, sep = "\r")
  feat <- !duplicated(key)
  labels <- .disambiguate_gene_labels(as.character(tab$gene)[feat],
                                      if (is.null(gid)) NULL else gid[feat])
  lab_of <- stats::setNames(labels, key[feat])
  data.frame(gene = factor(unname(lab_of[key]), levels = labels),
             group_a = tab$group_a, group_b = tab$group_b,
             frac_pairs_ibd = tab$frac_pairs_ibd, stringsAsFactors = FALSE)
}

#' IBD "triangle" panels for genes or specific SNPs
#'
#' Draws a group-by-group IBD-sharing triangle per feature: a **gene** (per-SNP
#' pairwise IBD for SNPs falling strictly inside the gene interval, aggregated) or a
#' **specific SNP** (the sharing at that single position). One facet per feature -- use
#' it to ask whether a gene or locus is itself shared between groups.
#'
#' Genes come from the `genes` track on the [IbdResults] object and membership is
#' strict: a SNP belongs to a gene only when its position is inside the interval (no
#' flanking). Target individual loci with `snps` when a specific variant matters more
#' than a whole gene.
#'
#' @details **How a gene's sharing is measured.** When the object carries IBD blocks
#' (`ibd_results(blocks=, meta=)`) or a precomputed overlap table (`gene_overlap=`, from
#' `plasgenomicsutils ibd_gene_overlap`), a gene's cell is the fraction of pairs whose IBD
#' **block overlaps the gene interval** -- so a pair counts when it shares a segment
#' spanning the gene even with no genotyped SNP inside it (see [gene_ibd_overlap()]).
#' Without blocks it falls back to aggregating the pairwise IBD of SNPs strictly inside the
#' gene. `snps=` always uses the per-SNP path (a single locus is a point, not an interval).
#'
#' @param x An [IbdResults] object. Gene triangles use its IBD `blocks` / `gene_overlap`
#'   table if present, else the `pairwise_group` per-SNP table.
#' @param genes Gene names to include from the track (default: all; case-insensitive).
#'   Ignored if only `snps` is given.
#' @param snps Specific SNPs to draw, as `"chr:pos"` ids (e.g. `"Pf3D7_07_v3:403222"`)
#'   or a data frame with `chr`, `pos` (and optional `name`). Each is one facet.
#' @param group For block-based overlap, the metadata column defining the groups
#'   (default: the first non-`sample` column of the object's `meta`).
#' @param within For block-based overlap, pad each gene interval by this many bp on both
#'   sides (default `0`).
#' @param agg For the SNP fallback, how a gene's value aggregates its in-gene SNPs:
#'   `"mean"` (default), `"median"`, or `"max"`.
#' @param individual If `TRUE`, return a **named list** of one-triangle plots (one per
#'   feature, e.g. to write a multi-page PDF) with the legend tucked into the empty
#'   upper triangle; if `FALSE` (default), one faceted grid plot.
#' @param label Draw the value in each tile.
#' @param digits Decimal places for the tile labels.
#' @param ncol Facet columns for the grid (default: ggplot2 chooses).
#' @param trans Fill-scale transform, e.g. `"identity"` (default), `"log2"`, `"sqrt"`.
#' @param colors Optional colour ramp for the fill (default: the pairwise-sharing ramp).
#' @param limits Fill limits. `NULL` (default) lets each plot scale to its own values;
#'   `"shared"` pins every feature to the range across all of them, which makes the pages
#'   from `individual = TRUE` colour-comparable; or give `c(lo, hi)` explicitly (values
#'   outside are squished). A faceted grid already shares one scale across its panels.
#' @param fill_scale Optional ggplot2 fill scale that fully overrides the above.
#' @return A ggplot object (grid), or a named list of ggplot objects when `individual`.
#' @export
plot_pairwise_ibd_for_genes <- function(x, genes = NULL, snps = NULL, group = NULL,
                                     within = 0, agg = c("mean", "median", "max"),
                                     individual = FALSE,
                                     label = TRUE, digits = 2, ncol = NULL,
                                     trans = "identity", colors = NULL, limits = NULL,
                                     fill_scale = NULL) {
  .need_package("ggplot2", "plot_pairwise_ibd_for_genes()")
  agg <- match.arg(agg)
  agg_fn <- switch(agg, mean = function(v) mean(v, na.rm = TRUE),
                   median = function(v) stats::median(v, na.rm = TRUE),
                   max = function(v) max(v, na.rm = TRUE))

  # gene triangles use block overlap when available (correct), SNPs always per-SNP
  use_blocks <- is.null(snps) && (!is.null(x$get_gene_overlap()) || !is.null(x$get_blocks()))
  if (use_blocks) {
    if (!is.null(x$get_gene_overlap())) {
      df <- .precomputed_gene_df(x$get_gene_overlap(), genes)
    } else {
      ov <- gene_ibd_overlap(x, genes = genes, group = group, within = within)
      df <- ov[, c("gene", "group_a", "group_b", "frac_pairs_ibd")]
      df$gene <- factor(as.character(df$gene), levels = levels(ov$gene))
    }
  } else {
    df <- .snp_gene_triangle_df(x, genes, snps, agg_fn)
  }

  # "shared" pins every feature to one scale, so the pages from individual = TRUE are
  # colour-comparable; taken over the values as drawn, so a log transform's dropped
  # non-positives do not pull the lower limit to -Inf
  if (is.character(limits)) {
    limits <- match.arg(limits, "shared")
    limits <- suppressWarnings(range(.na_nonpositive(df$frac_pairs_ibd, trans), na.rm = TRUE))
    if (!all(is.finite(limits)) || diff(limits) == 0) limits <- NULL
  }

  groups <- sort(unique(c(as.character(df$group_a), as.character(df$group_b))))
  fs <- .ibd_fill_scale("pairs IBD", trans = trans, colors = colors, limits = limits,
                        fill_scale = fill_scale)
  # colours used for the label-contrast luminance (fall back to the default ramp when a
  # full fill_scale override hides them)
  txt_cols <- if (is.null(colors)) .IBD_FILL_DEFAULT else colors
  feats_present <- levels(df$gene)

  if (individual) {
    plots <- lapply(feats_present, function(g) {
      pg <- .triangle_gg(df[df$gene == g, , drop = FALSE], groups, fs, label, digits,
                         title = g, legend_inside = TRUE,
                         colours = txt_cols, limits = limits, trans = trans)
      attr(pg, "plasgenomics_dims") <- .dims_triangles(1L, length(groups), 1L)
      pg
    })
    names(plots) <- feats_present
    return(plots)
  }

  p <- .triangle_gg(df, groups, fs, label, digits,
                    colours = txt_cols, limits = limits, trans = trans) +
    ggplot2::facet_wrap(~ .data$gene, ncol = ncol) +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"))
  attr(p, "plasgenomics_dims") <- .dims_triangles(length(feats_present), length(groups), ncol)
  p
}

# Readable label colour (white / dark) per value, from the *actual fill colour's*
# luminance rather than a value threshold -- so light tiles always get dark text even
# when their value is high relative to the (mostly-near-zero) rest of the data. Mirrors
# the fill scale's ramp / limits / transform so the luminance matches what's drawn.
.readable_text_colour <- function(values, colours, limits, trans = "identity") {
  .need_package("scales", "triangle label contrast")
  tfun <- switch(trans, sqrt = sqrt, log2 = log2, log10 = log10, log = log,
                 function(z) z)
  tv <- tfun(values); tl <- tfun(limits)
  from <- range(c(tv, tl)[is.finite(c(tv, tl))])
  r <- scales::rescale(tv, to = c(0, 1), from = from)
  r[!is.finite(r)] <- 0
  cols <- scales::gradient_n_pal(colours)(pmin(pmax(r, 0), 1))
  cols[is.na(cols)] <- colours[[1]]
  m <- grDevices::col2rgb(cols) / 255                       # relative luminance (WCAG)
  lin <- function(c) ifelse(c <= 0.03928, c / 12.92, ((c + 0.055) / 1.055)^2.4)
  L <- 0.2126 * lin(m[1, ]) + 0.7152 * lin(m[2, ]) + 0.0722 * lin(m[3, ])
  ifelse(L < 0.4, "white", "grey15")
}

# one group x group triangle for a single feature's aggregated data
.triangle_gg <- function(df, groups, fill_scale, label, digits, title = NULL,
                         legend_inside = FALSE, colours = .IBD_FILL_DEFAULT,
                         limits = NULL, trans = "identity") {
  df$group_a <- factor(df$group_a, levels = groups)
  df$group_b <- factor(df$group_b, levels = rev(groups))
  # fill maps a log-safe copy (0 -> NA), while the label keeps the true value
  df$.fill <- .na_nonpositive(df$frac_pairs_ibd, trans)
  if (label) {
    lim <- if (!is.null(limits)) limits else range(df$frac_pairs_ibd, na.rm = TRUE)
    df$.txt <- .readable_text_colour(df$frac_pairs_ibd, colours, lim, trans)
  }
  p <- ggplot2::ggplot(df, ggplot2::aes(.data$group_a, .data$group_b,
                                        fill = .data[[".fill"]])) +
    ggplot2::geom_tile(colour = "grey90") +
    fill_scale +
    ggplot2::coord_fixed() +
    ggplot2::labs(x = NULL, y = NULL, title = title) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5))
  if (legend_inside) p <- p + .legend_upper_triangle()
  if (label) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = formatC(.data$frac_pairs_ibd, format = "f", digits = digits),
                   colour = .data[[".txt"]]),
      size = 2.8, show.legend = FALSE) +
      ggplot2::scale_colour_identity()
  }
  p
}

# place the fill legend inside the empty upper-right triangle (ggplot2-version-safe)
.legend_upper_triangle <- function() {
  base <- ggplot2::theme(
    legend.background = ggplot2::element_rect(fill = "white", colour = NA),
    legend.key.size = grid::unit(0.9, "lines"),
    legend.title = ggplot2::element_text(size = 8),
    legend.text = ggplot2::element_text(size = 7))
  pos <- c(0.82, 0.80)
  if (utils::packageVersion("ggplot2") >= "3.5.0") {
    base + ggplot2::theme(legend.position = "inside", legend.position.inside = pos)
  } else {
    base + ggplot2::theme(legend.position = pos)
  }
}

.safe_scale <- function(v, mx) {
  if (!is.finite(mx) || mx <= 0) return(rep(0, length(v)))
  out <- v / mx
  out[!is.finite(out)] <- 0
  out
}

.resolve_groups <- function(group, sel, ibd) {
  if (!is.null(group)) return(group)
  regs <- unique(c(if ("group" %in% names(sel)) sel$group,
                   if ("group" %in% names(ibd)) ibd$group))
  if (length(regs)) regs else NULL
}
