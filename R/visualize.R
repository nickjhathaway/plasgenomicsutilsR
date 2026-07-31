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

# gene-name labels just above the panel (drawn in the top facet only when
# `facet_col`/`top_level` are given, so names aren't repeated down every panel).
# Every label sits at its TRUE genomic x (the gene mid-point, matching its reference
# line); labels closer than `min_gap` are lifted onto higher vertical tiers so they
# clear each other without ever being moved sideways (which would misrepresent the
# position). Needs `coord_cartesian(clip = "off")` + a top margin to sit outside.
.gene_label_layer <- function(genes, facet_col = NULL, top_level = NULL, min_gap = 0) {
  if (is.null(genes) || !nrow(genes) || !"name" %in% names(genes)) return(NULL)
  g <- genes[order(genes$cum_mid), , drop = FALSE]
  g$.tier <- 0L
  if (min_gap > 0 && nrow(g) >= 2) {
    last <- numeric(0)                       # last x placed on each tier
    for (i in seq_len(nrow(g))) {
      t <- 0L
      while (t < length(last) && g$cum_mid[i] - last[t + 1] < min_gap) t <- t + 1L
      g$.tier[i] <- t
      last[t + 1] <- g$cum_mid[i]
    }
  }
  if (!is.null(facet_col) && !is.null(top_level)) g[[facet_col]] <- top_level
  lapply(sort(unique(g$.tier)), function(t) {
    ggplot2::geom_text(
      data = g[g$.tier == t, , drop = FALSE],
      ggplot2::aes(x = .data$cum_mid, y = Inf, label = .data$name),
      inherit.aes = FALSE, vjust = -(0.4 + t * 1.3), hjust = 0, angle = 45,
      size = 2.7, colour = "grey20")
  })
}

# room above the panel for outside gene labels (only when labelling)
.gene_label_space <- function(label_genes) {
  if (!label_genes) return(NULL)
  list(ggplot2::coord_cartesian(clip = "off"),
       ggplot2::theme(plot.margin = ggplot2::margin(t = 42, r = 12, b = 6, l = 6)))
}

.manhattan_theme <- function() {
  ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
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

.filter_region <- function(df, regions) {
  if (is.null(regions) || is.null(df) || !"region" %in% names(df)) return(df)
  df[df$region %in% regions, , drop = FALSE]
}

.attach_band <- function(df, layout) {
  df$.band <- layout$band[match(df$chr, layout$chr)]
  df
}

# white -> pink -> magenta sequential ramp (RdPu-like), the default for the pairwise
# IBD-sharing tile plots (region x region heatmap, drug-gene triangles)
.IBD_FILL_DEFAULT <- c("white", "#fde0dd", "#fa9fb5", "#c51b8a", "#7a0177")

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
#' `region` column the plot is faceted one row per region.
#'
#' @param x An [IbdResults] object.
#' @param regions Optional character vector to keep (needs a `region` column).
#' @param chroms Optional chromosomes to keep (any spelling); others are dropped
#'   and the remaining ones re-laid-out contiguously.
#' @param skip_chr Optional chromosomes to drop (complement of `chroms`).
#' @param highlight_genes Optional gene names (from the object's `genes` track) to
#'   draw as reference lines; default all genes in the track.
#' @param label_genes Label the highlighted genes with their names.
#' @param point_size,point_alpha Point aesthetics.
#' @param colours Optional length-2 colour vector for the alternating chromosome bands.
#' @return A ggplot object.
#' @export
plot_ibd_manhattan <- function(x, regions = NULL, chroms = NULL, skip_chr = NULL,
                               highlight_genes = NULL, label_genes = FALSE,
                               point_size = 0.5, point_alpha = 0.6, colours = NULL) {
  .need_package("ggplot2", "plot_ibd_manhattan()")
  .need_package("scales", "plot_ibd_manhattan()")
  df <- x$get_per_snp_region()
  if (is.null(df)) stop("this IbdResults has no per_snp_region table", call. = FALSE)
  layout <- .select_layout(x$chrom_layout(), chroms, skip_chr)
  df <- .attach_band(.recum(.filter_region(df, regions), layout), layout)
  df <- df[is.finite(df$frac_pairs_ibd), , drop = FALSE]
  genes <- .genes_for_layout(x$get_genes(), layout, highlight_genes)
  faceted <- "region" %in% names(df)
  top_level <- if (faceted) sort(unique(df$region))[1] else NULL
  label_gap <- diff(range(c(layout$xmin, layout$xmax))) * 0.03

  p <- ggplot2::ggplot(df, ggplot2::aes(.data$cum_pos, .data$frac_pairs_ibd)) +
    .chr_band_layer(layout) +
    .gene_line_layer(genes) +
    (if (label_genes) .gene_label_layer(genes, if (faceted) "region" else NULL, top_level, label_gap)) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$.band), size = point_size,
                        alpha = point_alpha, stroke = 0) +
    ggplot2::scale_colour_manual(values = .band_colours(colours)) +
    .chr_axis(layout) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(x = "chromosome", y = "pairs IBD") +
    .manhattan_theme() +
    .gene_label_space(label_genes)
  if (faceted) {
    p <- p + ggplot2::facet_wrap(~ .data$region, ncol = 1, strip.position = "right")
  }
  attr(p, "plasgenomics_dims") <- .dims_genome(
    if (faceted) length(unique(df$region)) else 1L,
    sum(layout$len) / sum(x$chrom_layout()$len), label_genes = label_genes)
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
#' @param regions Optional character vector to keep (needs a `region` column).
#' @param draw_threshold Draw the significance threshold line (only for `neg_log10_p`).
#' @param point_size,point_alpha Point aesthetics.
#' @param colours Optional length-2 colour vector for the alternating chromosome bands.
#' @return A ggplot object.
#' @export
plot_selection_manhattan <- function(x, metric = c("neg_log10_p", "chi2_stat", "z_score"),
                                     regions = NULL, chroms = NULL, skip_chr = NULL,
                                     highlight_genes = NULL, label_genes = FALSE,
                                     draw_threshold = TRUE,
                                     point_size = 0.5, point_alpha = 0.6, colours = NULL) {
  .need_package("ggplot2", "plot_selection_manhattan()")
  metric <- match.arg(metric)
  df <- x$get_selection()
  if (is.null(df)) stop("this IbdResults has no selection table", call. = FALSE)
  if (!metric %in% names(df)) stop(sprintf("selection table has no '%s' column", metric),
                                    call. = FALSE)
  layout <- .select_layout(x$chrom_layout(), chroms, skip_chr)
  df <- .attach_band(.recum(.filter_region(df, regions), layout), layout)
  df$.y <- df[[metric]]
  df <- df[is.finite(df$.y), , drop = FALSE]
  genes <- .genes_for_layout(x$get_genes(), layout, highlight_genes)
  faceted <- "region" %in% names(df)
  top_level <- if (faceted) sort(unique(df$region))[1] else NULL
  label_gap <- diff(range(c(layout$xmin, layout$xmax))) * 0.03

  p <- ggplot2::ggplot(df, ggplot2::aes(.data$cum_pos, .data$.y)) +
    .chr_band_layer(layout) +
    .gene_line_layer(genes) +
    (if (label_genes) .gene_label_layer(genes, if (faceted) "region" else NULL, top_level, label_gap)) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$.band), size = point_size,
                        alpha = point_alpha, stroke = 0) +
    ggplot2::scale_colour_manual(values = .band_colours(colours)) +
    .chr_axis(layout) +
    ggplot2::labs(x = "chromosome", y = metric) +
    .manhattan_theme() +
    .gene_label_space(label_genes)

  thr <- x$get_thresholds()
  if (draw_threshold && metric == "neg_log10_p" && !is.null(thr)) {
    if (!is.na(thr$region[1]) && "region" %in% names(df)) {
      thr <- thr[thr$region %in% df$region, , drop = FALSE]
      p <- p + ggplot2::geom_hline(data = thr,
        ggplot2::aes(yintercept = .data$threshold), colour = "firebrick",
        linetype = "dashed", linewidth = 0.4)
    } else {
      p <- p + ggplot2::geom_hline(yintercept = thr$threshold[1], colour = "firebrick",
        linetype = "dashed", linewidth = 0.4)
    }
  }
  if ("region" %in% names(df)) {
    p <- p + ggplot2::facet_wrap(~ .data$region, ncol = 1, strip.position = "right")
  }
  attr(p, "plasgenomics_dims") <- .dims_genome(
    if (faceted) length(unique(df$region)) else 1L,
    sum(layout$len) / sum(x$chrom_layout()$len), label_genes = label_genes)
  p
}

# ---- region x region IBD heatmap along the genome --------------------------

#' Region-by-region IBD heatmap along the genome
#'
#' Per-SNP IBD sharing between region pairs as tiles along the genome, one facet
#' per anchor region. The stored upper-triangle (`region_a <= region_b`) is
#' mirrored so every anchor shows all partner regions.
#'
#' @param x An [IbdResults] object.
#' @param anchor Optional single region to show (one panel) instead of all.
#' @param trans Fill-scale transform, e.g. `"identity"` (default), `"log2"`, `"sqrt"`.
#' @param colors Optional colour ramp for the fill (defaults to a single-hue
#'   light-to-dark sequential scale that stays readable when most values are near 0).
#' @param limits Optional `c(lo, hi)` fill limits; values outside are squished into
#'   range, so a few extremes near 1 don't crush the rest of the scale.
#' @param fill_scale Optional ggplot2 fill scale that fully overrides the above.
#' @param highlight_genes Optional gene names from the `genes` track to mark with lines.
#' @param label_genes Label the highlighted genes (top panel only, outside the plot).
#' @return A ggplot object.
#' @export
plot_ibd_region_heatmap <- function(x, anchor = NULL, chroms = NULL, skip_chr = NULL,
                                    trans = "identity", colors = NULL, limits = NULL,
                                    fill_scale = NULL, highlight_genes = NULL,
                                    label_genes = FALSE) {
  .need_package("ggplot2", "plot_ibd_region_heatmap()")
  df <- x$get_pairwise_region()
  if (is.null(df)) stop("this IbdResults has no pairwise_region table", call. = FALSE)
  layout <- .select_layout(x$chrom_layout(), chroms, skip_chr)

  # mirror the upper triangle so each anchor region sees every partner
  swapped <- df
  swapped$region_a <- df$region_b
  swapped$region_b <- df$region_a
  full <- unique(rbind(df, swapped[swapped$region_a != swapped$region_b, , drop = FALSE]))
  if (!is.null(anchor)) full <- full[full$region_a %in% anchor, , drop = FALSE]
  full <- .recum(full, layout)
  genes <- .genes_for_layout(x$get_genes(), layout, highlight_genes)
  top_a <- sort(unique(full$region_a))[1]
  label_gap <- diff(range(c(layout$xmin, layout$xmax))) * 0.03

  spacing <- stats::median(diff(sort(unique(full$cum_pos))))
  tile_w <- if (is.finite(spacing) && spacing > 0) spacing * 1.5 else 1

  p <- ggplot2::ggplot(full, ggplot2::aes(.data$cum_pos, .data$region_b,
                                          fill = .data$frac_pairs_ibd)) +
    ggplot2::geom_tile(width = tile_w, height = 0.9) +
    .gene_line_layer(genes) +
    (if (label_genes) .gene_label_layer(genes, "region_a", top_a, label_gap)) +
    .ibd_fill_scale("pairs IBD", trans = trans, colors = colors, limits = limits,
                    fill_scale = fill_scale) +
    .chr_axis(layout) +
    ggplot2::labs(x = "chromosome", y = "region") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   strip.background = ggplot2::element_rect(fill = "grey95", colour = NA)) +
    .gene_label_space(label_genes)
  if (is.null(anchor)) {
    p <- p + ggplot2::facet_wrap(~ .data$region_a, ncol = 1, strip.position = "right")
  }
  attr(p, "plasgenomics_dims") <- .dims_heatmap(
    length(unique(full$region_a)), length(unique(c(full$region_a, full$region_b))),
    sum(layout$len) / sum(x$chrom_layout()$len), label_genes = label_genes)
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
#' @param x An [IbdResults] object (needs both selection and per_snp_region tables).
#' @param region Region(s) to plot. `NULL` (default) plots every region, faceted one
#'   per row; a single region gives one panel; a vector facets those regions.
#' @param metric Selection metric for the top track (default `"neg_log10_p"`).
#' @param scale `"common"` (default) scales every panel to a shared maximum so they
#'   are directly comparable; `"free"` scales each region to its own maximum for more
#'   per-region detail (the axis then reads as within-region percentages).
#' @param chroms Optional chromosomes to keep (any spelling); others are dropped and
#'   the rest re-laid-out contiguously.
#' @param skip_chr Optional chromosomes to drop (complement of `chroms`).
#' @param highlight_genes Optional gene names from the `genes` track to mark with lines.
#' @param label_genes Label the highlighted genes (top panel only, outside the plot).
#' @param draw_threshold Draw the per-region significance threshold (only for
#'   `neg_log10_p` with `scale = "common"`, when thresholds are available).
#' @param selection_colour,ibd_colour Bar and axis colours for the two tracks.
#' @return A ggplot object.
#' @export
plot_ibd_tugofwar <- function(x, region = NULL, metric = "neg_log10_p",
                              scale = c("common", "free"), chroms = NULL, skip_chr = NULL,
                              highlight_genes = NULL, label_genes = FALSE,
                              draw_threshold = TRUE,
                              selection_colour = "#fd8d3c", ibd_colour = "#2166ac") {
  .need_package("ggplot2", "plot_ibd_tugofwar()")
  .need_package("scales", "plot_ibd_tugofwar()")
  scale <- match.arg(scale)
  sel <- x$get_selection()
  ibd <- x$get_per_snp_region()
  if (is.null(sel) || is.null(ibd)) {
    stop("plot_ibd_tugofwar() needs both selection and per_snp_region tables", call. = FALSE)
  }
  if (!metric %in% names(sel)) stop(sprintf("selection table has no '%s' column", metric),
                                    call. = FALSE)
  regs <- .resolve_regions(region, sel, ibd)
  if (!is.null(regs)) {
    if ("region" %in% names(sel)) sel <- sel[sel$region %in% regs, , drop = FALSE]
    if ("region" %in% names(ibd)) ibd <- ibd[ibd$region %in% regs, , drop = FALSE]
  }
  layout <- .select_layout(x$chrom_layout(), chroms, skip_chr)
  sel <- .recum(sel[is.finite(sel[[metric]]), , drop = FALSE], layout)
  ibd <- .recum(ibd[is.finite(ibd$frac_pairs_ibd), , drop = FALSE], layout)
  genes <- .genes_for_layout(x$get_genes(), layout, highlight_genes)
  label_gap <- diff(range(c(layout$xmin, layout$xmax))) * 0.03

  present <- sort(unique(c(sel$region, ibd$region)))
  faceted <- ("region" %in% names(sel)) && length(present) > 1
  top_level <- if (faceted) present[1] else NULL

  # normalise each track so its tallest bar reaches the centre (y = 0). common =
  # one shared max (panels comparable); free = each region to its own max.
  normalized <- scale == "free" && "region" %in% names(sel)
  if (normalized) {
    sm <- tapply(sel[[metric]], sel$region, max, na.rm = TRUE)
    im <- tapply(ibd$frac_pairs_ibd, ibd$region, max, na.rm = TRUE)
    sm[!is.finite(sm) | sm <= 0] <- 1
    im[!is.finite(im) | im <= 0] <- 1
    sel$.tip <- 1 - sel[[metric]] / sm[as.character(sel$region)]
    ibd$.tip <- -1 + ibd$frac_pairs_ibd / im[as.character(ibd$region)]
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
    if (!is.null(thr) && !is.na(thr$region[1]) && "region" %in% names(sel)) {
      thr <- thr[thr$region %in% present, , drop = FALSE]
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
    (if (label_genes) .gene_label_layer(genes, if (faceted) "region" else NULL, top_level, label_gap)) +
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
    p <- p + ggplot2::facet_wrap(~ region, ncol = 1, strip.position = "right")
  } else if (!is.null(regs)) {
    p <- p + ggplot2::labs(title = paste0("tug-of-war: ", paste(regs, collapse = ", ")))
  }
  attr(p, "plasgenomics_dims") <- .dims_genome(
    if (faceted) length(present) else 1L,
    sum(layout$len) / sum(x$chrom_layout()$len), label_genes = label_genes)
  p
}

# ---- drug-gene region x region triangles -----------------------------------

# SNPs from a pairwise-region table overlapping one feature interval
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

#' IBD "triangle" panels for genes or specific SNPs
#'
#' Draws a region-by-region IBD-sharing triangle per feature: a **gene** (per-SNP
#' pairwise IBD for SNPs falling strictly inside the gene interval, aggregated) or a
#' **specific SNP** (the sharing at that single position). One facet per feature — use
#' it to ask whether a gene or locus is itself shared between regions.
#'
#' Genes come from the `genes` track on the [IbdResults] object and membership is
#' strict: a SNP belongs to a gene only when its position is inside the interval (no
#' flanking). Target individual loci with `snps` when a specific variant matters more
#' than a whole gene.
#'
#' @param x An [IbdResults] object (needs `pairwise_region`; a `genes` track for `genes`).
#' @param genes Gene names to include from the track (default: all in the track;
#'   case-insensitive). Ignored if only `snps` is given.
#' @param snps Specific SNPs to draw, as `"chr:pos"` ids (e.g. `"Pf3D7_07_v3:403222"`)
#'   or a data frame with `chr`, `pos` (and optional `name`). Each is one facet.
#' @param agg How a gene's value in each region-pair cell is computed from **all**
#'   its in-gene SNPs (no SNP is picked or dropped): `"mean"` (default), `"median"`,
#'   or `"max"` (the strongest sharing at that pair). A single-SNP feature is unaffected.
#' @param individual If `TRUE`, return a **named list** of one-triangle plots (one per
#'   feature, e.g. to write a multi-page PDF) with the legend tucked into the empty
#'   upper triangle; if `FALSE` (default), one faceted grid plot.
#' @param label Draw the value in each tile.
#' @param digits Decimal places for the tile labels.
#' @param ncol Facet columns for the grid (default: ggplot2 chooses).
#' @param trans Fill-scale transform, e.g. `"identity"` (default), `"log2"`, `"sqrt"`.
#' @param colors Optional colour ramp for the fill (default: the pairwise-sharing ramp).
#' @param limits Optional `c(lo, hi)` fill limits (values outside are squished).
#' @param fill_scale Optional ggplot2 fill scale that fully overrides the above.
#' @return A ggplot object (grid), or a named list of ggplot objects when `individual`.
#' @export
plot_drug_gene_triangles <- function(x, genes = NULL, snps = NULL,
                                     agg = c("mean", "median", "max"),
                                     individual = FALSE,
                                     label = TRUE, digits = 2, ncol = NULL,
                                     trans = "identity", colors = NULL, limits = NULL,
                                     fill_scale = NULL) {
  .need_package("ggplot2", "plot_drug_gene_triangles()")
  agg <- match.arg(agg)
  agg_fn <- switch(agg, mean = function(v) mean(v, na.rm = TRUE),
                   median = function(v) stats::median(v, na.rm = TRUE),
                   max = function(v) max(v, na.rm = TRUE))
  pw <- x$get_pairwise_region()
  if (is.null(pw)) stop("this IbdResults has no pairwise_region table", call. = FALSE)

  # assemble the features to draw: gene intervals (from the track) and/or specific
  # SNPs. Genes are used by default, or when explicitly named; passing only `snps`
  # skips the gene track so it doesn't warn about genes with no SNPs in the panel.
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
    ag <- stats::aggregate(frac_pairs_ibd ~ region_a + region_b, data = sub, FUN = agg_fn)
    ag$gene <- f$name
    parts[[length(parts) + 1]] <- ag
  }
  if (!length(parts)) stop("none of the requested features had overlapping SNPs", call. = FALSE)
  df <- do.call(rbind, parts)
  df$gene <- factor(df$gene, levels = feats$name[feats$name %in% df$gene])

  regions <- sort(unique(c(as.character(df$region_a), as.character(df$region_b))))
  fs <- .ibd_fill_scale("pairs IBD", trans = trans, colors = colors, limits = limits,
                        fill_scale = fill_scale)
  feats_present <- levels(df$gene)

  if (individual) {
    plots <- lapply(feats_present, function(g) {
      pg <- .triangle_gg(df[df$gene == g, , drop = FALSE], regions, fs, label, digits,
                         title = g, legend_inside = TRUE)
      attr(pg, "plasgenomics_dims") <- .dims_triangles(1L, length(regions), 1L)
      pg
    })
    names(plots) <- feats_present
    return(plots)
  }

  p <- .triangle_gg(df, regions, fs, label, digits) +
    ggplot2::facet_wrap(~ .data$gene, ncol = ncol) +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"))
  attr(p, "plasgenomics_dims") <- .dims_triangles(length(feats_present), length(regions), ncol)
  p
}

# one region x region triangle for a single feature's aggregated data
.triangle_gg <- function(df, regions, fill_scale, label, digits, title = NULL,
                         legend_inside = FALSE) {
  df$region_a <- factor(df$region_a, levels = regions)
  df$region_b <- factor(df$region_b, levels = rev(regions))
  mid <- stats::median(df$frac_pairs_ibd, na.rm = TRUE)
  p <- ggplot2::ggplot(df, ggplot2::aes(.data$region_a, .data$region_b,
                                        fill = .data$frac_pairs_ibd)) +
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
                   colour = .data$frac_pairs_ibd > mid),
      size = 2.8, show.legend = FALSE) +
      ggplot2::scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15"))
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

.resolve_regions <- function(region, sel, ibd) {
  if (!is.null(region)) return(region)
  regs <- unique(c(if ("region" %in% names(sel)) sel$region,
                   if ("region" %in% names(ibd)) ibd$region))
  if (length(regs)) regs else NULL
}
