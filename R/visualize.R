# Genome-wide IBD / selection plots over an IbdResults object. Each returns a
# ggplot object (composable / saveable by the caller); ggplot2 and scales are
# optional (Suggests), guarded at call time.

# ---- shared plotting pieces ------------------------------------------------

# Alternating chromosome bands behind a genome-wide track. By default the second band is
# the panel itself (`NA` = draw nothing), which is what a track of points wants.
#
# Two greys instead, for a track whose *fill* carries the meaning: the IBD scales start at
# white, so a zero tile drawn over a white band vanishes while the same tile over a grey
# band is plainly there -- the plot would then show sequenced-but-not-IBD positions on
# every other chromosome only. Both bands grey keeps a white tile visible throughout.
.CHR_BAND_TILE <- c("#d9d9d9", "#f0f0f0")

# One layer per band rather than one layer with a vector of fills: a faceted plot
# replicates each layer's data across panels, so a per-row fill would not line up.
.chr_band_layer <- function(layout, fills = c("#ebebeb", NA)) {
  out <- lapply(seq_along(c("a", "b")), function(i) {
    d <- layout[layout$band == c("a", "b")[i], , drop = FALSE]
    if (is.na(fills[i]) || !nrow(d)) return(NULL)
    ggplot2::geom_rect(
      data = d,
      ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = -Inf, ymax = Inf),
      fill = fills[i], inherit.aes = FALSE
    )
  })
  out <- Filter(Negate(is.null), out)
  if (!length(out)) NULL else out
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

# The facet value the gene labels are pinned to, as a FACTOR carrying the full level set.
# Passing a bare character here is the recurring bug in this package: a secondary layer whose
# facet column is character makes ggplot merge the layers' facet values and re-sort them
# alphabetically, which both scrambles the panel order and sends the labels to the wrong
# panel. `levels(v)[1]` returns character, so never use it directly for this.
.first_level <- function(v) {
  v <- .as_group_factor(v)
  factor(levels(v)[1], levels = levels(v))
}

# room above the panel for outside gene labels (only when labelling)
.gene_label_space <- function(label_genes) {
  if (!label_genes) return(NULL)
  list(ggplot2::coord_cartesian(clip = "off"),
       ggplot2::theme(plot.margin = ggplot2::margin(t = 42, r = 12, b = 6, l = 6)))
}

# `draw_threshold` takes a kind by name, or TRUE/FALSE for the Bonferroni line / no line.
# Resolving it in one place keeps the Manhattan and the tug-of-war from drifting apart.
.resolve_threshold <- function(x) {
  if (isTRUE(x)) return("bonferroni")
  if (isFALSE(x) || is.null(x)) return("none")
  if (!is.character(x))
    stop("`draw_threshold` must be TRUE/FALSE or one of \"bonferroni\", \"fdr\", ",
         "\"permutation\", \"empirical\", \"all\", \"both\", \"none\"", call. = FALSE)
  match.arg(x, c("bonferroni", "fdr", "permutation", "empirical", "all", "both", "none"))
}

.threshold_kinds <- function(which_thr) {
  switch(which_thr,
         none = character(0),
         both = c("bonferroni", "fdr"),
         all = c("bonferroni", "fdr", "permutation", "empirical"),
         which_thr)
}

# How a top-track metric is written on the shared axis. A column the plot knows about gets
# the name and the units the IBD half already uses, so the two halves read as one scheme
# rather than as a statistic mirrored against a percentage; anything else keeps its own
# column name, which is what a caller who passed an arbitrary table expects to see.
.TOP_METRIC_STYLE <- list(
  frac_extreme = list(label = "Extreme Fraction", percent = TRUE))

# colour / linetype per threshold kind, so one line always means one thing
# The two chi2(1)-based lines are warm, the two permutation-based ones cool-to-green, so
# which family a line belongs to reads off the plot before the legend does.
.THRESHOLD_STYLE <- list(
  bonferroni  = list(colour = "firebrick", linetype = "dashed"),
  fdr         = list(colour = "#E69F00",  linetype = "dotdash"),
  permutation = list(colour = "#117733",  linetype = "solid"),
  empirical   = list(colour = "#2271B2",  linetype = "longdash")
)

# A quantile of the top track itself is not one of the kinds above -- it is not a test, it
# is the height the plotted statistic reaches on its own -- so it gets its own colour
# rather than borrowing one that already means "threshold from the object's run".
.TOP_QUANTILE_STYLE <- list(colour = "#762A83", linetype = "dotted")

# A washed-out version of a bar colour, for the stretch of a track that is below its bar.
# Mixing towards white rather than dropping chroma alone: the panel carries grey chromosome
# bands behind the bars, and a colour that only loses saturation goes muddy against them
# instead of receding. `amount` is how far towards white, so 0 is the colour itself.
.wash <- function(col, amount = 0.55) {
  # 0 hands back exactly what it was given, rather than a re-spelled equal colour: at that
  # setting `shade` is meant to be a no-op, and it should be one in the data too
  if (amount <= 0) return(col)
  m <- grDevices::col2rgb(col) / 255
  grDevices::rgb(t(m + (1 - m) * amount))
}

# Runs of a track that sit above its bar, as merged genomic intervals, one set per group.
#
# Merged, because the unit a reader judges is a peak and not a window: a 100 kb window
# stepping 10 kb puts about twenty overlapping rows over one peak, and the two halves of the
# mirror rarely summit on exactly the same row. Classifying row by row therefore paints a
# fringe of "only this track" around the shoulders of every peak both tracks agree on, which
# is the one thing this figure must not invent. Rows within `gap` of each other join.
.peak_intervals <- function(df, value, bar, gap) {
  hi <- is.finite(df[[value]]) & df[[value]] > bar[as.character(df$group)]
  d <- df[hi, , drop = FALSE]
  if (!nrow(d)) return(NULL)
  # a windowed track carries its own extent; a per-SNP track is a point until `gap` joins it
  lo <- if (all(c("start", "end") %in% names(d))) as.numeric(d$start) else as.numeric(d$pos)
  hi2 <- if (all(c("start", "end") %in% names(d))) as.numeric(d$end) else as.numeric(d$pos)
  key <- paste(as.character(d$group), as.character(d$chr), sep = "\r")
  out <- lapply(split(data.frame(lo = lo, hi = hi2), key), function(b) {
    b <- b[order(b$lo), , drop = FALSE]
    # a new interval starts where this row begins beyond everything seen so far, plus gap
    brk <- c(TRUE, b$lo[-1] > cummax(b$hi)[-nrow(b)] + gap)
    g <- cumsum(brk)
    data.frame(start = tapply(b$lo, g, min), end = tapply(b$hi, g, max))
  })
  res <- do.call(rbind, lapply(names(out), function(k) {
    p <- strsplit(k, "\r", fixed = TRUE)[[1]]
    data.frame(group = p[1], chr = p[2], start = out[[k]]$start, end = out[[k]]$end,
               stringsAsFactors = FALSE)
  }))
  rownames(res) <- NULL
  res
}

# The peaks either half calls, and which halves called each one.
.tug_peaks <- function(top_iv, bot_iv, gap) {
  all_iv <- rbind(
    if (!is.null(top_iv)) cbind(top_iv, .side = "top") else NULL,
    if (!is.null(bot_iv)) cbind(bot_iv, .side = "bottom") else NULL)
  if (is.null(all_iv) || !nrow(all_iv)) return(NULL)
  key <- paste(all_iv$group, all_iv$chr, sep = "\r")
  res <- lapply(split(all_iv, key), function(b) {
    b <- b[order(b$start), , drop = FALSE]
    g <- cumsum(c(TRUE, b$start[-1] > cummax(b$end)[-nrow(b)] + gap))
    do.call(rbind, lapply(split(b, g), function(k) data.frame(
      group = k$group[1], chr = k$chr[1], start = min(k$start), end = max(k$end),
      top = any(k$.side == "top"), bottom = any(k$.side == "bottom"),
      stringsAsFactors = FALSE)))
  })
  out <- do.call(rbind, res)
  rownames(out) <- NULL
  out$cell <- ifelse(out$top & out$bottom, "both", ifelse(out$top, "top", "bottom"))
  out
}

# One half's bar, as a value per group name, for the shading and the ribbon.
#
# Taken on the table before any crop, for the same reason the drawn quantile line is: a
# figure that drops the quiet chromosomes must not thereby raise the bar the surviving peaks
# are judged against. A track with no `group` column gets one bar repeated across the panels.
.tug_bar <- function(df, value, p, groups) {
  v <- abs(df[[value]])
  ok <- is.finite(v)
  if (!any(ok)) return(stats::setNames(rep(NA_real_, length(groups)), groups))
  if (!"group" %in% names(df))
    return(stats::setNames(rep(unname(stats::quantile(v[ok], p)), length(groups)), groups))
  q <- tapply(v[ok], as.character(df$group)[ok], stats::quantile, probs = p)
  stats::setNames(as.numeric(q)[match(groups, names(q))], groups)
}

# The `p`-quantile of a top track's plotted magnitude, per group when it has one.
#
# Computed on the table as passed in, BEFORE any `chroms` / `skip_chr` / `zoom` crop, so
# the reference is genome-wide and cropping to the chromosomes that carry signal cannot
# raise the bar the signal is then judged against. Per group for the same reason a facet
# gets its own panel: a region's windows are the distribution its own peaks stand out of.
.top_quantile_lines <- function(df, metric, p, use_abs) {
  v <- df[[metric]]
  if (isTRUE(use_abs)) v <- abs(v)
  ok <- is.finite(v)
  if (!any(ok)) return(NULL)
  v <- v[ok]
  if (!"group" %in% names(df)) {
    return(data.frame(group = NA_character_, threshold = unname(stats::quantile(v, p)),
                      stringsAsFactors = FALSE))
  }
  g <- as.character(df$group)[ok]
  q <- tapply(v, g, stats::quantile, probs = p)
  data.frame(group = names(q), threshold = as.numeric(unlist(q)),
             stringsAsFactors = FALSE)
}

# One kind of threshold as a group/threshold frame, or NULL when it is not available.
# `strict` says the caller named this kind, so a run that lacks it is an error rather
# than something to skip; under "all" the kinds are implied and missing ones are skipped.
.threshold_line <- function(thr, kind, strict = TRUE) {
  if (is.null(thr)) return(NULL)
  col <- switch(kind, fdr = "neg_log10_p_fdr_threshold",
                permutation = "neg_log10_p_perm_threshold",
                empirical = "neg_log10_p_emp_fdr_threshold", "threshold")
  if (!col %in% names(thr)) {
    if (strict && kind == "fdr")
      stop("this selection run has no FDR threshold; regenerate with a current ",
           "`ibd_selection_statistic`", call. = FALSE)
    if (strict && kind %in% c("permutation", "empirical"))
      stop("this selection run has no ", kind, " threshold; rerun ",
           "`ibd_selection_statistic --permute 200`", call. = FALSE)
    return(NULL)
  }
  out <- thr[c(intersect("group", names(thr)), col)]
  names(out)[ncol(out)] <- "threshold"
  out <- out[is.finite(out$threshold), , drop = FALSE]
  if (!nrow(out)) NULL else out
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
    colours = colors, trans = trans, limits = limits,
    name = .fill_scale_name(name, trans),
    na.value = colors[1],
    # a function, so this works off the *trained* range when no limits were given --
    # otherwise a per-page plot keeps the transform's default breaks and their padding
    breaks = function(lims) {
      b <- .fill_breaks(lims, trans)
      if (is.null(b)) ggplot2::waiver() else b
    },
    labels = .fill_labels,
    oob = if (is.null(limits)) scales::censor else scales::squish
  )
}

# The transform's own round breaks, plus the maximum.
#
# A transform picks pleasant round numbers -- powers of two under log2 -- but they stop
# wherever they stop, which under a log scale is routinely a long way below the top: the
# bar then runs on past its last label and the strongest colour has no number attached.
# So keep the round breaks and add the limit itself, dropping any round break that would
# print on top of it.
.fill_breaks <- function(limits, trans) {
  if (is.null(limits) || length(limits) != 2 || !all(is.finite(limits))) return(NULL)
  lo <- limits[1]
  hi <- limits[2]
  if (hi <= lo) return(NULL)
  if (.log_trans(trans) && lo <= 0) return(NULL)

  tf <- .transform_fun(trans)
  nice <- tryCatch(.transform_breaks(trans, limits), error = function(e) NULL)
  nice <- nice[is.finite(nice) & nice >= lo & nice <= hi]
  if (!length(nice)) return(sort(unique(c(lo, hi))))

  # "too close to the top to label separately": within a tenth of the bar's height
  span <- tf(hi) - tf(lo)
  keep <- if (is.finite(span) && span > 0) abs(tf(hi) - tf(nice)) > 0.1 * span else TRUE
  sort(unique(c(nice[keep], hi)))
}

# scales renamed as.trans() to as.transform(); accept whichever this install has.
.as_transform <- function(trans) {
  f <- asNamespace("scales")$as.transform
  if (is.null(f)) f <- asNamespace("scales")$as.trans
  f(trans)
}

# The break function belonging to a transform, so log2 gives powers of two.
.transform_breaks <- function(trans, limits) .as_transform(trans)$breaks(limits)

.transform_fun <- function(trans) {
  if (identical(trans, "identity")) return(identity)
  .as_transform(trans)$transform
}

# Name the transform in the legend title: a "pairs IBD" bar running 0.0039 -> 0.52 in
# even visual steps is a log scale, and saying so is cheaper than expecting the reader to
# infer it from the spacing.
.fill_scale_name <- function(name, trans) {
  label <- if (is.character(trans)) trans
           else tryCatch(.as_transform(trans)$name, error = function(e) NULL)
  if (is.null(label) || !nzchar(label) || identical(label, "identity")) return(name)
  paste0(name, " (", label, ")")
}

# Compact labels: significant digits, not decimal places, and no trailing zeros.
#
# A log scale's default formatter pads every break to one fixed width, which is where
# labels like "0.0039062500" come from. Rounding to fixed decimals instead is worse: on a
# log scale the smallest break is often a few ten-thousandths, and two decimal places
# render it "0", which reads as nothing at all. Significant digits keep every break
# informative regardless of its magnitude, and the loop only adds precision when two
# breaks would otherwise print the same.
.fill_labels <- function(breaks) {
  fmt <- function(v, sig) {
    vapply(v, function(x) {
      if (!is.finite(x)) return(NA_character_)
      if (x == 0) return("0")
      format(signif(x, sig), scientific = FALSE, trim = TRUE, drop0trailing = TRUE)
    }, character(1))
  }
  for (sig in 2:6) {
    out <- fmt(breaks, sig)
    if (!anyDuplicated(out)) return(out)
  }
  fmt(breaks, 6)
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
#' @param colours,colors Optional length-2 colour vector for the alternating chromosome bands.
#' @param zoom Optional single interval to crop to, keeping the same data and the same
#'   coordinates as the genome-wide plot: a chromosome (`"7"`), a range
#'   (`"7:728,081-988,719"`), a gene name from the object's track, or a one-row data frame
#'   with chr/start/end. Every gene in the window is drawn and named unless
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
#' @return A ggplot object.
#' @examples
#' plot_ibd_sharing_manhattan(example_ibd_results())
#' @export
plot_ibd_sharing_manhattan <- function(x, groups = NULL, chroms = NULL, skip_chr = NULL,
                                       zoom = NULL, zoom_pad = 0.05,
                                       genes_for_track = NULL, gene_label_angle = 0,
                               highlight_genes = NULL, label_genes = NULL,
                               point_size = 0.5, point_alpha = 0.6, colours = NULL,
                                       colors = NULL) {
  colours <- .alias_arg("colours", "colors")
  .need_package("ggplot2", "plot_ibd_sharing_manhattan()")
  .need_package("scales", "plot_ibd_sharing_manhattan()")
  df <- x$get_per_snp_group()
  if (is.null(df)) stop("this IbdResults has no per_snp_group table", call. = FALSE)
  .check_gene_request(x$get_genes(), highlight_genes)
  label_arg <- label_genes
  layout <- .select_layout(x$chrom_layout(), chroms, skip_chr)
  z <- .zoom_setup(zoom, x$get_genes(), layout, label_arg, zoom_pad,
                   x$reference_id(), genes_for_track = genes_for_track)
  # zoomed, the gene names go in the track stacked underneath, not above the panel
  label_genes <- if (!is.null(z)) FALSE
    else if (is.null(label_arg)) !is.null(highlight_genes) else isTRUE(label_arg)
  df <- .attach_band(.recum(.filter_group(df, groups), layout), layout)
  df <- df[is.finite(df$frac_pairs_ibd), , drop = FALSE]
  df <- .crop_to_window(df, z, "IBD SNPs")
  genes <- if (is.null(z)) .genes_for_layout(x$get_genes(), layout, highlight_genes)
           else z$genes
  faceted <- "group" %in% names(df)
  top_level <- if (faceted) .first_level(df$group) else NULL
  genome_frac <- sum(layout$len) / sum(x$chrom_layout()$len)

  p <- ggplot2::ggplot(df, ggplot2::aes(.data$cum_pos, .data$frac_pairs_ibd)) +
    .chr_band_layer(layout) +
    .gene_line_layer(genes) +
    (if (label_genes) .gene_label_layer(genes, if (faceted) "group" else NULL, top_level,
       footprint = .label_footprint_bp(nchar(genes$name), .fp_layout(z, layout),
                                       .fp_frac(z, genome_frac)))) +
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
  n_panels <- if (faceted) length(unique(df$group)) else 1L
  p <- .apply_zoom(p, z, layout, n_panels, gene_label_angle)
  attr(p, "plasgenomics_dims") <- if (is.null(z))
    .dims_genome(n_panels, genome_frac, label_genes = label_genes)
    else .dims_zoom(n_panels, attr(p, "plasgenomics_track_in"))
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
#' @param draw_threshold Which significance line(s) to draw (only for `neg_log10_p`).
#'   The two read off a chi-square(1) are warm, the two built by permutation cool:
#'   `TRUE` / `"bonferroni"` (red dashed), `"fdr"` (orange dot-dash), `"permutation"`
#'   (family-wise, green solid), `"empirical"` (FDR over the permutation's own p-values,
#'   blue long-dash). Also `"both"` for the first two, `"all"` for every kind present, or
#'   `FALSE` for none.
#'   Prefer a permutation line where it disagrees with a parametric one: it is built by
#'   re-shuffling the IBD segments themselves, so it needs no chi-square(1) reference and
#'   it accounts for one segment spanning many SNPs, and it can land well above
#'   Bonferroni's. `lambda_gc` says how far that reference is from fitting;
#'   the further from 1, the less the parametric lines mean.
#'   Naming a kind the run did not write is an error; `"all"` draws whichever kinds the
#'   threshold table carries.
#' @param point_size,point_alpha Point aesthetics.
#' @param colours,colors Optional length-2 colour vector for the alternating chromosome bands.
#' @param zoom Optional single interval to crop to, keeping the same data and the same
#'   coordinates as the genome-wide plot: a chromosome (`"7"`), a range
#'   (`"7:728,081-988,719"`), a gene name from the object's track, or a one-row data frame
#'   with chr/start/end. Every gene in the window is drawn and named unless
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
#' @return A ggplot object.
#' @examples
#' plot_selection_manhattan(example_ibd_results(), genes = PF_EXAMPLE_DRUG_GENES)
#' @export
plot_selection_manhattan <- function(x, metric = c("neg_log10_p", "chi2_stat", "z_score"),
                                     groups = NULL, chroms = NULL, skip_chr = NULL,
                                     zoom = NULL, zoom_pad = 0.05,
                                     genes_for_track = NULL, gene_label_angle = 0,
                                     highlight_genes = NULL, label_genes = NULL,
                                     draw_threshold = TRUE,
                                     point_size = 0.5, point_alpha = 0.6, colours = NULL,
                                     colors = NULL) {
  colours <- .alias_arg("colours", "colors")
  .need_package("ggplot2", "plot_selection_manhattan()")
  metric <- match.arg(metric)
  df <- x$get_selection()
  if (is.null(df)) stop("this IbdResults has no selection table", call. = FALSE)
  if (!metric %in% names(df)) stop(sprintf("selection table has no '%s' column", metric),
                                    call. = FALSE)
  .check_gene_request(x$get_genes(), highlight_genes)
  label_arg <- label_genes
  layout <- .select_layout(x$chrom_layout(), chroms, skip_chr)
  z <- .zoom_setup(zoom, x$get_genes(), layout, label_arg, zoom_pad,
                   x$reference_id(), genes_for_track = genes_for_track)
  # zoomed, the gene names go in the track stacked underneath, not above the panel
  label_genes <- if (!is.null(z)) FALSE
    else if (is.null(label_arg)) !is.null(highlight_genes) else isTRUE(label_arg)
  df <- .attach_band(.recum(.filter_group(df, groups), layout), layout)
  df$.y <- df[[metric]]
  df <- df[is.finite(df$.y), , drop = FALSE]
  df <- .crop_to_window(df, z, "SNPs")
  genes <- if (is.null(z)) .genes_for_layout(x$get_genes(), layout, highlight_genes)
           else z$genes
  faceted <- "group" %in% names(df)
  top_level <- if (faceted) .first_level(df$group) else NULL
  genome_frac <- sum(layout$len) / sum(x$chrom_layout()$len)

  p <- ggplot2::ggplot(df, ggplot2::aes(.data$cum_pos, .data$.y)) +
    .chr_band_layer(layout) +
    .gene_line_layer(genes) +
    (if (label_genes) .gene_label_layer(genes, if (faceted) "group" else NULL, top_level,
       footprint = .label_footprint_bp(nchar(genes$name), .fp_layout(z, layout),
                                       .fp_frac(z, genome_frac)))) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$.band), size = point_size,
                        alpha = point_alpha, stroke = 0) +
    ggplot2::scale_colour_manual(values = .band_colours(colours)) +
    .chr_axis(layout) +
    ggplot2::labs(x = "chromosome", y = metric) +
    .manhattan_theme() +
    .gene_label_space(label_genes)

  thr <- x$get_thresholds()
  which_thr <- .resolve_threshold(draw_threshold)
  kinds <- .threshold_kinds(which_thr)
  if (metric == "neg_log10_p") for (kind in kinds) {
    # `.threshold_line()` drops non-finite thresholds: a group with no valid SNPs gets a
    # NaN threshold (written as an empty cell by `ibd_selection_statistic`), and feeding
    # NA to geom_hline draws nothing but warns ("Removed N rows ...").
    line <- .threshold_line(thr, kind, strict = which_thr == kind)
    if (is.null(line)) next
    sty <- .THRESHOLD_STYLE[[kind]]
    has_g <- "group" %in% names(line) && any(!is.na(line$group)) && "group" %in% names(df)
    if (has_g) line <- line[line$group %in% df$group, , drop = FALSE]
    if (!nrow(line)) next
    p <- p + if (has_g)
      ggplot2::geom_hline(data = line,
        ggplot2::aes(yintercept = .data$threshold), colour = sty$colour,
        linetype = sty$linetype, linewidth = 0.4)
    else
      ggplot2::geom_hline(yintercept = line$threshold[1], colour = sty$colour,
        linetype = sty$linetype, linewidth = 0.4)
  }
  if ("group" %in% names(df)) {
    p <- p + ggplot2::facet_wrap(~ .data$group, ncol = 1, strip.position = "right")
  }
  n_panels <- if (faceted) length(unique(df$group)) else 1L
  p <- .apply_zoom(p, z, layout, n_panels, gene_label_angle)
  attr(p, "plasgenomics_dims") <- if (is.null(z))
    .dims_genome(n_panels, genome_frac, label_genes = label_genes)
    else .dims_zoom(n_panels, attr(p, "plasgenomics_track_in"))
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
#' @param colors,colours Optional colour ramp for the fill (defaults to a single-hue
#'   light-to-dark sequential scale that stays readable when most values are near 0).
#' @param limits Optional `c(lo, hi)` fill limits; values outside are squished into
#'   range, so a few extremes near 1 don't crush the rest of the scale.
#' @param fill_scale Optional ggplot2 fill scale that fully overrides the above.
#' @param highlight_genes Optional gene names from the `genes` track to mark with lines;
#'   requesting a name not in the track is an error.
#' @param label_genes Label the genes above the top panel. `NULL` (default) labels them
#'   only when `highlight_genes` is given; `TRUE`/`FALSE` forces it.
#' @return A ggplot object.
#' @examples
#' plot_ibd_pairwise_group_heatmap(example_ibd_results())
#' @export
plot_ibd_pairwise_group_heatmap <- function(x, anchor = NULL, chroms = NULL, skip_chr = NULL,
                                    trans = "identity", colors = NULL, limits = NULL,
                                    fill_scale = NULL, highlight_genes = NULL,
                                    label_genes = NULL,
                                            colours = NULL) {
  colors <- .alias_arg("colors", "colours")
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
    .chr_band_layer(layout, .CHR_BAND_TILE) +          # two greys: see .chr_band_layer()
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

# The top half of the mirror is either the object's own IBD selection statistic or an
# external per-SNP scan (iHS, beta, ...), which is the same shape: chr, pos, a value column
# and optionally the group the IBD half is split by. Keeping the two interchangeable is
# what lets one plot ask whether the IBD signal coincides with a haplotype signal, not only
# with the IBD statistic derived from the same segments.
.top_track <- function(x, top, metric, label = NULL) {
  if (is.null(top)) {
    df <- x$get_selection()
    if (is.null(df)) stop("this IbdResults has no selection table", call. = FALSE)
    return(list(df = df, label = label %||% "selection", own = TRUE))
  }
  df <- as.data.frame(top)
  if (!all(c("chr", "pos") %in% names(df)))
    stop("`top` needs chr and pos columns", call. = FALSE)
  if (!"group" %in% names(df) && "pair" %in% names(df))
    stop("a two-population scan (Rsb / XP-EHH) is indexed by population pair, not by the ",
         "groups the IBD half is split by; use a run_ihs() scan for the top track",
         call. = FALSE)
  df$chr <- normalise_chr(df$chr)
  list(df = df, label = label %||% "scan", own = FALSE)
}


#' IBD / selection "tug-of-war" mirror plot
#'
#' The selection statistic hangs from the top (bars descending) while the per-SNP
#' IBD fraction rises from the bottom (bars ascending), sharing one centred axis so
#' peaks of each line up. A single left axis carries both halves, its tick labels
#' tinted to their track (selection on top, IBD on the bottom).
#'
#' @param x An [IbdResults] object (needs a per_snp_group table, and a selection table
#'   unless `top` supplies the upper track).
#' @param top Optional per-SNP scan to hang from the top instead of the object's own IBD
#'   selection statistic: a [run_ihs()] result, a [beta_score()] table, or any table with
#'   chr, pos, the `metric` column and (to face the IBD groups) a matching `group` column.
#'   Two-population scans ([run_rsb()], [run_xpehh()]) are indexed by population pair
#'   rather than group and cannot be mirrored against per-group IBD.
#' @param top_label Name for the upper track on the shared axis; defaults to
#'   `"selection"` for the object's statistic and `"scan"` for a `top` table.
#' @param group Group(s) to plot. `NULL` (default) plots every group, faceted one
#'   per row; a single group gives one panel; a vector facets those groups.
#' @param metric Selection metric for the top track (default `"neg_log10_p"`).
#' @param scan_abs Plot the magnitude of a signed `metric` (iHS, Rsb, a z-score). `NULL`
#'   (default) does so whenever the metric has negative values, labelling the axis
#'   `|metric|`. The mirror cannot show the sign -- its top half only spans zero to the
#'   maximum -- so `FALSE` on signed values is an error rather than a silently clipped plot;
#'   use [plot_ibd_locus()] when the sign matters.
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
#' @param draw_threshold Which significance line(s) to draw (only for `neg_log10_p`, and
#'   only when not `normalized`): `TRUE` / `"bonferroni"`, `"fdr"`, `"permutation"`,
#'   `"empirical"`, `"both"`, `"all"` or `FALSE` -- the same set
#'   [plot_selection_manhattan()] takes, in the same colours, resolved by the same helper.
#'   Each line is mapped through the same transform as the mirrored selection half, so it
#'   lands where the data does.
#' @param top_quantile Draw a reference line at this quantile of the top track's own
#'   distribution -- `0.99` for the 99th percentile, `NULL` (default) for no line. This is
#'   an empirical reference, not a test: it says how high the statistic reaches across the
#'   genome, so a peak can be read against the rest of the scan rather than against a
#'   nominal per-SNP tail. It is what a windowed fraction wants, because the fraction of
#'   SNPs over `threshold` in [ihs_windows()] has no significance line of its own -- the
#'   cutoff inside it defines the tail, and the evidence is the genome-wide distribution of
#'   the fraction. The quantile is taken on the table as passed in, before `chroms`,
#'   `skip_chr` and `zoom` crop it, so dropping the quiet chromosomes cannot raise the bar
#'   its own peaks are then judged against; and per group when the track has a group
#'   column, so each panel's line is that region's genome-wide quantile whether or not the
#'   other regions are drawn. Stacks with `draw_threshold`, in its own colour.
#' @param shade Draw the stretch of each track that sits below its bar in a washed-out
#'   version of its own colour, so the peaks that clear the bar carry the full one
#'   (default `FALSE`). It is one hue at two strengths per track, not a second hue, so it
#'   survives every dichromacy and reproduces in grey. The bar is `peak_quantile`.
#' @param shade_wash How far towards white `shade` takes the below-bar stretch, from 0 (no
#'   change) to 1 (invisible); default 0.55. Worth turning down for print, where a wash that
#'   reads on a screen can drop out altogether.
#' @param ribbon Draw a strip along the centre line saying, for each peak, whether both
#'   halves called it or only one (default `FALSE`). This is the comparison the mirror is
#'   for, and reading it off two bars either side of a gap is exactly what the eye is bad
#'   at. Both halves are cut at `peak_quantile`.
#'
#'   The unit is a peak, not a window, and that is the whole reason this is worth a function
#'   rather than a colour mapped per row. A sliding window puts many overlapping rows over
#'   one peak and the two halves seldom summit on the same row, so a row-by-row rule paints
#'   a fringe of "only this half" around the shoulders of every peak they agree on. Runs
#'   above the bar are merged into peaks first (joining across gaps up to `peak_gap`), the
#'   two halves' peaks are then unioned, and each resulting peak is labelled by which halves
#'   reach into it.
#'
#'   Read it knowing what the halves can and cannot see. A peak only the IBD half calls is
#'   as likely to be a sweep near fixation, where a statistic contrasting two alleles has
#'   nothing left to contrast, as it is to be nothing. A peak only the top half calls is
#'   consistent with several origins of the same allele, and equally with low recombination
#'   or unmodelled structure -- it is a screen, not a test, and counting haplotype
#'   backgrounds among the carriers is what settles it.
#' @param peak_quantile Quantile of each half's own genome-wide distribution that `shade`
#'   and `ribbon` treat as the bar (default `0.99`). Taken per group and before `chroms`,
#'   `skip_chr` and `zoom` crop anything, as [plot_ibd_tugofwar()]'s `top_quantile` is.
#'   Matched quantiles are the defensible choice here because neither half then borrows the
#'   other's stringency; giving one half a nominal cutoff and the other an empirical one can
#'   reverse which half looks the more sensitive, and the ribbon would report that as
#'   biology. It cuts both halves at the same place by construction, so the count of
#'   top-only and bottom-only peaks is close to matched -- it is *which* loci fall in each
#'   that carries the information, not how many.
#' @param peak_gap Distance in base pairs within which two runs above the bar are treated as
#'   one peak (default 50 kb). Also the distance at which the two halves' peaks are taken to
#'   be the same peak.
#' @param ribbon_colours,ribbon_colors Fills for the ribbon, named `both`, `top` and
#'   `bottom`. The default gives each single-half cell its own track's colour, so two of the
#'   three need no legend lookup, and `both` the darkest of the three.
#' @param metric_label Name for the top metric on the shared axis. `NULL` (default) uses
#'   the column name, except for a column the plot knows how to write -- `frac_extreme`
#'   from [ihs_windows()] reads `Extreme Fraction`, matching `IBD Fraction` below it.
#' @param top_percent Write the top track's tick labels as percentages. `NULL` (default)
#'   does so when the metric is a fraction (and always under `scale = "free"`, where both
#'   halves are percentages of their own maximum), so a fraction on top is not shown in
#'   different units from the fraction underneath it.
#' @param centre_gap Fraction of each half-axis left empty at the centre line, so the
#'   tallest bar of each track -- and its largest tick label -- stops short of the middle
#'   instead of meeting the other track's there. `0` restores the two halves meeting.
#' @param selection_colour,ibd_colour Bar and axis colours for the two tracks.
#' @param zoom Optional single interval to crop to, keeping the same data and the same
#'   coordinates as the genome-wide plot: a chromosome (`"7"`), a range
#'   (`"7:728,081-988,719"`), a gene name from the object's track, or a one-row data frame
#'   with chr/start/end. Every gene in the window is drawn and named unless
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
#' @return A ggplot object.
#' @examples
#' # IBD sharing above, the selection statistic mirrored below
#' plot_ibd_tugofwar(example_ibd_results())
#' @export
plot_ibd_tugofwar <- function(x, group = NULL, top = NULL, top_label = NULL,
                              metric = "neg_log10_p", scan_abs = NULL,
                              scale = c("common", "free"), chroms = NULL, skip_chr = NULL,
                              zoom = NULL, zoom_pad = 0.05,
                              genes_for_track = NULL, gene_label_angle = 0,
                              highlight_genes = NULL, label_genes = NULL,
                              draw_threshold = TRUE, top_quantile = NULL,
                              shade = FALSE, shade_wash = 0.55,
                              ribbon = FALSE, peak_quantile = 0.99,
                              peak_gap = 50000, ribbon_colours = NULL,
                              metric_label = NULL,
                              top_percent = NULL, centre_gap = 0.08,
                              selection_colour = "#fd8d3c", ibd_colour = "#2166ac",
                              ribbon_colors = NULL) {
  .need_package("ggplot2", "plot_ibd_tugofwar()")
  .need_package("scales", "plot_ibd_tugofwar()")
  scale <- match.arg(scale)
  sty <- .TOP_METRIC_STYLE[[metric]]
  metric_label <- metric_label %||% sty$label %||% metric
  top_percent <- top_percent %||% isTRUE(sty$percent)
  if (!is.numeric(centre_gap) || length(centre_gap) != 1L || centre_gap < 0 ||
      centre_gap >= 1)
    stop("`centre_gap` must be a single fraction of the half-axis, in [0, 1)", call. = FALSE)
  ribbon_colours <- .alias_arg("ribbon_colours", "ribbon_colors")
  if (!is.null(top_quantile) &&
      (!is.numeric(top_quantile) || length(top_quantile) != 1L ||
       !is.finite(top_quantile) || top_quantile <= 0 || top_quantile >= 1))
    stop("`top_quantile` must be a single probability in (0, 1), e.g. 0.99", call. = FALSE)
  want_peaks <- isTRUE(shade) || isTRUE(ribbon)
  if (want_peaks) {
    if (!is.numeric(peak_quantile) || length(peak_quantile) != 1L ||
        !is.finite(peak_quantile) || peak_quantile <= 0 || peak_quantile >= 1)
      stop("`peak_quantile` must be a single probability in (0, 1), e.g. 0.99",
           call. = FALSE)
    if (!is.numeric(peak_gap) || length(peak_gap) != 1L || !is.finite(peak_gap) ||
        peak_gap < 0)
      stop("`peak_gap` must be one distance in base pairs", call. = FALSE)
    if (!is.numeric(shade_wash) || length(shade_wash) != 1L || !is.finite(shade_wash) ||
        shade_wash < 0 || shade_wash > 1)
      stop("`shade_wash` must be a single fraction in [0, 1]", call. = FALSE)
    # the drawn line and the shading would otherwise be two different bars in one figure,
    # and the reader has no way to tell which the colours answered to
    if (!is.null(top_quantile) && !isTRUE(all.equal(top_quantile, peak_quantile)))
      message("`top_quantile` (", top_quantile, ") draws the line and `peak_quantile` (",
              peak_quantile, ") cuts the shading and ribbon, so the line will not sit ",
              "where the colour changes; set them equal unless that is deliberate")
  }
  tt <- .top_track(x, top, metric, top_label)
  # kept before any crop: the quantile line is a genome-wide reference, so it must not be
  # recomputed on whatever subset of the genome the figure ends up showing
  full_top <- tt$df
  sel <- tt$df
  ibd <- x$get_per_snp_group()
  if (is.null(ibd)) {
    stop("plot_ibd_tugofwar() needs a per_snp_group table", call. = FALSE)
  }
  full_ibd <- ibd     # the bottom half's bar is genome-wide too, so keep it before the crop
  if (!metric %in% names(sel))
    stop(sprintf("the top track has no '%s' column", metric), call. = FALSE)
  if (!tt$own && all(c("group") %in% names(sel)) && "group" %in% names(ibd) &&
      !any(as.character(sel$group) %in% as.character(ibd$group)))
    stop("`top` and the IBD table share no group: top has ",
         paste(unique(as.character(sel$group)), collapse = ", "), "; IBD has ",
         paste(unique(as.character(ibd$group)), collapse = ", "), call. = FALSE)
  .check_gene_request(x$get_genes(), highlight_genes)
  label_arg <- label_genes
  regs <- .resolve_groups(group, sel, ibd)
  if (!is.null(regs)) {
    if ("group" %in% names(sel)) sel <- sel[sel$group %in% regs, , drop = FALSE]
    if ("group" %in% names(ibd)) ibd <- ibd[ibd$group %in% regs, , drop = FALSE]
  }
  layout <- .select_layout(x$chrom_layout(), chroms, skip_chr)
  z <- .zoom_setup(zoom, x$get_genes(), layout, label_arg, zoom_pad,
                   x$reference_id(), genes_for_track = genes_for_track)
  # zoomed, the gene names go in the track stacked underneath, not above the panel
  label_genes <- if (!is.null(z)) FALSE
    else if (is.null(label_arg)) !is.null(highlight_genes) else isTRUE(label_arg)
  sel <- .crop_to_window(.recum(sel[is.finite(sel[[metric]]), , drop = FALSE], layout),
                         z, "SNPs in the top track")
  ibd <- .crop_to_window(.recum(ibd[is.finite(ibd$frac_pairs_ibd), , drop = FALSE], layout),
                         z, "IBD SNPs")
  genes <- if (is.null(z)) .genes_for_layout(x$get_genes(), layout, highlight_genes)
           else z$genes
  genome_frac <- sum(layout$len) / sum(x$chrom_layout()$len)

  # Only test for the column, never reach through it: a table written before the grouping
  # column was standardised names it after the grouping instead ("region"), and `$group`
  # there warns and silently unfacets. `.combine_groups()` keeps the declared group order
  # rather than the alphabet -- and `top_level` stays a factor, which matters because a
  # secondary layer carrying the facet column as bare character makes ggplot merge the
  # layers' facet values alphabetically.
  grp <- .combine_groups(if ("group" %in% names(sel)) sel$group,
                         if ("group" %in% names(ibd)) ibd$group)
  present <- levels(droplevels(grp))
  faceted <- ("group" %in% names(sel)) && length(present) > 1
  top_level <- if (faceted) factor(present[1], levels = present) else NULL
  # Both halves must carry the SAME factor. The two tracks are separate layers, and when
  # one holds `group` as a factor and the other as character, ggplot merges their facet
  # values into character and re-sorts them alphabetically, scrambling the panel order.
  if ("group" %in% names(sel)) sel$group <- factor(as.character(sel$group), levels = present)
  if ("group" %in% names(ibd)) ibd$group <- factor(as.character(ibd$group), levels = present)

  # A signed statistic (iHS, Rsb, a z-score) cannot be mirrored with its sign: the top half
  # of the mirror only has room for [0, max], so a negative value maps above the top and is
  # clipped away silently. Take the magnitude -- both tails of these statistics mean
  # selection anyway -- and say so on the axis.
  signed <- any(sel[[metric]] < 0, na.rm = TRUE)
  if (is.null(scan_abs)) scan_abs <- signed
  if (signed && !isTRUE(scan_abs))
    stop("`", metric, "` has negative values and the mirror has no room below zero on the ",
         "top half; leave `scan_abs` alone to plot the magnitude, or use plot_ibd_locus() ",
         "to keep the sign", call. = FALSE)
  if (isTRUE(scan_abs) && signed) {
    sel[[metric]] <- abs(sel[[metric]])
    metric_lab <- paste0("|", metric_label, "|")
  } else {
    metric_lab <- metric_label
  }

  # normalise each track so its tallest bar reaches the centre (y = 0). common =
  # one shared max (panels comparable); free = each group to its own max.
  # Both halves are scaled into [centre_gap, 1] rather than [0, 1]. Each track's tallest
  # bar otherwise reaches the centre line exactly, and so does its largest tick label --
  # putting the two tracks' most informative labels within a few points of each other.
  span <- 1 - centre_gap
  normalized <- scale == "free" && "group" %in% names(sel) && "group" %in% names(ibd)
  sm <- NULL
  if (normalized) {
    sm <- tapply(sel[[metric]], sel$group, max, na.rm = TRUE)
    im <- tapply(ibd$frac_pairs_ibd, ibd$group, max, na.rm = TRUE)
    sm[!is.finite(sm) | sm <= 0] <- 1
    im[!is.finite(im) | im <= 0] <- 1
    sel$.tip <- 1 - span * sel[[metric]] / sm[as.character(sel$group)]
    ibd$.tip <- -1 + span * ibd$frac_pairs_ibd / im[as.character(ibd$group)]
    sel_max <- 1; ibd_max <- 1
  } else {
    sel_max <- max(sel[[metric]], na.rm = TRUE)
    ibd_max <- max(ibd$frac_pairs_ibd, na.rm = TRUE)
    if (!is.finite(sel_max) || sel_max <= 0) sel_max <- 1
    if (!is.finite(ibd_max) || ibd_max <= 0) ibd_max <- 1
    sel$.tip <- 1 - span * .safe_scale(sel[[metric]], sel_max)
    ibd$.tip <- -1 + span * .safe_scale(ibd$frac_pairs_ibd, ibd_max)
  }
  sel_y <- function(v) 1 - span * .safe_scale(v, sel_max)
  ibd_y <- function(v) -1 + span * .safe_scale(v, ibd_max)

  # ---- peaks: the bar each half is judged against, the shading, and the centre ribbon ----
  sel$.col <- selection_colour
  ibd$.col <- ibd_colour
  ribbon_df <- NULL
  if (want_peaks) {
    panels <- if (length(present)) present else "all"
    withg <- function(d) {
      if (!"group" %in% names(d)) d$group <- panels[1]
      d$group <- as.character(d$group)
      d
    }
    top_bar <- .tug_bar(full_top, metric, peak_quantile, panels)
    bot_bar <- .tug_bar(full_ibd, "frac_pairs_ibd", peak_quantile, panels)
    st <- withg(sel); sb <- withg(ibd)
    if (isTRUE(shade)) {
      # the bar is on the raw statistic, so the shading says the same thing under
      # `scale = "free"` as under "common" even though the bar heights differ per panel
      above_t <- is.finite(st[[metric]]) & st[[metric]] > top_bar[st$group]
      above_b <- is.finite(sb$frac_pairs_ibd) & sb$frac_pairs_ibd > bot_bar[sb$group]
      above_t[is.na(above_t)] <- FALSE
      above_b[is.na(above_b)] <- FALSE
      sel$.col <- ifelse(above_t, selection_colour, .wash(selection_colour, shade_wash))
      ibd$.col <- ifelse(above_b, ibd_colour, .wash(ibd_colour, shade_wash))
    }
    if (isTRUE(ribbon)) {
      pk <- .tug_peaks(.peak_intervals(st, metric, top_bar, peak_gap),
                       .peak_intervals(sb, "frac_pairs_ibd", bot_bar, peak_gap),
                       peak_gap)
      if (is.null(pk)) {
        message("neither half of the mirror has a run above its ", peak_quantile,
                " bar inside the region drawn, so no ribbon is shown")
      } else {
        off <- layout$offset[match(pk$chr, layout$chr)]
        pk$xmin <- pk$start + off
        pk$xmax <- pk$end + off
        # a peak found on a single marker has no width, and a zero-width rectangle draws
        # nothing at all -- widen those about their midpoint until they are visible, at a
        # size taken from the span actually on screen so a zoomed panel does not lose them
        min_w <- diff(range(c(sel$cum_pos, ibd$cum_pos))) * 0.0015
        narrow <- (pk$xmax - pk$xmin) < min_w
        if (any(narrow)) {
          mid <- (pk$xmin[narrow] + pk$xmax[narrow]) / 2
          pk$xmin[narrow] <- mid - min_w / 2
          pk$xmax[narrow] <- mid + min_w / 2
        }
        cols <- c(both = "#333333", top = selection_colour, bottom = ibd_colour)
        if (!is.null(ribbon_colours)) {
          bad <- setdiff(names(ribbon_colours), names(cols))
          if (length(bad))
            stop("`ribbon_colours` takes the names both/top/bottom, not ",
                 paste(bad, collapse = ", "), call. = FALSE)
          cols[names(ribbon_colours)] <- unname(ribbon_colours)
        }
        lbl <- c(both = "both", top = paste0(tt$label, " only"), bottom = "IBD only")
        pk$cell <- factor(lbl[pk$cell], levels = unname(lbl))
        if (faceted) pk$group <- factor(as.character(pk$group), levels = present)
        else pk$group <- NULL
        names(cols) <- unname(lbl[names(cols)])
        ribbon_df <- list(d = pk, cols = cols)
      }
    }
  }

  # single left axis: selection breaks in the top half, IBD in the bottom half,
  # each tick label tinted to its track.
  sel_vals <- pretty(c(0, sel_max), 4); sel_vals <- sel_vals[sel_vals >= 0 & sel_vals <= sel_max]
  ibd_vals <- pretty(c(0, ibd_max), 4); ibd_vals <- ibd_vals[ibd_vals >= 0 & ibd_vals <= ibd_max]
  sel_lab <- if (normalized || top_percent) scales::percent(sel_vals, accuracy = 1)
             else as.character(sel_vals)
  yax <- data.frame(
    y   = c(sel_y(sel_vals), ibd_y(ibd_vals)),
    lab = c(sel_lab, scales::percent(ibd_vals, accuracy = 1)),
    col = c(rep(selection_colour, length(sel_vals)), rep(ibd_colour, length(ibd_vals))),
    stringsAsFactors = FALSE)
  yax <- yax[order(yax$y), ]

  use_ggtext <- requireNamespace("ggtext", quietly = TRUE)
  if (use_ggtext) {
    # colour each tick label via markdown so we avoid vectorised element_text()
    yax$lab <- paste0("<span style='color:", yax$col, ";'>", yax$lab, "</span>")
    ytitle <- paste0(
      "<span style='color:", selection_colour, ";'>", tt$label, " (", metric_lab, ", top)</span>",
      " / <span style='color:", ibd_colour, ";'>IBD Fraction (bottom)</span>")
    ytext_elem  <- ggtext::element_markdown()
    ytitle_elem <- ggtext::element_markdown(angle = 90)
  } else {
    ytitle <- paste0(tt$label, " ", metric_lab, " (top)  /  IBD Fraction (bottom)")
    ytext_elem  <- ggplot2::element_text()
    ytitle_elem <- ggplot2::element_text()
  }

  thr_layer <- NULL
  # The threshold table describes the object's own statistic, so it says nothing about an
  # external scan: that is read at a fixed tail (the same 1% the scan Manhattans draw) or
  # at whatever height the caller names.
  which_thr <- if (tt$own) .resolve_threshold(draw_threshold) else "none"
  if (!tt$own && !normalized) {
    lvl <- if (isTRUE(draw_threshold)) -log10(0.01)
           else if (is.numeric(draw_threshold)) draw_threshold else NULL
    if (!is.null(lvl) && lvl > sel_max) {
      # the top half only spans [0, sel_max]; a higher threshold would be drawn down in the
      # IBD half, where it would read as an IBD line
      message("the ", tt$label, " threshold (", signif(lvl, 3), ") is above every value in ",
              "the top track (max ", signif(sel_max, 3), "), so no line is drawn")
    } else if (!is.null(lvl)) {
      thr_layer <- ggplot2::geom_hline(yintercept = sel_y(lvl), colour = "firebrick",
                                       linetype = "dashed", linewidth = 0.4)
    }
  }
  if (which_thr != "none" && metric == "neg_log10_p" && !normalized) {
    all_thr <- x$get_thresholds()
    thr_layer <- list()
    for (kind in .threshold_kinds(which_thr)) {
      thr <- .threshold_line(all_thr, kind, strict = which_thr == kind)
      if (is.null(thr)) next
      if ("group" %in% names(thr) && any(!is.na(thr$group)) && "group" %in% names(sel))
        thr <- thr[thr$group %in% present, , drop = FALSE]
      if (!nrow(thr)) next
      # the selection half of the mirror is drawn upside down in [0, 1], so the line has
      # to be mapped through the same transform as the data -- and a threshold above the
      # tallest bar would land in the IBD half, where it would read as an IBD line
      thr$.y <- sel_y(thr$threshold)
      thr <- thr[thr$.y >= 0, , drop = FALSE]
      if (!nrow(thr)) next
      sty <- .THRESHOLD_STYLE[[kind]]
      thr_layer[[kind]] <- ggplot2::geom_hline(
        data = thr, ggplot2::aes(yintercept = .data$.y), inherit.aes = FALSE,
        colour = sty$colour, linetype = sty$linetype, linewidth = 0.4)
    }
    if (!length(thr_layer)) thr_layer <- NULL
  }

  # The empirical reference, drawn from the top track itself and so independent of whether
  # the object has a threshold table -- it stacks with whatever `draw_threshold` produced.
  q_layer <- NULL
  if (!is.null(top_quantile)) {
    # always the magnitude: a signed track that was not mirrored as `|metric|` has already
    # errored out above, so everything the top half draws is non-negative by this point --
    # and the uncropped table can carry a sign the cropped one happens not to
    ql <- .top_quantile_lines(full_top, metric, top_quantile, use_abs = TRUE)
    if (is.null(ql)) {
      message("the top track has no finite `", metric, "` values, so no quantile line ",
              "is drawn")
    } else {
      if (!is.na(ql$group[1]) && "group" %in% names(sel))
        ql <- ql[ql$group %in% present, , drop = FALSE]
      # each panel is scaled to its own maximum under `scale = "free"`, so the line has to
      # go through that group's transform rather than the shared one
      ql$.y <- if (normalized && !is.na(ql$group[1])) {
        # one max per group, so the division is vectorised here rather than through
        # `.safe_scale()`, which takes a single denominator
        mx <- unname(sm[as.character(ql$group)])
        r <- ql$threshold / mx
        r[!is.finite(r) | !is.finite(mx) | mx <= 0] <- 0
        1 - span * r
      } else sel_y(ql$threshold)
      n_in <- nrow(ql)
      # the top half only spans [0, sel_max]; a line above it would be drawn down in the
      # IBD half, where it would read as an IBD line
      ql <- ql[is.finite(ql$.y) & ql$.y >= 0, , drop = FALSE]
      if (n_in > nrow(ql))
        message(n_in - nrow(ql), " of ", n_in, " quantile line(s) sit above every value ",
                "drawn in the top track, so they are not shown")
      if (nrow(ql)) {
        if (!is.na(ql$group[1]) && "group" %in% names(sel))
          ql$group <- factor(as.character(ql$group), levels = present)
        else ql$group <- NULL
        q_layer <- ggplot2::geom_hline(
          data = ql, ggplot2::aes(yintercept = .data$.y), inherit.aes = FALSE,
          colour = .TOP_QUANTILE_STYLE$colour, linetype = .TOP_QUANTILE_STYLE$linetype,
          linewidth = 0.4)
      }
    }
  }

  # The ribbon sits in the empty band `centre_gap` leaves between the two tracks, so it
  # costs the figure no height and lands exactly where the halves are already compared.
  ribbon_layer <- NULL
  if (!is.null(ribbon_df)) {
    h <- max(0.012, min(0.03, centre_gap * 0.4))
    ribbon_layer <- list(
      ggplot2::geom_rect(
        data = ribbon_df$d,
        ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = -h, ymax = h,
                     fill = .data$cell),
        inherit.aes = FALSE),
      ggplot2::scale_fill_manual(values = ribbon_df$cols, name = "peak called by",
                                 drop = FALSE))
  }

  p <- ggplot2::ggplot() +
    .chr_band_layer(layout) +
    .gene_line_layer(genes) +
    (if (label_genes) .gene_label_layer(genes, if (faceted) "group" else NULL, top_level,
       footprint = .label_footprint_bp(nchar(genes$name), .fp_layout(z, layout),
                                       .fp_frac(z, genome_frac)))) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
    ggplot2::geom_segment(data = sel,
      ggplot2::aes(x = .data$cum_pos, xend = .data$cum_pos, y = 1, yend = .data$.tip,
                   colour = .data$.col), linewidth = 0.15) +
    ggplot2::geom_segment(data = ibd,
      ggplot2::aes(x = .data$cum_pos, xend = .data$cum_pos, y = -1, yend = .data$.tip,
                   colour = .data$.col), linewidth = 0.15) +
    # the bar colours are literal, so no legend is drawn for them: two strengths of one hue
    # explain themselves, and a four-key legend beside a two-track figure would not
    ggplot2::scale_colour_identity() +
    thr_layer +
    q_layer +
    ribbon_layer +
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
    # the Manhattan theme draws no legend, which is right for the bars but not for the
    # ribbon: three fills nobody can decode from the figure alone. Underneath rather than
    # beside, so a genome-wide panel keeps its width.
    (if (!is.null(ribbon_layer))
       ggplot2::theme(legend.position = "bottom", legend.direction = "horizontal")) +
    .gene_label_space(label_genes)
  if (faceted) {
    p <- p + ggplot2::facet_wrap(~ group, ncol = 1, strip.position = "right")
  } else if (!is.null(regs)) {
    p <- p + ggplot2::labs(title = paste0("tug-of-war: ", paste(regs, collapse = ", ")))
  }
  n_panels <- if (faceted) length(present) else 1L
  p <- .apply_zoom(p, z, layout, n_panels, gene_label_angle)
  attr(p, "plasgenomics_dims") <- if (is.null(z))
    .dims_genome(n_panels, genome_frac, label_genes = label_genes)
    else .dims_zoom(n_panels, attr(p, "plasgenomics_track_in"))
  p
}

# ---- drug-gene group x group triangles -----------------------------------

# SNPs from a pairwise-group table overlapping one feature interval; `pos` and the
# interval are both 0-based, the interval half-open.
.snps_in_gene <- function(pw, chr, start, end) {
  pw[pw$chr == chr & pw$pos >= start & pw$pos < end, , drop = FALSE]
}

# Turn user SNP identifiers into single-position features: either "chr:pos" ids or a
# data frame with chr + pos (+ optional name). `pos` is 0-based, so the feature it becomes
# is the half-open interval covering that one base, [pos, pos + 1).
.parse_snp_features <- function(snps) {
  if (is.data.frame(snps)) {
    if (!all(c("chr", "pos") %in% names(snps)))
      stop("`snps` data frame needs 'chr' and 'pos' columns", call. = FALSE)
    nm <- if ("name" %in% names(snps)) as.character(snps$name) else paste0(snps$chr, ":", snps$pos)
    return(data.frame(name = nm, chr = normalise_chr(snps$chr),
                      start = as.numeric(snps$pos), end = as.numeric(snps$pos) + 1,
                      stringsAsFactors = FALSE))
  }
  ids <- as.character(snps)
  parts <- strsplit(ids, ":", fixed = TRUE)
  chr <- vapply(parts, function(z) paste(z[-length(z)], collapse = ":"), "")
  pos <- suppressWarnings(as.numeric(vapply(parts, function(z) z[length(z)], "")))
  if (any(is.na(pos)))
    stop("`snps` must be 'chr:pos' ids (e.g. 'Pf3D7_07_v3:403222') or a chr/pos data frame",
         call. = FALSE)
  data.frame(name = ids, chr = normalise_chr(chr), start = pos, end = pos + 1,
             stringsAsFactors = FALSE)
}

# per-gene / per-SNP triangle values from the pairwise per-SNP table (the SNP-in-feature
# path: a gene's value aggregates the pairwise IBD of SNPs strictly inside it). Used for
# `snps=` and as the fallback when no IBD blocks / overlap table are available.
.snp_gene_triangle_df <- function(x, genes, snps, agg_fn, within = 0) {
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
      pad = within, stringsAsFactors = FALSE)
  }
  # an explicitly requested SNP is never padded -- that would silently fold in its
  # neighbours, which is not what naming a single variant asks for
  if (!is.null(snps)) feats[["snps"]] <- cbind(.parse_snp_features(snps), pad = 0)
  feats <- do.call(rbind, feats)
  if (is.null(feats) || !nrow(feats)) {
    stop("nothing to plot: give a genes track (ibd_results(genes=)) and/or snps=", call. = FALSE)
  }
  pw$chr <- normalise_chr(pw$chr)
  parts <- list()
  for (i in seq_len(nrow(feats))) {
    f <- feats[i, ]
    sub <- .snps_in_gene(pw, f$chr, f$start - f$pad, f$end + f$pad)
    if (!nrow(sub)) {
      warning(sprintf("feature '%s' has no overlapping SNPs%s; skipped", f$name,
                      if (f$pad > 0) sprintf(" within %g bp", f$pad) else ""),
              call. = FALSE)
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
#' @param within Pad each gene interval by this many bp on both sides (default `0`), on
#'   either path: it widens which IBD blocks count as overlapping the gene, and on the SNP
#'   fallback it widens which SNPs are taken as the gene's. Raising it is what makes the
#'   SNP fallback usable on a sparse panel, where a short gene may contain no genotyped
#'   SNP at all. Features named through `snps` are never padded -- naming one variant asks
#'   for that variant, not its neighbourhood.
#' @param agg For the SNP fallback, how a gene's value aggregates its in-gene SNPs:
#'   `"mean"` (default), `"median"`, or `"max"`.
#' @param individual If `TRUE`, return a **named list** of one-triangle plots (one per
#'   feature, e.g. to write a multi-page PDF) with the legend tucked into the empty
#'   upper triangle; if `FALSE` (default), one faceted grid plot.
#' @param label Draw the value in each tile.
#' @param digits Decimal places for the tile labels.
#' @param ncol Facet columns for the grid (default: ggplot2 chooses).
#' @param trans Fill-scale transform, e.g. `"identity"` (default), `"log2"`, `"sqrt"`.
#' @param colors,colours Optional colour ramp for the fill (default: the pairwise-sharing ramp).
#' @param limits Fill limits. `NULL` (default) lets each plot scale to its own values;
#'   `"shared"` pins every feature to the range across all of them, which makes the pages
#'   from `individual = TRUE` colour-comparable; or give `c(lo, hi)` explicitly (values
#'   outside are squished). A faceted grid already shares one scale across its panels.
#' @param fill_scale Optional ggplot2 fill scale that fully overrides the above.
#' @return A ggplot object (grid), or a named list of ggplot objects when `individual`.
#' @examples
#' ibd <- example_ibd_results()
#' # `within` widens the selection: the example panel is sparse, and a real gene often
#' # carries no SNP of its own, so nearby ones stand in for it
#' plot_pairwise_ibd_for_genes(ibd, genes = c("pfcrt", "pfdhps"), within = 20000)
#' @export
plot_pairwise_ibd_for_genes <- function(x, genes = NULL, snps = NULL, group = NULL,
                                     within = 0, agg = c("mean", "median", "max"),
                                     individual = FALSE,
                                     label = TRUE, digits = 2, ncol = NULL,
                                     trans = "identity", colors = NULL, limits = NULL,
                                     fill_scale = NULL,
                                        colours = NULL) {
  colors <- .alias_arg("colors", "colours")
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
    df <- .snp_gene_triangle_df(x, genes, snps, agg_fn, within = within)
  }

  # "shared" pins every feature to one scale, so the pages from individual = TRUE are
  # colour-comparable; taken over the values as drawn, so a log transform's dropped
  # non-positives do not pull the lower limit to -Inf
  if (is.character(limits)) {
    limits <- match.arg(limits, "shared")
    limits <- suppressWarnings(range(.na_nonpositive(df$frac_pairs_ibd, trans), na.rm = TRUE))
    if (!all(is.finite(limits)) || diff(limits) == 0) limits <- NULL
  }

  # axis order: the group columns' own levels (set by group_col_in_meta / set_group_order),
  # else a natural sort
  groups <- unique(c(.levels_of(df$group_a), .levels_of(df$group_b)))
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
                         limits = NULL, trans = "identity",
                         value = "frac_pairs_ibd", sublabel = NULL, na_label = NULL) {
  # the value column is named so a second statistic can reuse this geometry rather than
  # copy it; everything below works off `value` alone
  df$frac_pairs_ibd <- df[[value]]
  # place every pair on one side of the diagonal by the axis (`groups`) order, regardless
  # of how the input canonicalised its pairs -- so cells never straddle the diagonal when
  # the group order is not alphabetical
  ga <- as.character(df$group_a); gb <- as.character(df$group_b)
  ia <- match(ga, groups); ib <- match(gb, groups)
  swap <- !is.na(ia) & !is.na(ib) & ia > ib
  df$group_a <- ifelse(swap, gb, ga)
  df$group_b <- ifelse(swap, ga, gb)
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
    txt <- formatC(df$frac_pairs_ibd, format = "f", digits = digits)
    # a second line carries the count the value rests on, so a striking cell built from a
    # handful of pairs cannot be read as if it were built from hundreds
    if (!is.null(sublabel) && sublabel %in% names(df))
      txt <- paste0(txt, "\n(", df[[sublabel]], ")")
    # an untested cell gets the whole label replaced, count included -- printing "(NA)" in
    # an empty tile is worse than printing nothing
    if (!is.null(na_label)) txt[is.na(df$frac_pairs_ibd)] <- na_label
    df$.lab <- txt
    p <- p + ggplot2::geom_text(
      data = df, ggplot2::aes(label = .data$.lab, colour = .data[[".txt"]]),
      size = 2.5, lineheight = 0.9, show.legend = FALSE) +
      ggplot2::scale_colour_identity()
  }
  p
}

# place the legend(s) inside the empty upper-right triangle (ggplot2-version-safe).
# Anchored by the box's top-left, so several stacked legends grow *down* into the empty
# space rather than up over whatever sits above the panel (a dendrogram / annotation strip).
.legend_upper_triangle <- function(pos = c(0.60, 0.97)) {
  base <- ggplot2::theme(
    legend.background = ggplot2::element_rect(fill = "white", colour = NA),
    legend.key.size = grid::unit(0.9, "lines"),
    legend.title = ggplot2::element_text(size = 8),
    legend.text = ggplot2::element_text(size = 7))
  if (utils::packageVersion("ggplot2") >= "3.5.0") {
    base + ggplot2::theme(legend.position = "inside",
                          legend.position.inside = pos,
                          legend.justification.inside = c(0, 1))
  } else {
    base + ggplot2::theme(legend.position = pos, legend.justification = c(0, 1))
  }
}

.safe_scale <- function(v, mx) {
  if (!is.finite(mx) || mx <= 0) return(rep(0, length(v)))
  out <- v / mx
  out[!is.finite(out)] <- 0
  out
}

# Groups from two tables as one factor, declared order first. `c()` is not usable here:
# combining a factor with a character vector drops the factor to its integer codes, so the
# groups come back as "1", "2", ... whenever the two tables disagree about the column type.
.combine_groups <- function(a, b) {
  lv <- unique(c(if (!is.null(a)) levels(.as_group_factor(a)),
                 if (!is.null(b)) levels(.as_group_factor(b))))
  if (!length(lv)) return(NULL)
  factor(c(as.character(a), as.character(b)), levels = lv)
}

.resolve_groups <- function(group, sel, ibd) {
  if (!is.null(group)) return(group)
  regs <- .combine_groups(if ("group" %in% names(sel)) sel$group,
                          if ("group" %in% names(ibd)) ibd$group)
  if (is.null(regs)) NULL else levels(droplevels(regs))
}
