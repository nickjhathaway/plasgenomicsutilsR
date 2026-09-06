# Colour-blind checking: simulate a palette as the three kinds of dichromacy see it, and
# measure how close its nearest pair gets. The point of the plot is to make an accessibility
# claim inspectable -- reviewers see the swatches, not just a number.

# Simulations available, in the order they are drawn, with how common each condition is.
.cvd_kinds <- c(normal = "normal vision", deutan = "deuteranopia",
                protan = "protanopia", tritan = "tritanopia")
.cvd_prevalence <- c(normal = "", deutan = "~6% of men", protan = "~2% of men",
                     tritan = "~1 in 10,000")

# The simulator for one kind; `normal` is the identity so it can be treated uniformly.
.cvd_sim <- function(kind) {
  if (identical(kind, "normal")) return(function(x) x)
  .need_package("colorspace", "colour-blindness simulation")
  switch(kind,
         deutan = colorspace::deutan,
         protan = colorspace::protan,
         tritan = colorspace::tritan,
         stop("unknown colour-blindness type: ", kind, call. = FALSE))
}

#' Smallest colour difference within (or between) palettes, under colour blindness
#'
#' The worst-case pair: for one palette, the two colours a viewer would find hardest to
#' tell apart; for two, the closest a colour in one gets to a colour in the other. Reported
#' as CIEDE2000 distance -- roughly, below ~5 is a confusion risk, above ~10 is comfortable.
#'
#' @param x A vector of colours.
#' @param y Optional second vector. When given, distances are measured *between* `x` and
#'   `y` (e.g. region colours against ancestry-component colours) rather than within `x`.
#' @param vision Which vision types to report, from `"normal"`, `"deutan"`, `"protan"`,
#'   `"tritan"`.
#' @return A named numeric vector, one worst-case distance per vision type.
#' @examples
#' colour_blind_distance(color_palette(6))
#' colour_blind_distance(RColorBrewer::brewer.pal(6, "Paired"))  # deutan pairs collapse
#' @export
colour_blind_distance <- function(x, y = NULL,
                                  vision = c("normal", "deutan", "protan", "tritan")) {
  .need_package("farver", "colour_blind_distance()")
  vision <- match.arg(vision, several.ok = TRUE)
  vapply(vision, function(v) {
    f <- .cvd_sim(v)
    a <- farver::decode_colour(f(unname(x)))
    if (is.null(y)) {
      m <- farver::compare_colour(a, a, from_space = "rgb", method = "cie2000")
      if (nrow(m) < 2L) return(NA_real_)
      min(m[lower.tri(m)])
    } else {
      b <- farver::decode_colour(f(unname(y)))
      min(farver::compare_colour(a, b, from_space = "rgb", method = "cie2000"))
    }
  }, numeric(1))
}

#' @rdname colour_blind_distance
#' @export
color_blind_distance <- colour_blind_distance

# One block of swatches: the palette drawn once per vision type, worst-case distance in the
# row label so the claim and the evidence sit together.
.cvd_panel <- function(cols, title, subtitle, vision, base_size, show_distance, label_angle) {
  .need_package("ggplot2", "plot_colour_blind_check()")
  if (is.null(names(cols)) || any(!nzchar(names(cols))))
    names(cols) <- seq_along(cols)
  de <- if (show_distance) colour_blind_distance(cols, vision = vision) else NULL
  rows <- lapply(vision, function(v) {
    data.frame(vision = v, i = seq_along(cols), fill = .cvd_sim(v)(unname(cols)),
               stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows)
  lab <- vapply(vision, function(v) {
    p <- .cvd_prevalence[[v]]
    out <- .cvd_kinds[[v]]
    if (nzchar(p)) out <- paste0(out, "  (", p, ")")
    if (show_distance) out <- paste0(out, "\nΔE ", sprintf("%.1f", de[[v]]))
    out
  }, character(1))
  df$vision <- factor(df$vision, levels = rev(vision), labels = rev(lab))
  ggplot2::ggplot(df, ggplot2::aes(.data$i, .data$vision, fill = .data$fill)) +
    ggplot2::geom_tile(colour = "black", linewidth = 0.3, width = 0.96, height = 0.9) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_x_continuous(breaks = seq_along(cols), labels = names(cols),
                                expand = ggplot2::expansion(add = 0.05)) +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = label_angle,
                                          hjust = if (label_angle == 0) 0.5 else 1,
                                          size = base_size * 0.8),
      axis.text.y = ggplot2::element_text(hjust = 1, size = base_size * 0.8,
                                          lineheight = 1.05),
      plot.title = ggplot2::element_text(face = "bold", size = base_size * 1.2),
      plot.subtitle = ggplot2::element_text(size = base_size * 0.8, colour = "grey30"),
      plot.margin = ggplot2::margin(4, 14, 4, 4))
}

#' Show how a set of palettes reads under colour blindness
#'
#' Draws each palette once per vision type -- normal, deuteranope, protanope, tritanope --
#' with the worst-case CIEDE2000 distance for that row, so an accessibility claim about a
#' figure can be shown rather than asserted. Pass the colour maps a figure actually uses
#' (one per legend) and the panels come out in the order given.
#'
#' Names carry through: a named colour vector labels its swatches, and the names of
#' `palettes` become the panel titles.
#'
#' @param palettes A named list of colour vectors -- one entry per legend in the figure
#'   being checked, e.g. `list("Region colours" = region_cols, "Components" = k_cols)`.
#'   A bare colour vector is treated as a single unnamed palette.
#' @param subtitles Optional character vector, one per palette, describing where each is
#'   used (e.g. `"UMAP points and the strip over each admixture panel"`).
#' @param title Overall figure title; `NULL` for none.
#' @param caption Explanatory caption under the figure: `TRUE` (default) for the standard
#'   note on method and how to read the numbers, `FALSE` for none, or your own string.
#' @param vision Vision types to draw, from `"normal"`, `"deutan"`, `"protan"`, `"tritan"`.
#' @param show_distance Put the worst-case distance in each row label (default `TRUE`).
#' @param base_size Base font size.
#' @param label_angle Angle for the swatch labels (default `40`; use `0` for short names).
#' @return A [patchwork::patchwork] object; `ggplot2::ggsave()` or [save_plot()] it.
#' @examples
#' plot_colour_blind_check(list("A 6-colour palette" = color_palette(6)))
#' @export
plot_colour_blind_check <- function(palettes, subtitles = NULL, title = NULL,
                                    caption = TRUE,
                                    vision = c("normal", "deutan", "protan", "tritan"),
                                    show_distance = TRUE, base_size = 11,
                                    label_angle = 40) {
  .need_package("ggplot2", "plot_colour_blind_check()")
  .need_package("patchwork", "plot_colour_blind_check()")
  vision <- match.arg(vision, several.ok = TRUE)
  if (!is.list(palettes)) palettes <- list(palettes)
  if (!length(palettes)) stop("`palettes` is empty", call. = FALSE)
  nms <- names(palettes)
  if (is.null(nms)) nms <- rep("", length(palettes))
  if (!is.null(subtitles) && length(subtitles) != length(palettes))
    stop("`subtitles` must have one entry per palette", call. = FALSE)

  panels <- lapply(seq_along(palettes), function(i)
    .cvd_panel(palettes[[i]], if (nzchar(nms[i])) nms[i] else NULL,
               if (is.null(subtitles)) NULL else subtitles[i],
               vision, base_size, show_distance, label_angle))

  fig <- patchwork::wrap_plots(panels, ncol = 1)
  cap <- if (isTRUE(caption)) paste(
    "Colour blindness simulated with colorspace::deutan / protan / tritan",
    "(Machado et al. 2009). ΔE is the smallest CIEDE2000 distance between any two",
    "swatches in that row - higher is harder to confuse; above ~5 is distinguishable,",
    "above ~10 comfortably so.") else if (isFALSE(caption)) NULL else caption
  if (!is.null(cap) && requireNamespace("scales", quietly = TRUE))
    cap <- paste(strwrap(cap, width = 118), collapse = "\n")
  if (!is.null(title) || !is.null(cap))
    fig <- fig + patchwork::plot_annotation(
      title = title, caption = cap,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = base_size * 1.4),
        plot.caption = ggplot2::element_text(hjust = 0, size = base_size * 0.72,
                                             colour = "grey30")))
  fig
}

#' @rdname plot_colour_blind_check
#' @export
plot_color_blind_check <- plot_colour_blind_check
