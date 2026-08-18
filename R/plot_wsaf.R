# The per-sample distribution of heterozygous-site allele fractions, which is what
# distinguishes a dominant clone with minor companions from a mixture of comparable strains.

.WSAF_CLASS_COLOURS <- c(monoclonal = "#2271B2", dominant_clone = "#359B73",
                         mixed = "#D55E00", undetermined = "grey60")

#' Per-sample heterozygous allele-fraction distributions
#'
#' One small histogram per sample of the within-sample **minor** allele fraction at
#' heterozygous sites, which is where the strain proportions show. Reads the per-site table
#' from `plasgenomicsutils wsaf_profile --sites-out`.
#'
#' @section Reading it: At a heterozygous site the read fractions follow the strain
#'   proportions, so for `K` strains the sites fall into bands (Zhu et al. 2019). In the
#'   minor half of that picture:
#'
#'   * **mass squeezed against zero** -- one dominant strain, with a minor companion or
#'     just sequencing noise. Filtering minor alleles leaves the dominant haplotype intact,
#'     so such a sample can be treated as monoclonal.
#'   * **a band near 0.5** -- two strains of comparable size. No minor-allele filter can
#'     remove it, because a filter has to stay below 0.5 to keep the dominant call at all;
#'     forcing it would delete a real strain and leave a chimera.
#'   * **more than one band** -- more than two strains.
#'
#'   The dashed line marks `comparable_at`: a minor-allele filter has to stay below 0.5 to
#'   keep the dominant call, so mass past that line is what no filter can remove.
#'
#'   Judge the panels by how much *area* sits past the line, not by whether anything does.
#'   A handful of sites near 0.5 in an otherwise empty panel is a few repetitive or
#'   mismapped loci; a second strain shows up across the genome. That is why `wsaf_profile`
#'   decides on rates over covered sites, and reports `min_freq_needed` -- the filter that
#'   clears a given sample -- rather than a verdict alone.
#'
#' @param sites Per-site table from `wsaf_profile --sites-out` (path or data frame):
#'   `sample`, `minor_frac`, and optionally `wsmaf`.
#' @param profile Optional per-sample summary from `wsaf_profile` (path or data frame). When
#'   given, panels are ordered and coloured by its `class`, and the panel strip carries the
#'   class -- which is the quickest way to read a cohort.
#' @param samples Optional subset of samples, in the order given.
#' @param value Which fraction to draw: `"minor_frac"` (default, per-sample minor, `[0,
#'   0.5]`) or `"wsmaf"` (the population-level minor allele's within-sample frequency,
#'   `[0, 1]`, as used by the COI literature).
#' @param comparable_at Where to draw the comparable-clone line (default `0.35`, matching
#'   `wsaf_profile`).
#' @param bins Histogram bins (default `50`).
#' @param scales Panel y axes: `"free_y"` (default) gives each sample its own, which shows
#'   the *shape* of every distribution however few calls it has; `"fixed"` puts them on one
#'   axis, so panel heights compare and a sample with a handful of heterozygous sites reads
#'   as nearly empty next to one with thousands. Worth switching when the question is how
#'   *much* of the genome is heterozygous rather than where the minor mass sits -- the two
#'   are separate things, and `wsaf_profile` reports the first as `het_rate`.
#' @param ncol Panel columns (default chosen from the sample count).
#' @param colours,colors Named colours per class.
#' @return A ggplot object.
#' @seealso [plasgenomicsutilsR-package], and `plasgenomicsutils wsaf_profile`.
#' @references
#' Zhu, S. J. et al. (2019) The origins and relatedness structure of mixed infections vary
#' with local prevalence of *P. falciparum* malaria. \emph{eLife} 8, e40845.
#' \doi{10.7554/eLife.40845}
#'
#' Paschalidis, A. et al. (2023) coiaf: directly estimating complexity of infection with
#' allele frequencies. \emph{PLOS Computational Biology} 19, e1010247.
#' \doi{10.1371/journal.pcbi.1010247}
#' @examples
#' # a dominant clone with a 5% companion, and an even two-strain mixture
#' set.seed(1)
#' sites <- rbind(
#'   data.frame(sample = "dominant", minor_frac = stats::rbeta(400, 2, 40)),
#'   data.frame(sample = "even_mix", minor_frac = stats::rbeta(400, 40, 42) / 2))
#' plot_wsaf(sites)
#'
#' # one y axis across panels, to compare how many het calls each sample has
#' plot_wsaf(sites, scales = "fixed")
#' @export
plot_wsaf <- function(sites, profile = NULL, samples = NULL,
                      value = c("minor_frac", "wsmaf"), comparable_at = 0.35,
                      bins = 50, scales = c("free_y", "fixed"), ncol = NULL,
                      colours = NULL, colors = NULL) {
  colours <- .alias_arg("colours", "colors")
  .need_package("ggplot2", "plot_wsaf()")
  value <- match.arg(value)
  scales <- match.arg(scales)

  df <- as.data.frame(.read_maybe(sites, "the per-site WSAF table"),
                      stringsAsFactors = FALSE)
  if (!"sample" %in% names(df)) stop("`sites` needs a `sample` column", call. = FALSE)
  if (!value %in% names(df))
    stop("`sites` has no `", value, "` column. Available: ",
         paste(names(df), collapse = ", "), call. = FALSE)
  df$.v <- suppressWarnings(as.numeric(df[[value]]))
  df <- df[is.finite(df$.v), , drop = FALSE]
  if (!nrow(df)) stop("no finite `", value, "` values to plot", call. = FALSE)

  prof <- if (is.null(profile)) NULL else
    as.data.frame(.read_maybe(profile, "the wsaf_profile summary"), stringsAsFactors = FALSE)
  if (!is.null(prof) && !all(c("sample", "class") %in% names(prof)))
    stop("`profile` needs `sample` and `class` columns", call. = FALSE)

  if (!is.null(samples)) {
    df <- df[df$sample %in% samples, , drop = FALSE]
    if (!nrow(df)) stop("none of those samples are in `sites`", call. = FALSE)
    ord <- as.character(samples)[as.character(samples) %in% df$sample]
  } else if (!is.null(prof)) {
    # the mixtures last, so a cohort reads worst-to-best down the page
    lev <- names(.WSAF_CLASS_COLOURS)
    p <- prof[prof$sample %in% df$sample, , drop = FALSE]
    # within a class, by the filter each sample needs, so the hardest cases sit together
    key <- if ("min_freq_needed" %in% names(p)) p$min_freq_needed else
           if ("minor_q95" %in% names(p)) p$minor_q95 else 0
    p <- p[order(match(p$class, lev), key), , drop = FALSE]
    ord <- p$sample
  } else {
    ord <- sort(unique(df$sample))
  }

  lab <- ord
  fill <- rep("grey35", length(ord))
  if (!is.null(prof)) {
    cls <- prof$class[match(ord, prof$sample)]
    lab <- ifelse(is.na(cls), ord, paste0(ord, "\n", cls))
    pal <- .WSAF_CLASS_COLOURS
    if (!is.null(colours)) pal[names(colours)] <- unname(colours)
    fill <- unname(pal[cls])
    fill[is.na(fill)] <- "grey35"
  }
  df$.panel <- factor(lab[match(df$sample, ord)], levels = lab)
  df$.fill <- fill[match(df$sample, ord)]

  hi <- if (value == "wsmaf") 1 else 0.5
  if (is.null(ncol)) ncol <- max(1, min(5, ceiling(sqrt(length(ord)))))

  p <- ggplot2::ggplot(df, ggplot2::aes(.data$.v, fill = .data$.fill)) +
    # bins anchored at 0 and the range set by coord, so nothing is clipped into an NA bin
    ggplot2::geom_histogram(bins = bins, boundary = 0, colour = NA) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.01)) +
    ggplot2::coord_cartesian(xlim = c(0, hi)) +
    ggplot2::facet_wrap(~ .data$.panel, ncol = ncol, scales = scales) +
    ggplot2::labs(
      x = if (value == "minor_frac") "within-sample minor allele fraction at het sites"
          else "WSMAF (population-minor allele's within-sample frequency)",
      y = "het sites") +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(size = 7, lineheight = 1.1))

  if (is.finite(comparable_at) && comparable_at < hi) {
    p <- p + ggplot2::geom_vline(xintercept = comparable_at, linetype = "dashed",
                                 colour = "grey30", linewidth = 0.3)
    if (value == "wsmaf")
      p <- p + ggplot2::geom_vline(xintercept = 1 - comparable_at, linetype = "dashed",
                                   colour = "grey30", linewidth = 0.3)
  }
  n_rows <- ceiling(length(ord) / ncol)
  attr(p, "plasgenomics_dims") <- c(width = min(12, 2.4 * ncol + 0.6),
                                    height = min(14, 1.5 * n_rows + 0.8))
  p
}
