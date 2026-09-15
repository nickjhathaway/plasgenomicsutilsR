# The block-extension statistic as a group x group triangle, the length counterpart of the
# sharing triangles that plot_pairwise_ibd_for_genes() draws.
#
# The diagonal is the within-group extension that ibd_block_extension_test() reports. The
# off-diagonal is the question the diagonal cannot answer: whether the long-segment haplotype
# is shared *across* groups, or stops at each group's border.

# Diverging fill on log2, centred at 1, because the value is a ratio: 2 and 0.5 are the same
# distance from "no extension" and a sequential ramp would hide that. Limits are made
# symmetric in log space for the same reason.
.extension_fill_scale <- function(limits, colors = NULL, fill_scale = NULL,
                                  name = "block extension") {
  if (!is.null(fill_scale)) return(fill_scale)
  .need_package("scales", "the block-extension fill scale")
  if (is.null(colors)) colors <- c("#3A5FCD", "#9DBEEA", "grey93", "#F3A582", "#B2182B")
  ggplot2::scale_fill_gradientn(
    colours = colors, transform = "log2", limits = limits,
    name = paste0(name, "\n(ratio, log2)"),
    # an untested cell is absent, not average: grey it out rather than letting na.value fall
    # on a colour that reads as a real number
    na.value = "grey97",
    breaks = function(lims) {
      b <- 2 ^ seq(floor(log2(lims[1])), ceiling(log2(lims[2])))
      b[b >= lims[1] & b <= lims[2]]
    },
    labels = function(b) formatC(b, format = "fg", digits = 2))
}

#' Block extension between every pair of groups, as a triangle per locus
#'
#' The length counterpart of [plot_pairwise_ibd_for_genes()]. That plot fills each group-pair
#' cell with the *fraction of pairs sharing* a gene; this one fills it with how much **longer**
#' those shared segments are than the sharing pairs' own genome-wide background, the
#' `paired_ratio` of [ibd_block_extension_test()].
#'
#' The diagonal is the within-group result. The off-diagonal is what the diagonal cannot say:
#' whether the long segments continue across a group boundary, which is the difference between
#' a haplotype that is spreading between regions and one that is expanding inside each.
#'
#' @section Reading the tiles:
#' The fill diverges around **1** on a log2 scale, so 2 and 0.5 sit the same distance from the
#' middle. Note that 1 is *not* the null -- every locus is inflated by length-biased sampling,
#' so a whole triangle sitting above 1 is the expected picture, not a finding. Use
#' [ibd_block_extension_scan()] to say what a given ratio is worth in a given group.
#'
#' Cells are grey where the group pair had fewer than `min_pairs` sharing pairs, which is
#' common and is not the same as a ratio near 1. Every drawn cell carries its pair count under
#' the ratio, because a cell resting on five pairs and one resting on five hundred look
#' identical otherwise and routinely differ by more than the colour does.
#'
#' @param x An [IbdResults] built with `blocks =` and `meta =`, or a data frame already
#'   returned by `ibd_block_extension_test(..., pairs = "all")`.
#' @param loci Loci to draw, as [ibd_block_extension_test()] takes them. Ignored when `x` is
#'   already a result table.
#' @param min_pairs Group pairs with fewer sharing pairs are left grey (default `5`).
#' @param individual Return a named list of one plot per locus instead of a faceted grid.
#' @param label Draw the ratio and pair count on each tile (default `TRUE`).
#' @param digits Decimal places for the ratio (default `2`).
#' @param ncol Columns in the facet grid.
#' @param limits Fill limits, or `"shared"` (default) to pin every locus to one symmetric
#'   scale so the pages are comparable.
#' @param colors,colours Fill ramp, low to high through the midpoint.
#' @param fill_scale A complete `ggplot2` fill scale, replacing the built-in one.
#' @param group,within,sharing,min_ref_blocks,meta Passed to [ibd_block_extension_test()].
#' @return A `ggplot`, or a named list of them when `individual = TRUE`.
#' @seealso [ibd_block_extension_test()], [ibd_block_extension_scan()],
#'   [plot_pairwise_ibd_for_genes()] for the sharing-fraction counterpart.
#' @examples
#' \dontrun{
#' ibd <- ibd_results(blocks = "hmm.txt", meta = meta, group_col_in_meta = "region")
#' plot_pairwise_block_extension(ibd, c("pfcrt", "pfgch1"))
#' }
#' @export
plot_pairwise_block_extension <- function(x, loci = NULL, min_pairs = 5L, group = NULL,
                                          within = 0, sharing = c("overlap", "complete"),
                                          min_ref_blocks = 1L, individual = FALSE,
                                          label = TRUE, digits = 2, ncol = NULL,
                                          limits = "shared", colors = NULL,
                                          fill_scale = NULL, meta = NULL, colours = NULL) {
  colors <- .alias_arg("colors", "colours")
  .need_package("ggplot2", "plot_pairwise_block_extension()")
  sharing <- match.arg(sharing)
  meta <- .normalise_meta(meta)

  res <- if (inherits(x, "IbdResults")) {
    ibd_block_extension_test(x, loci = loci, group = group, within = within,
                             sharing = sharing, pairs = "all",
                             min_ref_blocks = min_ref_blocks, min_pairs = min_pairs,
                             adjust = "none", keep_pair_ratios = FALSE, meta = meta)
  } else {
    as.data.frame(x, stringsAsFactors = FALSE)
  }
  .require_cols(res, c("locus", "group_a", "group_b", "paired_ratio", "n_pairs"),
                "block-extension result")
  if (!nrow(res)) stop("no group pair had enough sharing pairs to draw", call. = FALSE)

  df <- data.frame(gene = factor(as.character(res$locus),
                                 levels = if (is.factor(res$locus)) levels(res$locus)
                                          else unique(as.character(res$locus))),
                   group_a = as.character(res$group_a),
                   group_b = as.character(res$group_b),
                   paired_ratio = as.numeric(res$paired_ratio),
                   n_pairs = as.integer(res$n_pairs), stringsAsFactors = FALSE)

  # every group that appears on either end, in the object's own order where it has one
  groups <- unique(c(.levels_of(res$group_a), .levels_of(res$group_b)))
  if (inherits(x, "IbdResults") && !is.null(x$get_group_order())) {
    ord <- x$get_group_order()
    groups <- c(intersect(ord, groups), setdiff(groups, ord))
  }

  if (is.character(limits)) {
    limits <- match.arg(limits, "shared")
    v <- df$paired_ratio[is.finite(df$paired_ratio) & df$paired_ratio > 0]
    # symmetric in log space, so the midpoint of the ramp really is 1
    limits <- if (!length(v)) NULL else {
      h <- max(abs(log2(v)))
      if (!is.finite(h) || h == 0) NULL else c(2^-h, 2^h)
    }
  }
  fs <- .extension_fill_scale(limits, colors = colors, fill_scale = fill_scale)
  txt_cols <- if (is.null(colors)) c("#3A5FCD", "#9DBEEA", "grey93", "#F3A582", "#B2182B")
              else colors
  feats <- levels(df$gene)

  one <- function(d, title = NULL, legend_inside = FALSE) {
    # fill in the missing cells so an untested group pair is drawn grey rather than dropped,
    # which would leave a hole the eye reads as the axis being wrong
    full <- expand.grid(group_a = groups, group_b = groups, stringsAsFactors = FALSE)
    full <- full[match(full$group_a, groups) <= match(full$group_b, groups), , drop = FALSE]
    key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "\r")
    full$paired_ratio <- d$paired_ratio[match(key(full$group_a, full$group_b),
                                              key(d$group_a, d$group_b))]
    full$n_pairs <- d$n_pairs[match(key(full$group_a, full$group_b),
                                    key(d$group_a, d$group_b))]
    full$gene <- d$gene[1]
    .triangle_gg(full, groups, fs, label, digits, title = title,
                 legend_inside = legend_inside, colours = txt_cols, limits = limits,
                 trans = "log2", value = "paired_ratio",
                 sublabel = if (label) "n_pairs" else NULL, na_label = "")
  }

  if (individual) {
    plots <- lapply(feats, function(g) {
      pg <- one(df[df$gene == g, , drop = FALSE], title = g, legend_inside = TRUE)
      attr(pg, "plasgenomics_dims") <- .dims_triangles(1L, length(groups), 1L)
      pg
    })
    names(plots) <- feats
    return(plots)
  }
  d <- do.call(rbind, lapply(feats, function(g) {
    dd <- df[df$gene == g, , drop = FALSE]
    full <- expand.grid(group_a = groups, group_b = groups, stringsAsFactors = FALSE)
    full <- full[match(full$group_a, groups) <= match(full$group_b, groups), , drop = FALSE]
    key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "\r")
    k <- match(key(full$group_a, full$group_b), key(dd$group_a, dd$group_b))
    full$paired_ratio <- dd$paired_ratio[k]
    full$n_pairs <- dd$n_pairs[k]
    full$gene <- factor(g, levels = feats)
    full
  }))
  p <- .triangle_gg(d, groups, fs, label, digits, colours = txt_cols, limits = limits,
                    trans = "log2", value = "paired_ratio",
                    sublabel = if (label) "n_pairs" else NULL, na_label = "") +
    ggplot2::facet_wrap(~ .data$gene, ncol = ncol) +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"))
  attr(p, "plasgenomics_dims") <- .dims_triangles(length(feats), length(groups), ncol)
  p
}
