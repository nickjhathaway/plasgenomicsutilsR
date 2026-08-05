# Population differentiation between metadata groups, per SNP and per pair: Jost's D,
# Nei's Gst, Hedrick's standardized G'st, and Hudson's Fst. Plus a group x group summary
# (mean / median / top-percentile / max across SNPs), a triangle heatmap with an optional
# dendrogram and metadata annotation, and a helper to pick the most differentiating SNPs.
#
# D / G'st use the Nei & Chesser (1983) bias-corrected within/total gene diversities
# Hs, Ht (genotypes are alt-allele dosages, so gene copies = 2 x called). Fst is the
# Hudson estimator (Bhatia et al. 2013), a robust "ratio of averages" pairwise Fst.

STATISTIC_LABELS <- c(jost_d = "Jost's D", gst_hedrick = "Hedrick's G'st",
                      fst = "Hudson's Fst")

# coerce a genotype override (a run_ld_prune() list or a matrix) to a labelled matrix
.coerce_geno <- function(g) {
  if (is.list(g) && !is.null(g$genotype)) {
    m <- as.matrix(g$genotype)
    if (is.null(rownames(m)) && !is.null(g$sample.id)) rownames(m) <- g$sample.id
    return(m)
  }
  as.matrix(g)
}

# per-group ALT-allele frequency and gene-copy count, one row per group
.group_freqs <- function(G, grp, levs) {
  P <- matrix(NA_real_, length(levs), ncol(G), dimnames = list(levs, colnames(G)))
  N <- P
  for (l in levs) {
    Gi <- G[which(grp == l), , drop = FALSE]
    called <- colSums(!is.na(Gi))
    P[l, ] <- colSums(Gi, na.rm = TRUE) / (2 * called)
    N[l, ] <- 2 * called                        # gene copies (0 where nothing called)
  }
  list(P = P, N = N)
}

# one differentiation statistic for a pair of groups, vectorised over SNPs
.pair_diff <- function(pa, pb, Na, Nb, statistic, clamp) {
  ok <- Na > 1 & Nb > 1 & is.finite(pa) & is.finite(pb)
  if (statistic == "fst") {
    num <- (pa - pb)^2 - pa * (1 - pa) / (Na - 1) - pb * (1 - pb) / (Nb - 1)
    den <- pa * (1 - pb) + pb * (1 - pa)
    v <- num / den
    v[!is.finite(den) | den == 0] <- 0
  } else {
    ha <- 1 - (pa^2 + (1 - pa)^2)                # within-group gene diversity (biallelic)
    hb <- 1 - (pb^2 + (1 - pb)^2)
    Nh <- 2 / (1 / Na + 1 / Nb)                  # harmonic mean gene copies
    Hs <- (Nh / (Nh - 1)) * (ha + hb) / 2        # bias-corrected within
    pbar <- (pa + pb) / 2
    Ht <- 1 - (pbar^2 + (1 - pbar)^2) + Hs / (Nh * 2)   # bias-corrected total (k = 2)
    v <- switch(statistic,
                jost_d      = (Ht - Hs) / (1 - Hs) * 2,
                gst_hedrick = ((Ht - Hs) / Ht) * (1 + Hs) / (1 - Hs))
    v[!is.finite(v)] <- NA_real_                 # e.g. monomorphic (Ht = 0)
  }
  v[!ok] <- NA_real_
  if (clamp) v[is.finite(v) & v < 0] <- 0
  v
}

#' Population differentiation between metadata groups, per SNP
#'
#' Computes a per-SNP differentiation statistic for every pair of `group` levels, from
#' alt-allele dosages. Run it on the **full, unpruned** genotypes -- LD-pruning removes
#' the differentiating SNPs you want here. Feed the result to [pop_diff_matrix()] /
#' [plot_diff_heatmap()] for a group summary, or [top_differentiating_snps()] for markers.
#'
#' @param x A [PopStructure] (uses its full genotype matrix for the active samples) or a
#'   genotype matrix (samples x SNPs, 0/1/2, `NA`).
#' @param group Metadata column (for a `PopStructure`) or a per-sample grouping vector
#'   (for a matrix).
#' @param statistic `"jost_d"` (Jost's D, default), `"gst_hedrick"` (Hedrick's
#'   standardized G'st), or `"fst"` (Hudson's Fst); all in \[0, 1]. Plain Nei's Gst is not
#'   offered -- it is strongly deflated when within-group diversity is high (typical
#'   genome-wide in *P. falciparum*); its standardized form G'st is provided instead.
#' @param meta When `x` is a matrix, an optional data frame with a `sample` column plus
#'   `group`; otherwise `group` is a vector aligned to the matrix rows.
#' @param clamp Clamp small negative estimates to 0 (default `TRUE`).
#' @param genotype Optional genotype matrix (samples x SNPs) or [run_ld_prune()] list to
#'   use **instead** of a `PopStructure`'s stored matrix -- pass the **full, unpruned**
#'   set here (e.g. `run_ld_prune(vcf, prune = FALSE)`) so differentiation is measured on
#'   every SNP while PCA/UMAP keep using the pruned matrix. Ignored when `x` is a matrix.
#' @return A `pop_diff` object: a list with `D` (a SNP x pair matrix of the statistic),
#'   `snp`, `groups`, `pairs`, `statistic`, and the group `freqs`.
#' @references
#' Jost, L. (2008) G_ST and its relatives do not measure differentiation.
#' \emph{Molecular Ecology} 17, 4015-4026. \doi{10.1111/j.1365-294X.2008.03887.x}
#'
#' Nei, M. & Chesser, R. K. (1983) Estimation of fixation indices and gene diversities.
#' \emph{Annals of Human Genetics} 47, 253-259. \doi{10.1111/j.1469-1809.1983.tb00993.x}
#'
#' Hedrick, P. W. (2005) A standardized genetic differentiation measure.
#' \emph{Evolution} 59, 1633-1638. \doi{10.1111/j.0014-3820.2005.tb01814.x}
#'
#' Hudson, R. R., Slatkin, M. & Maddison, W. P. (1992) Estimation of levels of gene flow
#' from DNA sequence data. \emph{Genetics} 132, 583-589. \doi{10.1093/genetics/132.2.583}
#'
#' Bhatia, G., Patterson, N., Sankararaman, S. & Price, A. L. (2013) Estimating and
#' interpreting F_ST: the impact of rare variants. \emph{Genome Research} 23, 1514-1521.
#' \doi{10.1101/gr.154831.113}
#' @export
pop_diff <- function(x, group = NULL,
                     statistic = c("jost_d", "gst_hedrick", "fst"),
                     meta = NULL, clamp = TRUE, genotype = NULL) {
  statistic <- match.arg(statistic)
  if (inherits(x, "PopStructure")) {
    G <- if (is.null(genotype)) x$genotype() else .coerce_geno(genotype)
    meta <- x$get_meta()
    if (is.null(group)) group <- setdiff(names(meta), "sample")[1]
    grp <- as.character(meta[[group]])[match(rownames(G), meta$sample)]
    levs <- .levels_of(meta[[group]])
  } else {
    G <- as.matrix(x)
    if (!is.null(meta) && length(group) == 1 && group %in% names(meta)) {
      grp <- as.character(meta[[group]])[match(rownames(G), meta$sample)]
      levs <- .levels_of(meta[[group]])
    } else {
      grp <- as.character(group)
      if (length(grp) != nrow(G)) stop("`group` must align to the genotype rows", call. = FALSE)
      levs <- sort(unique(grp[!is.na(grp)]))
    }
  }
  if (length(levs) < 2) stop("need at least two groups", call. = FALSE)

  gf <- .group_freqs(G, grp, levs)
  pairs <- utils::combn(levs, 2)
  D <- matrix(NA_real_, ncol(G), ncol(pairs),
              dimnames = list(colnames(G), apply(pairs, 2, paste, collapse = " vs ")))
  for (k in seq_len(ncol(pairs))) {
    a <- pairs[1, k]; b <- pairs[2, k]
    D[, k] <- .pair_diff(gf$P[a, ], gf$P[b, ], gf$N[a, ], gf$N[b, ], statistic, clamp)
  }
  structure(list(D = D, snp = colnames(G), groups = levs,
                 pairs = data.frame(a = pairs[1, ], b = pairs[2, ], stringsAsFactors = FALSE),
                 statistic = statistic,
                 group_col = if (length(group) == 1 && is.character(group)) group else NULL,
                 freqs = gf),
            class = c(statistic, "pop_diff"))
}

#' Per-SNP Jost's D between metadata groups
#'
#' Convenience wrapper for `pop_diff(..., statistic = "jost_d")`. See [pop_diff()].
#'
#' @inheritParams pop_diff
#' @return A `pop_diff` object (of Jost's D).
#' @export
jost_d <- function(x, group = NULL, meta = NULL, clamp = TRUE) {
  pop_diff(x, group = group, statistic = "jost_d", meta = meta, clamp = clamp)
}

#' @export
print.pop_diff <- function(x, ...) {
  cat("<pop_diff>", STATISTIC_LABELS[[x$statistic]], "--", length(x$snp), "SNPs,",
      length(x$groups), "groups,", ncol(x$D), "pairwise comparisons\n")
  cat("  mean per pair:", paste(sprintf("%.3f", colMeans(x$D, na.rm = TRUE)),
                                collapse = " "), "\n")
  invisible(x)
}

#' Group x group differentiation summary matrix
#'
#' Collapses the per-SNP pairwise values into a symmetric group-by-group matrix. Because
#' most of the genome is barely differentiated in *P. falciparum*, a genome-wide `"mean"`
#' looks near-zero; `"top_mean"` (mean of the top `top` fraction of SNPs per pair) or
#' `"max"` surface the differentiation that lives in a minority of loci.
#'
#' @param pd A [pop_diff()] / [jost_d()] result.
#' @param stat `"mean"` (default), `"median"`, `"top_mean"`, or `"max"`.
#' @param top Fraction of highest-value SNPs to average for `stat = "top_mean"`
#'   (default `0.05`, i.e. the top 5%).
#' @return A symmetric numeric matrix (0 on the diagonal).
#' @export
pop_diff_matrix <- function(pd, stat = c("mean", "median", "top_mean", "max"), top = 0.05) {
  stat <- match.arg(stat)
  summ <- function(v) {
    v <- v[is.finite(v)]
    if (!length(v)) return(NA_real_)
    switch(stat,
           mean     = mean(v),
           median   = stats::median(v),
           max      = max(v),
           top_mean = { k <- max(1L, ceiling(top * length(v)))
                        mean(sort(v, decreasing = TRUE)[seq_len(k)]) })
  }
  g <- pd$groups
  M <- matrix(0, length(g), length(g), dimnames = list(g, g))
  for (k in seq_len(nrow(pd$pairs))) {
    v <- summ(pd$D[, k])
    M[pd$pairs$a[k], pd$pairs$b[k]] <- v
    M[pd$pairs$b[k], pd$pairs$a[k]] <- v
  }
  M
}

#' @rdname pop_diff_matrix
#' @export
jost_d_matrix <- function(pd, stat = "mean", top = 0.05) pop_diff_matrix(pd, stat, top)

#' Group-pair differentiation summary across statistics
#'
#' One tidy table of every pairwise comparison with, for each statistic and summary, its
#' value -- so Jost's D, Nei's Gst, Hedrick's G'st and Hudson's Fst can be read together.
#'
#' @inheritParams pop_diff
#' @param statistics Measures to include (any of `"jost_d"`, `"gst_hedrick"`, `"fst"`;
#'   default all three). Nei's Gst is not offered -- see [pop_diff()].
#' @param stats Per-pair summaries to include (`"mean"`, `"median"`, `"top_mean"`,
#'   `"max"`; default mean + top_mean + max).
#' @param top Fraction of top SNPs for `"top_mean"` (default `0.05`).
#' @param genotype Optional full/unpruned genotype override (see [pop_diff()]).
#' @return A data frame: `a`, `b`, `n_snps`, then one `statistic_stat` column each.
#' @examples
#' \dontrun{
#' ps <- example_pop_structure("africa")
#' pop_diff_table(ps, group = "site")
#' }
#' @export
pop_diff_table <- function(x, group = NULL,
                           statistics = c("jost_d", "gst_hedrick", "fst"),
                           stats = c("mean", "top_mean", "max"), top = 0.05,
                           meta = NULL, clamp = TRUE, genotype = NULL) {
  statistics <- match.arg(statistics, c("jost_d", "gst_hedrick", "fst"),
                          several.ok = TRUE)
  stats <- match.arg(stats, c("mean", "median", "top_mean", "max"), several.ok = TRUE)
  summ <- function(v, st) {
    v <- v[is.finite(v)]
    if (!length(v)) return(NA_real_)
    switch(st, mean = mean(v), median = stats::median(v), max = max(v),
           top_mean = { k <- max(1L, ceiling(top * length(v)))
                        mean(sort(v, decreasing = TRUE)[seq_len(k)]) })
  }
  out <- NULL
  for (s in statistics) {
    pd <- pop_diff(x, group = group, statistic = s, meta = meta, clamp = clamp,
                   genotype = genotype)
    if (is.null(out))
      out <- data.frame(a = pd$pairs$a, b = pd$pairs$b,
                        n_snps = colSums(is.finite(pd$D)),
                        stringsAsFactors = FALSE, row.names = NULL)
    for (st in stats) out[[paste0(s, "_", st)]] <- round(apply(pd$D, 2, summ, st = st), 4)
  }
  out
}

# resolve `annotate` into a named list of value-vectors aligned to the ordered groups
.resolve_annotations <- function(annotate, ord, group_col, meta) {
  if (is.null(annotate)) return(NULL)
  specs <- if (is.character(annotate)) stats::setNames(as.list(annotate), annotate) else annotate
  out <- list()
  for (i in seq_along(specs)) {
    a <- specs[[i]]
    nm <- names(specs)[i]
    if (is.character(a) && length(a) == 1) {                 # a metadata column name
      if (is.null(meta) || is.null(group_col) || !all(c(a, group_col) %in% names(meta)))
        stop("supply `meta` containing the grouping column and the `annotate` column(s)",
             call. = FALSE)
      if (is.null(nm) || !nzchar(nm)) nm <- a
      out[[nm]] <- as.character(meta[[a]])[match(ord, as.character(meta[[group_col]]))]
    } else {                                                 # a named group -> value vector
      if (is.null(nm) || !nzchar(nm)) nm <- paste0("annotation", i)
      out[[nm]] <- unname(a[ord])
    }
  }
  out
}

# one annotation colour strip aligned above the heatmap columns
.annotation_panel <- function(name, vals, ord, cols, base_size) {
  adf <- data.frame(col = factor(ord, levels = ord), y = 1,
                    value = factor(vals, levels = unique(vals)), stringsAsFactors = FALSE)
  if (is.null(cols)) cols <- meta_colors(data.frame(value = adf$value))$value
  ggplot2::ggplot(adf, ggplot2::aes(.data$col, .data$y, fill = .data$value)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::scale_fill_manual(values = cols, name = name) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(legend.position = "right")
}

# horizontal-dendrogram segment coordinates from an hclust (leaves at x = 1..n in order)
.dendro_segments <- function(hc) {
  m <- hc$merge; h <- hc$height; n <- nrow(m) + 1L
  xleaf <- match(seq_len(n), hc$order)           # leaf i -> x position
  nodex <- numeric(nrow(m))
  getx <- function(e) if (e < 0) xleaf[-e] else nodex[e]
  gety <- function(e) if (e < 0) 0 else h[e]
  segs <- vector("list", nrow(m))
  for (k in seq_len(nrow(m))) {
    a <- m[k, 1]; b <- m[k, 2]
    xa <- getx(a); xb <- getx(b); ya <- gety(a); yb <- gety(b)
    nodex[k] <- (xa + xb) / 2
    segs[[k]] <- data.frame(x = c(xa, xb, xa), xend = c(xa, xb, xb),
                            y = c(ya, yb, h[k]), yend = c(h[k], h[k], h[k]))
  }
  do.call(rbind, segs)
}

#' Per-SNP differentiation in long form
#'
#' Unpacks a [pop_diff()] result (a SNP x pair matrix) into a tidy long table, parsing the
#' `chr:pos` SNP ids into genomic coordinates -- so you can look at differentiation SNP by
#' SNP or feed it to [plot_diff_manhattan()].
#'
#' @param pd A [pop_diff()] / [jost_d()] result (its SNPs must be `chr:pos` ids, as from
#'   [run_ld_prune()]).
#' @return A tibble with `snp`, `chr`, `pos`, `a`, `b`, `pair`, `statistic`, and `value`
#'   (one row per SNP x pair).
#' @export
pop_diff_snps <- function(pd) {
  ids <- pd$snp
  chr <- sub(":[^:]*$", "", ids)
  pos <- suppressWarnings(as.numeric(sub(".*:", "", ids)))
  np  <- ncol(pd$D)
  long <- data.frame(
    snp   = rep(ids, times = np),
    chr   = rep(chr, times = np),
    pos   = rep(pos, times = np),
    a     = rep(pd$pairs$a, each = nrow(pd$D)),
    b     = rep(pd$pairs$b, each = nrow(pd$D)),
    value = as.vector(pd$D),
    stringsAsFactors = FALSE)
  long$pair <- paste(long$a, "vs", long$b)
  long$statistic <- pd$statistic
  tibble::as_tibble(long[, c("snp", "chr", "pos", "a", "b", "pair", "statistic", "value")])
}

#' Genome-wide differentiation Manhattan plot
#'
#' The per-SNP differentiation statistic along the genome, styled like the IBD Manhattan
#' plots. By default it collapses the pairwise comparisons to one value per SNP
#' (`combine = "max"`, the strongest differentiation at that SNP in any pair); pass
#' `pair` to plot a single group pair instead.
#'
#' @param pd A [pop_diff()] / [jost_d()] result with `chr:pos` SNP ids.
#' @param pair Optional single group pair to plot: `c("A", "B")` or `"A vs B"`. Default
#'   combines all pairs per SNP.
#' @param combine How to collapse pairs per SNP when `pair` is `NULL`: `"max"` (default)
#'   or `"mean"`.
#' @param reference Reference id for the chromosome layout (default `"pf3d7"`).
#' @param chroms,skip_chr Optional chromosomes to keep / drop (the rest re-laid-out).
#' @param point_size,point_alpha Point aesthetics.
#' @param colours Optional length-2 colour vector for the alternating chromosome bands.
#' @return A ggplot object.
#' @export
plot_diff_manhattan <- function(pd, pair = NULL, combine = c("max", "mean"),
                                reference = "pf3d7", chroms = NULL, skip_chr = NULL,
                                point_size = 0.6, point_alpha = 0.6, colours = NULL) {
  .need_package("ggplot2", "plot_diff_manhattan()")
  long <- pop_diff_snps(pd)
  long <- long[is.finite(long$value) & is.finite(long$pos), , drop = FALSE]
  if (!is.null(pair)) {
    keys <- if (length(pair) == 2) paste(pair, collapse = " vs ") else pair
    keys <- c(keys, if (length(pair) == 2) paste(rev(pair), collapse = " vs "))
    df <- long[long$pair %in% keys, c("chr", "pos", "value"), drop = FALSE]
    if (!nrow(df)) {
      stop("no such pair; available: ", paste(unique(long$pair), collapse = ", "), call. = FALSE)
    }
    title <- keys[1]
  } else {
    combine <- match.arg(combine)
    fn <- if (combine == "max") max else mean
    ag <- stats::aggregate(value ~ snp + chr + pos, data = long, FUN = fn)
    df <- ag[, c("chr", "pos", "value")]
    title <- sprintf("%s across %d group pairs", combine, ncol(pd$D))
  }
  df$chr <- normalise_chr(df$chr)
  layout <- .select_layout(.chrom_layout(reference), chroms, skip_chr)
  df <- .attach_band(.recum(df, layout), layout)
  p <- ggplot2::ggplot(df, ggplot2::aes(.data$cum_pos, .data$value)) +
    .chr_band_layer(layout) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$.band), size = point_size,
                        alpha = point_alpha, stroke = 0) +
    ggplot2::scale_colour_manual(values = .band_colours(colours)) +
    .chr_axis(layout) +
    ggplot2::labs(x = "chromosome", y = STATISTIC_LABELS[[pd$statistic]], title = title) +
    .manhattan_theme()
  attr(p, "plasgenomics_dims") <- .dims_genome(
    1L, sum(layout$len) / sum(.chrom_layout(reference)$len))
  p
}

#' Triangle heatmap of group x group differentiation
#'
#' A lower-triangle heatmap of a group summary (see [pop_diff_matrix()]), styled like the
#' drug-gene triangles: the fill legend sits in the empty upper triangle. Optionally adds
#' a clustering dendrogram across the top and a metadata annotation strip (e.g. colour
#' each site by its country).
#'
#' @param pd A [pop_diff()] / [jost_d()] result.
#' @param stat,top Summary passed to [pop_diff_matrix()] (`"mean"`, `"median"`,
#'   `"top_mean"`, `"max"`; `top` for the top-percentile mean).
#' @param cluster Order groups by hierarchical clustering of the matrix.
#' @param dendrogram Draw the clustering dendrogram across the top (needs `cluster`).
#' @param triangle Show only the lower triangle (with the legend in the empty corner).
#' @param annotate One or more annotations drawn as colour strips above the columns:
#'   a metadata column name (looked up in `meta`), a vector of several column names, a
#'   named `group -> value` vector, or a named list mixing these (the names label the
#'   strips/legends). Each strip aligns to the clustered column order.
#' @param meta Metadata data frame for resolving `annotate` column names.
#' @param annotate_colours Named list `annotation -> (value -> colour)` giving custom
#'   colours per annotation (unlisted annotations get an automatic palette).
#' @param group_colours Optional named `group -> colour` vector used to colour the
#'   dendrogram leaf tips (so the tree carries the same colours as the UMAP / admixture);
#'   defaults to an automatic palette. A `PopStructure` passes its shared map here.
#' @param label Print the value in each cell (default `FALSE`; the text can distract).
#' @param digits Cell-label digits.
#' @param colors Fill ramp (low -> high differentiation).
#' @param trans Fill transform, e.g. `"identity"` (default) or `"sqrt"` to lift a scale
#'   dominated by near-zero values.
#' @param base_size Base font size.
#' @param ... For `plot_jost_d_heatmap()`, arguments forwarded to `plot_diff_heatmap()`.
#' @return A ggplot (or \pkg{patchwork}) object.
#' @export
plot_diff_heatmap <- function(pd, stat = c("mean", "median", "top_mean", "max"), top = 0.05,
                              cluster = TRUE, dendrogram = TRUE, triangle = TRUE,
                              annotate = NULL, meta = NULL, annotate_colours = NULL,
                              group_colours = NULL, label = FALSE, digits = 2,
                              colors = c("white", "#fde0dd", "#fa9fb5", "#c51b8a", "#7a0177"),
                              trans = "identity", base_size = 11) {
  .need_package("ggplot2", "plot_diff_heatmap()")
  stat <- match.arg(stat)
  M <- pop_diff_matrix(pd, stat, top)
  n <- nrow(M)
  hc <- if (cluster && n > 2) stats::hclust(stats::as.dist(max(M, na.rm = TRUE) - M)) else NULL
  ord <- if (is.null(hc)) rownames(M) else rownames(M)[hc$order]

  df <- expand.grid(row = ord, col = ord, stringsAsFactors = FALSE)
  df$D <- M[cbind(df$row, df$col)]
  df$xi <- match(df$col, ord); df$yi <- match(df$row, ord)
  if (triangle) df <- df[df$yi > df$xi, , drop = FALSE]        # strictly lower triangle
  df$col <- factor(df$col, levels = ord)
  df$row <- factor(df$row, levels = rev(ord))
  stat_lab <- c(mean = "mean", median = "median", max = "max",
                top_mean = sprintf("top %g%% mean", top * 100))[[stat]]
  lab <- sprintf("%s\n(%s)", STATISTIC_LABELS[[pd$statistic]], stat_lab)

  fill_args <- list(colours = colors, name = lab, limits = c(0, max(df$D, na.rm = TRUE)))
  if (!identical(trans, "identity")) {
    tn <- if (utils::packageVersion("ggplot2") >= "3.5.0") "transform" else "trans"
    fill_args[[tn]] <- trans
  }
  hm <- ggplot2::ggplot(df, ggplot2::aes(.data$col, .data$row, fill = .data$D)) +
    ggplot2::geom_tile(colour = "grey90") +
    do.call(ggplot2::scale_fill_gradientn, fill_args) +
    ggplot2::scale_x_discrete(drop = FALSE) + ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   panel.grid = ggplot2::element_blank())
  if (triangle) hm <- hm + .legend_upper_triangle()
  if (label) {
    mid <- stats::median(df$D, na.rm = TRUE)
    hm <- hm + ggplot2::geom_text(
      ggplot2::aes(label = formatC(.data$D, format = "f", digits = digits),
                   colour = .data$D > mid), size = base_size / 3.4, show.legend = FALSE) +
      ggplot2::scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15"))
  }

  # optional metadata annotation strips (one per annotation, aligned above the columns)
  ann_list <- .resolve_annotations(annotate, ord, pd$group_col, meta)
  ann_panels <- NULL
  if (!is.null(ann_list)) {
    ann_panels <- lapply(names(ann_list), function(nm)
      .annotation_panel(nm, ann_list[[nm]], ord, annotate_colours[[nm]], base_size))
  }

  dnd_panel <- NULL
  if (dendrogram && !is.null(hc)) {
    seg <- .dendro_segments(hc)
    # colour the leaf tips by group so the tree matches the UMAP / admixture colours
    tips <- data.frame(x = seq_along(ord), y = 0, grp = factor(ord, levels = ord),
                       stringsAsFactors = FALSE)
    gcols <- if (is.null(group_colours)) meta_colors(data.frame(grp = ord))[["grp"]]
             else group_colours
    dnd_panel <- ggplot2::ggplot(seg) +
      ggplot2::geom_segment(ggplot2::aes(.data$x, .data$y, xend = .data$xend, yend = .data$yend),
                            linewidth = 0.3, colour = "grey35") +
      ggplot2::geom_point(data = tips, ggplot2::aes(.data$x, .data$y, colour = .data$grp),
                          size = 2.6, show.legend = FALSE) +
      ggplot2::scale_colour_manual(values = gcols) +
      ggplot2::scale_x_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.04, 0.05))) +
      ggplot2::theme_void()
  }

  if (is.null(ann_panels) && is.null(dnd_panel)) return(hm)
  .need_package("patchwork", "the dendrogram / annotation in plot_diff_heatmap()")
  panels <- c(if (!is.null(dnd_panel)) list(dnd_panel), ann_panels, list(hm))
  heights <- c(if (!is.null(dnd_panel)) 0.16, rep(0.06, length(ann_panels)), 1)
  # collect the annotation-strip legends into one tidy column (the fill legend stays
  # inside the triangle); patchwork leaves in-panel "inside" legends where they are.
  patchwork::wrap_plots(panels, ncol = 1, heights = heights) +
    patchwork::plot_layout(guides = "collect")
}

#' @rdname plot_diff_heatmap
#' @export
plot_jost_d_heatmap <- function(pd, stat = "mean", ...) plot_diff_heatmap(pd, stat = stat, ...)

#' The SNPs that most differentiate groups
#'
#' Picks markers from a [pop_diff()] / [jost_d()] result. `"roundrobin"` (default) walks
#' the pairwise comparisons taking each one's next-highest SNP in turn until `n` unique
#' SNPs are collected (so every pair contributes its top differentiators); `"max"` ranks
#' SNPs by their largest value across all pairs.
#'
#' @param jd A [pop_diff()] / [jost_d()] result.
#' @param n How many SNPs to return.
#' @param method `"roundrobin"` or `"max"`.
#' @return A character vector of SNP ids (a subset of the result's SNPs).
#' @export
top_differentiating_snps <- function(jd, n, method = c("roundrobin", "max")) {
  method <- match.arg(method)
  n <- min(n, length(jd$snp))
  if (method == "max") {
    score <- apply(jd$D, 1, max, na.rm = TRUE)
    return(jd$snp[order(score, decreasing = TRUE)][seq_len(n)])
  }
  orders <- lapply(seq_len(ncol(jd$D)), function(k) jd$snp[order(jd$D[, k], decreasing = TRUE)])
  chosen <- character(0)
  seen <- logical(0)
  r <- 1L
  while (length(chosen) < n && r <= length(jd$snp)) {
    for (o in orders) {
      s <- o[r]
      if (!is.na(s) && is.na(seen[s])) {
        chosen <- c(chosen, s); seen[s] <- TRUE
        if (length(chosen) >= n) break
      }
    }
    r <- r + 1L
  }
  chosen
}
