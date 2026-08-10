# Genotypes over one region, sample by SNP, with the samples clustered.

.GENO_LEVELS <- c("reference", "mixed", "alternate")
.HAP_ROWS_IN <- 9        # inches the rows get when they are not individually labelled
.HAP_ANN_IN <- 0.22      # inches per annotation strip
# Every call gets a colour of its own, from the package's colour-blind-friendly set, and only
# missing data is grey: a near-white fill for one of the calls is hard to tell from both the
# panel and the grey of a missing call, which is the one distinction that must stay obvious.
.GENO_FILL <- c(reference = "#2271B2", mixed = "#359B73", alternate = "#D55E00")

# What a dosage means depends on which allele it counts, and the two codings are
# indistinguishable from the matrix alone -- 2 is homozygous alternate under alt dosage and
# homozygous reference under ref dosage. Getting it backwards silently mislabels the whole
# plot, so the object is asked rather than assumed.
.geno_calls <- function(v, allele) {
  idx <- if (identical(allele, "ref")) 3L - v else v + 1L
  factor(.GENO_LEVELS[idx], levels = .GENO_LEVELS)
}

.resolve_allele <- function(x, allele) {
  if (!is.null(allele)) return(match.arg(allele, c("alt", "ref")))
  recorded <- if (is.function(x$allele)) x$allele() else NULL
  if (!is.null(recorded)) return(recorded)
  message("this object does not record which allele its dosages count; assuming `alt` -- ",
          "pass `allele = \"ref\"` if it counts reference alleles instead")
  "alt"
}

# Samples in clustered order within each split block. Clustering per block rather than over
# everything is the point of splitting: the blocks are fixed by the annotation and the
# ordering inside each one is still learned from the genotypes, which is what
# ComplexHeatmap's row_split does.
.cluster_within <- function(G, blocks, cluster = TRUE) {
  ord <- lapply(levels(blocks), function(b) {
    ids <- rownames(G)[blocks == b]
    if (!cluster || length(ids) < 3) return(list(ids = ids, hc = NULL))
    sub <- G[ids, , drop = FALSE]
    # a block whose samples are identical (or all-missing) has nothing to cluster on
    d <- try(stats::dist(sub), silent = TRUE)
    if (inherits(d, "try-error") || !any(is.finite(d)) || max(d, na.rm = TRUE) == 0)
      return(list(ids = ids, hc = NULL))
    d[!is.finite(d)] <- max(d, na.rm = TRUE)
    hc <- stats::hclust(d, method = "ward.D2")
    list(ids = ids[hc$order], hc = hc)
  })
  names(ord) <- levels(blocks)
  ord
}

# Dendrogram segments for every block, in the heatmap's row coordinates: leaves sit on the
# row they label, so the two panels line up facet by facet.
.block_dendro <- function(ord, rows) {
  segs <- lapply(names(ord), function(b) {
    hc <- ord[[b]]$hc
    if (is.null(hc)) return(NULL)
    s <- .dendro_segments(hc)
    # .dendro_segments() puts leaf k of the *clustered order* at x = k; the heatmap numbers
    # its rows the same way inside the block, so the leaf index IS the row index
    base <- min(rows$.row[rows$.split == b]) - 1
    data.frame(y = s$x + base, yend = s$xend + base, x = s$y, xend = s$yend,
               .split = b, stringsAsFactors = FALSE)
  })
  segs <- do.call(rbind, segs[!vapply(segs, is.null, logical(1))])
  if (is.null(segs) || !nrow(segs)) return(NULL)
  segs$.split <- factor(segs$.split, levels = levels(rows$.split))
  segs
}

# Which SNPs to mark: `chr:pos` ids, bare positions, or a gene name resolved in `genes`.
.resolve_marks <- function(marks, loci, genes, reference) {
  if (is.null(marks) || !length(marks)) return(numeric(0))
  if (is.numeric(marks)) return(marks)
  out <- lapply(marks, function(m) {
    hit <- which(loci$snp_id == m | loci$id_norm == m)
    if (length(hit)) return(loci$pos[hit])
    iv <- try(.resolve_region(m, genes, reference), silent = TRUE)
    if (inherits(iv, "try-error"))
      stop("`mark_snps = \"", m, "\"` is not a SNP in the window, a position, or a gene ",
           "name in `genes`", call. = FALSE)
    hit <- loci$pos[loci$pos >= iv$start & loci$pos < iv$end]
    # a gene can resolve fine and still hold no genotyped SNP, which would otherwise just
    # draw nothing and look like the argument was ignored
    if (!length(hit)) message("no genotyped SNP inside ", m, ", so nothing is marked there")
    hit
  })
  unlist(out, use.names = FALSE)
}

#' Genotypes over one region, clustered by sample
#'
#' A genotype heatmap for a single interval: one row per sample, one column per SNP, with the
#' samples clustered and a dendrogram beside them, and a gene track underneath. Splitting the
#' rows by a metadata column fixes the blocks and clusters *within* each one, so a haplotype
#' shared across a group shows up as a solid band rather than being scattered by a
#' genome-wide ordering.
#'
#' The genotypes should be the **full** panel, not an LD-pruned one. Pruning keeps one SNP out
#' of each correlated run, and a correlated run is what a shared haplotype is -- so a pruned
#' panel shows fewer SNPs, chosen to be as uncorrelated as possible, and the blocks come out
#' thinner than they are. Build the object with `load_genotypes(prune = FALSE)` for this plot
#' and keep the pruned one for PCA / UMAP / admixture, where pruning is what you want; the
#' plot says so when the object records that it was pruned.
#'
#' `spacing` decides what the horizontal axis means, and the two answers show different
#' things. `"even"` gives every SNP the same width, which is how the haplotype structure is
#' easiest to read but says nothing about distance. `"genomic"` puts each SNP at its real
#' coordinate, so a dense cluster of SNPs looks dense -- correct about position, but sparse
#' stretches become wide empty bands.
#'
#' @param x A [PopStructure] object (its genotypes supply the calls, its metadata the split).
#' @param region The interval to draw: a gene name from `genes`, a range
#'   (`"13:1,720,000-1,730,000"`), a whole chromosome, or a one-row data frame with
#'   chr/start/end.
#' @param split Optional metadata column whose levels block the rows. Samples are clustered
#'   inside each block, and the blocks keep the column's level order.
#' @param genotypes Optional alternative calls to draw: a genotype matrix (samples x SNPs
#'   with `chr:pos` column names), a [load_genotypes()] list, or another [PopStructure]. `NULL`
#'   (default) uses `x`'s own matrix. Metadata, grouping and the active sample set always come
#'   from `x`, so one pruned object can supply the annotations while the full panel supplies
#'   the calls -- which is the combination this plot wants.
#' @param annotations Optional metadata columns to draw as coloured strips down the right,
#'   one column each, sharing the object's colour maps (see [meta_colors()]) so a level keeps
#'   the colour it has in the other plots. Each gets its own legend.
#' @param border Outline every call (default `TRUE`), which is what makes single SNPs
#'   readable as cells rather than a wash of colour.
#' @param border_colour Colour of that outline.
#' @param allele Which allele the dosages count, `"alt"` or `"ref"`. `NULL` (default) asks
#'   the object; if it does not record it, alt is assumed and a message says so. The two are
#'   indistinguishable from the matrix, and getting it backwards mislabels every call.
#' @param samples Optional sample ids to keep.
#' @param spacing `"even"` (default) gives every SNP equal width; `"genomic"` places each at
#'   its real coordinate.
#' @param cluster Cluster the samples (default `TRUE`). `FALSE` keeps them in the order they
#'   arrive, which is worth doing when the metadata order is the point.
#' @param dendrogram Draw the dendrogram beside the rows (needs `cluster`).
#' @param dend_width Width of the dendrogram panel, as a fraction of the heatmap's.
#' @param mark_snps Optional SNPs to mark with a vertical line: `chr:pos` ids, bare
#'   positions, or a gene name from `genes` (marking every SNP in it).
#' @param mark_colour Colour for those lines.
#' @param genes Gene table for the track and for resolving names (e.g. [PF3D7_GENES]).
#' @param gene_track Draw the gene track under the heatmap (default `TRUE` when `genes` is
#'   given). Under `"even"` spacing the boxes are mapped onto the SNP columns they cover, so
#'   the track still says which columns sit in which gene.
#' @param gene_label_angle Rotation for the gene names, in degrees; `45` or `90` keeps long
#'   systematic ids from colliding.
#' @param pad,min_span Context around `region`, as in [plot_ibd_locus()]: one value pads both
#'   sides, two the left and the right (named `left` / `right` if you like).
#' @param max_snps Refuse to draw more than this many columns (default 2000). A window with
#'   thousands of SNPs is unreadable as tiles; narrow it rather than have it silently thinned.
#' @param snp_width Width of each mark under `"genomic"` spacing, in base pairs. `NULL`
#'   (default) uses 0.5% of the window, wide enough to see and narrow enough to leave the
#'   gaps between SNPs visible.
#' @param colours Named fill colours for `reference` / `mixed` / `alternate`.
#' @param na_colour Fill for missing calls.
#' @param show_sample_names Label the rows. `NULL` (default) labels them when there are at
#'   most 40 samples.
#' @param reference Reference id, used when `region` names a whole chromosome.
#' @return A patchwork of the dendrogram, heatmap and gene track, or a plain ggplot when
#'   neither of those two is drawn.
#' @examples
#' ps <- example_pop_structure(umap = FALSE)
#' plot_region_haplotypes(ps, "7", split = "country", genes = PF_EXAMPLE_DRUG_GENES)
#' @export
plot_region_haplotypes <- function(x, region, split = NULL, annotations = NULL,
                                   genotypes = NULL, samples = NULL,
                                   spacing = c("even", "genomic"),
                                   cluster = TRUE, dendrogram = TRUE, dend_width = 0.15,
                                   border = TRUE, border_colour = "grey45",
                                   allele = NULL,
                                   mark_snps = NULL, mark_colour = "#B2182B",
                                   genes = NULL, gene_track = NULL, gene_label_angle = 0,
                                   pad = 0, min_span = 0, max_snps = 2000, snp_width = NULL,
                                   colours = NULL, na_colour = "grey85",
                                   show_sample_names = NULL,
                                   reference = DEFAULT_REFERENCE) {
  .need_package("ggplot2", "plot_region_haplotypes()")
  spacing <- match.arg(spacing)
  if (is.null(gene_track)) gene_track <- !is.null(genes)

  gt <- .haplotype_genotypes(x, genotypes)
  G <- gt$G
  if (!is.null(samples)) {
    keep <- rownames(G) %in% samples
    if (!any(keep)) stop("none of `samples` are in the genotypes", call. = FALSE)
    G <- G[keep, , drop = FALSE]
  }

  # ---- the window and the SNPs in it ---------------------------------------
  loci <- .parse_snp_ids(colnames(G))
  loci$snp_id <- colnames(G)
  loci$chr <- normalise_chr(loci$chr)
  loci$id_norm <- paste0(loci$chr, ":", loci$pos)
  iv <- .pad_region(.resolve_region(region, genes, reference), pad,
                    .chrom_layout(reference), min_span)
  sel <- loci[loci$chr == iv$chr & loci$pos >= iv$start & loci$pos <= iv$end, , drop = FALSE]
  if (!nrow(sel))
    stop("no genotyped SNPs in ", iv$chr, ":", format(iv$start, scientific = FALSE), "-",
         format(iv$end, scientific = FALSE), "; widen the window with `pad`", call. = FALSE)
  if (nrow(sel) > max_snps)
    stop(nrow(sel), " SNPs in this window is more than `max_snps` (", max_snps, "); narrow ",
         "the window, or raise the limit if you really want that many columns",
         call. = FALSE)
  sel <- sel[order(sel$pos), , drop = FALSE]
  G <- G[, sel$snp_id, drop = FALSE]

  # ---- rows: blocked by the metadata, clustered inside each block ----------
  blocks <- .row_blocks(x, split, rownames(G))
  if (!is.null(blocks$dropped))
    message("dropped ", blocks$dropped, " sample(s) with no ", split)
  G <- G[names(blocks$f), , drop = FALSE]
  ord <- .cluster_within(G, blocks$f, cluster)
  row_ids <- unlist(lapply(ord, `[[`, "ids"), use.names = FALSE)
  rows <- data.frame(
    sample = row_ids,
    .split = factor(as.character(blocks$f[row_ids]), levels = levels(blocks$f)),
    stringsAsFactors = FALSE)
  # rows are numbered top to bottom within the whole plot; each facet then shows its own slice
  rows$.row <- seq_len(nrow(rows))

  # ---- long form, with the x axis the chosen spacing means -----------------
  long <- data.frame(
    sample = rep(row_ids, times = nrow(sel)),
    snp_id = rep(sel$snp_id, each = length(row_ids)),
    pos = rep(sel$pos, each = length(row_ids)),
    col = rep(seq_len(nrow(sel)), each = length(row_ids)),
    value = as.vector(G[row_ids, , drop = FALSE]),
    stringsAsFactors = FALSE)
  long$.row <- rows$.row[match(long$sample, rows$sample)]
  long$.split <- rows$.split[match(long$sample, rows$sample)]
  long$call <- .geno_calls(long$value, allele %||% gt$allele %||% .resolve_allele(x, NULL))
  # LD pruning keeps one SNP out of each correlated run, which is exactly what a shared
  # haplotype is made of, so a pruned panel understates the very structure this plot is for.
  if (isTRUE(gt$pruned))
    message("these genotypes are LD-pruned, so this window shows only the SNPs that survived ",
            "pruning and the haplotype blocks will look thinner than they are; rebuild with ",
            "`load_genotypes(..., prune = FALSE)` for a haplotype view")

  tiles <- .snp_tile_x(sel, spacing, snp_width)
  long$xmin <- tiles$xmin[long$col]
  long$xmax <- tiles$xmax[long$col]
  # under genomic spacing the marks no longer reach the window edges, but the axis should
  xlim <- if (spacing == "genomic") c(min(iv$start, min(tiles$xmin)),
                                      max(iv$end, max(tiles$xmax)))
          else c(min(tiles$xmin), max(tiles$xmax))

  # ---- the heatmap ---------------------------------------------------------
  fills <- .GENO_FILL
  if (!is.null(colours)) fills[names(colours)] <- unname(colours)
  labels <- if (is.null(show_sample_names)) nrow(rows) <= 40 else isTRUE(show_sample_names)

  p <- ggplot2::ggplot(long) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                   ymin = .data$.row - 0.5, ymax = .data$.row + 0.5, fill = .data$call),
      colour = if (border) border_colour else NA,
      linewidth = if (border) 0.06 else 0) +
    ggplot2::scale_fill_manual(values = fills, na.value = na_colour, drop = FALSE,
                               name = "call") +
    ggplot2::scale_y_reverse(
      breaks = if (labels) rows$.row, labels = if (labels) rows$sample,
      expand = ggplot2::expansion(0)) +
    .haplotype_x_scale(sel, spacing, iv, xlim) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   panel.spacing.y = grid::unit(2, "pt"),
                   axis.ticks.y = if (labels) ggplot2::element_line() else
                     ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(),
                   strip.background = ggplot2::element_rect(fill = "grey95", colour = NA),
                   strip.text.y.right = ggplot2::element_text(angle = 0))

  marks <- .resolve_marks(mark_snps, sel, genes, reference)
  if (length(marks)) {
    mx <- .marks_to_x(marks, sel, tiles, spacing)
    if (length(mx))
      p <- p + ggplot2::geom_vline(xintercept = mx, colour = mark_colour, linewidth = 0.35)
  }
  faceted <- nlevels(rows$.split) > 1
  if (faceted)
    p <- p + ggplot2::facet_grid(rows = ggplot2::vars(.data$.split), scales = "free_y",
                                 space = "free_y", switch = NULL)

  # ---- dendrogram beside the rows, gene track underneath -------------------
  dend <- if (dendrogram && cluster) .dendro_panel(.block_dendro(ord, rows), rows,
                                                   faceted) else NULL
  ann <- .hap_annotation_panel(x, annotations, rows, faceted, border, border_colour)
  # only the rightmost panel names the blocks, or every block is labelled twice
  if (!is.null(ann) && faceted)
    p <- p + ggplot2::theme(strip.text.y = ggplot2::element_blank(),
                            strip.text.y.right = ggplot2::element_blank(),
                            strip.background = ggplot2::element_blank())
  track <- if (gene_track) {
    boxes <- .gene_boxes_in_x(genes, iv, sel, tiles, spacing)
    .gene_track_panel(boxes, xlim, angle = gene_label_angle, width_in = .ZOOM_WIDTH_IN)
  } else NULL

  out <- .assemble_haplotype_panels(p, dend, ann, track, dend_width)
  # A labelled row has to be tall enough for its text; an unlabelled one only has to be
  # visible, so hundreds of samples compress into a readable page instead of a 30-inch one.
  row_in <- if (labels) 0.16 else max(0.015, min(0.09, .HAP_ROWS_IN / nrow(rows)))
  attr(out, "plasgenomics_dims") <- c(
    width = .ZOOM_WIDTH_IN + (if (labels) 1 else 0) +
      (if (is.null(ann)) 0 else attr(ann, "ann_in") + 0.6),
    height = max(3, row_in * nrow(rows) + 0.3 * max(1L, nlevels(rows$.split)) + 1) +
      if (is.null(track)) 0 else attr(track, "track_in"))
  out
}

# Coloured strips down the right, one column per annotation. Each needs its own fill scale --
# a level of "region" and a level of "year" are unrelated and must not share a palette -- which
# is what ggnewscale is for: without it a single panel can only carry one fill mapping.
# Colours come from the object's own maps, so a level keeps the colour it has in the UMAP or
# the admixture bars.
.hap_annotation_panel <- function(x, cols, rows, faceted, border, border_colour) {
  if (is.null(cols) || !length(cols)) return(NULL)
  .need_package("ggnewscale", "annotation strips in plot_region_haplotypes()")
  meta <- x$get_meta()
  if (is.null(meta)) stop("`annotations` needs metadata on the object", call. = FALSE)
  missing <- setdiff(cols, names(meta))
  if (length(missing))
    stop("not a metadata column: ", paste(missing, collapse = ", "), call. = FALSE)

  key <- if ("sample" %in% names(meta)) meta$sample else rownames(meta)
  maps <- x$get_colors()
  if (is.null(maps)) maps <- list()
  need <- setdiff(cols, names(maps))
  if (length(need)) maps[need] <- meta_colors(meta, cols = need)[need]

  p <- ggplot2::ggplot() +
    ggplot2::scale_x_continuous(breaks = seq_along(cols), labels = cols,
                                expand = ggplot2::expansion(0), position = "top") +
    ggplot2::scale_y_reverse(expand = ggplot2::expansion(0))
  for (k in seq_along(cols)) {
    cc <- cols[k]
    d <- data.frame(
      .row = rows$.row, .split = rows$.split, x = k,
      value = .as_group_factor(meta[[cc]][match(rows$sample, key)]))
    p <- p +
      ggplot2::geom_rect(
        data = d,
        ggplot2::aes(xmin = .data$x - 0.5, xmax = .data$x + 0.5,
                     ymin = .data$.row - 0.5, ymax = .data$.row + 0.5, fill = .data$value),
        colour = if (border) border_colour else NA,
        linewidth = if (border) 0.06 else 0) +
      ggplot2::scale_fill_manual(values = maps[[cc]], na.value = "grey85", drop = FALSE,
                                 name = cc) +
      ggnewscale::new_scale_fill()
  }
  p <- p +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(), panel.spacing.y = grid::unit(2, "pt"),
      axis.title = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text.x.top = ggplot2::element_text(angle = 90, hjust = 0, vjust = 0.5, size = 7),
      plot.margin = ggplot2::margin(l = 2, r = 2, t = 0, b = 0),
      strip.background = ggplot2::element_rect(fill = "grey95", colour = NA),
      strip.text.y.right = ggplot2::element_text(angle = 0))
  if (faceted)
    p <- p + ggplot2::facet_grid(rows = ggplot2::vars(.data$.split), scales = "free_y",
                                 space = "free_y")
  attr(p, "ann_in") <- .HAP_ANN_IN * length(cols)
  p
}

# The calls to draw. By default the object's own matrix, but a haplotype view usually wants
# the unpruned panel while PCA / admixture want the pruned one, and holding both is cheaper
# than rebuilding: `genotypes` takes a matrix, a load_genotypes() list, or another PopStructure,
# and metadata still comes from `x`. Whatever it is, it also carries whichever of the allele
# coding and the pruning flag it knows, since those are what the plot cannot infer.
.haplotype_genotypes <- function(x, genotypes) {
  src <- list(allele = NULL, pruned = NULL)
  if (is.null(genotypes)) {
    # the full panel when the object has one: pruning drops the correlated SNPs that make a
    # shared haplotype a solid band, which is the whole point of this plot
    panel <- if ("full" %in% x$panels()) "full" else NULL
    G <- x$genotype(prefer = "full")
    src$allele <- x$allele(panel)
    src$pruned <- x$pruned(panel)
  } else if (inherits(genotypes, "PopStructure")) {
    G <- genotypes$genotype()
    src$allele <- genotypes$allele(); src$pruned <- genotypes$pruned()
  } else if (is.list(genotypes) && !is.null(genotypes$genotype)) {
    G <- genotypes$genotype
    if (is.null(rownames(G)) && !is.null(genotypes$sample.id))
      rownames(G) <- genotypes$sample.id
    src$allele <- genotypes$allele; src$pruned <- genotypes$pruned
  } else {
    G <- as.matrix(genotypes)
  }
  if (is.null(G) || !nrow(G) || !ncol(G))
    stop("no genotypes to draw", call. = FALSE)
  if (is.null(rownames(G)) || is.null(colnames(G)))
    stop("`genotypes` needs sample row names and `chr:pos` column names", call. = FALSE)

  # keep the object's active samples, so restrict() / subset() still decide who is drawn
  if (!is.null(genotypes)) {
    want <- x$get_samples()
    keep <- intersect(want, rownames(G))
    if (!length(keep))
      stop("`genotypes` has none of this object's samples", call. = FALSE)
    if (length(keep) < length(want))
      message("`genotypes` is missing ", length(want) - length(keep), " of the object's ",
              length(want), " samples; drawing the ", length(keep), " it has")
    G <- G[keep, , drop = FALSE]
  }
  list(G = G, allele = src$allele, pruned = src$pruned)
}

# The metadata column that blocks the rows, as a factor over the samples being drawn.
.row_blocks <- function(x, split, ids) {
  if (is.null(split)) {
    f <- factor(rep("all", length(ids)), levels = "all")
    names(f) <- ids
    return(list(f = f, dropped = NULL))
  }
  meta <- x$get_meta()
  if (is.null(meta) || !split %in% names(meta))
    stop("`split = \"", split, "\"` is not a metadata column", call. = FALSE)
  key <- if ("sample" %in% names(meta)) meta$sample else rownames(meta)
  v <- meta[[split]][match(ids, key)]
  f <- .as_group_factor(v)
  names(f) <- ids
  keep <- !is.na(f)
  list(f = droplevels(f[keep]), dropped = if (all(keep)) NULL else sum(!keep))
}

# Tile edges per SNP. Even spacing gives each column the same width. So does genomic spacing,
# but centred on the SNP's real coordinate: equal marks at true positions are what let you see
# how far apart the SNPs are. Stretching each tile to meet its neighbours would fill the gaps
# back in and hide exactly that -- and would make an isolated SNP a wide block purely because
# nothing was called near it.
.snp_tile_x <- function(sel, spacing, snp_width = NULL) {
  n <- nrow(sel)
  if (spacing == "even" || n == 1)
    return(list(xmin = seq_len(n) - 0.5, xmax = seq_len(n) + 0.5))
  pos <- sel$pos
  span <- max(diff(range(pos)), 1)
  w <- if (!is.null(snp_width)) snp_width else max(span * 0.005, 1)
  list(xmin = pos - w / 2, xmax = pos + w / 2)
}

# Under even spacing the axis counts SNPs, so it is labelled with the coordinates of a few of
# them rather than pretending the numbers are positions.
.haplotype_x_scale <- function(sel, spacing, iv, xlim) {
  if (spacing == "genomic") {
    unit <- if (diff(xlim) >= 1e4) 1e3 else 1
    return(ggplot2::scale_x_continuous(
      name = paste0("chromosome ", iv$chr, " position (", if (unit == 1) "bp" else "kb", ")"),
      labels = function(v) format(round(v / unit, if (unit == 1) 0 else 1), big.mark = ",",
                                  trim = TRUE, scientific = FALSE),
      expand = ggplot2::expansion(0)))
  }
  at <- unique(round(seq(1, nrow(sel), length.out = min(nrow(sel), 6))))
  ggplot2::scale_x_continuous(
    name = paste0("chromosome ", iv$chr, ": ", nrow(sel), " SNPs, evenly spaced"),
    breaks = at,
    labels = format(round(sel$pos[at] / 1000, 1), big.mark = ",", trim = TRUE),
    expand = ggplot2::expansion(0))
}

# Positions -> the plot's x units (a column index under even spacing).
.marks_to_x <- function(marks, sel, tiles, spacing) {
  if (spacing == "genomic") return(marks)
  i <- match(marks, sel$pos)
  # a marked position with no genotyped SNP has no column to sit on
  (seq_len(nrow(sel)))[i[!is.na(i)]]
}

# Gene boxes in the plot's x units. Under even spacing a gene spans the columns it covers,
# interpolated so a gene between two SNPs still lands between their columns.
.gene_boxes_in_x <- function(genes, iv, sel, tiles, spacing) {
  if (is.null(genes)) return(NULL)
  g <- .gene_track(genes)
  g <- g[g$chr == iv$chr & as.numeric(g$end) >= iv$start & as.numeric(g$start) <= iv$end, ,
         drop = FALSE]
  if (!nrow(g)) return(NULL)
  if (spacing == "genomic") {
    g$.gene_xmin <- as.numeric(g$start)
    g$.gene_xmax <- as.numeric(g$end)
    return(g)
  }
  to_col <- function(v) {
    if (nrow(sel) == 1) return(rep(1, length(v)))
    stats::approx(sel$pos, seq_len(nrow(sel)), xout = v, rule = 2)$y
  }
  g$.gene_xmin <- to_col(as.numeric(g$start))
  g$.gene_xmax <- to_col(as.numeric(g$end))
  g
}

# The dendrogram panel: rows on the y axis so it lines up with the heatmap, distance growing
# leftwards away from it.
.dendro_panel <- function(segs, rows, faceted) {
  if (is.null(segs)) return(NULL)
  # Anchor each facet to the same y range the heatmap uses (its rows, plus the half row the
  # tiles extend past the first and last). Without them the dendrogram facet takes the range
  # of whatever segments it holds, and the leaves drift off the rows they label; a shared
  # global range is just as wrong, since `space = "free_y"` then gives every block the same
  # height regardless of how many samples it has.
  lim <- lapply(split(rows$.row, rows$.split), function(r) c(min(r) - 0.5, max(r) + 0.5))
  anchor <- data.frame(
    y = unlist(lim, use.names = FALSE), x = 0,
    .split = factor(rep(names(lim), each = 2), levels = levels(rows$.split)))

  p <- ggplot2::ggplot(segs) +
    ggplot2::geom_segment(ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y,
                                       yend = .data$yend), linewidth = 0.25,
                          colour = "grey25") +
    ggplot2::geom_blank(data = anchor, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::scale_x_reverse(expand = ggplot2::expansion(mult = c(0.04, 0))) +
    ggplot2::scale_y_reverse(expand = ggplot2::expansion(0)) +
    ggplot2::theme_void() +
    # theme_void() leaves the facet strips, and the blocks are already named once on the far
    # right -- naming them here as well repeats every label and pushes the dendrogram away
    # from the genotypes it is describing
    ggplot2::theme(plot.margin = ggplot2::margin(l = 2, r = 0, t = 0, b = 0),
                   panel.spacing.y = grid::unit(2, "pt"),
                   strip.text = ggplot2::element_blank(),
                   strip.background = ggplot2::element_blank())
  if (faceted)
    p <- p + ggplot2::facet_grid(rows = ggplot2::vars(.data$.split), scales = "free_y",
                                 space = "free_y")
  p
}

# dendrogram | heatmap | annotations, with the gene track under the heatmap. The track and the
# annotation strips keep the absolute sizes they measured for themselves, so their labels
# cannot be clipped however tall or wide the page ends up; the heatmap takes what is left.
#
# patchwork's `design` letters are positional -- A is the FIRST plot passed, B the second --
# so the letters have to be generated in the order the plots go in, not chosen to be mnemonic.
.assemble_haplotype_panels <- function(p, dend, ann, track, dend_width) {
  if (is.null(dend) && is.null(ann) && is.null(track)) return(p)
  .need_package("patchwork", "the panels around plot_region_haplotypes()")

  plots <- list(); widths <- list(); role <- character(0)
  add <- function(pl, w, what) {
    plots[[length(plots) + 1L]] <<- pl
    widths[[length(widths) + 1L]] <<- w
    role <<- c(role, what)
  }
  if (!is.null(dend)) add(dend, grid::unit(dend_width, "null"), "dend")
  add(p, grid::unit(1, "null"), "heat")                       # the one elastic column
  if (!is.null(ann)) add(ann, grid::unit(attr(ann, "ann_in"), "in"), "ann")

  design <- paste(LETTERS[seq_along(plots)], collapse = "")
  heights <- grid::unit(1, "null")
  if (!is.null(track)) {
    # the track spans only the heatmap column: it shares that x axis and nothing else does
    tl <- LETTERS[length(plots) + 1L]
    plots[[length(plots) + 1L]] <- track
    design <- paste(design, paste(ifelse(role == "heat", tl, "#"), collapse = ""), sep = "\n")
    heights <- grid::unit.c(heights, grid::unit(attr(track, "track_in"), "in"))
  }
  out <- patchwork::wrap_plots(plots, design = design) +
    patchwork::plot_layout(widths = do.call(grid::unit.c, widths), heights = heights,
                           guides = "collect")
  attr(out, "plasgenomics_track_in") <- if (is.null(track)) 0 else attr(track, "track_in")
  out
}
