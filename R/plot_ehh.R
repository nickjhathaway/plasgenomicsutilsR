# Extended haplotype homozygosity around one SNP, allele by allele.

.EHH_LEVELS <- c("reference", "alternate")
.EHH_FILL <- c(reference = "#2271B2", alternate = "#D55E00")

# The SNP the decay is measured from. A `chr:pos` id or a bare position names one outright; a
# gene usually holds several, and then the one with the most balanced alleles is the only
# defensible automatic choice -- EHH from a singleton is a line at 1 that says nothing. Which
# one was used is reported, since for a drug-resistance locus the caller often means a
# particular mutation and should name it.
.focal_marker <- function(focal, map, genes, reference) {
  if (is.numeric(focal) && length(focal) == 1) {
    hit <- which(map$pos == focal)
    if (!length(hit)) stop("no SNP at position ", format(focal, scientific = FALSE),
                           call. = FALSE)
    return(hit[1])
  }
  if (!is.character(focal) || length(focal) != 1)
    stop("`focal` must be a `chr:pos` id, a position, or a gene name", call. = FALSE)

  hit <- which(map$snp_id == focal | paste0(normalise_chr(map$chr), ":", map$pos) == focal)
  if (length(hit)) return(hit[1])

  iv <- try(.resolve_region(focal, genes, reference), silent = TRUE)
  if (inherits(iv, "try-error"))
    stop("`focal = \"", focal, "\"` is not a SNP in the haplotypes, a position, or a gene ",
         "name in `genes`", call. = FALSE)
  inside <- which(normalise_chr(map$chr) == iv$chr & map$pos >= iv$start & map$pos < iv$end)
  if (!length(inside))
    stop("no genotyped SNP inside ", focal, call. = FALSE)
  inside
}

# EHH either side of one marker, for one set of haplotypes, as a long table. Returns a string
# instead when there is nothing to draw, so the caller can say which reason it was.
# `chr` is the map's own chromosome name: .haplohh_list() keys its objects by that, not by the
# normalised form used for intervals and gene tracks.
.ehh_curve <- function(hap, rows, mrk_pos, chr, polarized, limehh) {
  objs <- .haplohh_list(hap, rows)
  o <- objs[[as.character(chr)]]
  if (is.null(o)) return("too few polymorphic SNPs on that chromosome")
  mrk <- which(o@positions == mrk_pos)
  # .haplohh_list() drops columns that are monomorphic within the group, focal SNP included
  if (!length(mrk)) return("the focal SNP is not variable in it")
  e <- try(rehh::calc_ehh(o, mrk = mrk[1], polarized = polarized, limehh = limehh,
                          include_zero_values = TRUE, phased = TRUE), silent = TRUE)
  if (inherits(e, "try-error") || is.null(e$ehh) || !nrow(e$ehh))
    return("rehh returned no EHH values there")

  d <- as.data.frame(e$ehh)
  cols <- intersect(c("EHH_MAJ", "EHH_MIN", "EHH_A", "EHH_D", "EHH"), names(d))
  if (!length(cols)) return("rehh returned no EHH columns")
  # unpolarized: MAJ / MIN are the two alleles as coded, 0 then 1, so reference then alternate
  lab <- if (length(cols) == 1) .EHH_LEVELS[2] else .EHH_LEVELS[seq_along(cols)]
  out <- do.call(rbind, lapply(seq_along(cols), function(k) data.frame(
    pos = d$POSITION, ehh = d[[cols[k]]],
    allele = factor(lab[k], levels = .EHH_LEVELS), stringsAsFactors = FALSE)))
  freq <- e$freq
  attr(out, "freq") <- stats::setNames(as.numeric(freq), lab[seq_along(freq)])
  attr(out, "n") <- length(rows)
  out
}

# Candidate focal SNPs and the balance metric the focal is chosen on. Shared so that
# ehh_candidates() ranks by what plot_ehh() actually acts on -- a table that merely resembled
# the choice would be no use for overriding it.
.focal_candidates <- function(x, focal, genes, reference) {
  idx <- .focal_marker(focal, x$map, genes, reference)
  maf <- vapply(idx, function(i) {
    p <- mean(x$hap[, i], na.rm = TRUE)
    min(p, 1 - p)
  }, numeric(1))
  list(idx = idx, maf = maf, best = idx[which.max(maf)])
}

#' The SNPs an EHH plot had to choose between
#'
#' [plot_ehh()] measures decay from a single focal SNP, and a gene usually holds many. It takes
#' the one with the most balanced alleles and reports which -- this returns the whole shortlist
#' it chose from, chosen SNP first, so the decision is visible and can be overridden by passing
#' a `chr:pos` back to `focal`.
#'
#' `maf` is the metric that decides it: minor-allele frequency across the haplotypes, which
#' [plot_ehh()] maximises because EHH measured from a near-singleton is a flat line at 1 that
#' says nothing about a sweep.
#'
#' With `group`, each group gets a `maf_<group>` column and `n_groups_variable` counts the groups
#' the SNP actually varies in. That is the other way a panel comes back empty: a SNP can be well
#' balanced overall and monomorphic inside one group, which drops that group's curve with "the
#' focal SNP is not variable in it". Sorting on `n_groups_variable` finds a SNP that works
#' everywhere, when one exists.
#'
#' @param x A [parasite_haplotypes()] object, or a [PopStructure].
#' @param focal As in [plot_ehh()]: a `chr:pos` id, a bare position, or a gene name in `genes`.
#' @param group Optional metadata column, for the per-group columns.
#' @param genes Gene table used to resolve `focal` (e.g. [PF3D7_GENES]).
#' @param min_haplotypes Groups smaller than this are left out, as in [plot_ehh()].
#' @param reference Reference id, used when `focal` names a whole chromosome.
#' @return A tibble of `snp_id`, `chr`, `pos`, `maf`, `n_hap` and `chosen`, plus the per-group
#'   columns when `group` is given. The chosen SNP is the first row; the rest follow by
#'   descending `maf`.
#' @seealso [plot_ehh()]
#' @examples
#' ps <- example_pop_structure(umap = FALSE)
#' hap <- parasite_haplotypes(ps, maf = 0.05)
#' ehh_candidates(hap, "pfcrt", genes = PF_EXAMPLE_DRUG_GENES)
#' @export
ehh_candidates <- function(x, focal, group = NULL, genes = NULL, min_haplotypes = 10,
                           reference = DEFAULT_REFERENCE) {
  if (inherits(x, "PopStructure")) x <- parasite_haplotypes(x)
  if (is.null(x$hap) || is.null(x$map))
    stop("`x` must be a parasite_haplotypes() object", call. = FALSE)

  cand <- .focal_candidates(x, focal, genes, reference)
  idx <- cand$idx
  out <- tibble::tibble(
    snp_id = x$map$snp_id[idx],
    chr = as.character(x$map$chr[idx]),
    pos = x$map$pos[idx],
    maf = round(cand$maf, 4),
    n_hap = vapply(idx, function(i) sum(!is.na(x$hap[, i])), integer(1)),
    chosen = idx == cand$best)

  if (!is.null(group)) {
    rows <- .ihs_rows(x, group, x$meta, min_haplotypes)
    for (g in names(rows)) {
      out[[paste0("maf_", g)]] <- round(vapply(idx, function(i) {
        p <- mean(x$hap[rows[[g]], i], na.rm = TRUE)
        if (is.nan(p)) NA_real_ else min(p, 1 - p)
      }, numeric(1)), 4)
    }
    gcols <- paste0("maf_", names(rows))
    m <- as.matrix(out[, gcols, drop = FALSE])
    out$n_groups_variable <- as.integer(rowSums(!is.na(m) & m > 0))
  }
  out[order(!out$chosen, -out$maf), , drop = FALSE]
}

#' EHH decay around one SNP
#'
#' Extended haplotype homozygosity either side of a focal SNP, one curve per allele: how far
#' the haplotype carrying each allele stays identical as you walk away from it. A sweep shows
#' as one allele holding EHH near 1 far past the point where the other has decayed -- the
#' picture behind a single point on an [run_ihs()] scan.
#'
#' The alleles are the two states at the focal SNP itself, so this is the mutant-versus-
#' reference comparison without needing the SNPs annotated: `reference` is the allele coded 0
#' and `alternate` the one coded 1. `group` adds a panel per metadata group; without it every
#' haplotype is pooled, which is usually what you want first, since EHH knows nothing about
#' population structure and a group with few carriers gives a ragged curve.
#'
#' @param x A [parasite_haplotypes()] object, or a [PopStructure] (haplotypes are then built
#'   with `parasite_haplotypes()` defaults, which is worth doing yourself when the Fws or MAF
#'   cutoffs matter).
#' @param focal The SNP to measure from: a `chr:pos` id, a bare position, or a gene name from
#'   `genes`. A gene holding several SNPs resolves to the one with the most balanced alleles,
#'   reported in a message -- name a `chr:pos` to pick a particular mutation.
#' @param group Optional metadata column; one panel per level. `NULL` (default) pools every
#'   haplotype.
#' @param span How far either side of the focal SNP to draw, in base pairs (default 50 kb).
#'   One value is symmetric, two are the left and the right (named `left` / `right` if you
#'   like), as elsewhere in the package.
#' @param min_haplotypes Skip a group with fewer haplotypes than this (default 10). EHH from a
#'   handful of haplotypes is mostly noise.
#' @param polarized Treat the alleles as ancestral / derived (default `FALSE`, matching
#'   [run_ihs()] on unpolarized calls, where they are simply the two states).
#' @param limehh Stop each curve once EHH falls below this (rehh's `limehh`, default 0.05).
#' @param genes Gene table for the track and for resolving `focal` (e.g. [PF3D7_GENES]).
#' @param gene_track Draw the gene track underneath (default `FALSE`). `genes` is usually
#'   supplied only to resolve `focal`, and an EHH window is wide enough that a full annotation
#'   would crowd a hundred names under it, so this is opt-in.
#' @param gene_label_angle Rotation for the gene names, in degrees.
#' @param colours,colors Named colours for `reference` / `alternate`.
#' @param show_freq Note each panel's haplotype count and allele frequencies inside it
#'   (default `TRUE`); `FALSE` leaves the panel clean.
#' @param freq_position Which corner that note sits in: `"topleft"` (default), `"topright"`,
#'   `"bottomleft"` or `"bottomright"`. The top corners are usually clear, since EHH is 1 at
#'   the focal SNP and both curves have flattened along the bottom by the window's edges.
#' @param reference Reference id, used when `focal` names a whole chromosome.
#' @return A patchwork of the curves over the gene track, or a plain ggplot without one.
#' @examples
#' ps <- example_pop_structure(umap = FALSE)
#' hap <- parasite_haplotypes(ps, maf = 0.05)
#' plot_ehh(hap, "pfcrt", genes = PF_EXAMPLE_DRUG_GENES, span = 30000)
#' @export
plot_ehh <- function(x, focal, group = NULL, span = 50000, min_haplotypes = 10,
                     polarized = FALSE, limehh = 0.05, genes = NULL, gene_track = NULL,
                     gene_label_angle = 0, colours = NULL, show_freq = TRUE,
                     freq_position = c("topleft", "topright", "bottomleft", "bottomright"),
                     reference = DEFAULT_REFERENCE,
                     colors = NULL) {
  colours <- .alias_arg("colours", "colors")
  .need_package("ggplot2", "plot_ehh()")
  .need_package("rehh", "plot_ehh()")
  freq_position <- match.arg(freq_position)
  if (inherits(x, "PopStructure")) {
    message("building haplotypes with parasite_haplotypes() defaults; pass a ",
            "parasite_haplotypes() object to control the Fws / MAF filtering")
    x <- parasite_haplotypes(x)
  }
  if (is.null(x$hap) || is.null(x$map))
    stop("`x` must be a parasite_haplotypes() object", call. = FALSE)
  # `genes` is usually passed just to resolve `focal`, and an EHH window is wide enough that
  # a full annotation would put a hundred names under it, so the track is opt-in here.
  gene_track <- isTRUE(gene_track)

  map <- x$map
  # most balanced alleles: EHH from a near-singleton is a flat line that says nothing
  cc <- .focal_candidates(x, focal, genes, reference)
  cand <- cc$best
  if (length(cc$idx) > 1) {
    message("`", focal, "` holds ", length(cc$idx), " SNPs; measuring from ",
            map$snp_id[cand], " (minor allele ", format(round(max(cc$maf), 3)), ") -- name a ",
            "`chr:pos` to pick another, or see ehh_candidates() for the shortlist")
  }
  mrk_pos <- map$pos[cand]
  mrk_chr <- as.character(map$chr[cand])          # as the haplotypes name it
  mrk_chr_norm <- normalise_chr(mrk_chr)          # as intervals and gene tracks name it

  rows <- if (is.null(group)) list(all = seq_len(nrow(x$hap)))
          else .ihs_rows(x, group, x$meta, min_haplotypes)
  curves <- lapply(names(rows), function(g) {
    if (length(rows[[g]]) < min_haplotypes) {
      message("skipping ", g, ": ", length(rows[[g]]), " haplotype(s), fewer than ",
              min_haplotypes)
      return(NULL)
    }
    d <- .ehh_curve(x, rows[[g]], mrk_pos, mrk_chr, polarized, limehh)
    if (is.character(d)) {
      message("skipping ", g, ": ", d)
      return(NULL)
    }
    d$group <- g
    attr(d, "freq") <- attr(d, "freq")
    d
  })
  names(curves) <- names(rows)
  keep <- !vapply(curves, is.null, logical(1))
  if (!any(keep)) stop("no group has a usable EHH curve at this SNP", call. = FALSE)
  freqs <- lapply(curves[keep], attr, "freq")
  df <- do.call(rbind, curves[keep])

  iv <- .pad_region(list(chr = mrk_chr_norm, start = mrk_pos, end = mrk_pos), span,
                    .chrom_layout(reference))
  df <- df[df$pos >= iv$start & df$pos <= iv$end, , drop = FALSE]
  if (!nrow(df)) stop("no EHH values within `span` of the focal SNP", call. = FALSE)
  faceted <- !is.null(group) && length(unique(df$group)) > 1
  if (faceted) df$group <- factor(df$group, levels = names(rows)[keep])

  fills <- .EHH_FILL
  if (!is.null(colours)) fills[names(colours)] <- unname(colours)
  xlim <- c(iv$start, iv$end)

  p <- ggplot2::ggplot(df, ggplot2::aes(.data$pos, .data$ehh, colour = .data$allele)) +
    ggplot2::geom_vline(xintercept = mrk_pos, colour = "grey55", linewidth = 0.3) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::geom_point(size = 0.5) +
    ggplot2::scale_colour_manual(values = fills, drop = FALSE, name = "allele at focal SNP") +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = ggplot2::expansion(mult = 0.02)) +
    .region_axis(iv, 0) +
    ggplot2::labs(y = "EHH", title = paste0("EHH around ", map$snp_id[cand])) +
    ggplot2::coord_cartesian(xlim = xlim) +
    .manhattan_theme() +
    ggplot2::theme(legend.position = "right")
  if (faceted)
    p <- p + ggplot2::facet_wrap(~ .data$group, ncol = 1, strip.position = "right")

  if (isTRUE(show_freq)) {
    lab <- vapply(names(freqs), function(g) paste0(
      "n = ", attr(curves[[g]], "n"), "; ",
      paste(sprintf("%s %.0f%%", names(freqs[[g]]), 100 * freqs[[g]]), collapse = ", ")),
      character(1))
    ann <- data.frame(group = names(freqs), label = unname(lab), stringsAsFactors = FALSE)
    if (faceted) ann$group <- factor(ann$group, levels = levels(df$group))
    # A corner the curves are least likely to occupy: EHH is 1 at the focal SNP and decays
    # outwards, so the far edges are low and the top corners stay clear. The default is the top
    # left because the bottom is where both curves flatten out once they have decayed.
    at <- switch(freq_position,
                 topleft     = list(x = xlim[1], y = 0.99, h = 0, v = 1),
                 topright    = list(x = xlim[2], y = 0.99, h = 1, v = 1),
                 bottomleft  = list(x = xlim[1], y = 0.01, h = 0, v = 0),
                 bottomright = list(x = xlim[2], y = 0.01, h = 1, v = 0))
    p <- p + ggplot2::geom_text(
      data = ann, ggplot2::aes(x = at$x, y = at$y, label = .data$label),
      inherit.aes = FALSE, hjust = at$h, vjust = at$v, size = 2.6, colour = "grey30")
  }

  n_panels <- if (faceted) nlevels(df$group) else 1L
  if (gene_track) {
    z <- list(interval = iv, offset = 0, xlim = xlim,
              track = .genes_in_span(genes, iv), label = TRUE)
    if (!is.null(z$track))
      p <- .stack_gene_track(p, z, data.frame(chr = mrk_chr_norm, offset = 0), n_panels,
                             gene_label_angle)
  }
  # a curve needs more height than a Manhattan row to be readable
  attr(p, "plasgenomics_dims") <- c(
    width = .ZOOM_WIDTH_IN,
    height = max(3.4, 1.4 + 2 * max(1L, n_panels)) +
      (attr(p, "plasgenomics_track_in") %||% 0))
  p
}

# Genes overlapping an interval, in the track's own columns.
.genes_in_span <- function(genes, iv) {
  if (is.null(genes)) return(NULL)
  g <- .gene_track(genes)
  g <- g[g$chr == iv$chr & as.numeric(g$end) >= iv$start & as.numeric(g$start) <= iv$end, ,
         drop = FALSE]
  if (!nrow(g)) NULL else g
}
