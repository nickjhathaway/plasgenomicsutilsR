# Extended haplotype homozygosity around one SNP, allele by allele.

.EHH_LEVELS <- c("reference", "alternate")
.EHH_FILL <- c(reference = "#2271B2", alternate = "#D55E00")

# rehh names a marker's allele columns EHH_MAJ, EHH_MIN1, EHH_MIN2, ... and those names are
# **positional**, not frequency ranks: EHH_MAJ is allele 0 whether or not it is the commonest
# one (on a 3-allele marker coded 10/30/20, FREQ_MAJ comes back as the rarest, 0.167). So the
# k-th column is allele k-1, and that is the only thing the labels can safely be built from.
# Two alleles keep the reference / alternate wording the plot has always used; more get
# numbered, because "the alternate" stops meaning anything once there are several.
# Name an allele by *which allele it is*, never by which column of a scan it landed in.
# rehh's columns are dense over the alleles present in the haplotypes it was handed, so a
# group that happens to lack allele 2 gets two columns and would otherwise call its allele 1
# "alternate" while the group beside it calls the same allele "alternate 1" -- one allele,
# two names, in one faceted plot.
.ehh_allele_label <- function(a, n_alleles) {
  alt <- if (n_alleles <= 2) rep("alternate", length(a)) else paste("alternate", a)
  ifelse(a == 0L, "reference", alt)
}

.ehh_allele_levels <- function(k) {
  if (k <= 2) return(.EHH_LEVELS[seq_len(max(k, 1))])
  c("reference", paste("alternate", seq_len(k - 1)))
}

.ehh_allele_fill <- function(levels) {
  if (all(levels %in% names(.EHH_FILL))) return(.EHH_FILL[levels])
  stats::setNames(.pick_palette(length(levels)), levels)
}

# The focal SNP's row in an iHS scan, one value per panel of the plot.
#
# iHS is standardised against the whole genome, so it cannot be recomputed from the window
# this plot draws -- the number has to come from a genome-wide scan, either one handed in or
# one run here. Matching is by chromosome and position, and then by panel name: a scan
# carrying `group` is read per group, one without a `group` column applies to every panel.
# A panel with no row in the scan gets `NA` rather than a neighbour's score.
.ehh_ihs_at <- function(scan, chr, pos, panels, contrast = NULL) {
  df <- as.data.frame(scan)
  for (nm in c("chr", "pos", "ihs"))
    if (!nm %in% names(df))
      stop("`add_ihs` needs a run_ihs() result; this table has no '", nm, "' column",
           call. = FALSE)
  # one focal marker can be asked more than one question -- reference against each of a
  # multiallelic site's alternates, say -- so a `contrast` column splits the table into
  # named sets and each is looked up on its own.
  if (!is.null(contrast)) df <- df[as.character(df$contrast) == contrast, , drop = FALSE]
  # the scan may name chromosomes as the haplotypes do or in the short form, so compare
  # both through the same normalisation rather than trusting either spelling
  hit <- df[normalise_chr(as.character(df$chr)) == normalise_chr(chr) &
              as.numeric(df$pos) == as.numeric(pos), , drop = FALSE]
  out <- stats::setNames(rep(NA_real_, length(panels)), panels)
  freq <- out
  if (!nrow(hit)) return(list(ihs = out, freq_minor = freq, grouped = "group" %in% names(df)))
  if ("group" %in% names(df)) {
    g <- as.character(hit$group)
    m <- match(panels, g)
    out[] <- hit$ihs[m]
    if ("freq_minor" %in% names(hit)) freq[] <- hit$freq_minor[m]
  } else {
    # an ungrouped scan describes one population, so it labels every panel it is shown with
    out[] <- hit$ihs[1]
    if ("freq_minor" %in% names(hit)) freq[] <- hit$freq_minor[1]
  }
  list(ihs = out, freq_minor = freq, grouped = "group" %in% names(df))
}

# How one panel's iHS is written. Unpolarized, the sign is not the derived allele's -- it is
# just which of major/minor carried the longer haplotype -- so only the magnitude is shown,
# the same restraint run_ihs() documents.
#
# `abs(iHS)` rather than the conventional `|iHS|`: at this text size the vertical bars are
# the same stroke as a lowercase l in every sans face the device is likely to pick, so
# "|iHS|" renders as "liHSl" and is genuinely misread -- checked at 600 dpi, and spacing the
# bars out does not fix it. The wording follows run_ihs()'s own docs, which say to read
# `abs(ihs)`.
.ehh_ihs_label <- function(v, polarized) {
  if (!is.finite(v)) return(NA_character_)
  if (isTRUE(polarized)) sprintf("iHS %+.2f", v) else sprintf("abs(iHS) %.2f", abs(v))
}

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

  # The alleles at this marker across the whole object, so a label means the same allele in
  # every facet, and the alleles this group actually carries, which may be a subset.
  gcol <- which(hap$map$chr == chr & hap$map$pos == mrk_pos)[1]
  all_alleles <- sort(unique(as.integer(hap$hap[, gcol])))
  alleles <- o@haplo[, mrk[1]]
  present <- sort(unique(as.integer(alleles[!is.na(alleles)])))

  # rehh numbers its columns densely over the alleles it is shown, and truncates rather than
  # reporting a leading absent one (given alleles 1 and 2 it returns FREQ_MAJ = 0 and drops
  # allele 2 outright). Hand it a dense 0..k-1 coding, which is a bijection on the focal
  # column and so leaves the haplotype partition untouched, and keep the map back.
  if (!identical(present, seq_along(present) - 1L))
    o@haplo[, mrk[1]] <- match(alleles, present) - 1L

  e <- try(rehh::calc_ehh(o, mrk = mrk[1], polarized = polarized, limehh = limehh,
                          include_zero_values = TRUE, phased = TRUE), silent = TRUE)
  if (inherits(e, "try-error") || is.null(e$ehh) || !nrow(e$ehh))
    return("rehh returned no EHH values there")

  d <- as.data.frame(e$ehh)
  # Take whatever allele columns rehh returned rather than a fixed list of names: a marker
  # with three alleles comes back as EHH_MAJ / EHH_MIN1 / EHH_MIN2, and a fixed list that
  # knows only EHH_MIN would keep one curve of three and label it as the other one.
  cols <- grep("^EHH(_|$)", names(d), value = TRUE)
  if (!length(cols)) return("rehh returned no EHH columns")
  freq <- e$freq
  if (length(cols) != length(present) || (length(freq) && length(freq) != length(cols)))
    return(sprintf("rehh returned %d curves and %d frequencies for the %d allele(s) here",
                   length(cols), length(freq), length(present)))
  lab <- .ehh_allele_label(present, length(all_alleles))
  levs <- .ehh_allele_label(all_alleles, length(all_alleles))
  out <- do.call(rbind, lapply(seq_along(cols), function(k) data.frame(
    pos = d$POSITION, ehh = d[[cols[k]]],
    allele = factor(lab[k], levels = levs), stringsAsFactors = FALSE)))
  attr(out, "freq") <- stats::setNames(as.numeric(freq), lab[seq_along(freq)])
  attr(out, "n") <- length(rows)
  attr(out, "levels") <- levs
  # Count the alleles rather than recovering them from the frequency: rounding a share back
  # into a count is right until it is not, and the haplotypes are right here. Counted against
  # the real allele values, not the dense recoding rehh was handed.
  attr(out, "count") <- stats::setNames(
    vapply(present, function(a) sum(alleles == a, na.rm = TRUE), integer(1)), lab)
  out
}

# Candidate focal SNPs and the balance metric the focal is chosen on. Shared so that
# ehh_candidates() ranks by what plot_ehh() actually acts on -- a table that merely resembled
# the choice would be no use for overriding it.
.focal_candidates <- function(x, focal, genes, reference) {
  idx <- .focal_marker(focal, x$map, genes, reference)
  # `.minor_af()` rather than `min(mean(v), 1 - mean(v))`: the latter is a dosage formula
  # and returns a negative number on an allele-index column, which `which.max()` below can
  # never choose -- so a multiallelic focal was silently passed over for a biallelic
  # neighbour. The two agree exactly on 0/1 data.
  maf <- vapply(idx, function(i) .minor_af(x$hap[, i]), numeric(1))
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
        .minor_af(x$hap[rows[[g]], i])
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
#' @param colours,colors Named colours for the focal alleles. A biallelic marker has
#'   `reference` and `alternate`; one with more alleles has `reference`, `alternate 1`,
#'   `alternate 2`, ... and takes its default colours from the shared palette.
#' @param show_freq Note each panel's haplotype count and allele frequencies inside it
#'   (default `TRUE`); `FALSE` leaves the panel clean.
#' @param freq_position Which corner that note sits in: `"topleft"` (default), `"topright"`,
#'   `"bottomleft"` or `"bottomright"`. The top corners are usually clear, since EHH is 1 at
#'   the focal SNP and both curves have flattened along the bottom by the window's edges.
#' @param add_ihs Add the focal SNP's iHS to that corner note. A [run_ihs()] result is read
#'   for the focal SNP -- the cheap path, and the one to prefer, since it reuses a scan you
#'   already have and so the number in the corner is the same one the genome-wide figures
#'   were drawn from. `TRUE` runs [run_ihs()] here instead, on `x`, with this plot's `group`
#'   and `polarized` and anything in `ihs_args`; that is a whole-genome scan per call, so it
#'   is slow and worth doing once into a variable rather than once per plot. `NULL` (default)
#'   or `FALSE` adds nothing. iHS is standardised against the whole genome and cannot be
#'   recovered from the window drawn here, which is why there is no third option. Read
#'   against the same caution [run_ihs()] carries: unpolarized, only the magnitude is shown,
#'   because the sign is major-versus-minor and not ancestral-versus-derived.
#'
#'   A `contrast` column in the table asks the same focal marker more than one question --
#'   the reference against each alternate of a multiallelic codon, say. Each named set gets
#'   its own line in the corner, prefixed by its name, and a panel with no value for one of
#'   them simply omits that line rather than printing a blank. Without the column the table
#'   is read as a single unnamed contrast, as before.
#' @param ihs_args Extra arguments for the [run_ihs()] call made by `add_ihs = TRUE`, as a
#'   named list -- `maxgap`, `maf_bands`, `min_maf` and the rest. `group` and `polarized` come
#'   from this plot so the two halves cannot disagree, and naming either here is an error.
#'   Ignored, with a warning, when `add_ihs` is a scan you computed yourself.
#' @param reference Reference id, used when `focal` names a whole chromosome.
#' @param title Plot title: `NULL` (default) uses `"EHH around <snp>"`, a string sets a custom
#'   one, and `NA`/`FALSE` draws none. Set it here rather than adding `labs(title = )` to the
#'   result -- with `gene_track = TRUE` the result is a patchwork, and `+ labs()` lands on the
#'   gene track at the bottom, putting the title in the middle of the figure.
#' @param subtitle Line under the title; `NULL` (default) draws none.
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
                     reference = DEFAULT_REFERENCE, title = NULL, subtitle = NULL,
                     add_ihs = NULL, ihs_args = list(), colors = NULL) {
  colours <- .alias_arg("colours", "colors")
  .need_package("ggplot2", "plot_ehh()")
  .need_package("rehh", "plot_ehh()")
  freq_position <- match.arg(freq_position)
  # checked before anything expensive runs: `add_ihs = TRUE` scans the whole genome, and a
  # mistyped `ihs_args` should not be found out on the far side of that
  want_ihs <- !is.null(add_ihs) && !isFALSE(add_ihs)
  if (want_ihs) {
    if (isTRUE(add_ihs)) {
      if (!is.list(ihs_args) || (length(ihs_args) && is.null(names(ihs_args))))
        stop("`ihs_args` must be a named list of arguments for run_ihs()", call. = FALSE)
      # these two describe the same thing as the plot's own arguments; letting them be set
      # twice is how the corner note ends up describing a different scan from the curves
      clash <- intersect(names(ihs_args), c("hap", "group", "polarized"))
      if (length(clash))
        stop("`ihs_args` must not set ", paste(sprintf("`%s`", clash), collapse = " or "),
             ": it is taken from plot_ehh()'s own argument so the scan and the curves ",
             "cannot disagree", call. = FALSE)
    } else {
      if (!is.data.frame(add_ihs))
        stop("`add_ihs` must be TRUE, FALSE, or a run_ihs() result", call. = FALSE)
      if (length(ihs_args))
        warning("`ihs_args` is ignored when `add_ihs` is a scan you computed yourself; ",
                "those settings belong in the run_ihs() call that made it", call. = FALSE)
    }
  }
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

  # the levels the curves actually carry, which is two only when the marker is biallelic
  lv <- attr(curves[keep][[1]], "levels")
  fills <- .ehh_allele_fill(if (is.null(lv)) .EHH_LEVELS else lv)
  if (!is.null(colours)) fills[names(colours)] <- unname(colours)
  xlim <- c(iv$start, iv$end)

  p <- ggplot2::ggplot(df, ggplot2::aes(.data$pos, .data$ehh, colour = .data$allele)) +
    ggplot2::geom_vline(xintercept = mrk_pos, colour = "grey55", linewidth = 0.3) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::geom_point(size = 0.5) +
    ggplot2::scale_colour_manual(values = fills, drop = FALSE, name = "allele at focal SNP") +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = ggplot2::expansion(mult = 0.02)) +
    .region_axis(iv, 0) +
    ggplot2::labs(y = "EHH",
                  title = .ehh_title(title, map$snp_id[cand]),
                  subtitle = if (is.character(subtitle)) subtitle else NULL) +
    ggplot2::coord_cartesian(xlim = xlim) +
    .manhattan_theme() +
    ggplot2::theme(legend.position = "right")
  if (faceted)
    p <- p + ggplot2::facet_wrap(~ .data$group, ncol = 1, strip.position = "right")

  ihs_lab <- NULL
  if (want_ihs) {
    scan <- if (isTRUE(add_ihs))
      do.call(run_ihs, c(list(hap = x, group = group, polarized = polarized), ihs_args))
    else add_ihs
    # a `contrast` column asks the same focal marker more than one question. Each named set
    # becomes its own line in the corner, prefixed by its name, because two bare numbers
    # stacked in a corner say nothing about which is which.
    contrasts <- if ("contrast" %in% names(as.data.frame(scan)))
      unique(as.character(as.data.frame(scan)$contrast)) else NULL
    if (!is.null(contrasts)) {
      per <- lapply(contrasts, function(ct) {
        g <- .ehh_ihs_at(scan, mrk_chr, mrk_pos, names(freqs), contrast = ct)
        blank <- !is.finite(g$ihs)
        if (all(blank))
          message("no iHS for the ", ct, " contrast at ", map$snp_id[cand],
                  " in any panel drawn")
        else if (any(blank))
          message("no iHS for the ", ct, " contrast in ",
                  paste(names(freqs)[blank], collapse = ", "))
        stats::setNames(vapply(names(freqs), function(nm) {
          v <- .ehh_ihs_label(g$ihs[[nm]], polarized)
          if (is.na(v)) NA_character_ else paste(ct, v)
        }, character(1)), names(freqs))
      })
      ihs_lab <- stats::setNames(lapply(names(freqs), function(nm) {
        v <- vapply(per, function(p) p[[nm]], character(1))
        v[!is.na(v)]
      }), names(freqs))
      if (!length(unlist(ihs_lab))) ihs_lab <- NULL
    } else {
    got <- .ehh_ihs_at(scan, mrk_chr, mrk_pos, names(freqs))
    miss <- names(freqs)[!is.finite(got$ihs)]
    if (length(miss) == length(freqs)) {
      if (got$grouped && is.null(group))
        message("the scan is grouped and this plot pools every haplotype, so no panel ",
                "matches it; pass a run_ihs() result without `group`, or facet this plot ",
                "by the same column")
      else
        message("no iHS for ", map$snp_id[cand], " in the scan (it may be below the scan's ",
                "`min_maf`, or dropped at a `maxgap` border), so none is shown")
    } else if (length(miss)) {
      message("no iHS for ", paste(miss, collapse = ", "), " at ", map$snp_id[cand],
              "; those panels are labelled without one")
    }
    # A scan computed on a different set of haplotypes than the curves would put a number
    # in the corner that describes other samples. There is no sample list in a scan to
    # compare, but its minor-allele frequency at this SNP is one the plot also knows.
    if (any(is.finite(got$freq_minor))) {
      own <- vapply(names(freqs), function(g) {
        # only a biallelic marker has one minor allele to compare; and only a panel that
        # will actually show a score is worth warning about
        f <- freqs[[g]]
        if (length(f) == 2 && is.finite(got$ihs[[g]])) min(f) else NA_real_
      }, numeric(1))
      d <- abs(own - got$freq_minor)
      off <- names(freqs)[is.finite(d) & d > 0.02]
      if (length(off))
        message("the scan's minor-allele frequency at ", map$snp_id[cand], " differs from ",
                "this plot's for ", paste(off, collapse = ", "),
                " -- the scan looks to have been run on a different set of haplotypes, so ",
                "its iHS describes those samples, not the curves drawn here")
    }
    ihs_lab <- stats::setNames(lapply(names(freqs), function(g) {
      v <- .ehh_ihs_label(got$ihs[[g]], polarized)
      if (is.na(v)) character(0) else v
    }), names(freqs))
    if (!length(unlist(ihs_lab))) ihs_lab <- NULL
    }
  }

  if (isTRUE(show_freq) || !is.null(ihs_lab)) {
    # the count as well as the share: reading "8%" against "n = 60" to get 5 haplotypes is
    # arithmetic the reader should not have to do, and 8% of 60 reads very differently from
    # 8% of 600. The counts are the ones the haplotypes were counted into, not the share
    # multiplied back out.
    lab <- vapply(names(freqs), function(g) {
      lines <- character(0)
      if (isTRUE(show_freq)) {
        f <- freqs[[g]]
        k <- attr(curves[[g]], "count")
        shares <- if (is.null(k)) sprintf("%s %.0f%%", names(f), 100 * f)
                  else sprintf("%s %d (%.0f%%)", names(f), k[names(f)], 100 * f)
        lines <- paste0("n = ", attr(curves[[g]], "n"), "; ",
                        paste(shares, collapse = ", "))
      }
      # its own line rather than appended to the counts: the counts line is already long
      # on a faceted plot, and the score is the thing a reader is looking for
      if (!is.null(ihs_lab)) lines <- c(lines, ihs_lab[[g]])
      paste(lines, collapse = "\n")
    }, character(1))
    ann <- data.frame(group = names(freqs), label = unname(lab), stringsAsFactors = FALSE)
    ann <- ann[nzchar(ann$label), , drop = FALSE]
    if (faceted) ann$group <- factor(ann$group, levels = levels(df$group))
    # A corner the curves are least likely to occupy: EHH is 1 at the focal SNP and decays
    # outwards, so the far edges are low and the top corners stay clear. The default is the top
    # left because the bottom is where both curves flatten out once they have decayed.
    at <- switch(freq_position,
                 topleft     = list(x = xlim[1], y = 0.99, h = 0, v = 1),
                 topright    = list(x = xlim[2], y = 0.99, h = 1, v = 1),
                 bottomleft  = list(x = xlim[1], y = 0.01, h = 0, v = 0),
                 bottomright = list(x = xlim[2], y = 0.01, h = 1, v = 0))
    if (nrow(ann))
      p <- p + ggplot2::geom_text(
        data = ann, ggplot2::aes(x = at$x, y = at$y, label = .data$label),
        inherit.aes = FALSE, hjust = at$h, vjust = at$v, size = 2.6, colour = "grey30",
        lineheight = 0.95)
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

# `NULL` keeps the written title, a string replaces it, `NA`/`FALSE` drops it -- the same
# three cases the network plots take.
.ehh_title <- function(title, snp_id) {
  if (is.null(title)) return(paste0("EHH around ", snp_id))
  if (isFALSE(title) || (length(title) == 1 && is.na(title))) return(NULL)
  as.character(title)
}
