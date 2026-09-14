# The Schaffner control for a carrier-stratified locus result: are the carriers just more
# related to each other across the *rest* of the genome to begin with? If they are, a locus
# statistic normalised only at the population level counts that background as locus signal.
#
# This is the background that ibd_block_extension_by_allele()'s locus contrast has to survive.
# It is not that function -- that one measures the signal *at* the locus; this measures the
# confound everywhere *except* the locus's chromosome. It is also not pair_fraction_summary():
# that summarises the whole-genome per-pair fraction and is not carrier-split, and it reads a
# single genome-wide scalar per pair, which structurally cannot drop one chromosome. Excluding
# the focal chromosome is the whole point here (so the sweep cannot seed its own baseline), so
# the background is summed from the blocks instead.

# One row per carrier class within a group, from a frame of per-pair (class, mb).
.carrier_background_group <- function(d, off_focal_bp, min_pairs) {
  ref <- d$mb[d$class == "reference/reference"]
  ref_mean <- if (length(ref) >= min_pairs) mean(ref) else NA_real_
  parts <- split(d$mb, d$class)
  out <- lapply(names(parts), function(cl) {
    mb <- parts[[cl]]
    data.frame(
      class = cl,
      n_pairs = length(mb),
      mean_ibd_mb = mean(mb),
      median_ibd_mb = stats::median(mb),
      # a share of the callable genome off the focal chromosome, the scale the rest of the
      # document reads IBD on
      mean_ibd_frac = mean(mb) * 1e6 / off_focal_bp,
      # anchored on a reference stratum that itself clears min_pairs, so the fold is never
      # divided by a number too thin to report
      fold_vs_ref = mean(mb) / ref_mean,
      # rank-sum against that same group's reference/reference stratum
      p_vs_ref = if (identical(cl, "reference/reference") ||
                     length(mb) < min_pairs || length(ref) < min_pairs) NA_real_
                 else .rank_sum_p(mb, ref),
      stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

#' Are a variant's carriers more related genome-wide to begin with?
#'
#' The genome-wide background control for [ibd_block_extension_by_allele()]. Schaffner et
#' al.'s caution about reading IBD inside a sweep is that isolates carrying the swept
#' haplotype are typically more related to each other across the *rest* of the genome as
#' well, so a locus statistic normalised only at the population level counts that background
#' as locus signal. This measures exactly that background: within each group it splits the
#' pairs by their state at a core variant and reports each stratum's mean per-pair IBD off
#' the **focal chromosome**, which is dropped in full so the sweep cannot contribute to its
#' own baseline.
#'
#' Where `fold_vs_ref` is above 1 and `p_vs_ref` is small, the carriers are a more related
#' set genome-wide and a population-level normalisation at the locus would be reading part of
#' that as locus signal; [ibd_block_extension_test()]'s `paired_ratio` and everything
#' [ibd_block_extension_by_allele()] reports already divide each pair by its own baseline and
#' are immune to it, so the point of the table is to say how much that per-pair control is
#' removing, and where. Where the fold sits at 1, carriage says nothing about background
#' relatedness and the two normalisations agree.
#'
#' @section What counts as carrier and reference:
#' Unlike [ibd_block_extension_by_allele()], which contrasts exactly two named states,
#' `carrier` and `reference` are **sets** here, matched with `%in%`: at *pfkelch13* every
#' propeller allele can count as carrier against a single `"REF"` reference, and a mixed call
#' naming both a carrier and a reference allele lands in neither. A pair is `carrier/carrier`
#' when both ends are in `carrier`, `reference/reference` when both are in `reference`, and
#' `discordant` when one is a carrier and the other a reference. A pair with an end in neither
#' set (a third allele, or `NA`) is left out of every stratum.
#'
#' @param x An [IbdResults] built with `blocks =` and `meta =`.
#' @param allele The variant to split on: the name of a metadata column, or a named vector of
#'   `sample -> state` (as [allele_states()] returns). Samples with `NA` are left out.
#' @param focal_chr The chromosome to exclude from the background, as the locus sits on it.
#'   Any spelling `normalise_chr()` accepts (`"Pf3D7_07_v3"` or `"7"`).
#' @param carrier,reference Which states count as carrying and as reference; each may be a
#'   vector. With exactly two observed states either may be left out and the other is implied
#'   (`reference` the first level), with a message. With three or more states both must be
#'   named, since there is nothing to imply the second from.
#' @param group Metadata column defining the groups. Defaults to the object's declared group
#'   column, then to the first non-`sample` column of `meta`.
#' @param min_pairs A stratum needs this many pairs before it is reported, and the
#'   reference stratum needs it before a fold or a test is anchored on it (default `5`).
#' @param meta Sample metadata; taken from `x` when not given.
#' @return A tibble, one row per group x class, ordered by group then class:
#'   \describe{
#'     \item{`group`, `class`}{the group and one of `carrier/carrier`,
#'       `reference/reference`, `discordant`.}
#'     \item{`n_pairs`}{every within-group pair of the class, including the pairs that share
#'       nothing -- taking the mean over only the sharing pairs would compare two different
#'       denominators and hide the effect being tested.}
#'     \item{`mean_ibd_mb`, `median_ibd_mb`}{megabases of IBD per pair, summed over that
#'       pair's segments off the focal chromosome.}
#'     \item{`mean_ibd_frac`}{the mean as a share of the callable genome outside the focal
#'       chromosome (its core length), the scale the rest of an IBD document uses.}
#'     \item{`fold_vs_ref`}{`mean_ibd_mb` over the group's `reference/reference` mean.}
#'     \item{`p_vs_ref`}{two-sided rank-sum of the stratum against that same group's
#'       `reference/reference` pairs; `NA` on the reference row itself and where either side
#'       is under `min_pairs`.}
#'   }
#' @seealso [ibd_block_extension_by_allele()] for the locus contrast this backstops,
#'   [pair_fraction_summary()] for the whole-genome per-pair sharing by group,
#'   [ibd_block_extension_test()] for the per-pair normalisation that absorbs this background.
#' @examples
#' \dontrun{
#' ibd <- ibd_results(blocks = "hmm.txt", meta = meta, group_col_in_meta = "region")
#' carrier_genome_wide_relatedness(ibd, allele = "PIN_Haplotype", focal_chr = "Pf3D7_07_v3",
#'                                 carrier = "Present", reference = "Absent")
#' }
#' @export
carrier_genome_wide_relatedness <- function(x, allele, focal_chr,
                                            carrier = NULL, reference = NULL,
                                            group = NULL, min_pairs = 5L, meta = NULL) {
  if (!inherits(x, "IbdResults"))
    stop("`x` must be an IbdResults; build one with ibd_results(blocks = , meta = )",
         call. = FALSE)
  if (missing(allele) || is.null(allele)) stop("`allele` is required", call. = FALSE)
  if (missing(focal_chr) || is.null(focal_chr)) stop("`focal_chr` is required", call. = FALSE)
  blocks <- x$get_blocks()
  if (is.null(blocks) || !nrow(blocks))
    stop("this IbdResults has no IBD blocks; build it with ibd_results(blocks = , meta = )",
         call. = FALSE)
  if (is.null(meta)) meta <- x$get_meta()
  meta <- .normalise_meta(meta)
  if (is.null(meta) || !"sample" %in% names(meta))
    stop("splitting the pairs needs meta with a 'sample' column", call. = FALSE)
  if (is.null(group)) group <- x$get_group_col()
  if (is.null(group)) group <- setdiff(names(meta), "sample")[1]
  if (is.na(group) || !group %in% names(meta))
    stop("meta has no column '", group, "'.\n  columns available: ",
         paste(setdiff(names(meta), "sample"), collapse = ", "), call. = FALSE)

  fc <- normalise_chr(focal_chr)

  # per-pair total IBD bp off the focal chromosome, keyed by the unordered pair
  bl <- blocks
  bl$chr <- normalise_chr(bl$chr)
  bl <- bl[bl$chr != fc, , drop = FALSE]
  if (!nrow(bl))
    stop("every IBD segment is on the focal chromosome ", fc, "; no background to measure",
         call. = FALSE)
  key <- paste(pmin(bl$sample1, bl$sample2), pmax(bl$sample1, bl$sample2), sep = "\r")
  ibd_bp <- tapply(as.numeric(bl$end) - as.numeric(bl$start), key, sum)

  # per-sample state (reusing the allele reader) and group
  m <- meta[!duplicated(meta$sample), , drop = FALSE]
  ss <- .allele_states(allele, m, m$sample)
  st <- ss$state
  lv <- ss$levels
  if (is.null(carrier) && is.null(reference)) {
    if (length(lv) < 2)
      stop("`allele` has fewer than two observed states, so there is nothing to contrast",
           call. = FALSE)
    if (length(lv) > 2)
      stop("`allele` has ", length(lv), " states (", paste(lv, collapse = ", "),
           "); name `carrier =` and `reference =`", call. = FALSE)
    reference <- lv[1]; carrier <- lv[2]
    message("reading `", reference, "` as reference and `", carrier, "` as carrier; ",
            "pass carrier=/reference= to swap them")
  } else if (is.null(reference) || is.null(carrier)) {
    if (length(lv) > 2) {
      named <- if (is.null(reference)) "carrier" else "reference"
      missing <- if (is.null(reference)) "reference" else "carrier"
      stop("`allele` has ", length(lv), " states (", paste(lv, collapse = ", "),
           "); with `", named, " =` given, `", missing,
           " =` cannot be inferred -- name it too.", call. = FALSE)
    }
    if (is.null(reference)) reference <- setdiff(lv, carrier)
    else carrier <- setdiff(lv, reference)
  }
  carrier <- as.character(carrier); reference <- as.character(reference)
  if (length(intersect(carrier, reference)))
    stop("`carrier` and `reference` overlap (", paste(intersect(carrier, reference),
         collapse = ", "), "); a state cannot be both", call. = FALSE)

  gp <- stats::setNames(as.character(m[[group]]), m$sample)
  keep <- intersect(x$get_analyzed_samples(), names(st)[!is.na(st)])
  if (length(keep) < 2)
    stop("fewer than two analysed samples have a known `allele` state", call. = FALSE)

  per_group <- lapply(split(keep, gp[keep]), function(s) {
    if (length(s) < 2) return(NULL)
    cb <- utils::combn(sort(s), 2)
    a <- st[cb[1, ]]; b <- st[cb[2, ]]
    cls <- ifelse(a %in% carrier   & b %in% carrier,   "carrier/carrier",
           ifelse(a %in% reference & b %in% reference, "reference/reference",
           ifelse(a %in% c(carrier, reference) & b %in% c(carrier, reference),
                  "discordant", NA_character_)))
    v <- unname(ibd_bp[paste(cb[1, ], cb[2, ], sep = "\r")]); v[is.na(v)] <- 0
    data.frame(class = cls, mb = v / 1e6, stringsAsFactors = FALSE)[!is.na(cls), , drop = FALSE]
  })
  per_group <- per_group[!vapply(per_group, is.null, logical(1))]
  if (!length(per_group))
    stop("no group has two analysed samples with a known `allele` state", call. = FALSE)

  # callable span off the focal chromosome: the reference's core chromosome lengths, so the
  # denominator is the same accessible span the package's genome-wide coordinates use rather
  # than a per-run guess read off the blocks
  lens <- get_reference(x$reference_id())$core_chrom_lengths_bp
  off_focal_bp <- sum(as.numeric(lens[names(lens) != fc]))

  out <- do.call(rbind, lapply(names(per_group), function(g) {
    d <- .carrier_background_group(per_group[[g]], off_focal_bp, min_pairs)
    if (is.null(d) || !nrow(d)) return(NULL)
    cbind(group = g, d, stringsAsFactors = FALSE)
  }))
  out <- out[out$n_pairs >= min_pairs, , drop = FALSE]
  if (!nrow(out))
    stop("no group x class stratum has at least ", min_pairs, " pairs", call. = FALSE)

  # order groups by the object's declared order where it has one, classes in a fixed order
  gord <- x$get_group_order()
  glev <- if (!is.null(gord)) c(intersect(gord, out$group), setdiff(unique(out$group), gord))
          else unique(out$group)
  clev <- c("carrier/carrier", "discordant", "reference/reference")
  out <- out[order(match(out$group, glev), match(out$class, clev)), , drop = FALSE]
  out$group <- factor(out$group, levels = glev)
  out$class <- factor(out$class, levels = clev[clev %in% out$class])
  rownames(out) <- NULL
  res <- tibble::as_tibble(out)
  attr(res, "states") <- list(carrier = carrier, reference = reference)
  attr(res, "focal_chr") <- fc
  attr(res, "off_focal_bp") <- off_focal_bp
  res
}
