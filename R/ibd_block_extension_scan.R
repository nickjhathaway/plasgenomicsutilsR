# Genome-wide background for the locus block-extension statistic: the same test run over
# fixed-width windows, so a locus can be read against a distribution instead of against 1.
#
# Why windows rather than genes. A `paired_ratio` of 1 is not the null -- a locus is a fixed
# target and a segment covers it with probability rising in the segment's own length, so
# every locus is inflated whether or not anything happened there. The scan measures that
# inflation empirically. Windows beat genes as the background for two reasons: they are all
# the same width, which matters under `sharing = "complete"` where width is a floor under the
# numerator; and they do not depend on annotation, so an intergenic target is visible.

# Fixed-width windows tiling a set of regions, dropping any region too short to hold one.
.tile_regions <- function(regions, width, step) {
  r <- .as_gene_track(regions)
  starts <- lapply(seq_len(nrow(r)), function(i) {
    lo <- as.numeric(r$start[i])
    hi <- as.numeric(r$end[i]) - width
    if (!is.finite(lo) || !is.finite(hi) || hi < lo) return(NULL)   # region shorter than a window
    s <- seq(lo, hi, by = step)
    if (!length(s)) return(NULL)
    data.frame(name = paste0(r$chr[i], ":", format(s, scientific = FALSE, trim = TRUE)),
               chr = r$chr[i], start = s, end = s + width, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, starts[!vapply(starts, is.null, logical(1))])
  if (is.null(out) || !nrow(out))
    stop("no region is at least `width` (", width, " bp) long", call. = FALSE)
  out
}

# Where `value` falls in `ref`, as a percentage of `ref` at or below it.
.percentile_in <- function(value, ref) {
  ref <- ref[is.finite(ref)]
  if (!length(ref)) return(rep(NA_real_, length(value)))
  100 * stats::ecdf(ref)(value)
}

#' Genome-wide scan of IBD block extension, and loci placed against it
#'
#' Runs [ibd_block_extension_test()] over fixed-width windows tiling the genome, then reports
#' where each window -- and each locus you name -- sits in its own group's distribution. This
#' is the interface to reach for first: a single locus's `paired_ratio` has no interpretable
#' scale on its own, because the statistic is inflated at *every* locus by length-biased
#' sampling, and the scan is what measures that inflation for the cohort in hand.
#'
#' Windows rather than genes as the background, for two reasons. They are all the same width,
#' which matters under `sharing = "complete"`, where a segment cannot be shorter than the
#' locus is wide and so wider loci score higher for no biological reason. And they owe
#' nothing to annotation, so a target between genes still shows up.
#'
#' @section What the resolution is, and is not:
#' Do not read a window as localising selection to its own width. The value at a point is set
#' by the segments covering it, and those run tens of kb either side, so neighbouring windows
#' repeat each other almost exactly: on one real cohort the correlation of
#' `log2(paired_ratio)` between windows 250 bp apart was 0.997, still 0.973 at 2 kb, and only
#' fell to 0.65 at the median segment length of 33 kb. Every window inside *pfcrt*, *pfgch1*
#' and *pfkelch13* returned an identical ratio. The resolution floor is the segment-length
#' scale, so a sub-genic target is out of reach here -- use [run_ihs()] and [plot_ehh()],
#' which work at SNP resolution.
#'
#' That is also why `step` defaults to 2 kb rather than something finer: a denser grid costs
#' time and memory and returns the same numbers.
#'
#' @section Why there is no q-value:
#' Windows overlap the same segments, so they are not independent tests. The core genome over
#' the median segment length is on the order of a few hundred independent positions against
#' tens of thousands of windows, which makes a Benjamini-Hochberg q across windows both
#' invalid and wildly conservative. `percentile` is the quantity to quote; it does not care
#' how finely the genome was cut. Per-window p-values are still reported, unadjusted, so you
#' can see which windows are tested at all.
#'
#' @param x An [IbdResults] built with `blocks =` and `meta =`.
#' @param loci Optional loci to place against the scan: gene names from the object's track,
#'   or an interval data frame as [ibd_block_extension_test()] takes. Each comes back with
#'   its `percentile` in its group's window distribution. `NULL` (default) returns the
#'   windows alone.
#' @param regions Regions to tile (`chr`/`chrom`, `start`, `end`). Defaults to
#'   [PF3D7_CORE_REGIONS] for a `pf3d7` object, which keeps the hypervariable subtelomeres
#'   out of the background; pass your own for any other reference.
#' @param width Window width in bp (default `500`).
#' @param step Distance between window starts in bp (default `2000`). See the resolution
#'   note: finer buys nothing.
#' @param min_windows Groups with fewer than this many tested windows get `NA` percentiles
#'   and a warning, since a distribution over a handful of windows is not one (default `50`).
#' @param group,within,sharing,min_ref_blocks,min_pairs,meta Passed to
#'   [ibd_block_extension_test()].
#' @return A tibble with [ibd_block_extension_test()]'s columns plus:
#'   \describe{
#'     \item{`source`}{`"window"` for the background rows, `"locus"` for anything named in
#'       `loci`.}
#'     \item{`percentile`}{where this row's `paired_ratio` falls among the **windows** of the
#'       same group, as a percentage at or below it. The number to quote.}
#'     \item{`n_windows`}{tested windows behind that group's distribution.}
#'   }
#'   Loci come first, then windows by descending `paired_ratio`. The per-pair ratios are not
#'   retained; call [ibd_block_extension_test()] directly on a locus when you want them.
#' @seealso [ibd_block_extension_test()] for one locus and what the columns mean,
#'   [ibd_block_extension_by_allele()] for splitting the pairs by carriage of a variant.
#' @examples
#' \dontrun{
#' ibd <- ibd_results(blocks = "hmm.txt", meta = meta, group_col_in_meta = "region")
#' scan <- ibd_block_extension_scan(ibd, loci = c("pfcrt", "pfgch1"))
#' subset(scan, source == "locus", c(locus, group, paired_ratio, percentile))
#' }
#' @export
ibd_block_extension_scan <- function(x, loci = NULL, regions = NULL, width = 500,
                                     step = 2000, group = NULL, within = 0,
                                     sharing = c("overlap", "complete"),
                                     min_ref_blocks = 1L, min_pairs = 5L,
                                     min_windows = 50L, meta = NULL) {
  sharing <- match.arg(sharing)
  if (!inherits(x, "IbdResults"))
    stop("`x` must be an IbdResults; build one with ibd_results(blocks = , meta = )",
         call. = FALSE)
  meta <- .normalise_meta(meta)
  if (!is.numeric(width) || width < 1) stop("`width` must be a positive number of bp", call. = FALSE)
  if (!is.numeric(step) || step < 1) stop("`step` must be a positive number of bp", call. = FALSE)
  if (is.null(regions)) {
    ref <- x$reference_id()
    if (!identical(tolower(ref), "pf3d7"))
      stop("no default regions for reference '", ref, "'; pass `regions =` (chr, start, end)",
           call. = FALSE)
    regions <- PF3D7_CORE_REGIONS
  }

  win <- .tile_regions(regions, width, step)
  common <- list(x = x, group = group, within = within, sharing = sharing,
                 min_ref_blocks = min_ref_blocks, min_pairs = min_pairs,
                 keep_pair_ratios = FALSE, meta = meta)
  w <- do.call(ibd_block_extension_test, c(common, list(loci = win, adjust = "none")))
  if (!nrow(w)) stop("no window had enough sharing pairs to test; loosen `min_pairs`",
                     call. = FALSE)
  w$source <- "window"

  res <- if (is.null(loci)) w else {
    l <- do.call(ibd_block_extension_test, c(common, list(loci = loci, adjust = "none")))
    l$source <- "locus"
    # `locus` is a factor of each call's own levels, so bind on the labels
    w$locus <- as.character(w$locus)
    l$locus <- as.character(l$locus)
    dplyr::bind_rows(l, w)
  }

  n_win <- table(w$group)
  res$n_windows <- as.integer(n_win[as.character(res$group)])
  thin <- names(n_win)[n_win < min_windows]
  if (length(thin))
    warning(sprintf("%d group(s) have fewer than %d tested windows (%s); their percentiles ",
                    length(thin), min_windows, paste(thin, collapse = ", ")),
            "are NA -- a distribution over that few windows is not one", call. = FALSE)

  res$percentile <- NA_real_
  for (g in names(n_win)[n_win >= min_windows]) {
    ref_ratios <- w$paired_ratio[as.character(w$group) == g]
    i <- as.character(res$group) == g
    res$percentile[i] <- .percentile_in(res$paired_ratio[i], ref_ratios)
  }

  res <- res[order(factor(res$source, levels = c("locus", "window")),
                   -res$paired_ratio), , drop = FALSE]
  res$locus <- factor(as.character(res$locus), levels = unique(as.character(res$locus)))
  keep <- c("locus", "name", "gene_id", "chr", "start", "end", "span_bp", "source", "group",
            "n_pairs", "gw_median", "locus_median", "naive_ratio", "ref_median",
            "paired_ratio", "percentile", "n_windows", "z_paired", "p_paired", "p_two_sided")
  res <- res[, intersect(keep, names(res)), drop = FALSE]
  attr(res, "sharing") <- sharing
  attr(res, "window") <- c(width = width, step = step)
  res
}
