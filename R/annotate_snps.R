# Attach interval annotations (genes, or any BED-like table) to a per-SNP scan.

#' Which intervals each SNP falls in
#'
#' Takes any per-SNP table -- [run_ihs()], [run_rsb()], [run_xpehh()], [beta_score()],
#' `pop_diff_snps()`, an `IbdResults` selection table, anything with a `snp_id` or a
#' `chr`/`pos` pair -- and reports the intervals covering each SNP. The usual use is "which
#' gene is this hit in", which otherwise means writing the overlap join by hand every time.
#'
#' Positions are 0-based and intervals half-open `[start, end)`, the package convention
#' (`?"plasgenomicsutilsR-coordinates"`), so a SNP at the interval's `end` is *outside* it.
#' Chromosome names are normalised on both sides, so `Pf3D7_07_v3`, `chr7` and `7` match.
#'
#' A SNP in no interval is kept with `NA` (or dropped by `keep = "hits"`), and a SNP in
#' several -- overlapping gene spans, or nested intervals -- yields one row per interval
#' rather than being silently collapsed. Use `collapse = TRUE` for one row per SNP with the
#' names pasted together instead.
#'
#' @param scan A data frame with `snp_id` (`"chr:pos0"`) or with `chr` and `pos` columns.
#' @param intervals A data frame with `name`, `chr` (or `chrom`), `start` and `end` --
#'   [PF3D7_GENES], [PF_EXAMPLE_DRUG_GENES], [PF3D7_CORE_REGIONS] or your own BED.
#' @param within Widen every interval by this many bp on both sides before testing, for
#'   catching a hit just outside a short gene (default `0`).
#' @param keep `"all"` (default) keeps SNPs matching nothing, with `NA` annotations;
#'   `"hits"` keeps only SNPs inside an interval.
#' @param collapse Return one row per SNP, with multiple hits pasted into `name` and counted
#'   in `n_intervals`, instead of one row per SNP x interval (default `FALSE`).
#' @param prefix Prepend this to the added column names, to keep two annotations side by
#'   side (e.g. `prefix = "core_"`).
#' @param one_based_snps The positions in `scan` are 1-based (VCF `POS`), so shift them before
#'   testing. Genotype-matrix column names from [load_genotypes()] are 1-based, and joining
#'   those against these 0-based intervals without saying so moves every SNP one base and can
#'   turn a near-miss into a hit. Default `FALSE`, the package convention.
#' @return `scan` with `name`, `interval_start`, `interval_end`, `distance_to_midpoint` and
#'   any `gene_id` carried through -- plus `n_intervals` when `collapse = TRUE`. The
#'   original columns and their order are preserved.
#' @seealso [bed_intersect()] for interval-to-interval overlap; [selection_peaks()] to merge
#'   neighbouring significant SNPs into loci before annotating.
#' @examples
#' ihs <- data.frame(snp_id = c("Pf3D7_07_v3:403500", "Pf3D7_07_v3:1"),
#'                   ihs = c(4.2, 0.1))
#' annotate_snps(ihs, PF_EXAMPLE_DRUG_GENES)
#' annotate_snps(ihs, PF_EXAMPLE_DRUG_GENES, keep = "hits")
#' @export
annotate_snps <- function(scan, intervals, within = 0, keep = c("all", "hits"),
                          collapse = FALSE, prefix = "", one_based_snps = FALSE) {
  keep <- match.arg(keep)
  df <- as.data.frame(scan, stringsAsFactors = FALSE)
  if (!nrow(df)) return(scan)

  # position: explicit chr/pos if present, else split the snp_id
  if (all(c("chr", "pos") %in% names(df))) {
    chr <- normalise_chr(df$chr)
    pos <- as.numeric(df$pos)
  } else if ("snp_id" %in% names(df)) {
    id <- as.character(df$snp_id)
    chr <- normalise_chr(sub(":[^:]*$", "", id))
    pos <- suppressWarnings(as.numeric(sub("^.*:", "", id)))
    if (anyNA(pos))
      stop("could not read a position out of `snp_id`; expected \"chr:pos\"", call. = FALSE)
  } else {
    stop("`scan` needs a `snp_id` column, or `chr` and `pos` columns", call. = FALSE)
  }
  # A genotype matrix from load_genotypes() names its columns with 1-based VCF positions
  # (SNPRelate's convention), while intervals here are 0-based half-open. Joining the two
  # without saying so shifts every SNP one base and quietly turns near-misses into hits.
  if (isTRUE(one_based_snps)) pos <- pos - 1

  iv <- .gene_track(intervals)
  gid <- if ("gene_id" %in% names(intervals)) as.character(intervals$gene_id) else NULL
  # .gene_track drops extra columns, so re-attach gene_id by the row order it preserved
  by_chr <- split(seq_len(nrow(iv)), iv$chr)

  hit_snp <- integer(0); hit_iv <- integer(0)
  for (i in seq_len(nrow(df))) {
    cand <- by_chr[[chr[i]]]
    if (is.null(cand) || is.na(pos[i])) next
    # half-open, widened by `within`
    m <- cand[(iv$start[cand] - within) <= pos[i] & pos[i] < (iv$end[cand] + within)]
    if (!length(m)) next
    hit_snp <- c(hit_snp, rep(i, length(m)))
    hit_iv <- c(hit_iv, m)
  }

  cols <- function(n) paste0(prefix, n)
  if (collapse) {
    out <- df
    out[[cols("name")]] <- NA_character_
    out[[cols("n_intervals")]] <- 0L
    if (length(hit_snp)) {
      byi <- split(hit_iv, hit_snp)
      idx <- as.integer(names(byi))
      out[[cols("name")]][idx] <- vapply(byi, function(k)
        paste(iv$name[k[order(iv$start[k])]], collapse = ","), character(1))
      out[[cols("n_intervals")]][idx] <- lengths(byi)
    }
    if (keep == "hits") out <- out[out[[cols("n_intervals")]] > 0, , drop = FALSE]
    rownames(out) <- NULL
    return(tibble::as_tibble(out))
  }

  rows <- if (keep == "hits") hit_snp else {
    miss <- setdiff(seq_len(nrow(df)), hit_snp)
    c(hit_snp, miss)
  }
  ivrow <- if (keep == "hits") hit_iv else c(hit_iv, rep(NA_integer_, nrow(df) - length(unique(hit_snp))))
  ord <- order(rows)
  rows <- rows[ord]; ivrow <- ivrow[ord]

  out <- df[rows, , drop = FALSE]
  out[[cols("name")]] <- iv$name[ivrow]
  out[[cols("interval_start")]] <- iv$start[ivrow]
  out[[cols("interval_end")]] <- iv$end[ivrow]
  mid <- (iv$start[ivrow] + iv$end[ivrow]) / 2
  out[[cols("distance_to_midpoint")]] <- pos[rows] - mid
  if (!is.null(gid)) out[[cols("gene_id")]] <- gid[ivrow]
  rownames(out) <- NULL
  tibble::as_tibble(out)
}
