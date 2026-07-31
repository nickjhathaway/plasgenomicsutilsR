# IBD post-analysis results container. Ingests the tables produced by the Python
# `plasgenomicsutils ibd` tools and precomputes the cumulative-genome coordinate
# every genome-wide plot shares, so the plot_*() functions stay thin.

#' @importFrom R6 R6Class
NULL

# ---- internal helpers ------------------------------------------------------

.read_maybe <- function(x, what) {
  if (is.null(x)) return(NULL)
  if (is.character(x)) {
    .need_package("readr", "reading IBD tables from a path")
    return(tibble::as_tibble(readr::read_tsv(x, show_col_types = FALSE, progress = FALSE)))
  }
  if (is.data.frame(x)) return(tibble::as_tibble(x))
  stop("expected a file path, a data frame, or NULL, not ", class(x)[1], call. = FALSE)
}

.require_cols <- function(df, cols, what) {
  missing <- setdiff(cols, names(df))
  if (length(missing)) {
    stop(sprintf("%s is missing required column(s): %s", what,
                 paste(missing, collapse = ", ")), call. = FALSE)
  }
  invisible(df)
}

# Cumulative-genome layout (offsets, chromosome bands, axis mid-points) from a
# reference's core chromosome lengths, ordered numerically.
.chrom_layout <- function(reference) {
  lens <- get_reference(reference)$core_chrom_lengths_bp
  chr <- names(lens)
  ord <- order(suppressWarnings(as.integer(chr)))
  chr <- chr[ord]
  len <- as.numeric(lens[ord])
  offset <- cumsum(c(0, len))[seq_along(len)]
  tibble::tibble(
    chr = chr, len = len, offset = offset,
    xmin = offset, xmax = offset + len, mid = offset + len / 2,
    band = rep(c("a", "b"), length.out = length(chr))
  )
}

# Add a cumulative-genome x-coordinate (`cum_pos`) to a table carrying chr + pos.
.add_cum_pos <- function(df, layout, chr_col = "chr", pos_col = "pos") {
  if (is.null(df)) return(NULL)
  key <- normalise_chr(df[[chr_col]])
  df$cum_pos <- layout$offset[match(key, layout$chr)] + as.numeric(df[[pos_col]])
  df$chr <- key
  df
}

# Normalise the many shapes a per-region significance threshold arrives in
# (scalar, named vector, or a data frame of region + threshold) to a tibble.
.normalise_threshold <- function(threshold) {
  if (is.null(threshold)) return(NULL)
  if (is.data.frame(threshold)) {
    df <- tibble::as_tibble(threshold)
    tcol <- intersect(c("threshold", "neg_log10_p_threshold"), names(df))[1]
    if (is.na(tcol)) stop("threshold data frame needs a 'threshold' or ",
                          "'neg_log10_p_threshold' column", call. = FALSE)
    rcol <- intersect(c("region", "population"), names(df))[1]
    if (is.na(rcol)) return(tibble::tibble(region = NA_character_, threshold = df[[tcol]][1]))
    return(tibble::tibble(region = as.character(df[[rcol]]), threshold = as.numeric(df[[tcol]])))
  }
  if (is.character(threshold) && length(threshold) == 1L && file.exists(threshold)) {
    return(.normalise_threshold(as.numeric(readLines(threshold, warn = FALSE)[1])))
  }
  vals <- as.numeric(threshold)
  if (!is.null(names(threshold))) {
    return(tibble::tibble(region = names(threshold), threshold = vals))
  }
  tibble::tibble(region = NA_character_, threshold = vals[1])
}

# ---- the class -------------------------------------------------------------

#' IBD post-analysis results
#'
#' A container for the per-SNP, pairwise-region, and selection-statistic tables
#' produced by the Python `plasgenomicsutils ibd` tools, plus the shared
#' cumulative-genome coordinate the genome-wide plots use. Pass file paths or
#' data frames; each argument is optional so you can hold only what you plan to
#' plot. The `plot_*()` functions read from an object of this class.
#'
#' @details Expected columns (superset; extras are kept):
#'   * `per_snp_region`: `chr`, `pos`, `frac_pairs_ibd`, optionally `region`
#'     (output of `analyze_ibd_matrix` per-SNP / per-SNP-per-region).
#'   * `pairwise_region`: `chr`, `pos`, `region_a`, `region_b`, `frac_pairs_ibd`
#'     (per-SNP pairwise-region output).
#'   * `selection`: `chr`, `pos`, a metric column (`neg_log10_p`, `chi2_stat`, or
#'     `z_score`), optionally `region` and `significant`
#'     (output of `ibd_selection_statistic`).
#'
#' @export
IbdResults <- R6::R6Class(
  "IbdResults",
  public = list(
    #' @description Create an IbdResults object.
    #' @param per_snp_region Per-SNP (optionally per-region) IBD table: path or data frame.
    #' @param pairwise_region Per-SNP pairwise-region IBD table: path or data frame.
    #' @param selection IBD selection-statistic table: path or data frame.
    #' @param threshold Significance threshold(s): a scalar, a named vector or
    #'   `region`/`threshold` data frame (per region), or a path to a one-number file.
    #' @param genes Optional gene-annotation track (`name`, `chr`, `start`, `end`)
    #'   drawn as reference lines on the Manhattan plots.
    #' @param reference Reference id for chromosome lengths (default `"pf3d7"`).
    #' @return An `IbdResults` object (invisibly self).
    initialize = function(per_snp_region = NULL, pairwise_region = NULL,
                          selection = NULL, threshold = NULL, genes = NULL,
                          reference = "pf3d7") {
      private$reference <- reference
      private$layout <- .chrom_layout(reference)

      per_snp <- .read_maybe(per_snp_region, "per_snp_region")
      if (!is.null(per_snp)) {
        .require_cols(per_snp, c("chr", "pos", "frac_pairs_ibd"), "per_snp_region")
        per_snp <- .add_cum_pos(per_snp, private$layout)
      }
      private$per_snp <- per_snp

      pw <- .read_maybe(pairwise_region, "pairwise_region")
      if (!is.null(pw)) {
        .require_cols(pw, c("chr", "pos", "region_a", "region_b", "frac_pairs_ibd"),
                      "pairwise_region")
        pw <- .add_cum_pos(pw, private$layout)
      }
      private$pairwise <- pw

      sel <- .read_maybe(selection, "selection")
      if (!is.null(sel)) {
        .require_cols(sel, c("chr", "pos"), "selection")
        sel <- .add_cum_pos(sel, private$layout)
      }
      private$selection <- sel

      genes <- .read_maybe(genes, "genes")
      if (!is.null(genes)) {
        .require_cols(genes, c("chr", "start", "end"), "genes")
        genes$chr <- normalise_chr(genes$chr)
        if (!"name" %in% names(genes)) {
          genes$name <- paste0(genes$chr, ":", genes$start, "-", genes$end)
        }
        genes$cum_mid <- private$layout$offset[match(genes$chr, private$layout$chr)] +
          (as.numeric(genes$start) + as.numeric(genes$end)) / 2
      }
      private$genes <- genes
      private$threshold <- .normalise_threshold(threshold)
      invisible(self)
    },

    #' @description Per-SNP (optionally per-region) IBD table with `cum_pos`.
    get_per_snp_region = function() private$per_snp,
    #' @description Per-SNP pairwise-region IBD table with `cum_pos`.
    get_pairwise_region = function() private$pairwise,
    #' @description Selection-statistic table with `cum_pos`.
    get_selection = function() private$selection,
    #' @description Threshold tibble (`region`, `threshold`) or `NULL`.
    get_thresholds = function() private$threshold,
    #' @description Gene-annotation track with `cum_mid`, or `NULL`.
    get_genes = function() private$genes,
    #' @description Chromosome layout tibble (offsets, bands, axis mid-points).
    chrom_layout = function() private$layout,
    #' @description The reference id used for chromosome lengths.
    reference_id = function() private$reference,

    #' @description Compact summary of what the object holds.
    print = function(...) {
      cat("<IbdResults>  reference:", private$reference, "\n")
      shape <- function(df) if (is.null(df)) "-" else paste(nrow(df), "rows")
      cat("  per_snp_region :", shape(private$per_snp), "\n")
      cat("  pairwise_region:", shape(private$pairwise), "\n")
      cat("  selection      :", shape(private$selection), "\n")
      if (!is.null(private$threshold)) cat("  thresholds     :", nrow(private$threshold), "\n")
      if (!is.null(private$genes)) cat("  genes          :", nrow(private$genes), "\n")
      invisible(self)
    }
  ),
  private = list(
    reference = NULL, layout = NULL, per_snp = NULL, pairwise = NULL,
    selection = NULL, threshold = NULL, genes = NULL
  )
)

#' Create an [IbdResults] object
#'
#' Convenience wrapper for `IbdResults$new()`.
#'
#' @inheritParams IbdResults
#' @param ... Passed to [IbdResults]'s constructor.
#' @return An [IbdResults] object.
#' @export
ibd_results <- function(...) {
  IbdResults$new(...)
}
