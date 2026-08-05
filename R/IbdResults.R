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

# Normalise the many shapes a per-group significance threshold arrives in
# (scalar, named vector, or a data frame of group + threshold) to a tibble.
.normalise_threshold <- function(threshold) {
  if (is.null(threshold)) return(NULL)
  if (is.data.frame(threshold)) {
    df <- tibble::as_tibble(threshold)
    tcol <- intersect(c("threshold", "neg_log10_p_threshold"), names(df))[1]
    if (is.na(tcol)) stop("threshold data frame needs a 'threshold' or ",
                          "'neg_log10_p_threshold' column", call. = FALSE)
    rcol <- intersect(c("group", "population"), names(df))[1]
    if (is.na(rcol)) return(tibble::tibble(group = NA_character_, threshold = df[[tcol]][1]))
    return(tibble::tibble(group = as.character(df[[rcol]]), threshold = as.numeric(df[[tcol]])))
  }
  if (is.character(threshold) && length(threshold) == 1L && file.exists(threshold)) {
    first <- readLines(threshold, n = 1L, warn = FALSE)
    # a tabular threshold file (e.g. group, alpha, n_tests, neg_log10_p_threshold as written
    # by `ibd_selection_statistic`) vs a bare scalar in a text file
    if (grepl("\t", first) || grepl("threshold", first, ignore.case = TRUE)) {
      return(.normalise_threshold(utils::read.delim(threshold, stringsAsFactors = FALSE)))
    }
    return(.normalise_threshold(as.numeric(first)))
  }
  vals <- as.numeric(threshold)
  if (!is.null(names(threshold))) {
    return(tibble::tibble(group = names(threshold), threshold = vals))
  }
  tibble::tibble(group = NA_character_, threshold = vals[1])
}

# ---- the class -------------------------------------------------------------

#' IBD post-analysis results
#'
#' A container for the per-SNP, pairwise-group, and selection-statistic tables
#' produced by the Python `plasgenomicsutils ibd` tools, plus the shared
#' cumulative-genome coordinate the genome-wide plots use. Pass file paths or
#' data frames; each argument is optional so you can hold only what you plan to
#' plot. The `plot_*()` functions read from an object of this class.
#'
#' @details Expected columns (superset; extras are kept):
#'   * `per_snp_group`: `chr`, `pos`, `frac_pairs_ibd`, optionally `group`
#'     (output of `analyze_ibd_matrix` per-SNP / per-SNP-per-group).
#'   * `pairwise_group`: `chr`, `pos`, `group_a`, `group_b`, `frac_pairs_ibd`
#'     (per-SNP pairwise-group output).
#'   * `selection`: `chr`, `pos`, a metric column (`neg_log10_p`, `chi2_stat`, or
#'     `z_score`), optionally `group` and `significant`
#'     (output of `ibd_selection_statistic`).
#'
#' @export
IbdResults <- R6::R6Class(
  "IbdResults",
  public = list(
    #' @description Create an IbdResults object.
    #' @param per_snp_group Per-SNP (optionally per-group) IBD table: path or data frame.
    #' @param pairwise_group Per-SNP pairwise-group IBD table: path or data frame.
    #' @param selection IBD selection-statistic table: path or data frame.
    #' @param threshold Significance threshold(s): a scalar, a named vector or
    #'   `group`/`threshold` data frame (per group), or a path to a one-number file.
    #' @param genes Optional gene-annotation track (`name`, `chr`, `start`, `end`)
    #'   drawn as reference lines on the Manhattan plots.
    #' @param blocks Optional hmmibd-rs IBD segments (path or data frame:
    #'   `sample1`, `sample2`, `chr`, `start`, `end`, optional `different`). Enables
    #'   block-based gene triangles ([gene_ibd_overlap()] / [plot_pairwise_ibd_for_genes()]):
    #'   only IBD segments (`different == 0`) are kept, but the analyzed-sample set (the
    #'   denominator) is taken from every row first, so pairs that are never IBD still count.
    #' @param meta Optional sample metadata (path or data frame with a `sample` column plus
    #'   grouping columns), used to group `blocks` pairs.
    #' @param gene_overlap Optional precomputed per-gene per-group-pair block-overlap table
    #'   (from `plasgenomicsutils ibd_gene_overlap`: `gene`, `group_a`, `group_b`,
    #'   `frac_pairs_ibd`, ...), used directly by the gene triangles.
    #' @param reference Reference id for chromosome lengths (default `"pf3d7"`).
    #' @return An `IbdResults` object (invisibly self).
    initialize = function(per_snp_group = NULL, pairwise_group = NULL,
                          selection = NULL, threshold = NULL, genes = NULL,
                          blocks = NULL, meta = NULL, gene_overlap = NULL,
                          reference = "pf3d7") {
      private$reference <- reference
      private$layout <- .chrom_layout(reference)

      per_snp <- .read_maybe(per_snp_group, "per_snp_group")
      if (!is.null(per_snp)) {
        .require_cols(per_snp, c("chr", "pos", "frac_pairs_ibd"), "per_snp_group")
        per_snp <- .add_cum_pos(per_snp, private$layout)
      }
      private$per_snp <- per_snp

      pw <- .read_maybe(pairwise_group, "pairwise_group")
      if (!is.null(pw)) {
        .require_cols(pw, c("chr", "pos", "group_a", "group_b", "frac_pairs_ibd"),
                      "pairwise_group")
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
        # accept `chrom` (as in PF3D7_GENES / PF_EXAMPLE_DRUG_GENES) as an alias for `chr`
        if (!"chr" %in% names(genes) && "chrom" %in% names(genes)) genes$chr <- genes$chrom
        .require_cols(genes, c("chr", "start", "end"), "genes")
        genes$chr <- normalise_chr(genes$chr)
        if (!"name" %in% names(genes)) {
          genes$name <- paste0(genes$chr, ":", genes$start, "-", genes$end)
        }
        genes$cum_mid <- private$layout$offset[match(genes$chr, private$layout$chr)] +
          (as.numeric(genes$start) + as.numeric(genes$end)) / 2
      }
      private$genes <- genes

      bl <- .read_maybe(blocks, "blocks")
      if (!is.null(bl)) {
        .require_cols(bl, c("sample1", "sample2", "chr", "start", "end"), "blocks")
        private$analyzed_samples <- unique(c(bl$sample1, bl$sample2))   # before IBD filter
        if ("different" %in% names(bl)) bl <- bl[bl$different == 0, , drop = FALSE]
        bl$chr <- normalise_chr(bl$chr)
        # hmmibd-rs reports both endpoints as the 0-based position of the first and last
        # SNP in the segment; shift `end` to make the interval half-open [start, end) so it
        # matches the gene / region tables and every interval test in the package
        bl$end <- as.numeric(bl$end) + 1
        private$blocks <- bl[, c("sample1", "sample2", "chr", "start", "end"), drop = FALSE]
      }
      private$meta <- .read_maybe(meta, "meta")
      private$gene_overlap <- .read_maybe(gene_overlap, "gene_overlap")

      private$threshold <- .normalise_threshold(threshold)
      invisible(self)
    },

    #' @description Per-SNP (optionally per-group) IBD table with `cum_pos`.
    get_per_snp_group = function() private$per_snp,
    #' @description Per-SNP pairwise-group IBD table with `cum_pos`.
    get_pairwise_group = function() private$pairwise,
    #' @description Selection-statistic table with `cum_pos`.
    get_selection = function() private$selection,
    #' @description Threshold tibble (`group`, `threshold`) or `NULL`.
    get_thresholds = function() private$threshold,
    #' @description Gene-annotation track with `cum_mid`, or `NULL`.
    get_genes = function() private$genes,
    #' @description IBD segment table (`sample1`, `sample2`, `chr`, `start`, `end`), or `NULL`.
    get_blocks = function() private$blocks,
    #' @description Analyzed-sample ids (from every block row, pre IBD filter), or `NULL`.
    get_analyzed_samples = function() private$analyzed_samples,
    #' @description Sample metadata data frame, or `NULL`.
    get_meta = function() private$meta,
    #' @description Precomputed per-gene block-overlap table, or `NULL`.
    get_gene_overlap = function() private$gene_overlap,
    #' @description Chromosome layout tibble (offsets, bands, axis mid-points).
    chrom_layout = function() private$layout,
    #' @description The reference id used for chromosome lengths.
    reference_id = function() private$reference,

    #' @description Compact summary of what the object holds.
    print = function(...) {
      cat("<IbdResults>  reference:", private$reference, "\n")
      shape <- function(df) if (is.null(df)) "-" else paste(nrow(df), "rows")
      cat("  per_snp_group :", shape(private$per_snp), "\n")
      cat("  pairwise_group:", shape(private$pairwise), "\n")
      cat("  selection      :", shape(private$selection), "\n")
      if (!is.null(private$threshold)) cat("  thresholds     :", nrow(private$threshold), "\n")
      if (!is.null(private$genes)) cat("  genes          :", nrow(private$genes), "\n")
      if (!is.null(private$blocks)) cat("  IBD blocks     :", nrow(private$blocks), "rows,",
                                        length(private$analyzed_samples), "samples\n")
      if (!is.null(private$gene_overlap)) cat("  gene_overlap   :",
                                              nrow(private$gene_overlap), "rows\n")
      invisible(self)
    },

    # ---- plotting: thin methods over the plot_*() functions -----------------
    # Each forwards `self` to the same-named exported function, so `ibd$plot_*()`
    # and `plot_*(ibd)` are interchangeable (see those functions for arguments).

    #' @description Genome-wide per-SNP IBD-sharing Manhattan. See [plot_ibd_sharing_manhattan()].
    #' @param ... Passed to [plot_ibd_sharing_manhattan()].
    plot_ibd_sharing_manhattan = function(...) plot_ibd_sharing_manhattan(self, ...),

    #' @description IBD selection-statistic Manhattan. See [plot_selection_manhattan()].
    #' @param ... Passed to [plot_selection_manhattan()].
    plot_selection_manhattan = function(...) plot_selection_manhattan(self, ...),

    #' @description Selection/IBD "tug-of-war" mirror. See [plot_ibd_tugofwar()].
    #' @param ... Passed to [plot_ibd_tugofwar()].
    plot_ibd_tugofwar = function(...) plot_ibd_tugofwar(self, ...),

    #' @description Group x group IBD heatmap along the genome. See [plot_ibd_pairwise_group_heatmap()].
    #' @param ... Passed to [plot_ibd_pairwise_group_heatmap()].
    plot_ibd_pairwise_group_heatmap = function(...) plot_ibd_pairwise_group_heatmap(self, ...),

    #' @description Per-gene (or per-SNP) group x group IBD triangles. See [plot_pairwise_ibd_for_genes()].
    #' @param ... Passed to [plot_pairwise_ibd_for_genes()].
    plot_pairwise_ibd_for_genes = function(...) plot_pairwise_ibd_for_genes(self, ...),

    #' @description Per-gene IBD-block overlap between groups (see [gene_ibd_overlap()]).
    #' @param ... Passed to [gene_ibd_overlap()].
    gene_ibd_overlap = function(...) gene_ibd_overlap(self, ...),

    #' @description Sample-level IBD network at a gene / locus (see [plot_ibd_network()]).
    #' @param ... Passed to [plot_ibd_network()].
    plot_ibd_network = function(...) plot_ibd_network(self, ...),

    #' @description Genes overlapping (or within `within` bp of) above-threshold
    #'   selection SNPs. See [pos_selection_genes()].
    #' @param ... Passed to [pos_selection_genes()].
    pos_selection_genes = function(...) pos_selection_genes(self, ...)
  ),
  private = list(
    reference = NULL, layout = NULL, per_snp = NULL, pairwise = NULL,
    selection = NULL, threshold = NULL, genes = NULL,
    blocks = NULL, analyzed_samples = NULL, meta = NULL, gene_overlap = NULL
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
