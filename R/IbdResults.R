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

# Conventional thresholds for discarding short, SNP-poor IBD segments, which are commonly
# spurious. Applied to the IBD evidence only -- `analyzed_samples`, the denominator behind
# every fraction, is taken from the unfiltered blocks first.
IBD_MIN_BLOCK_SNP <- 15L
IBD_MIN_BLOCK_KB <- 15

.filter_ibd_blocks <- function(bl, min_snp = IBD_MIN_BLOCK_SNP, min_kb = IBD_MIN_BLOCK_KB) {
  if (!nrow(bl)) return(bl)
  keep <- rep(TRUE, nrow(bl))
  if (!is.null(min_kb) && is.finite(min_kb) && min_kb > 0) {
    keep <- keep & (as.numeric(bl$end) - as.numeric(bl$start)) >= min_kb * 1000
  }
  if (!is.null(min_snp) && is.finite(min_snp) && min_snp > 0) {
    if (!"Nsnp" %in% names(bl)) {
      warning("min_block_snp = ", min_snp, " was requested but the blocks have no 'Nsnp' ",
              "column; only the length filter was applied", call. = FALSE)
    } else {
      keep <- keep & as.numeric(bl$Nsnp) >= min_snp
    }
  }
  bl[keep, , drop = FALSE]
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
    # the FDR cutoff and the calibration diagnostic ride along when present, so a plot can
    # draw either line and the printed object can show why neither may be trustworthy
    extra <- intersect(c("neg_log10_p_fdr_threshold", "fdr_alpha", "alpha", "n_tests",
                         "n_significant", "n_significant_fdr", "lambda_gc"), names(df))
    out <- if (is.na(rcol))
      tibble::tibble(group = NA_character_, threshold = as.numeric(df[[tcol]])[1])
    else
      tibble::tibble(group = as.character(df[[rcol]]), threshold = as.numeric(df[[tcol]]))
    for (k in extra) out[[k]] <- if (is.na(rcol)) df[[k]][1] else df[[k]]
    return(out)
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
    #' @param min_block_snp,min_block_kb Discard IBD segments carrying fewer than
    #'   `min_block_snp` SNPs or shorter than `min_block_kb` kb (defaults `15` and `15`).
    #'   Short, SNP-poor segments are commonly spurious, and this filter is conventionally
    #'   applied before any IBD summary -- it is built in so the blocks need not be
    #'   pre-filtered. `0` disables either criterion. Only the IBD evidence is filtered:
    #'   the analyzed-sample set (the denominator behind every fraction) still comes from
    #'   every row of the blocks file, so a pair whose only segment is short still counts
    #'   as compared. The SNP criterion needs the `Nsnp` column hmmibd-rs writes.
    #' @param group_col_in_meta Name of the `meta` column that defines the grouping. It
    #'   becomes the default `group` for the block-based tools, and if the column is a
    #'   factor its levels set the group order for every loaded table (equivalent to
    #'   calling `$set_group_order(levels(meta[[group_col_in_meta]]))`).
    #' @param reference Reference id for chromosome lengths (default `"pf3d7"`).
    #' @return An `IbdResults` object (invisibly self).
    initialize = function(per_snp_group = NULL, pairwise_group = NULL,
                          selection = NULL, threshold = NULL, genes = NULL,
                          blocks = NULL, meta = NULL, gene_overlap = NULL,
                          group_col_in_meta = NULL,
                          min_block_snp = IBD_MIN_BLOCK_SNP,
                          min_block_kb = IBD_MIN_BLOCK_KB,
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
        n_before <- nrow(bl)
        bl <- .filter_ibd_blocks(bl, min_block_snp, min_block_kb)
        private$block_filter <- c(min_snp = min_block_snp %||% 0, min_kb = min_block_kb %||% 0,
                                  dropped = n_before - nrow(bl), kept = nrow(bl))
        private$blocks <- bl[, c("sample1", "sample2", "chr", "start", "end"), drop = FALSE]
      }
      private$meta <- .read_maybe(meta, "meta")
      private$gene_overlap <- .read_maybe(gene_overlap, "gene_overlap")

      private$threshold <- .normalise_threshold(threshold)

      if (!is.null(group_col_in_meta)) {
        if (is.null(private$meta) || !group_col_in_meta %in% names(private$meta)) {
          stop("group_col_in_meta = '", group_col_in_meta, "' is not a column of meta",
               call. = FALSE)
        }
        private$group_col <- group_col_in_meta
        # a factor column carries the intended order; anything else is natural-sorted
        private$adopt_group_order(.levels_of(private$meta[[group_col_in_meta]]),
                                  sprintf("meta$%s", group_col_in_meta))
      }
      invisible(self)
    },

    #' @description Set the order of the groups for every loaded table, so facets,
    #'   legends and axes follow it. Errors if a group present in the results is missing
    #'   from `levels` (it would silently become `NA`); warns about levels no result uses.
    #' @param levels Group names in the desired order.
    #' @return Invisibly self.
    set_group_order = function(levels) {
      private$adopt_group_order(levels, "the requested order")
      invisible(self)
    },

    #' @description The current group order (`NULL` when none has been set).
    get_group_order = function() private$group_order,

    #' @description The `meta` column naming the grouping, if one was declared.
    get_group_col = function() private$group_col,

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
      if (!is.null(private$blocks)) {
        cat("  IBD blocks     :", nrow(private$blocks), "rows,",
            length(private$analyzed_samples), "samples\n")
        bf <- private$block_filter
        # the filter changes every downstream fraction, so say what it removed
        if (!is.null(bf) && bf[["dropped"]] > 0) {
          cat(sprintf("                   dropped %d short/SNP-poor segment(s) (>= %g SNPs, >= %g kb)\n",
                      bf[["dropped"]], bf[["min_snp"]], bf[["min_kb"]]))
        }
      }
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

    #' @description Sample pairs sharing IBD over each gene (see [gene_ibd_pairs()]).
    #' @param ... Passed to [gene_ibd_pairs()].
    gene_ibd_pairs = function(...) gene_ibd_pairs(self, ...),

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
    blocks = NULL, analyzed_samples = NULL, meta = NULL, gene_overlap = NULL,
    group_order = NULL, group_col = NULL, block_filter = NULL,

    # every group label appearing anywhere in the loaded tables
    groups_present = function() {
      g <- character(0)
      one <- function(df, cols) {
        if (is.null(df)) return(character(0))
        unlist(lapply(cols[cols %in% names(df)], function(cc) as.character(df[[cc]])))
      }
      g <- c(g, one(private$per_snp, "group"), one(private$selection, "group"),
             one(private$pairwise, c("group_a", "group_b")),
             one(private$gene_overlap, c("group_a", "group_b")),
             one(private$threshold, "group"))
      unique(g[!is.na(g)])
    },

    # stamp the order onto every group column, so plots, facets and legends all follow it
    apply_group_order = function() {
      levs <- private$group_order
      if (is.null(levs)) return(invisible(NULL))
      set <- function(df, cols) {
        if (is.null(df)) return(df)
        for (cc in cols[cols %in% names(df)]) {
          df[[cc]] <- factor(as.character(df[[cc]]), levels = levs)
        }
        df
      }
      private$per_snp <- set(private$per_snp, "group")
      private$selection <- set(private$selection, "group")
      private$pairwise <- set(private$pairwise, c("group_a", "group_b"))
      private$gene_overlap <- set(private$gene_overlap, c("group_a", "group_b"))
      # the per-group threshold rows become a geom_hline layer of their own; leaving them
      # character makes ggplot merge the two layers' facet values alphabetically
      private$threshold <- set(private$threshold, "group")
      if (!is.null(private$meta) && !is.null(private$group_col) &&
          private$group_col %in% names(private$meta)) {
        private$meta[[private$group_col]] <-
          factor(as.character(private$meta[[private$group_col]]), levels = levs)
      }
      invisible(NULL)
    },

    # shared by set_group_order() and the constructor's group_col_in_meta
    adopt_group_order = function(levs, source) {
      levs <- as.character(levs)
      levs <- levs[!is.na(levs)]
      if (anyDuplicated(levs)) {
        stop("duplicated group(s) in the requested order: ",
             paste(unique(levs[duplicated(levs)]), collapse = ", "), call. = FALSE)
      }
      present <- private$groups_present()
      dropped <- setdiff(present, levs)
      if (length(dropped)) {
        stop("group(s) in the IBD results are missing from ", source, ": ",
             paste(sort(dropped), collapse = ", "),
             ". They would become NA; add them to keep the results intact.", call. = FALSE)
      }
      unused <- setdiff(levs, present)
      if (length(unused) && length(present)) {
        warning(source, " has group(s) that no IBD result uses: ",
                paste(unused, collapse = ", "), call. = FALSE)
      }
      private$group_order <- levs
      private$apply_group_order()
      invisible(NULL)
    }
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
