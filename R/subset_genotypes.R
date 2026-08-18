# Narrowing a genotype panel to a set of samples -- most often the monoclonals a haplotype
# set kept, so a diversity or differentiation run reads the same cohort the scans did.

# Anything that names a set of samples -> the ids.
.as_sample_ids <- function(x) {
  if (is.character(x)) return(unique(x))
  if (inherits(x, "parasite_haplotypes")) return(haplotype_samples(x))
  if (inherits(x, "PopStructure")) return(x$get_samples())
  if (is.matrix(x)) return(rownames(x))
  if (is.list(x) && !is.null(x$sample.id)) return(as.character(x$sample.id))
  if (is.list(x) && !is.null(x$genotype)) return(rownames(as.matrix(x$genotype)))
  if (is.data.frame(x) && "sample" %in% names(x)) return(as.character(x$sample))
  stop("cannot read sample ids from an object of class ", class(x)[1], call. = FALSE)
}

#' The samples a haplotype set kept
#'
#' [parasite_haplotypes()] drops samples twice over: the polyclonal ones the Fws gate
#' removes, and any left too incomplete by `max_sample_missing`. This is what came through
#' both -- the cohort the EHH and iHS scans actually ran on.
#'
#' It is the set to reuse when another analysis should describe the same samples.
#' `hap$filtering` says how many went to each cause, so it is visible whether the Fws gate
#' was the whole story.
#'
#' @param hap A [parasite_haplotypes()] object.
#' @return A character vector of sample ids, in the order the haplotype matrix holds them.
#' @seealso [subset_genotypes()] to narrow a panel to them, [parasite_haplotypes()].
#' @examples
#' ps <- example_pop_structure(umap = FALSE)
#' hap <- parasite_haplotypes(ps, maf = 0.05)
#' length(haplotype_samples(hap))
#' hap$filtering[c("n_dropped_polyclonal", "n_dropped_sample_missing")]
#' @export
haplotype_samples <- function(hap) {
  if (!inherits(hap, "parasite_haplotypes"))
    stop("`hap` is not a parasite_haplotypes object (got ", class(hap)[1], ")",
         call. = FALSE)
  rownames(hap$hap)
}

#' Restrict a genotype panel to a set of samples
#'
#' Narrows the rows and leaves everything else alone, returning the same kind of object it
#' was given -- a [load_genotypes()] list comes back a list, with `sample.id` narrowed to
#' match and `snp.id` / `allele` / `pruned` / `positions` / `variants` carried through.
#'
#' The use it was written for is running a diversity or differentiation analysis on the
#' monoclonals a haplotype set kept:
#'
#' ```
#' mono <- parasite_haplotypes(geno$genotype, fws = fws, min_fws = 0.92)
#' pop_diversity(subset_genotypes(geno, mono)$genotype, group = "region", meta = meta)
#' ```
#'
#' Note it subsets **samples**, not SNPs: the full SNP panel is kept, which is what a
#' diversity estimate wants. Handing `hap$hap` to such an analysis instead would silently
#' use the MAF-filtered, imputed 0/1 matrix the scans need, which is a different panel.
#'
#' @param x A [load_genotypes()] list, a genotype matrix (samples x SNPs), or a
#'   [PopStructure] (which is narrowed with its own `$subset()`).
#' @param samples The samples to keep: a character vector of ids, a
#'   [parasite_haplotypes()] object (its kept samples), a `PopStructure`, another
#'   `load_genotypes()` list, or a metadata data frame with a `sample` column.
#' @param strict Error when `samples` names ids the panel does not hold (default `FALSE`,
#'   which warns and keeps the intersection).
#' @return The same kind of object as `x`, holding only those samples, in the panel's own
#'   order.
#' @seealso [haplotype_samples()], [PopStructure]'s `$subset()`.
#' @examples
#' ps <- example_pop_structure(umap = FALSE)
#' hap <- parasite_haplotypes(ps, maf = 0.05)
#'
#' # a bare matrix
#' dim(subset_genotypes(ps$genotype(), hap))
#'
#' # a load_genotypes()-shaped list keeps its other slots
#' geno <- list(genotype = ps$genotype(), sample.id = ps$get_samples(),
#'              snp.id = colnames(ps$genotype()), pruned = TRUE, positions = "0-based")
#' str(subset_genotypes(geno, hap)[c("sample.id", "pruned", "positions")], max.level = 1)
#' @export
subset_genotypes <- function(x, samples, strict = FALSE) {
  ids <- .as_sample_ids(samples)
  if (!length(ids)) stop("`samples` names no samples", call. = FALSE)

  if (inherits(x, "PopStructure")) return(x$subset(samples = ids))

  is_list <- is.list(x) && !is.data.frame(x) && !is.null(x$genotype)
  G <- if (is_list) as.matrix(x$genotype) else as.matrix(x)
  have <- rownames(G)
  if (is.null(have)) {
    if (is_list && !is.null(x$sample.id)) have <- as.character(x$sample.id)
    else stop("the genotype matrix has no sample row names to match against",
              call. = FALSE)
  }

  missing <- setdiff(ids, have)
  if (length(missing)) {
    msg <- paste0(length(missing), " of ", length(ids),
                  " requested samples are not in the panel: ",
                  paste(utils::head(missing, 5), collapse = ", "),
                  if (length(missing) > 5) ", ..." else "")
    if (strict) stop(msg, call. = FALSE)
    warning(msg, call. = FALSE)
  }
  # the panel's own order, not the order they were asked for, so a genotype matrix and the
  # metadata joined to it stay in step
  keep <- have[have %in% ids]
  if (!length(keep)) stop("none of those samples are in the panel", call. = FALSE)
  if (length(keep) < length(have))
    message("keeping ", length(keep), " of ", length(have), " samples")

  idx <- match(keep, have)
  Gk <- G[idx, , drop = FALSE]
  rownames(Gk) <- keep
  if (!is_list) return(Gk)

  out <- x
  out$genotype <- Gk
  if (!is.null(out$sample.id)) out$sample.id <- keep
  out
}
