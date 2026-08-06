# The multilocus index of association (Ia / rbarD): a whole-genome measure of how far the
# population is from random mating, cheap enough to belong on this side.
#
# Its local counterpart, r-squared decay with physical distance, is the one genuinely
# quadratic statistic in the set and lives in the Python package
# (`plasgenomicsutils ld_decay`); [read_ld_decay()] brings the result back for plotting.
# The two answer different questions -- decay is local and physical, Ia is global and
# unlinked -- so they are reported separately rather than blended.

#: SNPs sampled per chromosome before a pairwise scan, to bound an O(n^2) problem.
LD_MAX_SNPS <- 3000L

# Thin a locus table to at most `max_snps` per chromosome, evenly along the chromosome.
.thin_loci <- function(loci, max_snps) {
  do.call(rbind, lapply(split(loci, loci$chr), function(sub) {
    if (nrow(sub) <= max_snps) return(sub)
    sub[round(seq(1, nrow(sub), length.out = max_snps)), , drop = FALSE]
  }))
}

.maf_ok <- function(Hg, cols, maf) {
  f <- .snp_freqs(Hg[, cols, drop = FALSE])
  m <- pmin(f$p, 1 - f$p)
  !is.na(m) & m >= maf
}

.ld_prepare <- function(x, group, meta, genotype, het) {
  if (inherits(x, "PopStructure")) {
    G <- if (is.null(genotype)) x$genotype() else .coerce_geno(genotype)
    if (is.null(meta)) meta <- x$get_meta()
  } else {
    G <- .coerce_geno(x)
  }
  if (is.null(colnames(G)))
    stop("genotypes need `chr:pos` column names", call. = FALSE)
  grp <- .diversity_group(group, meta, rownames(G), nrow(G))
  list(H = .haploid_calls(G, het), loci = .parse_snp_ids(colnames(G)),
       grp = grp, levs = .group_order(grp))
}

#' Multilocus linkage disequilibrium: the index of association
#'
#' `Ia` and its sample-size-standardised form `rbarD` measure whether alleles at
#' *unlinked* loci travel together -- the genome-wide signature of clonal propagation.
#' Both are 0 under free recombination and rise as reproduction becomes more clonal.
#' `rbarD` is the one to compare between datasets, since `Ia` grows with the number of
#' loci. Both are also upward-biased in small samples, so a group of ten will score above
#' a group of a hundred drawn from the same population -- read them beside `n_samples`.
#'
#' A multilocus genotype is only defined where every locus is called, so the loci
#' carrying a gap are dropped and every sample is kept -- over thousands of loci even
#' slight missingness would otherwise leave almost no complete sample, and it is the
#' samples the statistic is about. `n_loci` reports what survived. The estimator is the
#' pairwise-distance decomposition of Brown et al. (1980), standardised by Agapow & Burt
#' (2001).
#'
#' @param x A [PopStructure] or a genotype matrix (samples x SNPs, alt dosage, `chr:pos`
#'   column names).
#' @param group Metadata column name, or a vector aligned to the rows; `NULL` pools
#'   every sample.
#' @param maf Skip loci below this minor-allele frequency within the group.
#' @param max_snps Loci sampled before the pairwise distances, evenly along the genome.
#' @param het How a heterozygous call is read; see [pop_diversity()].
#' @param min_samples Skip groups smaller than this.
#' @param meta,genotype As in [pop_diversity()].
#' @return A tibble with `group`, `n_samples`, `n_loci`, `ia`, `rbar_d`.
#' @references
#' Brown, A. H. D., Feldman, M. W. & Nevo, E. (1980) Multilocus structure of natural
#' populations of \emph{Hordeum spontaneum}. \emph{Genetics} 96, 523-536.
#' \doi{10.1093/genetics/96.2.523}
#'
#' Agapow, P.-M. & Burt, A. (2001) Indices of multilocus linkage disequilibrium.
#' \emph{Molecular Ecology Notes} 1, 101-102. \doi{10.1046/j.1471-8278.2000.00014.x}
#' @examples
#' ps <- example_pop_structure(umap = FALSE)
#' ld_index(ps, group = "country")
#' @export
ld_index <- function(x, group = NULL, maf = 0, max_snps = LD_MAX_SNPS,
                     het = c("missing", "dosage"), min_samples = 4,
                     meta = NULL, genotype = NULL) {
  het <- match.arg(het)
  prep <- .ld_prepare(x, group, meta, genotype, het)
  loci <- .thin_loci(prep$loci, max_snps)

  rows <- list()
  for (l in prep$levs) {
    idx <- which(prep$grp == l)
    if (length(idx) < min_samples) next
    Hg <- prep$H[idx, loci$idx, drop = FALSE]
    if (maf > 0) Hg <- Hg[, .maf_ok(prep$H[idx, , drop = FALSE], loci$idx, maf), drop = FALSE]
    # A multilocus genotype needs every locus called. Drop the loci with a gap rather
    # than the samples carrying one: over thousands of loci even slight missingness
    # leaves almost no complete sample, and it is the samples the statistic is about.
    Hg <- Hg[, colSums(is.na(Hg)) == 0, drop = FALSE]
    Hg <- Hg[, apply(Hg, 2, function(v) length(unique(v)) > 1), drop = FALSE]
    stat <- .ia_rbard(Hg)
    rows[[length(rows) + 1L]] <- data.frame(
      group = l, n_samples = nrow(Hg), n_loci = ncol(Hg),
      ia = stat$ia, rbar_d = stat$rbar_d, stringsAsFactors = FALSE)
  }
  if (!length(rows)) {
    warning("no group had a usable complete-case block", call. = FALSE)
    return(tibble::tibble())
  }
  out <- do.call(rbind, rows)
  out$group <- factor(out$group, levels = prep$levs)
  tibble::as_tibble(out)
}

# Ia and rbarD from the per-locus pairwise mismatch vectors.
.ia_rbard <- function(H) {
  if (nrow(H) < 3 || ncol(H) < 2) return(list(ia = NA_real_, rbar_d = NA_real_))
  # d_j: pairwise mismatches at locus j; D: their sum over loci
  dj <- vapply(seq_len(ncol(H)), function(j) as.numeric(stats::dist(H[, j]) > 0),
               numeric(nrow(H) * (nrow(H) - 1) / 2))
  if (is.null(dim(dj))) dj <- matrix(dj, ncol = ncol(H))
  D <- rowSums(dj)
  # population variances, matching the pairwise-distance formulation
  pvar <- function(v) mean((v - mean(v))^2)
  vo <- pvar(D)
  vj <- apply(dj, 2, pvar)
  ve <- sum(vj)
  if (ve == 0) return(list(ia = NA_real_, rbar_d = NA_real_))
  cross <- sum(outer(sqrt(vj), sqrt(vj))[upper.tri(diag(length(vj)))])
  list(ia = vo / ve - 1,
       rbar_d = if (cross > 0) (vo - ve) / (2 * cross) else NA_real_)
}

#' Read an LD-decay table
#'
#' Reads what `plasgenomicsutils ld_decay` writes, ready for [plot_ld_decay()]. The
#' half-decay distance is attached as an `ld_half_decay` attribute -- recomputed from the
#' binned means if the separate table is not supplied -- and the scan's `max_dist`,
#' `max_snps` and `maf` are attached from the file header, so the thinning behind a curve
#' travels with it.
#'
#' The decay scan lives in the Python package because it is quadratic in the number of
#' SNPs inside each window; on a 249-sample, 28k-SNP callset it runs in about 7 seconds
#' there. `ld_index()` stays here, being linear and quick.
#'
#' @param path TSV(.gz) written by `plasgenomicsutils ld_decay`.
#' @param half_decay Optional TSV from `--half-decay-output`; recomputed when absent.
#' @param levels Optional group order; defaults to a natural sort.
#' @return A tibble of `group`, `bin_start`, `bin_end`, `bin_mid`, `n_pairs`, `mean_r2`,
#'   `median_r2`.
#' @seealso [ld_index()], [plot_ld_decay()]
#' @export
read_ld_decay <- function(path, half_decay = NULL, levels = NULL) {
  .need_package("readr", "read_ld_decay()")
  hdr <- readLines(path, n = 1L, warn = FALSE)
  df <- readr::read_tsv(path, comment = "#", show_col_types = FALSE, progress = FALSE,
                        col_types = readr::cols(group = readr::col_character()))
  need <- c("group", "bin_mid", "mean_r2")
  missing <- setdiff(need, names(df))
  if (length(missing))
    stop(sprintf("`%s` is missing column(s): %s", basename(path),
                 paste(missing, collapse = ", ")), call. = FALSE)
  df$group <- factor(df$group, levels = if (is.null(levels)) .levels_of(df$group) else levels)

  hd <- if (!is.null(half_decay)) {
    h <- readr::read_tsv(half_decay, show_col_types = FALSE, progress = FALSE,
                         col_types = readr::cols(group = readr::col_character()))
    as.data.frame(h)
  } else {
    .half_decay_from_bins(df)
  }
  attr(df, "ld_half_decay") <- hd
  for (k in c("max_dist", "max_snps", "maf")) {
    m <- regmatches(hdr, regexpr(paste0("#", k, "=[^\t]+"), hdr))
    if (length(m)) attr(df, paste0("ld_", k)) <- as.numeric(sub(".*=", "", m))
  }
  df
}

# Distance at which a group's binned mean first drops below half its closest-bin value.
.half_decay_from_bins <- function(df) {
  res <- lapply(split(df, df$group), function(s) {
    s <- s[!is.na(s$mean_r2), , drop = FALSE]
    s <- s[order(s$bin_mid), , drop = FALSE]
    if (nrow(s) < 2) return(NA_real_)
    target <- s$mean_r2[1] / 2
    below <- which(s$mean_r2 <= target)
    if (!length(below)) return(NA_real_)          # still above half at max_dist
    k <- below[1]
    if (k == 1) return(s$bin_mid[1])
    x0 <- s$bin_mid[k - 1]; x1 <- s$bin_mid[k]
    y0 <- s$mean_r2[k - 1]; y1 <- s$mean_r2[k]
    if (y0 == y1) return(x1)
    x0 + (target - y0) * (x1 - x0) / (y1 - y0)
  })
  data.frame(group = names(res), half_decay_bp = unlist(res, use.names = FALSE),
             row.names = NULL, stringsAsFactors = FALSE)
}
