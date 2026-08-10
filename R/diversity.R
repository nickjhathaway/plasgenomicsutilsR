# Within-population diversity: nucleotide diversity, expected heterozygosity, Watterson's
# theta, Tajima's D, and haplotype/multilocus-genotype diversity, summarised over metadata
# groups and reported per gene, per window, or genome-wide.
#
# Two things here are easy to get wrong and are handled explicitly:
#
#  * pi is *per base pair*, so its denominator is the number of accessible sites, not the
#    number of SNPs. Dividing by the SNP count gives a number that changes with SNP
#    density and missingness and is not comparable between windows -- so that quantity is
#    reported too, but under its own name (`he`), never as pi.
#  * the parasite is haploid, so a heterozygous call is a mixed infection at that site
#    rather than a diploid genotype. By default those calls are treated as missing (see
#    the `het` argument); every count of "gene copies" is then a count of samples.

#: Tajima's D needs a few segregating sites before it means anything.
DIVERSITY_MIN_SNPS <- 3L

# Alt-allele dosages (0/1/2, NA) -> haploid 0/1/NA calls.
.haploid_calls <- function(G, het = "missing") {
  H <- matrix(NA_real_, nrow(G), ncol(G), dimnames = dimnames(G))
  H[G == 0] <- 0
  H[G == 2] <- 1
  if (identical(het, "dosage")) H[G == 1] <- 0.5   # split a mixed call between alleles
  H
}

# Per-SNP allele frequency and called-sample count for one block of samples.
.snp_freqs <- function(H) {
  n <- colSums(!is.na(H))
  p <- ifelse(n > 0, colSums(H, na.rm = TRUE) / n, NA_real_)
  list(p = p, n = n)
}

# Unbiased per-site heterozygosity == per-site pi: 2p(1-p) * n/(n-1).
.site_het <- function(p, n) {
  h <- 2 * p * (1 - p) * n / (n - 1)
  h[!is.finite(h) | n < 2] <- NA_real_
  h
}

#' Tajima's D from segregating sites and per-site heterozygosity
#'
#' @param h Per-site heterozygosity at each usable site (see [pop_diversity()]).
#' @param n Number of sampled gene copies (haploid: the sample count). A mean over sites
#'   with unequal missingness is fine; it is rounded for the harmonic sums.
#' @return Tajima's D, or `NA` when there are too few segregating sites or samples.
#' @references Tajima, F. (1989) Statistical method for testing the neutral mutation
#'   hypothesis by DNA polymorphism. \emph{Genetics} 123, 585-595.
#'   \doi{10.1093/genetics/123.3.585}
#' @export
tajima_d <- function(h, n) {
  h <- h[is.finite(h)]
  S <- sum(h > 0)
  n <- round(n)
  if (n < 4 || S < 1) return(NA_real_)
  i <- seq_len(n - 1)
  a1 <- sum(1 / i)
  a2 <- sum(1 / i^2)
  b1 <- (n + 1) / (3 * (n - 1))
  b2 <- 2 * (n^2 + n + 3) / (9 * n * (n - 1))
  c1 <- b1 - 1 / a1
  c2 <- b2 - (n + 2) / (a1 * n) + a2 / a1^2
  e1 <- c1 / a1
  e2 <- c2 / (a1^2 + a2)
  denom <- sqrt(e1 * S + e2 * S * (S - 1))
  if (!is.finite(denom) || denom == 0) return(NA_real_)
  (sum(h) - S / a1) / denom
}

#' p-value for a Tajima's D
#'
#' Two-sided significance against the standard neutral model, either by Tajima's own beta
#' approximation (the default, and what \pkg{pegas} reports as `Pval.beta`) or by treating
#' D as standard normal (`Pval.normal`). Both agree with \pkg{pegas} to floating point.
#'
#' **Read it as conservative, and low-powered.** The variance term behind D assumes *no
#' recombination*, which is the most variance the statistic can have; a heavily
#' recombining organism like *P. falciparum* therefore under-rejects. Measured on a real
#' 249-sample cohort, per-gene D reached `p < 0.05` for only 1-2.5% of genes -- less than
#' chance would give -- while 72% of genes had negative D. Per gene there are usually only
#' a handful of segregating sites, which leaves almost no power.
#'
#' What the test asks is also not quite what a scan wants. It tests the standard neutral
#' model as a whole, so a population that has expanded gives negative D genome-wide and a
#' small p-value says "not a constant-size neutral population", not "this gene is under
#' selection". For shortlisting loci, rank on `tajima_d` itself, or use the
#' `tajima_percentile` that [pop_diversity()] reports, and keep the p-value as context.
#'
#' @param D Tajima's D.
#' @param n Number of sampled gene copies.
#' @param S Number of segregating sites.
#' @param method `"beta"` (Tajima's approximation over D's attainable range) or `"normal"`.
#' @return A two-sided p-value, or `NA` when D is undefined.
#' @references Tajima, F. (1989) Statistical method for testing the neutral mutation
#'   hypothesis by DNA polymorphism. \emph{Genetics} 123, 585-595.
#'   \doi{10.1093/genetics/123.3.585}
#' @examples
#' tajima_d_pvalue(-2.1, n = 40, S = 25)
#' @export
tajima_d_pvalue <- function(D, n, S, method = c("beta", "normal")) {
  method <- match.arg(method)
  n <- round(n)
  if (!is.finite(D) || !is.finite(n) || !is.finite(S) || n < 4 || S < 1) return(NA_real_)
  if (identical(method, "normal")) return(2 * stats::pnorm(-abs(D)))
  i <- seq_len(n - 1)
  a1 <- sum(1 / i); a2 <- sum(1 / i^2)
  b1 <- (n + 1) / (3 * (n - 1)); b2 <- 2 * (n^2 + n + 3) / (9 * n * (n - 1))
  c1 <- b1 - 1 / a1
  c2 <- b2 - (n + 2) / (a1 * n) + a2 / a1^2
  e2 <- c2 / (a1^2 + a2)
  # D is bounded; Tajima fits a beta over [Dmin, Dmax] with mean 0 and variance 1
  d_min <- (2 / n - 1 / a1) / sqrt(e2)
  d_max <- ((n + 1) / (2 * n) - 1 / a1) / sqrt(e2)
  if (!is.finite(d_min) || !is.finite(d_max) || d_max <= d_min) return(NA_real_)
  shape2 <- -(1 + d_min * d_max) * d_max / (d_max - d_min)
  shape1 <- (1 + d_min * d_max) * d_min / (d_max - d_min)
  p <- stats::pbeta((D - d_min) / (d_max - d_min), shape1, shape2)
  if (!is.finite(p)) return(NA_real_)
  2 * min(p, 1 - p)
}

# Nei's unbiased haplotype diversity plus the multilocus-genotype summaries poppr reports.
#
# A haplotype needs every one of its sites called, and something has to give when calls
# are missing: either drop the samples with a gap or drop the sites. Dropping samples is
# hopeless at scale -- at 0.05% missing over 27,000 SNPs a sample survives with
# probability ~1e-6, so a whole cohort collapses to a handful. Sites are the plentiful
# axis, so sites are what gets dropped, keeping every sample in the comparison. (This is
# also what published Pf analyses do: pick the SNPs with no missingness in any sample.)
.haplotype_stats <- function(H, max_missing = 1) {
  out <- list(n_hap_snps = 0L, n_hap_samples = 0L, n_hap = NA_integer_,
              hap_div = NA_real_, shannon_h = NA_real_, simpson_lambda = NA_real_,
              evenness = NA_real_)
  if (!ncol(H) || !nrow(H)) return(out)
  H <- H[, colSums(is.na(H)) == 0, drop = FALSE]
  n <- nrow(H)
  out$n_hap_snps <- ncol(H)
  out$n_hap_samples <- n
  if (n < 2 || !ncol(H)) return(out)
  keys <- apply(H, 1, paste0, collapse = "")
  tab <- table(keys)
  f <- as.numeric(tab) / n
  sum_sq <- sum(f^2)
  shannon <- -sum(f * log(f))
  G <- 1 / sum_sq                                  # Stoddart & Taylor's G
  out$n_hap <- length(tab)
  out$hap_div <- (1 - sum_sq) * n / (n - 1)        # Nei (1987) eq. 8.4
  out$shannon_h <- shannon
  out$simpson_lambda <- 1 - sum_sq                 # Simpson's index of diversity
  out$evenness <- if (exp(shannon) > 1) (G - 1) / (exp(shannon) - 1) else NA_real_
  out
}

# Accessible base pairs in [start, end) on `chr`, given a BED-like data frame.
.accessible_bp <- function(chr, start, end, accessible) {
  if (is.null(accessible)) return(as.numeric(end - start))
  a <- accessible[normalise_chr(accessible$chrom) == normalise_chr(chr), , drop = FALSE]
  if (!nrow(a)) return(0)
  lo <- pmax(a$start, start)
  hi <- pmin(a$end, end)
  sum(pmax(0, hi - lo))
}

# Normalise an accessible-region argument to chrom/start/end.
.as_regions <- function(x) {
  if (is.null(x)) return(NULL)
  df <- as.data.frame(x)
  if (!"chrom" %in% names(df) && "chr" %in% names(df)) df$chrom <- df$chr
  if (!"chrom" %in% names(df) && "Pf3D7_chrom" %in% names(df)) df$chrom <- df$Pf3D7_chrom
  need <- c("chrom", "start", "end")
  if (!all(need %in% names(df)))
    stop("`accessible` needs chrom, start and end columns", call. = FALSE)
  df[need]
}

# One row of statistics for a set of SNP columns and a set of samples.
.diversity_row <- function(H, cols, L, min_snps, max_missing) {
  Hs <- H[, cols, drop = FALSE]
  n_samples <- nrow(Hs)
  out <- list(n_samples = n_samples, n_snps = length(cols), n_sites = L,
              seg_sites = NA_integer_, he = NA_real_, pi = NA_real_,
              theta_w = NA_real_, tajima_d = NA_real_, tajima_p = NA_real_,
              n_taj_snps = 0L, n_taj_samples = NA_real_)
  out <- c(out, .haplotype_stats(Hs, max_missing))
  if (!length(cols) || n_samples < 2) return(out)

  f <- .snp_freqs(Hs)
  h <- .site_het(f$p, f$n)
  usable <- f$n >= 2 & is.finite(h)
  out$seg_sites <- as.integer(sum(f$p > 0 & f$p < 1, na.rm = TRUE))
  out$he <- if (any(usable)) mean(h[usable]) else NA_real_
  out$pi <- if (L > 0 && any(usable)) sum(h[usable]) / L else NA_real_

  # Tajima's D needs one sample count for its variance, but the sites have different
  # numbers of calls. Rather than throw away every sample with a gap, keep all the sites
  # and take the mean number of calls across them as n -- so missingness costs precision
  # in the variance term instead of costing most of the cohort.
  if (sum(usable) >= min_snps) {
    nt <- mean(f$n[usable])
    if (nt >= 4) {
      St <- sum(f$p[usable] > 0 & f$p[usable] < 1)
      a1 <- sum(1 / seq_len(floor(nt) - 1))
      out$theta_w <- if (L > 0) (St / a1) / L else NA_real_
      out$tajima_d <- tajima_d(h[usable], nt)
      out$tajima_p <- tajima_d_pvalue(out$tajima_d, nt, St)
      out$n_taj_snps <- sum(usable)
      out$n_taj_samples <- nt
    }
  }
  out
}

#' Within-population genetic diversity
#'
#' Nucleotide diversity, expected heterozygosity, Watterson's theta, Tajima's D and
#' haplotype diversity, computed for each metadata group and reported genome-wide, per
#' gene, or in sliding windows.
#'
#' **`pi` is per accessible base pair** -- the sum of per-site heterozygosity over the
#' number of callable sites (`n_sites`), so windows of different SNP density stay
#' comparable. Pass `accessible` to say which bases were callable; without it every base
#' of the unit is assumed callable, which inflates nothing but does assume your VCF was
#' called across the whole span. `he` is the same per-site quantity averaged over the
#' **SNPs** instead, which is what papers usually label expected heterozygosity (Hs/He);
#' it is not nucleotide diversity and the two are never comparable.
#'
#' The parasite is haploid: `het = "missing"` (the default) treats a heterozygous call as
#' a mixed infection and drops it at that site, so an allele count is a count of samples.
#' `het = "dosage"` instead splits the call evenly between the two alleles.
#'
#' Tajima's D needs a single sample count, so it is computed on a complete-case block --
#' SNPs missing in more than `max_missing` of the group are dropped, then any sample
#' still missing a call. `n_taj_snps` / `n_taj_samples` report what survived, and the
#' statistic is `NA` below `min_snps` SNPs or four samples.
#'
#' @param x A [PopStructure], or a genotype matrix (samples x SNPs, alt dosage 0/1/2,
#'   `NA` for missing) with `chr:pos` column names.
#' @param group Metadata column naming the grouping (for a `PopStructure`) or a vector
#'   aligned to the matrix rows. `NULL` treats every sample as one population.
#' @param by `"genome"` (default), `"gene"`, or `"window"`.
#' @param genes Gene table for `by = "gene"` (columns `name`, `chr`, `start`, `end`);
#'   defaults to [PF3D7_GENES]. Coordinates are 0-based half-open CDS spans.
#' @param window,step Window size and step in bp for `by = "window"`. `step` defaults to
#'   `window`, giving abutting windows; set it smaller for a **sliding** scan --
#'   `window = 5000, step = 2500` steps a 5 kb window along in 2.5 kb hops, so consecutive
#'   windows share half their span. Sliding is the usual way to scan Tajima's D: a fixed
#'   grid can split a signal across two windows and dilute it in both, and overlapping
#'   windows also make the track read more smoothly. The cost is that neighbouring windows
#'   are no longer independent -- do not count them as separate findings (see
#'   [selection_peaks()]).
#' @param accessible Callable regions as a data frame with `chrom`, `start`, `end` --
#'   e.g. [PF3D7_CORE_REGIONS]. Sets the denominator of `pi` and `theta_w`.
#' @param het How to read a heterozygous call: `"missing"` or `"dosage"`.
#' @param min_snps Fewest SNPs for Tajima's D (default `r DIVERSITY_MIN_SNPS`).
#' @param max_missing Per-SNP missingness allowed into the Tajima's D block.
#' @param min_samples Skip a group with fewer samples than this.
#' @param genotype Genotype matrix to use instead of a `PopStructure`'s stored one --
#'   pass the **full, unpruned** set, since LD-pruning removes the very sites diversity
#'   is measured over.
#' @param meta When `x` is a matrix, a data frame with `sample` plus `group`.
#' @return A tibble with one row per group x unit: `group`, the unit's identity, then
#'   `n_samples`, `n_snps`, `n_sites`, `seg_sites`, `he`, `pi`, `theta_w`, `tajima_d`,
#'   `n_hap`, `hap_div`, `shannon_h`, `simpson_lambda`, `evenness`. `tajima_p` is the
#'   two-sided test against the standard neutral model and `tajima_percentile` places each
#'   unit's D within its own group's distribution -- read [tajima_d_pvalue()] before using
#'   the former, which is conservative under recombination.
#' @references
#' Nei, M. (1987) \emph{Molecular Evolutionary Genetics}. Columbia University Press.
#'
#' Watterson, G. A. (1975) On the number of segregating sites in genetical models without
#' recombination. \emph{Theoretical Population Biology} 7, 256-276.
#' \doi{10.1016/0040-5809(75)90020-9}
#'
#' Korunes, K. L. & Samuk, K. (2021) pixy: Unbiased estimation of nucleotide diversity and
#' divergence in the presence of missing data. \emph{Molecular Ecology Resources} 21,
#' 1359-1368. \doi{10.1111/1755-0998.13326}
#' @seealso [pop_diff()] for between-group differentiation, [ld_index()] and
#'   [read_ld_decay()] for linkage disequilibrium, [run_ihs()] for selection.
#' @examples
#' ps <- example_pop_structure(umap = FALSE)
#' pop_diversity(ps, group = "country")
#' @export
pop_diversity <- function(x, group = NULL, by = c("genome", "gene", "window"),
                          genes = NULL, window = 10000, step = NULL,
                          accessible = NULL, het = c("missing", "dosage"),
                          min_snps = DIVERSITY_MIN_SNPS, max_missing = 0.1,
                          min_samples = 4, genotype = NULL, meta = NULL) {
  by <- match.arg(by)
  het <- match.arg(het)
  accessible <- .as_regions(accessible)

  if (inherits(x, "PopStructure")) {
    G <- .geno_for(x, genotype)
    if (is.null(meta)) meta <- x$get_meta()
    G <- G[rownames(G) %in% x$get_samples(), , drop = FALSE]
  } else {
    G <- .coerce_geno(x)
  }
  if (is.null(colnames(G)))
    stop("genotypes need `chr:pos` column names", call. = FALSE)

  grp <- .diversity_group(group, meta, rownames(G), nrow(G))
  levs <- .group_order(grp)

  loci <- .parse_snp_ids(colnames(G))
  H <- .haploid_calls(G, het)

  units <- switch(
    by,
    genome = .genome_units(loci, accessible),
    gene   = .gene_units(if (is.null(genes)) PF3D7_GENES else genes, loci, accessible),
    window = .window_units(loci, window, if (is.null(step)) window else step, accessible))

  rows <- list()
  for (l in levs) {
    idx <- which(grp == l)
    if (length(idx) < min_samples) next
    Hg <- H[idx, , drop = FALSE]
    for (u in seq_len(nrow(units))) {
      cols <- units$cols[[u]]
      stats <- .diversity_row(Hg, cols, units$n_sites[u], min_snps, max_missing)
      rows[[length(rows) + 1L]] <- c(
        list(group = l, chr = units$chr[u], start = units$start[u], end = units$end[u],
             unit = units$unit[u]), stats)
    }
  }
  if (!length(rows)) {
    warning("no group had at least `min_samples` samples", call. = FALSE)
    return(.empty_diversity())
  }
  out <- do.call(rbind, lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
  out$group <- factor(out$group, levels = levs)
  # where each unit's D sits in this group's own distribution of D. With the parametric
  # test conservative under recombination (see `tajima_d_pvalue`), the empirical rank is
  # usually the more useful way to shortlist loci.
  out$tajima_percentile <- NA_real_
  for (l in levs) {
    k <- which(out$group == l & is.finite(out$tajima_d))
    if (length(k) > 2)
      out$tajima_percentile[k] <- stats::ecdf(out$tajima_d[k])(out$tajima_d[k])
  }
  rownames(out) <- NULL
  tibble::as_tibble(out)
}

.empty_diversity <- function() {
  tibble::as_tibble(data.frame(
    group = character(), chr = character(), start = numeric(), end = numeric(),
    unit = character(), n_samples = integer(), n_snps = integer(), n_sites = numeric(),
    seg_sites = integer(), he = numeric(), pi = numeric(), theta_w = numeric(),
    tajima_d = numeric(), n_taj_snps = integer(), n_taj_samples = integer(),
    n_hap_samples = integer(), n_hap = integer(), hap_div = numeric(),
    shannon_h = numeric(), simpson_lambda = numeric(), evenness = numeric(),
    stringsAsFactors = FALSE))
}

.diversity_group <- function(group, meta, samples, n) {
  if (is.null(group)) return(factor(rep("all", n)))
  if (length(group) == 1L && is.character(group)) {
    if (is.null(meta) || !group %in% names(meta))
      stop(sprintf("`%s` is not a column of the metadata", group), call. = FALSE)
    return(.group_factor(meta, group, samples))
  }
  if (length(group) != n)
    stop("`group` must be one metadata column name or one value per sample", call. = FALSE)
  .as_group_factor(group)
}

# "chr:pos" column names -> chr / pos, 0-based
.parse_snp_ids <- function(ids) {
  parts <- strsplit(ids, ":", fixed = TRUE)
  bad <- lengths(parts) < 2
  if (any(bad))
    stop("genotype column names must look like `chr:pos`", call. = FALSE)
  data.frame(chr = vapply(parts, `[`, character(1), 1),
             pos = as.numeric(vapply(parts, `[`, character(1), 2)),
             idx = seq_along(ids), stringsAsFactors = FALSE)
}

.genome_units <- function(loci, accessible) {
  n_sites <- if (is.null(accessible)) {
    sum(PF3D7_CORE_CHROM_LENGTHS_BP)
  } else {
    sum(pmax(0, accessible$end - accessible$start))
  }
  out <- data.frame(chr = NA_character_, start = NA_real_, end = NA_real_,
                    unit = "genome", n_sites = n_sites, stringsAsFactors = FALSE)
  out$cols <- I(list(loci$idx))
  out
}

.gene_units <- function(genes, loci, accessible) {
  g <- as.data.frame(genes)
  if (!"chr" %in% names(g)) {
    for (alt in c("chrom", "Pf3D7_chrom")) if (alt %in% names(g)) { g$chr <- g[[alt]]; break }
  }
  need <- c("name", "chr", "start", "end")
  if (!all(need %in% names(g)))
    stop("`genes` needs name, chr (or chrom), start and end columns", call. = FALSE)
  key <- normalise_chr(loci$chr)
  cols <- vector("list", nrow(g))
  n_sites <- numeric(nrow(g))
  for (i in seq_len(nrow(g))) {
    hit <- key == normalise_chr(g$chr[i]) & loci$pos >= g$start[i] & loci$pos < g$end[i]
    cols[[i]] <- loci$idx[hit]
    n_sites[i] <- .accessible_bp(g$chr[i], g$start[i], g$end[i], accessible)
  }
  out <- data.frame(chr = as.character(g$chr), start = as.numeric(g$start),
                    end = as.numeric(g$end), unit = as.character(g$name),
                    n_sites = n_sites, stringsAsFactors = FALSE)
  out$cols <- I(cols)
  out
}

.window_units <- function(loci, window, step, accessible) {
  pieces <- list()
  for (chr in unique(loci$chr)) {
    sub <- loci[loci$chr == chr, , drop = FALSE]
    starts <- seq(0, max(sub$pos), by = step)
    cols <- lapply(starts, function(s) sub$idx[sub$pos >= s & sub$pos < s + window])
    keep <- lengths(cols) > 0
    if (!any(keep)) next
    starts <- starts[keep]
    cols <- cols[keep]
    df <- data.frame(chr = chr, start = starts, end = starts + window,
                     unit = sprintf("%s:%d-%d", chr, starts, starts + window),
                     n_sites = vapply(starts, function(s)
                       .accessible_bp(chr, s, s + window, accessible), numeric(1)),
                     stringsAsFactors = FALSE)
    df$cols <- I(cols)
    pieces[[length(pieces) + 1L]] <- df
  }
  if (!length(pieces)) stop("no window contained a SNP", call. = FALSE)
  do.call(rbind, pieces)
}
