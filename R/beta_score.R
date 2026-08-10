# Beta: the allele-frequency signature of long-term balancing selection.
#
# Under balancing selection a variant is held at intermediate frequency for a long time,
# and the mutations that arise on either side of it are dragged towards *the same*
# frequency. Beta scores that clustering: around each core SNP it compares an estimate of
# theta weighted by how close each neighbour's folded frequency sits to the core's
# (theta_beta) against Watterson's theta from the plain count of segregating sites. A
# large positive value means the neighbourhood's frequencies are piled up around the
# core's, which drift alone rarely produces.
#
# This is the natural counterpart to iHS for Plasmodium: sweeps at drug-resistance loci
# are directional and show up in extended haplotype homozygosity, while the antigens are
# under balancing selection and show up here.
#
# Follows Siewert & Voight (2017) and matches the folded Beta1 of their BetaScan
# reference implementation, including its defaults (a 1 kb window, p = 2).

#: Total width in bp of the window of neighbouring SNPs around each core SNP.
BETA_WINDOW <- 1000L
#: Sharpness of the frequency-similarity weighting; higher is more peaked on the core.
BETA_P <- 2

# Similarity of neighbour frequencies `freq` to core frequency `x`, folded.
.beta_d <- function(freq, x, p) {
  xf <- min(x, 1 - x)
  f <- pmin(freq, 1 - freq)
  maxdiff <- max(xf, 0.5 - xf)
  ((maxdiff - abs(xf - f)) / maxdiff)^p
}

# Folded Beta1 for one core SNP: theta_beta - theta_Watterson over its window.
.beta_folded <- function(win_freq, core_freq, n, p) {
  if (!length(win_freq) || n < 4) return(NA_real_)
  i <- seq_len(n - 1)
  a1 <- sum(1 / i)
  theta_w <- length(win_freq) / a1
  num <- sum(.beta_d(win_freq, core_freq, p))
  den <- sum((1 / i) * .beta_d(i / n, core_freq, p))
  if (!is.finite(den) || den == 0) return(NA_real_)
  num / den - theta_w
}

#' Beta: balancing selection from clustered allele frequencies
#'
#' For each core SNP, compares an estimate of theta weighted towards neighbours whose
#' folded allele frequency matches the core's against Watterson's theta from the same
#' window. Large positive values mark neighbourhoods where frequencies are piled up
#' around an intermediate-frequency variant -- the footprint of long-term balancing
#' selection, and in *P. falciparum* typically an antigen rather than a drug target.
#'
#' @param x A [PopStructure] or a genotype matrix (samples x SNPs, alt dosage, `chr:pos`
#'   column names).
#' @param group Metadata column name, or a vector aligned to the rows; `NULL` pools
#'   every sample.
#' @param window Total window width in bp; neighbours are taken within `window / 2` on
#'   each side (default `r BETA_WINDOW`).
#' @param p Sharpness of the frequency-similarity weighting (default `r BETA_P`).
#' @param min_freq Skip core SNPs whose folded frequency is at or below this.
#' @param min_window_snps Skip core SNPs with fewer neighbours than this in the window.
#' @param het How a heterozygous call is read; see [pop_diversity()].
#' @param min_samples Skip groups smaller than this.
#' @param meta,genotype As in [pop_diversity()].
#' @return A tibble with `group`, `chr`, `pos`, `snp_id`, `freq`, `n_called`,
#'   `n_window_snps` and `beta`, one row per core SNP per group.
#'
#'   **Reading `beta`.** This is Siewert & Voight's folded Beta1: how much more the SNPs
#'   around this one share its allele frequency than drift alone would produce. Long-term
#'   balancing selection holds a variant at intermediate frequency for many generations, and
#'   its neighbours hitch-hike to *similar* frequencies -- so a cluster of matched
#'   frequencies is the signature, which is what beta measures.
#'
#'   The scale is not a z-score and has no p-value here: it depends on the window, the SNP
#'   density and the sample size, so a value meaningful in one scan is not comparable to
#'   another. **Rank rather than threshold.** Larger is more evidence of balancing
#'   selection; around zero is what neutrality gives; negative means the neighbourhood is
#'   *less* frequency-matched than expected, which is a directional-sweep pattern, not
#'   balancing selection -- for that read [run_ihs()] instead.
#'
#'   In practice take the top tail --
#'   `selection_peaks(b, criterion = "top", top = 0.01, metric = "beta")` -- and annotate it
#'   with [annotate_snps()]. In *P. falciparum* the expected hits are the antigens under
#'   long-term frequency-dependent selection from host immunity -- on a real 249-sample
#'   cohort the top 1% is led by *pfama1* and *pfdblmsp*, which is the positive control to
#'   look for; a scan that surfaces none of them is a reason to check the input before
#'   believing anything else in it. Note the classic examples *var* and *rifin* will
#'   **not** appear: they are subtelomeric, so a core-genome filter removes them before the
#'   scan ever sees them.
#' @references Siewert, K. M. & Voight, B. F. (2017) Detecting long-term balancing
#'   selection using allele frequency correlation. \emph{Molecular Biology and Evolution}
#'   34, 2996-3005. \doi{10.1093/molbev/msx209}
#' @seealso [beta_genes()] to summarise per gene, [run_ihs()] for directional selection.
#' @examples
#' ps <- example_pop_structure(umap = FALSE)
#' beta_score(ps, group = "country", window = 200000, min_window_snps = 1)
#' @export
beta_score <- function(x, group = NULL, window = BETA_WINDOW, p = BETA_P,
                       min_freq = 0, min_window_snps = 2,
                       het = c("missing", "dosage"), min_samples = 4,
                       meta = NULL, genotype = NULL) {
  het <- match.arg(het)
  prep <- .ld_prepare(x, group, meta, genotype, het)
  half <- window / 2
  rows <- list()

  for (l in prep$levs) {
    idx <- which(prep$grp == l)
    if (length(idx) < min_samples) next
    Hg <- prep$H[idx, , drop = FALSE]
    f <- .snp_freqs(Hg)
    for (chr in unique(prep$loci$chr)) {
      s <- prep$loci[prep$loci$chr == chr, , drop = FALSE]
      s <- s[order(s$pos), , drop = FALSE]
      freq <- f$p[s$idx]
      n <- f$n[s$idx]
      lo <- findInterval(s$pos - half, s$pos, left.open = TRUE) + 1L
      hi <- findInterval(s$pos + half, s$pos)
      for (i in seq_len(nrow(s))) {
        fi <- freq[i]
        if (is.na(fi) || fi <= min_freq || fi >= 1 - min_freq || n[i] < 4) next
        j <- setdiff(seq.int(lo[i], hi[i]), i)          # the core is not its own neighbour
        wf <- freq[j]
        wf <- wf[is.finite(wf) & wf > 0 & wf < 1]
        if (length(wf) < min_window_snps) next
        rows[[length(rows) + 1L]] <- data.frame(
          group = l, chr = chr, pos = s$pos[i],
          snp_id = paste0(chr, ":", format(s$pos[i], scientific = FALSE, trim = TRUE)),
          freq = fi, n_called = n[i], n_window_snps = length(wf),
          beta = .beta_folded(wf, fi, n[i], p), stringsAsFactors = FALSE)
      }
    }
  }
  if (!length(rows)) {
    warning("no core SNP had enough neighbours in its window", call. = FALSE)
    return(tibble::tibble())
  }
  out <- do.call(rbind, rows)
  out$group <- factor(out$group, levels = prep$levs)
  tibble::as_tibble(out)
}

#' Summarise beta scores per gene
#'
#' Mean and maximum beta over the SNPs inside each gene, giving the per-gene view used to
#' rank candidates for balancing selection.
#'
#' @param b The tibble from [beta_score()].
#' @param genes Gene table (`name`, `chr` or `chrom`, `start`, `end`); defaults to
#'   [PF3D7_GENES]. Coordinates are 0-based half-open.
#' @param min_snps Genes with fewer scored SNPs than this are dropped.
#' @return A tibble with `group`, `gene`, `chr`, `start`, `end`, `n_snps`, `beta_mean`,
#'   `beta_max`, sorted by `beta_mean` within each group.
#' @examples
#' ps <- example_pop_structure(umap = FALSE)
#' b <- beta_score(ps, group = "country", window = 200000, min_window_snps = 1)
#' beta_genes(b, genes = PF_EXAMPLE_DRUG_GENES, min_snps = 1)
#' @export
beta_genes <- function(b, genes = NULL, min_snps = 3) {
  if (!nrow(b)) return(tibble::tibble())
  g <- as.data.frame(if (is.null(genes)) PF3D7_GENES else genes)
  if (!"chr" %in% names(g)) {
    for (alt in c("chrom", "Pf3D7_chrom")) if (alt %in% names(g)) { g$chr <- g[[alt]]; break }
  }
  need <- c("name", "chr", "start", "end")
  if (!all(need %in% names(g)))
    stop("`genes` needs name, chr (or chrom), start and end columns", call. = FALSE)

  key <- normalise_chr(b$chr)
  rows <- list()
  for (i in seq_len(nrow(g))) {
    hit <- key == normalise_chr(g$chr[i]) & b$pos >= g$start[i] & b$pos < g$end[i]
    if (!any(hit)) next
    sub <- b[hit, , drop = FALSE]
    for (l in unique(sub$group)) {
      v <- sub$beta[sub$group == l]
      v <- v[is.finite(v)]
      if (length(v) < min_snps) next
      rows[[length(rows) + 1L]] <- data.frame(
        group = l, gene = as.character(g$name[i]), chr = as.character(g$chr[i]),
        start = as.numeric(g$start[i]), end = as.numeric(g$end[i]),
        n_snps = length(v), beta_mean = mean(v), beta_max = max(v),
        stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(tibble::tibble())
  out <- do.call(rbind, rows)
  out$group <- factor(out$group, levels = levels(b$group))
  out <- out[order(out$group, -out$beta_mean), , drop = FALSE]
  rownames(out) <- NULL
  tibble::as_tibble(out)
}
