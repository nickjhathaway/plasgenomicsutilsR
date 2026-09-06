# Recent directional selection from extended haplotype homozygosity: iHS within a
# population, Rsb and XP-EHH between two.
#
# A recent sweep leaves a long haplotype at high frequency, because the favoured allele
# rose faster than recombination could break up the background it arose on. iHS scores
# that asymmetry at each SNP; Rsb and XP-EHH ask the same question across two populations,
# which finds sweeps that are complete (and so invisible to iHS) in one of them.
#
# The statistics come from \pkg{rehh}, which needs phased haplotypes with no missing
# calls. Neither holds for a *P. falciparum* VCF: infections can be polyclonal, and calls
# are missing. [parasite_haplotypes()] is the bridge -- gate to monoclonal infections,
# resolve or drop the remaining mixed calls, and impute what is left -- and it reports
# exactly what it removed, because how the haplotypes were made determines what the scan
# can honestly claim.

#: Fws at or above which an infection is treated as monoclonal.
IHS_MIN_FWS <- 0.95
#: Default minor-allele frequency floor for SNPs entering a haplotype scan.
IHS_MIN_MAF <- 0.03

#' Build phased haplotypes for a haplotype-homozygosity scan
#'
#' Turns a genotype matrix into the complete, unambiguous 0/1 haplotypes that
#' [run_ihs()] and its cross-population relatives need.
#'
#' The steps, in order, each reported in the result:
#'
#' 1. **Monoclonal gate.** With `fws` supplied, samples below `min_fws` are dropped: a
#'    polyclonal infection is a mixture of haplotypes, not one haplotype.
#' 2. **Mixed calls.** Whatever heterozygous calls remain are either resolved by drawing
#'    the allele at the population frequency (`het = "sample"`, the usual choice for
#'    *Plasmodium*) or set to missing.
#' 3. **Filtering.** SNPs below `maf` or above `max_snp_missing`, then samples above
#'    `max_sample_missing`.
#' 4. **Imputation.** Remaining gaps are filled by drawing at the SNP's allele frequency,
#'    so the matrix is complete.
#'
#' Both the allele draw and the imputation are random; `seed` makes a run reproducible.
#' Because the result depends on that draw, a signal worth reporting should survive
#' repeating the whole thing under a different seed.
#'
#' @param x A [PopStructure] or a genotype matrix (samples x SNPs, alt dosage 0/1/2,
#'   `NA` missing, `chr:pos` column names).
#' @param samples Restrict to these samples before anything else.
#' @param fws Per-sample Fws: a named numeric vector, or a data frame with `sample` and
#'   `fws` columns (e.g. read from `plasgenomicsutils calculate_fws`). `NULL` skips the
#'   monoclonal gate.
#' @param min_fws Fws floor for the gate (default `r IHS_MIN_FWS`).
#' @param het `"sample"` draws the allele at the population frequency; `"missing"` leaves
#'   the call to imputation.
#' @param maf Minor-allele frequency floor (default `r IHS_MIN_MAF`).
#' @param max_snp_missing,max_sample_missing Missingness ceilings for SNPs and samples.
#' @param alleles How the input matrix codes its calls. `"dosage"` (default) is ALT dosage
#'   0/1/2, which can only ever describe two alleles -- 2 is two copies of the one ALT, not
#'   a second ALT. `"index"` is one allele index per haplotype (0 = reference, then 1, 2,
#'   ...), which is what a marker with more than two alleles needs. An index matrix is
#'   already haploid, so nothing is made haploid and no heterozygote is resolved (`het` does
#'   not apply); the MAF floor becomes the share of everything that is not the commonest
#'   allele, and imputation draws from the alleles seen at that SNP. Only the EHH curves
#'   read more than two alleles -- [run_ihs()], [run_rsb()] and [run_xpehh()] rest on a
#'   ratio between exactly two, and \pkg{rehh}'s `scan_hh()` silently keeps the two
#'   commonest.
#' @param impute Fill remaining gaps by drawing at the SNP's allele frequency. `FALSE`
#'   instead drops every SNP that still has a gap.
#' @param seed Random seed.
#' @param meta,genotype As in [pop_diversity()].
#' @return A `parasite_haplotypes` object: `hap` (samples x SNPs, 0/1), `map`
#'   (`chr`, `pos`, `snp_id`), `meta`, and a `filtering` record of what was dropped.
#' @seealso [run_ihs()], [run_rsb()], [run_xpehh()]
#' @examples
#' ps <- example_pop_structure(umap = FALSE)
#' parasite_haplotypes(ps, maf = 0.05)
#' @export
parasite_haplotypes <- function(x, samples = NULL, fws = NULL, min_fws = IHS_MIN_FWS,
                                het = c("sample", "missing"), maf = IHS_MIN_MAF,
                                max_snp_missing = 0.1, max_sample_missing = 0.2,
                                impute = TRUE, seed = 42, meta = NULL, genotype = NULL,
                                alleles = c("dosage", "index")) {
  meta <- .normalise_meta(meta)
  het <- match.arg(het)
  alleles <- match.arg(alleles)
  if (inherits(x, "PopStructure")) {
    G <- .geno_for(x, genotype)
    if (is.null(meta)) meta <- x$get_meta()
  } else {
    G <- .coerce_geno(x)
  }
  if (is.null(colnames(G)))
    stop("genotypes need `chr:pos` column names", call. = FALSE)
  set.seed(seed)
  rec <- list(n_samples_in = nrow(G), n_snps_in = ncol(G))

  if (!is.null(samples)) {
    keep <- rownames(G) %in% samples
    if (!any(keep)) stop("none of `samples` are in the genotype matrix", call. = FALSE)
    G <- G[keep, , drop = FALSE]
  }
  rec$n_dropped_not_requested <- rec$n_samples_in - nrow(G)

  n_before <- nrow(G)
  if (!is.null(fws)) {
    f <- .fws_vector(fws)
    unknown <- setdiff(rownames(G), names(f))
    if (length(unknown))
      stop(sprintf("no Fws for %d sample(s), e.g. %s", length(unknown),
                   paste(utils::head(unknown, 3), collapse = ", ")), call. = FALSE)
    G <- G[f[rownames(G)] >= min_fws, , drop = FALSE]
    if (!nrow(G))
      stop(sprintf("no sample has Fws >= %.2f", min_fws), call. = FALSE)
  }
  rec$n_dropped_polyclonal <- n_before - nrow(G)

  if (identical(alleles, "index")) {
    # Already one allele index per haplotype, so there is nothing to make haploid and no
    # heterozygote to resolve: an index says which allele this haplotype carries, and a
    # monoclonal infection carries one. This is the route for a marker with more than two
    # alleles, which a dosage matrix cannot express -- dosage 2 is two copies of the one
    # ALT, not a second one.
    seen <- G[!is.na(G)]
    if (length(seen) && (any(seen < 0) || any(seen != round(seen))))
      stop("`alleles = \"index\"` needs non-negative whole-number allele indices",
           call. = FALSE)
    H <- G
    rec$n_het_calls <- 0L
  } else {
    H <- .haploid_calls(G, "missing")
    rec$n_het_calls <- sum(G == 1, na.rm = TRUE)
    if (identical(het, "sample") && rec$n_het_calls > 0) {
      p <- .snp_freqs(H)$p
      hit <- which(!is.na(G) & G == 1, arr.ind = TRUE)
      pr <- p[hit[, "col"]]
      pr[is.na(pr)] <- 0.5
      H[hit] <- stats::rbinom(nrow(hit), 1, pr)
    }
  }

  n_snp <- ncol(H)
  keep_snp <- colSums(is.na(H)) / nrow(H) <= max_snp_missing
  H <- H[, keep_snp, drop = FALSE]
  rec$n_dropped_snp_missing <- n_snp - ncol(H)

  n_samp <- nrow(H)
  keep_samp <- rowSums(is.na(H)) / max(1L, ncol(H)) <= max_sample_missing
  H <- H[keep_samp, , drop = FALSE]
  rec$n_dropped_sample_missing <- n_samp - nrow(H)
  if (nrow(H) < 4) stop("fewer than 4 samples survived filtering", call. = FALSE)

  n_snp <- ncol(H)
  # With two alleles the minor-allele frequency is min(p, 1 - p); with more it is whatever
  # is left over from the commonest one, which reduces to the same thing when there are two.
  minor_af <- if (identical(alleles, "index")) .index_minor_af(H) else {
    m <- .snp_freqs(H); pmin(m$p, 1 - m$p)
  }
  keep_maf <- !is.na(minor_af) & minor_af >= maf
  H <- H[, keep_maf, drop = FALSE]
  rec$n_dropped_maf <- n_snp - ncol(H)
  if (!ncol(H)) stop(sprintf("no SNP passed maf >= %.3f", maf), call. = FALSE)

  rec$n_imputed <- sum(is.na(H))
  if (rec$n_imputed > 0) {
    if (impute) {
      if (identical(alleles, "index")) {
        # draw from the alleles actually seen at that SNP -- the k-allele form of filling a
        # gap from the ALT frequency
        for (j in unique(which(is.na(H), arr.ind = TRUE)[, "col"])) {
          obs <- H[!is.na(H[, j]), j]
          gap <- which(is.na(H[, j]))
          H[gap, j] <- if (length(obs)) sample(obs, length(gap), replace = TRUE) else 0L
        }
      } else {
        p <- .snp_freqs(H)$p
        hit <- which(is.na(H), arr.ind = TRUE)
        H[hit] <- stats::rbinom(nrow(hit), 1, p[hit[, "col"]])
      }
    } else {
      n_snp <- ncol(H)
      H <- H[, colSums(is.na(H)) == 0, drop = FALSE]
      rec$n_dropped_incomplete <- n_snp - ncol(H)
      rec$n_imputed <- 0L
    }
  }

  map <- .parse_snp_ids(colnames(H))
  map$snp_id <- colnames(H)
  storage.mode(H) <- "integer"          # rehh requires an integer haplotype matrix
  structure(list(hap = H, map = map[c("chr", "pos", "snp_id")],
                 meta = meta, filtering = rec, seed = seed, het = het, alleles = alleles),
            class = "parasite_haplotypes")
}

# Minor-allele frequency of an allele-index column: everything not the commonest allele.
.index_minor_af <- function(H) {
  vapply(seq_len(ncol(H)), function(j) {
    v <- H[!is.na(H[, j]), j]
    if (!length(v)) return(NA_real_)
    1 - max(tabulate(as.integer(v) + 1L)) / length(v)
  }, numeric(1))
}


.fws_vector <- function(fws) {
  if (is.data.frame(fws)) {
    nm <- names(fws)
    col <- nm[tolower(nm) %in% c("fws", "f_ws")][1]
    if (is.na(col) || !"sample" %in% nm)
      stop("`fws` data frame needs `sample` and `fws` columns", call. = FALSE)
    return(stats::setNames(as.numeric(fws[[col]]), as.character(fws$sample)))
  }
  if (is.null(names(fws))) stop("`fws` needs sample names", call. = FALSE)
  stats::setNames(as.numeric(fws), names(fws))
}

#' @export
print.parasite_haplotypes <- function(x, ...) {
  r <- x$filtering
  cat("<parasite_haplotypes>", nrow(x$hap), "haplotypes x", ncol(x$hap), "SNPs\n")
  cat("  from            :", r$n_samples_in, "samples x", r$n_snps_in, "SNPs\n")
  if (r$n_dropped_polyclonal) cat("  polyclonal      :", r$n_dropped_polyclonal, "dropped\n")
  if (r$n_het_calls) cat("  mixed calls     :", r$n_het_calls,
                         if (identical(x$het, "sample")) "resolved by allele draw"
                         else "set to missing", "\n")
  cat("  SNPs dropped    :", r$n_dropped_snp_missing, "missing,", r$n_dropped_maf, "MAF")
  if (!is.null(r$n_dropped_incomplete)) cat(",", r$n_dropped_incomplete, "incomplete")
  cat("\n")
  cat("  samples dropped :", r$n_dropped_sample_missing, "missing\n")
  cat("  imputed calls   :", r$n_imputed, " seed:", x$seed, "\n")
  # a subset is easy to lose track of, and it changes what every scan off this object means
  if (!is.null(x$subset))
    cat("  subset          :", nrow(x$hap), "of", x$subset$from, "haplotypes",
        if (length(x$subset$by)) paste0("(", paste(x$subset$by, collapse = "; "), ")"), "\n")
  invisible(x)
}

# rehh needs one haplohh object per chromosome; build them in memory.
.haplohh_list <- function(hap, rows) {
  .need_package("rehh", "the haplotype scans")
  map <- hap$map
  out <- list()
  for (chr in unique(map$chr)) {
    k <- which(map$chr == chr)
    k <- k[order(map$pos[k])]
    if (length(k) < 2) next
    h <- hap$hap[rows, k, drop = FALSE]
    # rehh reads a monomorphic column as uninformative; drop them per group
    poly <- apply(h, 2, function(v) length(unique(v)) > 1)
    if (sum(poly) < 2) next
    k <- k[poly]
    h <- hap$hap[rows, k, drop = FALSE]
    storage.mode(h) <- "integer"
    out[[chr]] <- methods::new("haplohh",
                               haplo = h,
                               positions = as.numeric(map$pos[k]),
                               chr.name = as.character(chr))
  }
  out
}

.scan_group <- function(hap, rows, polarized, threads, maxgap = NA, scalegap = NA,
                        discard_at_border = NULL) {
  objs <- .haplohh_list(hap, rows)
  if (!length(objs)) return(NULL)
  discard <- .resolve_border(discard_at_border, maxgap)
  scans <- lapply(objs, function(o)
    rehh::scan_hh(o, polarized = polarized, maxgap = maxgap, scalegap = scalegap,
                  discard_integration_at_border = discard, threads = threads))
  do.call(rbind, scans)
}

# Discarding at the border returns NA for a marker whose EHH had not decayed before the
# integration ran out of data. A few of those at the chromosome ends are the point. Most of
# the genome means the markers are too sparse for EHH to decay within them at all, and the
# scan comes back as a table of NA -- which is worth saying out loud, since the alternative
# is a caller wondering why an argument that sounded conservative emptied their results.
.border_na_frac <- function(scan) {
  if (is.null(scan) || !nrow(scan)) return(NA_real_)
  ihh <- scan[, grep("^IHH_", names(scan)), drop = FALSE]
  mean(!stats::complete.cases(ihh))
}

# Once per scan, not once per group: the same sparse markers produce the same warning in
# every group, and a five-region run should not say it five times.
.warn_border_na <- function(fracs, maxgap, discard_at_border) {
  if (!.resolve_border(discard_at_border, maxgap)) return(invisible(NULL))
  frac <- mean(fracs, na.rm = TRUE)
  if (is.finite(frac) && frac > 0.5)
    warning(sprintf(paste("%.0f%% of markers have no iHH: the integration reached a border",
                          "(a chromosome end, or a gap wider than `maxgap`) before EHH",
                          "decayed. Widen `maxgap`, or set `discard_at_border = FALSE` to",
                          "keep the truncated integrals."), 100 * frac), call. = FALSE)
  invisible(NULL)
}

# Without `maxgap` the only border is the end of a chromosome, and discarding there costs
# the markers nearest the telomeres for nothing. With it, the border is a hole in the data:
# a truncated integral over one describes the hole rather than the haplotypes, so it is
# better returned as NA than as a number.
.resolve_border <- function(discard, maxgap) {
  if (!is.null(discard)) return(isTRUE(discard))
  !is.na(maxgap)
}

.ihs_rows <- function(hap, group, meta, min_samples) {
  meta <- if (is.null(meta)) hap$meta else meta
  grp <- .diversity_group(group, meta, rownames(hap$hap), nrow(hap$hap))
  # A haplotype whose sample is absent from `meta` has no group. Say so rather than
  # letting it quietly shrink every group -- dropping samples from the metadata is a
  # normal thing to do (clones, QC failures) and the haplotypes usually still hold them.
  n_ungrouped <- sum(is.na(grp))
  if (n_ungrouped)
    warning(sprintf("%d haplotype(s) have no `%s` in the metadata and are excluded",
                    n_ungrouped, if (is.character(group) && length(group) == 1) group
                                 else "group"), call. = FALSE)
  levs <- .group_order(grp)
  keep <- levs[vapply(levs, function(l) sum(grp == l, na.rm = TRUE) >= min_samples,
                      logical(1))]
  if (!length(keep))
    stop(sprintf("no group has at least %d haplotypes", min_samples), call. = FALSE)
  stats::setNames(lapply(keep, function(l) which(grp == l)), keep)
}

# rehh standardises the log iHH ratio within bins of the frequency it was computed against.
# For an unpolarized scan that frequency is FREQ_MIN, which cannot exceed 0.5, while rehh
# cuts its bins across [min_maf, 1 - min_maf] -- a grid built for derived-allele
# frequencies, which do span the whole range. Half those bins are therefore empty by
# construction and each one warns, which is why a single bin is the default here.
#
# A single bin has its own cost: the spread of the ratio depends on how common the minor
# allele is (a rarer allele is carried by fewer haplotypes, so its iHH is noisier), and one
# global mean and sd cannot take that out. Bands cut at quantiles of the observed minor
# allele frequency do -- they cover only the range the data occupies, and each holds a
# similar number of markers by construction, so the sparse-bin problem cannot arise.
.band_standardise <- function(raw, freq, n_bands) {
  u <- raw$UNIHS
  br <- unique(stats::quantile(freq, seq(0, 1, length.out = n_bands + 1L), na.rm = TRUE))
  if (length(br) < 2) {
    bands <- factor(rep("all", length(u)))
  } else {
    bands <- cut(freq, br, include.lowest = TRUE)
  }
  n <- tapply(u, bands, function(x) sum(is.finite(x)))
  thin <- names(n)[!is.na(n) & n < 10]
  if (length(thin))
    warning(sprintf("%d minor-allele-frequency band(s) hold fewer than 10 markers (%s); ",
                    length(thin), paste(thin, collapse = ", ")),
            "lower `maf_bands`", call. = FALSE)
  m <- tapply(u, bands, mean, na.rm = TRUE)
  sd <- tapply(u, bands, stats::sd, na.rm = TRUE)
  z <- (u - m[bands]) / sd[bands]
  lp <- -log10(2 * stats::pnorm(-abs(z)))
  # a p-value that underflows to zero would plot as Inf; rehh puts those one above the
  # largest finite value, and the two paths should agree on what the axis means
  if (any(is.finite(lp))) lp[is.infinite(lp)] <- max(lp[is.finite(lp)]) + 1
  data.frame(CHR = raw$CHR, POSITION = raw$POSITION, IHS = as.numeric(z),
             LOGPVALUE = as.numeric(lp), stringsAsFactors = FALSE)
}

# One group's scan -> its standardised iHS, by rehh's own binning or by our frequency bands.
.standardise_ihs <- function(scan, freqbin, min_maf, maf_bands) {
  if (is.null(maf_bands))
    return(rehh::ihh2ihs(scan, freqbin = freqbin, min_maf = min_maf,
                         verbose = FALSE)$ihs)
  # standardize = FALSE returns the raw log ratio, and skips the binning entirely
  raw <- rehh::ihh2ihs(scan, freqbin = 1, min_maf = min_maf, standardize = FALSE,
                       include_freq = TRUE, verbose = FALSE)$ihs
  if (is.null(raw) || !nrow(raw)) return(raw)
  # the frequency the ratio is oriented against: FREQ_MIN unpolarized, FREQ_D polarized --
  # the same column rehh would have binned on, and always the last of the two it returns
  freq <- raw[[grep("^FREQ_", names(raw))[2]]]
  .band_standardise(raw, freq, maf_bands)
}

# Frequency binning only means something when there IS a derived allele to bin by, so an
# unpolarized scan gets a single bin. Warn rather than silently override a deliberate choice.
.resolve_freqbin <- function(freqbin, polarized) {
  if (is.null(freqbin)) return(if (polarized) 0.05 else 1)
  if (!polarized && freqbin < 1) {
    warning("freqbin = ", freqbin, " bins an unpolarized scan by major-allele frequency, ",
            "which is not the derived-allele frequency the standardisation is for. ",
            "freqbin = 1 (the default here) is what rehh recommends. The two can rank ",
            "loci quite differently, so this is not a cosmetic choice.", call. = FALSE)
  }
  freqbin
}

#' Integrated haplotype score (iHS)
#'
#' Scans each group for recent positive directional selection, standardising the
#' integrated EHH ratio within allele-frequency bins so scores are comparable along the
#' genome. Large `abs(ihs)` marks a SNP whose haplotype background is unusually long for
#' its frequency.
#'
#' Without an outgroup there is no ancestral state to polarise by, so `polarized = FALSE`
#' by default and the comparison is major versus minor allele rather than ancestral
#' versus derived. That is the standard treatment for *P. falciparum* and it means the
#' *sign* of `ihs` should not be read as "selection on the derived allele" -- use
#' `abs(ihs)` and `neg_log10_p`.
#'
#' @param hap A [parasite_haplotypes()] object.
#' @param group Metadata column naming the grouping, a vector aligned to the haplotype
#'   rows, or `NULL` to scan every sample as one population.
#' @param meta Metadata (defaults to the one carried by `hap`).
#' @param polarized Treat allele 1 as derived (needs a real ancestral state).
#' @param freqbin Width of the allele-frequency bins iHS is standardised within, or `NULL`
#'   (default) to pick one from `polarized`: **1** (a single bin) when unpolarized, `0.05`
#'   when polarized. The binning exists to control for *derived* allele frequency, and an
#'   unpolarized scan has no ancestral state -- only `FREQ_MAJ`/`FREQ_MIN` -- so major/minor
#'   is not derived/ancestral and binning by it controls nothing. rehh warns about this and
#'   about the resulting sparse bins above 0.5; a single bin silences both because it is the
#'   right answer, not because the warning was noise. Not cosmetic: on a real 249-sample
#'   cohort the two settings share only 25 of their top 50 |iHS| hits, and 4-12% of SNPs
#'   change sign (most in the 0.05-0.1 and >0.3 MAF bands).
#' @param min_maf Minor-allele frequency floor applied at standardisation.
#' @param maf_bands Standardise within this many minor-allele-frequency bands, cut at
#'   quantiles of the observed frequencies so each band holds a similar number of markers.
#'   `NULL` (default) leaves the standardisation to \pkg{rehh} and `freqbin`. The spread of
#'   the log iHH ratio depends on how common the minor allele is -- a rarer allele is
#'   carried by fewer haplotypes, so its integrals are noisier -- and a single bin cannot
#'   remove that, which leaves rarer SNPs over-represented in the tail. Bands do remove it,
#'   without running into the empty bins that make \pkg{rehh}'s own binning unusable on an
#'   unpolarized scan. Around 10 is reasonable; it replaces `freqbin` rather than combining
#'   with it, and it changes every score, so a result computed with it is not comparable
#'   to one computed without.
#' @param min_samples Skip groups smaller than this.
#' @param maxgap Largest gap between consecutive SNPs, in base pairs, that the EHH
#'   integration may cross; `NA` (the default, and \pkg{rehh}'s) lets it cross any gap.
#'   This matters more than its default suggests. A region with no SNPs -- a centromere, a
#'   masked hypervariable block -- has nothing to break the haplotype, so EHH runs flat
#'   across it and the integral accumulates `EHH x gap length`. The SNPs flanking such a
#'   hole then score on the width of the hole rather than on their haplotypes, in either
#'   direction: the ratio is diluted towards zero when both alleles carry EHH into the gap,
#'   and inflated when only one does. Pick a value from the data's own spacing (several
#'   dozen times the median gap leaves ordinary density untouched while stopping at a real
#'   hole) rather than from a round number.
#' @param scalegap Gaps wider than this are counted as being exactly this wide, rather than
#'   stopping the integration outright; `NA` (default) does not rescale. A softer form of
#'   `maxgap` -- it caps a hole's contribution instead of refusing to cross it.
#' @param discard_at_border Return `NA` instead of a truncated integral when the
#'   integration runs into the end of a chromosome or a gap wider than `maxgap`. `NULL`
#'   (default) ties it to `maxgap`: off when no `maxgap` is set (so the markers nearest the
#'   telomeres are still scored), on when one is. On sparse markers this can empty the scan
#'   -- if EHH never decays before the data runs out, every marker is at a border -- so a
#'   scan that comes back mostly `NA` says so.
#' @param threads Threads for \pkg{rehh}.
#' @return A tibble with `group`, `chr`, `pos`, `snp_id`, `freq_minor`, `ihs` and
#'   `neg_log10_p`.
#' @references
#' Voight, B. F., Kudaravalli, S., Wen, X. & Pritchard, J. K. (2006) A map of recent
#' positive selection in the human genome. \emph{PLoS Biology} 4, e72.
#' \doi{10.1371/journal.pbio.0040072}
#'
#' Gautier, M., Klassmann, A. & Vitalis, R. (2017) rehh 2.0: a reimplementation of the R
#' package rehh to detect positive selection from haplotype structure.
#' \emph{Molecular Ecology Resources} 17, 78-90. \doi{10.1111/1755-0998.12634}
#' @seealso [ihs_windows()], [ihs_genes()], [plot_ihs()], [run_rsb()], [beta_score()]
#' @examples
#' ps <- example_pop_structure(umap = FALSE)
#' hap <- parasite_haplotypes(ps, maf = 0.05)
#' run_ihs(hap, group = "country")
#' @export
run_ihs <- function(hap, group = NULL, meta = NULL, polarized = FALSE, freqbin = NULL,
                    min_maf = 0.05, maf_bands = NULL, min_samples = 4, maxgap = NA,
                    scalegap = NA, discard_at_border = NULL, threads = 1) {
  meta <- .normalise_meta(meta)
  .need_package("rehh", "run_ihs()")
  stopifnot(inherits(hap, "parasite_haplotypes"))
  freqbin <- .resolve_freqbin(freqbin, polarized)
  rows <- .ihs_rows(hap, group, meta, min_samples)

  out <- list()
  na_frac <- numeric(0)
  for (l in names(rows)) {
    scan <- .scan_group(hap, rows[[l]], polarized, threads, maxgap, scalegap,
                        discard_at_border)
    if (is.null(scan)) next
    na_frac <- c(na_frac, .border_na_frac(scan))
    res <- .standardise_ihs(scan, freqbin, min_maf, maf_bands)
    if (is.null(res) || !nrow(res)) next
    freq <- if ("FREQ_MIN" %in% names(scan))
      scan$FREQ_MIN[match(paste(res$CHR, res$POSITION), paste(scan$CHR, scan$POSITION))]
    else NA_real_
    out[[length(out) + 1L]] <- data.frame(
      group = l, chr = as.character(res$CHR), pos = as.numeric(res$POSITION),
      snp_id = paste0(res$CHR, ":", format(res$POSITION, scientific = FALSE, trim = TRUE)),
      freq_minor = freq, ihs = res$IHS, neg_log10_p = res$LOGPVALUE,
      stringsAsFactors = FALSE)
  }
  .warn_border_na(na_frac, maxgap, discard_at_border)
  if (!length(out)) {
    warning("no group produced an iHS scan", call. = FALSE)
    return(tibble::tibble())
  }
  df <- do.call(rbind, out)
  df$group <- factor(df$group, levels = names(rows))
  tibble::as_tibble(df)
}

# One block (a group on a chromosome) -> its windows. Positions are sorted, so a window is
# a contiguous index range and findInterval() finds both ends without scanning the block
# once per window.
.window_rows <- function(d, window, step, threshold, metric) {
  d <- d[order(d$pos), , drop = FALSE]
  v <- abs(d[[metric]])
  pos <- d$pos
  # anchor the grid at a multiple of `step` rather than at the first SNP, so every group
  # and every chromosome is cut on the same lines and their windows are comparable
  starts <- seq(floor(min(pos) / step) * step, max(pos), by = step)
  ends <- starts + window
  lo <- findInterval(starts, pos, left.open = TRUE) + 1L   # first SNP with pos >= start
  hi <- findInterval(ends, pos, left.open = TRUE)          # last SNP with pos < end
  n <- hi - lo + 1L
  cs <- c(0, cumsum(v > threshold))
  mx <- vapply(seq_along(starts),
               function(i) if (n[i] > 0) max(v[lo[i]:hi[i]]) else NA_real_, numeric(1))
  data.frame(
    chr = d$chr[1], start = starts, end = ends,
    # the last window reaches past the last SNP; keep its midpoint inside the data so the
    # genome layout cannot place it in the next chromosome's span
    pos = pmin((starts + ends) / 2, max(pos)),
    n_snps = n, n_extreme = cs[hi + 1L] - cs[lo], max_abs = mx,
    stringsAsFactors = FALSE)
}

#' Windowed summary of an iHS scan
#'
#' The fraction of SNPs in each window whose `abs(ihs)` exceeds `threshold` -- the summary
#' iHS is normally read through, rather than SNP by SNP.
#'
#' Per-SNP iHS is a high-variance statistic: each score rests on one focal SNP's EHH decay,
#' so neighbouring SNPs in the same haplotype disagree freely and a genome-wide plot of
#' them looks like grass. A sweep does not raise one SNP, it raises a *run* of them, and
#' counting the run is both what the original description of the statistic proposed and
#' what makes the result comparable to any other windowed track -- an IBD fraction, say.
#'
#' It also sidesteps the per-SNP p-value, which \pkg{rehh} computes from a normal
#' approximation and does not correct for multiple testing: a nominal `p < 0.01` line is
#' crossed by 1% of a neutral genome, so on a whole-genome scan it marks a tail, not a
#' finding.
#'
#' @param scan A [run_ihs()] result, or any table with `chr`, `pos`, the `metric` column
#'   and optionally `group`.
#' @param window Window width in base pairs.
#' @param step Distance between window starts; defaults to `window` (windows that tile
#'   without overlapping). A smaller `step` slides the window and smooths the track, at the
#'   cost of neighbouring windows sharing SNPs.
#' @param threshold The `abs(metric)` a SNP must exceed to be counted as extreme. The
#'   conventional 2 for iHS.
#' @param min_snps Drop windows holding fewer scored SNPs than this. A window with three
#'   SNPs can read 100% and mean nothing; this is what keeps the sparse edges of the data
#'   from becoming the tallest peaks.
#' @param metric Column to summarise (default `"ihs"`); its magnitude is used, so an
#'   unpolarized scan needs no other handling.
#' @return A tibble with `group` (when the scan has one), `chr`, `start`, `end`, `pos` (the
#'   window midpoint), `n_snps`, `n_extreme`, `frac_extreme` and `max_abs`.
#' @references
#' Voight, B. F., Kudaravalli, S., Wen, X. & Pritchard, J. K. (2006) A map of recent
#' positive selection in the human genome. \emph{PLoS Biology} 4, e72.
#' \doi{10.1371/journal.pbio.0040072}
#' @seealso [run_ihs()], [plot_ihs()], [plot_ibd_tugofwar()]
#' @examples
#' ps <- example_pop_structure(umap = FALSE)
#' hap <- parasite_haplotypes(ps, maf = 0.05)
#' ihs_windows(run_ihs(hap, group = "country"), window = 1e5, min_snps = 3)
#' @export
ihs_windows <- function(scan, window = 50000, step = NULL, threshold = 2,
                        min_snps = 10, metric = "ihs") {
  df <- as.data.frame(scan)
  for (nm in c("chr", "pos", metric))
    if (!nm %in% names(df)) stop(sprintf("`scan` has no '%s' column", nm), call. = FALSE)
  pos_num <- function(x, what) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0)
      stop("`", what, "` must be one positive width in base pairs", call. = FALSE)
    x
  }
  window <- pos_num(window, "window")
  step <- if (is.null(step)) window else pos_num(step, "step")

  df <- df[is.finite(df[[metric]]), , drop = FALSE]
  if (!nrow(df)) stop(sprintf("no finite `%s` values to summarise", metric), call. = FALSE)
  df$chr <- as.character(df$chr)
  grp <- if ("group" %in% names(df)) as.character(df$group) else NULL

  blocks <- split(df, if (is.null(grp)) list(df$chr) else list(grp, df$chr), drop = TRUE)
  rows <- lapply(blocks, function(d) {
    w <- .window_rows(d, window, step, threshold, metric)
    if (!is.null(grp)) w <- cbind(group = as.character(d$group[1]), w,
                                  stringsAsFactors = FALSE)
    w
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out <- out[out$n_snps >= min_snps, , drop = FALSE]
  if (!nrow(out))
    stop(sprintf("no window holds %d scored SNP(s); widen `window` or lower `min_snps`",
                 min_snps), call. = FALSE)
  out$frac_extreme <- out$n_extreme / out$n_snps
  # keep the scan's own group order rather than the alphabetical one split() left behind
  if (!is.null(grp) && is.factor(scan$group))
    out$group <- factor(out$group, levels = levels(scan$group))
  out <- out[, c(if (!is.null(grp)) "group", "chr", "start", "end", "pos", "n_snps",
                 "n_extreme", "frac_extreme", "max_abs")]
  tibble::as_tibble(out[order(out$chr, out$start), , drop = FALSE])
}

.cross_pop <- function(hap, group, meta, pairs, polarized, min_samples, threads, fn,
                       value_col, label, maxgap = NA, scalegap = NA,
                       discard_at_border = NULL) {
  .need_package("rehh", label)
  stopifnot(inherits(hap, "parasite_haplotypes"))
  rows <- .ihs_rows(hap, group, meta, min_samples)
  levs <- names(rows)
  if (length(levs) < 2) stop("need at least two groups", call. = FALSE)
  if (is.null(pairs)) pairs <- utils::combn(levs, 2, simplify = FALSE)
  scans <- lapply(rows, function(r) .scan_group(hap, r, polarized, threads, maxgap,
                                               scalegap, discard_at_border))
  .warn_border_na(vapply(scans, .border_na_frac, numeric(1)), maxgap, discard_at_border)

  out <- list()
  for (pr in pairs) {
    a <- pr[1]; b <- pr[2]
    if (!a %in% levs || !b %in% levs)
      stop(sprintf("unknown group in `pairs`: %s vs %s", a, b), call. = FALSE)
    if (is.null(scans[[a]]) || is.null(scans[[b]])) next
    res <- fn(scans[[a]], scans[[b]], popname1 = a, popname2 = b, verbose = FALSE)
    if (is.null(res) || !nrow(res)) next
    val <- res[[grep(value_col, names(res))[1]]]
    lp <- res[[grep("LOGPVALUE", names(res))[1]]]
    out[[length(out) + 1L]] <- data.frame(
      pair = paste(a, "vs", b), pop1 = a, pop2 = b,
      chr = as.character(res$CHR), pos = as.numeric(res$POSITION),
      snp_id = paste0(res$CHR, ":", format(res$POSITION, scientific = FALSE, trim = TRUE)),
      value = val, neg_log10_p = lp, stringsAsFactors = FALSE)
  }
  if (!length(out)) {
    warning("no group pair produced a scan", call. = FALSE)
    return(tibble::tibble())
  }
  tibble::as_tibble(do.call(rbind, out))
}

#' Cross-population extended haplotype homozygosity (Rsb)
#'
#' Compares the site-specific integrated EHH of two groups. Unlike [run_ihs()], Rsb finds
#' sweeps that have gone to completion in one population, because the comparison is
#' against the other population rather than against the other allele.
#'
#' @inheritParams run_ihs
#' @param group Metadata column naming the grouping (required -- there must be at least
#'   two groups to compare).
#' @param pairs Optional list of `c(group_a, group_b)` pairs; defaults to all pairs.
#' @return A tibble with `pair`, `pop1`, `pop2`, `chr`, `pos`, `snp_id`, `value` (Rsb)
#'   and `neg_log10_p`.
#'
#'   **Reading `value`.** Rsb is a log ratio of site-specific extended haplotype
#'   homozygosity between the two populations, standardised to roughly a standard normal
#'   under neutrality. So it is a z-score, and its **sign says which population**:
#'   positive means haplotype homozygosity extends further in `pop1` than `pop2`, i.e. the
#'   sweep is in `pop1`; negative points at `pop2`. Swapping the pair flips the sign.
#'
#'   Magnitude reads like any z: |Rsb| above ~2 is unremarkable in a genome scan, above ~4
#'   is worth a look, and real sweeps in *P. falciparum* run higher still. Rather than
#'   picking a cutoff by eye, use `neg_log10_p` -- the two-sided normal p-value rehh derives
#'   from `value` -- and correct it: at ~20k SNPs a Bonferroni line sits near 5.6, which is
#'   the convention in the literature. [selection_peaks()] will merge the survivors into
#'   loci, since one sweep spans many SNPs, and [annotate_snps()] says which genes they are
#'   in.
#'
#'   Two cautions. The standardisation assumes most of the genome is neutral, so a
#'   genome-wide p-value is relative to *this* comparison and not comparable across pairs
#'   with different sample sizes. And Rsb contrasts two populations, so a high score means
#'   they differ -- it cannot by itself tell a sweep in one from a loss of variation in the
#'   other; that is what the EHH decay curves are for.
#' @references Tang, K., Thornton, K. R. & Stoneking, M. (2007) A new approach for using
#'   genome scans to detect recent positive selection in the human genome.
#'   \emph{PLoS Biology} 5, e171. \doi{10.1371/journal.pbio.0050171}
#' @examples
#' \dontrun{
#' # Rsb contrasts two groups' haplotype homozygosity at the same SNP
#' hap <- parasite_haplotypes(ps, fws = fws, min_fws = 0.95)
#' rsb <- run_rsb(hap, group = "region", pop1 = "north", pop2 = "south")
#' plot_ihs(rsb)
#' }
#' @export
run_rsb <- function(hap, group, meta = NULL, pairs = NULL, polarized = FALSE,
                    min_samples = 4, maxgap = NA, scalegap = NA,
                    discard_at_border = NULL, threads = 1) {
  meta <- .normalise_meta(meta)
  .cross_pop(hap, group, meta, pairs, polarized, min_samples, threads,
             rehh::ines2rsb, "^RSB", "run_rsb()",
             maxgap = maxgap, scalegap = scalegap, discard_at_border = discard_at_border)
}

#' Cross-population extended haplotype homozygosity (XP-EHH)
#'
#' The allele-aware sibling of [run_rsb()]: it contrasts the integrated EHH of the same
#' allele between two populations.
#'
#' @inheritParams run_rsb
#' @return A tibble shaped like [run_rsb()]'s, with `value` holding XP-EHH. Read it exactly
#'   as Rsb -- a standardised log ratio, positive when the extended haplotype is in `pop1` --
#'   the difference being that XP-EHH integrates to a fixed distance while Rsb uses the
#'   site-specific EHH, so XP-EHH is the more sensitive of the two to a sweep that has gone
#'   nearly to fixation, where within-population statistics like iHS lose power.
#' @references Sabeti, P. C. et al. (2007) Genome-wide detection and characterization of
#'   positive selection in human populations. \emph{Nature} 449, 913-918.
#'   \doi{10.1038/nature06250}
#' @examples
#' \dontrun{
#' hap <- parasite_haplotypes(ps, fws = fws, min_fws = 0.95)
#' xp <- run_xpehh(hap, group = "region", pop1 = "north", pop2 = "south")
#' }
#' @export
run_xpehh <- function(hap, group, meta = NULL, pairs = NULL, polarized = FALSE,
                      min_samples = 4, maxgap = NA, scalegap = NA,
                      discard_at_border = NULL, threads = 1) {
  meta <- .normalise_meta(meta)
  .cross_pop(hap, group, meta, pairs, polarized, min_samples, threads,
             rehh::ies2xpehh, "^XPEHH", "run_xpehh()",
             maxgap = maxgap, scalegap = scalegap, discard_at_border = discard_at_border)
}

#' Summarise a haplotype scan per gene
#'
#' The strongest signal inside each gene, which is how selection scans are usually
#' reported: one row per gene rather than per SNP.
#'
#' @param scan The tibble from [run_ihs()], [run_rsb()] or [run_xpehh()].
#' @param genes Gene table (`name`, `chr` or `chrom`, `start`, `end`); defaults to
#'   [PF3D7_GENES]. Coordinates are 0-based half-open.
#' @param within Widen each gene by this many bp on both sides before matching SNPs.
#'   Note what this does to a table of "top genes": one strong SNP is then credited to
#'   every gene within `within` bp of it, so a single sweep can fill several rows that
#'   share a `peak_pos`. `peak_in_gene` says whether the peak is actually inside the
#'   gene or was pulled in from the flank, and `n_snps` counts the widened window.
#' @param min_snps Genes with fewer scored SNPs than this are dropped.
#' @return A tibble with the grouping column, `gene`, `chr`, `start`, `end`, `n_snps`,
#'   `max_neg_log10_p`, `max_abs_value`, `peak_pos` and `peak_in_gene`.
#' @examples
#' ps <- example_pop_structure(umap = FALSE)
#' hap <- parasite_haplotypes(ps, maf = 0.05)
#' ihs_genes(run_ihs(hap, group = "country"), genes = PF_EXAMPLE_DRUG_GENES, min_snps = 1)
#' @export
ihs_genes <- function(scan, genes = NULL, within = 0, min_snps = 1) {
  if (!nrow(scan)) return(tibble::tibble())
  by <- if ("group" %in% names(scan)) "group" else "pair"
  value <- if ("ihs" %in% names(scan)) "ihs" else "value"
  g <- as.data.frame(if (is.null(genes)) PF3D7_GENES else genes)
  if (!"chr" %in% names(g)) {
    for (alt in c("chrom", "Pf3D7_chrom")) if (alt %in% names(g)) { g$chr <- g[[alt]]; break }
  }
  need <- c("name", "chr", "start", "end")
  if (!all(need %in% names(g)))
    stop("`genes` needs name, chr (or chrom), start and end columns", call. = FALSE)

  key <- normalise_chr(scan$chr)
  rows <- list()
  for (i in seq_len(nrow(g))) {
    hit <- key == normalise_chr(g$chr[i]) &
      scan$pos >= g$start[i] - within & scan$pos < g$end[i] + within
    if (!any(hit)) next
    sub <- scan[hit, , drop = FALSE]
    for (l in unique(sub[[by]])) {
      s <- sub[sub[[by]] == l, , drop = FALSE]
      s <- s[is.finite(s$neg_log10_p), , drop = FALSE]
      if (nrow(s) < min_snps) next
      k <- which.max(s$neg_log10_p)
      rows[[length(rows) + 1L]] <- data.frame(
        by = as.character(l), gene = as.character(g$name[i]),
        chr = as.character(g$chr[i]), start = as.numeric(g$start[i]),
        end = as.numeric(g$end[i]), n_snps = nrow(s),
        max_neg_log10_p = s$neg_log10_p[k],
        max_abs_value = max(abs(s[[value]]), na.rm = TRUE),
        peak_pos = s$pos[k],
        peak_in_gene = s$pos[k] >= g$start[i] && s$pos[k] < g$end[i],
        stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(tibble::tibble())
  out <- do.call(rbind, rows)
  names(out)[1] <- by
  out[[by]] <- factor(out[[by]], levels = levels(scan[[by]]) %||% unique(out[[by]]))
  out <- out[order(out[[by]], -out$max_neg_log10_p), , drop = FALSE]
  rownames(out) <- NULL
  tibble::as_tibble(out)
}

#' Keep only some of the haplotypes
#'
#' A [parasite_haplotypes()] object restricted to certain samples or metadata groups, for when a
#' scan or an EHH curve is only worth reading within one part of the cohort -- a mutation that
#' segregates in a single region, say, where pooling everything buries it.
#'
#' Metadata columns are matched the way [PopStructure]'s `$subset()` does, `column = values`,
#' and a column takes several values:
#'
#' ```
#' subset_haplotypes(hap, region = c("North", "Southwest"))
#' ```
#'
#' The **SNP panel is left alone**: the same columns, the same MAF and missingness filtering that
#' built them. That is deliberate, so two subsets stay comparable to each other and to the whole
#' -- and it costs nothing for haplotype work, since a SNP that is monomorphic within the subset
#' is dropped when the scan or curve is computed anyway. Rebuild with `parasite_haplotypes()` on
#' a subset [PopStructure] instead when you want the filtering itself redone within the group.
#'
#' @param x A [parasite_haplotypes()] object.
#' @param samples Optional sample ids to keep.
#' @param meta Metadata to match `...` against; defaults to the object's own.
#' @param ... Metadata filters as `column = values`, e.g. `region = "Southwest"`. A value that
#'   no sample has is an error rather than a silently empty result.
#' @return `x` with fewer haplotypes; `print()` reports the restriction.
#' @seealso [parasite_haplotypes()], [plot_ehh()], [run_ihs()]
#' @examples
#' ps <- example_pop_structure(umap = FALSE)
#' hap <- parasite_haplotypes(ps, maf = 0.05)
#' subset_haplotypes(hap, country = "Ghana")
#' @export
subset_haplotypes <- function(x, samples = NULL, meta = NULL, ...) {
  meta <- .normalise_meta(meta)
  if (!inherits(x, "parasite_haplotypes"))
    stop("`x` must be a parasite_haplotypes() object", call. = FALSE)
  ids <- rownames(x$hap)
  keep <- ids
  note <- character(0)

  if (!is.null(samples)) {
    samples <- as.character(samples)
    unknown <- setdiff(samples, ids)
    if (length(unknown) == length(samples))
      stop("none of `samples` are in these haplotypes", call. = FALSE)
    if (length(unknown))
      warning(length(unknown), " of ", length(samples), " `samples` are not in these ",
              "haplotypes and are ignored", call. = FALSE)
    keep <- intersect(keep, samples)
    note <- c(note, paste0("samples: ", length(samples)))
  }

  filt <- list(...)
  if (length(filt)) {
    m <- meta %||% x$meta
    if (is.null(m)) stop("no metadata on these haplotypes to match `...` against; pass `meta`",
                         call. = FALSE)
    m <- as.data.frame(m)
    if (!"sample" %in% names(m)) stop("`meta` needs a `sample` column", call. = FALSE)
    bad <- setdiff(names(filt), names(m))
    if (length(bad))
      stop("not a metadata column: ", paste(bad, collapse = ", "), call. = FALSE)
    ok <- rep(TRUE, nrow(m))
    for (nm in names(filt)) {
      want <- as.character(filt[[nm]])
      have <- as.character(m[[nm]])
      missing <- setdiff(want, unique(have))
      # a typo'd level would otherwise just return nothing, which reads as "no data"
      if (length(missing))
        stop("no sample has ", nm, " = ", paste(missing, collapse = ", "), "; available: ",
             paste(sort(unique(have)), collapse = ", "), call. = FALSE)
      ok <- ok & have %in% want
      note <- c(note, paste0(nm, ": ", paste(want, collapse = ", ")))
    }
    keep <- intersect(keep, as.character(m$sample[ok]))
  }

  if (!length(keep))
    stop("no haplotype satisfies every filter at once", call. = FALSE)
  out <- x
  out$hap <- x$hap[keep, , drop = FALSE]
  if (!is.null(out$meta) && "sample" %in% names(out$meta))
    out$meta <- out$meta[as.character(out$meta$sample) %in% keep, , drop = FALSE]
  out$subset <- list(from = nrow(x$hap), by = note)
  out
}

# One marker's per-sample allele index, read straight from the callset rather than through a
# dosage matrix. SNPRelate cannot express a third allele -- `copy.num.of.ref` counts
# reference copies, so ALT1 and ALT2 both come back as "not reference" -- which is exactly
# the distinction a multiallelic EHH is about.
.read_allele_indices <- function(vcf, het = "missing", seed = 42) {
  if (!nzchar(Sys.which("bcftools")))
    stop("reading allele indices needs bcftools on PATH", call. = FALSE)
  if (!file.exists(vcf)) stop("no such file: ", vcf, call. = FALSE)
  samples <- system2("bcftools", c("query", "-l", shQuote(vcf)), stdout = TRUE,
                     stderr = FALSE)
  lines <- system2("bcftools",
                   c("query", "-f", shQuote("%CHROM\\t%POS[\\t%GT]\\n"), shQuote(vcf)),
                   stdout = TRUE, stderr = FALSE)
  if (!length(lines)) stop("no records in ", basename(vcf), call. = FALSE)
  parts <- strsplit(lines, "\t", fixed = TRUE)
  ids <- vapply(parts, function(p) paste0(normalise_chr(p[1]), ":",
                                          as.integer(p[2]) - 1L), character(1))
  H <- matrix(NA_integer_, nrow = length(samples), ncol = length(parts),
              dimnames = list(samples, ids))
  n_het <- 0L
  set.seed(seed)
  for (j in seq_along(parts)) {
    gt <- parts[[j]][-(1:2)]
    a <- strsplit(gt, "[/|]")
    idx <- vapply(a, function(v) {
      v <- v[v != "."]
      if (!length(v)) return(NA_integer_)
      u <- unique(as.integer(v))
      if (length(u) == 1L) return(u)
      -1L                                  # heterozygous: two alleles, one haplotype wanted
    }, integer(1))
    n_het <- n_het + sum(idx == -1L, na.rm = TRUE)
    mixed <- which(!is.na(idx) & idx == -1L)
    if (length(mixed)) {
      if (identical(het, "missing")) {
        idx[mixed] <- NA_integer_
      } else {
        # draw between the alleles that sample actually carries, weighted by how common
        # each is among the samples that are unambiguous -- the k-allele form of the draw
        # `parasite_haplotypes()` makes at the ALT frequency
        clean <- idx[!is.na(idx) & idx >= 0]
        for (i in mixed) {
          have <- unique(as.integer(a[[i]][a[[i]] != "."]))
          w <- vapply(have, function(v) sum(clean == v), numeric(1))
          idx[i] <- if (sum(w) > 0) sample(have, 1L, prob = w) else sample(have, 1L)
        }
      }
    }
    H[, j] <- idx
  }
  attr(H, "n_het") <- n_het
  H
}

#' Add markers to a haplotype set from a callset of their own
#'
#' Splices one or more markers into a [parasite_haplotypes()] object, reading them as
#' **allele indices** so a site with more than two alleles keeps them apart.
#'
#' This exists because the dosage route cannot carry the information. A genotype matrix
#' counts copies of one allele, so at a triallelic site every non-reference call collapses
#' to the same number no matter which alternate it carries -- `load_genotypes()` says as
#' much when asked for `variants = "all"`. A marker that was dropped from the main callset
#' *for being multiallelic*, then called separately to get it back, therefore has to enter
#' at the haplotype layer, where an allele is an index rather than a count.
#'
#' Samples are matched by name, not position. The marker is inserted in coordinate order, so
#' the haplotypes still read along the genome, and a position already present is an error
#' rather than a silent replacement.
#'
#' @param hap A [parasite_haplotypes()] object.
#' @param vcf A VCF/BCF holding the marker(s) to add. Needs `bcftools` on `PATH`.
#' @param het What to do with a call carrying two different alleles, which cannot be one
#'   haplotype: `"missing"` (default) sets it aside for the object's imputation, `"draw"`
#'   picks between the alleles that sample carries, weighted by how common each is among the
#'   unambiguous samples.
#' @param impute Fill the resulting gaps by drawing at the marker's allele frequencies
#'   (default `TRUE`). \pkg{rehh} cannot read a missing call, so a marker left with gaps
#'   would drop those haplotypes from every curve.
#' @param seed Seed for the draws.
#' @return The `parasite_haplotypes` object with the markers added.
#' @examples
#' \dontrun{
#' hap <- parasite_haplotypes(ps, fws = fws)
#' hap <- add_haplotype_markers(hap, "PF3D7_1343700.1-AA469.bcf", het = "draw")
#' plot_ehh(hap, "Pf3D7_13_v3:1725591")      # three curves, one per allele
#' }
#' @seealso [parasite_haplotypes()], [plot_ehh()]
#' @export
add_haplotype_markers <- function(hap, vcf, het = c("missing", "draw"), impute = TRUE,
                                  seed = 42) {
  stopifnot(inherits(hap, "parasite_haplotypes"))
  het <- match.arg(het)
  new <- .read_allele_indices(vcf, het = het, seed = seed)

  # The two callsets need not spell chromosomes the same way: SNPRelate keeps `Pf3D7_13_v3`
  # from one file and reduces `chr13` to `13` in another, depending on what it recognised. A
  # marker that keeps its own spelling lands on a chromosome of its own, with no neighbours
  # to decay against -- which surfaces as "too few polymorphic SNPs on that chromosome"
  # rather than as a naming problem. So match on the normalised name and adopt the spelling
  # the haplotypes already use.
  have <- unique(hap$map$chr)
  lookup <- stats::setNames(have, normalise_chr(have))
  as_have <- lookup[normalise_chr(sub(":.*", "", colnames(new)))]
  colnames(new) <- ifelse(is.na(as_have), colnames(new),
                          paste0(as_have, ":", sub(".*:", "", colnames(new))))

  miss <- setdiff(rownames(hap$hap), rownames(new))
  if (length(miss))
    stop(sprintf("%d of the haplotype set's samples are not in %s (e.g. %s)",
                 length(miss), basename(vcf), paste(utils::head(miss, 3), collapse = ", ")),
         call. = FALSE)
  # A left join on the haplotypes: the marker is usually called on the whole cohort while the
  # haplotypes are a monoclonal subset, so extras are expected -- but say how many, or a name
  # mismatch that silently drops most of the callset looks like a clean merge.
  spare <- setdiff(rownames(new), rownames(hap$hap))
  if (length(spare))
    message(sprintf("%s has %d sample(s) not in the haplotypes; they are left out",
                    basename(vcf), length(spare)))
  new <- new[rownames(hap$hap), , drop = FALSE]

  clash <- intersect(colnames(hap$hap), colnames(new))
  if (length(clash))
    stop(sprintf("%s is already in the haplotypes; drop it first if the new call replaces it",
                 paste(clash, collapse = ", ")), call. = FALSE)

  n_het <- attr(new, "n_het") %||% 0L
  if (impute) {
    for (j in seq_len(ncol(new))) {
      gap <- which(is.na(new[, j]))
      obs <- new[!is.na(new[, j]), j]
      if (length(gap) && length(obs)) new[gap, j] <- sample(obs, length(gap), replace = TRUE)
    }
  }
  H <- cbind(hap$hap, new)
  loc <- .parse_snp_ids(colnames(H))
  chr_num <- suppressWarnings(as.numeric(loc$chr))
  H <- H[, order(is.na(chr_num), chr_num, loc$chr, loc$pos), drop = FALSE]
  storage.mode(H) <- "integer"

  map <- .parse_snp_ids(colnames(H))
  map$snp_id <- colnames(H)
  hap$hap <- H
  hap$map <- map[c("chr", "pos", "snp_id")]
  hap$filtering$n_added_markers <- (hap$filtering$n_added_markers %||% 0L) + ncol(new)
  message(sprintf("added %d marker(s) from %s%s", ncol(new), basename(vcf),
                  if (n_het) sprintf("; %d mixed call(s) %s", n_het,
                                     if (identical(het, "draw")) "resolved by an allele draw"
                                     else "set to missing") else ""))
  hap
}
