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
                                alleles = c("auto", "dosage", "index")) {
  meta <- .normalise_meta(meta)
  het <- match.arg(het)
  alleles <- match.arg(alleles)
  if (inherits(x, "PopStructure")) {
    # "auto" takes whatever the object holds. Haplotype work is the case that *can* use a
    # multiallelic panel, so an object built from allele indices should reach it without
    # being asked twice -- and `alleles = "dosage"` is the one-argument way back to the
    # biallelic reading, which the object derives on demand.
    if (identical(alleles, "auto")) alleles <- if (identical(x$encoding(), "allele_index"))
      "index" else "dosage"
    G <- .geno_for(x, genotype, what = "parasite_haplotypes()",
                   needs = if (identical(alleles, "index")) "allele_index" else "dosage")
    if (is.null(meta)) meta <- x$get_meta()
  } else {
    # A raw `load_genotypes()` list still carries its `encoding`, so honour it: "auto" reads
    # the list rather than assuming dosage, and an index list reaches the haplotype path
    # without being routed through a PopStructure first. Otherwise a list built with
    # `encoding = "allele_index"` -- the whole point of which is to keep the alternates
    # apart -- was refused here as "not dosages".
    if (identical(alleles, "auto"))
      alleles <- if (is.list(x) && identical(x$encoding %||% "dosage", "allele_index"))
        "index" else "dosage"
    G <- .coerce_geno(x, "parasite_haplotypes()",
                      needs = if (identical(alleles, "index")) "allele_index" else "dosage")
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

# Minor-allele frequency of one allele-index column: everything not the commonest allele.
#
# This is the k-allele form, and it is a strict generalisation rather than an alternative:
# on a 0/1 column `1 - max(n_0, n_1)/n` is `min(n_0, n_1)/n`, which is exactly
# `min(p, 1 - p)`. So it can be used wherever the dosage formula was without moving a
# biallelic number. On allele *indices* the dosage formula is not a frequency at all --
# `mean()` of indices 0/1/2 is a mean index, and `min(p, 1 - p)` on it goes negative as
# soon as the mean exceeds one.
.minor_af <- function(v) {
  v <- v[!is.na(v)]
  if (!length(v)) return(NA_real_)
  1 - max(tabulate(as.integer(v) + 1L)) / length(v)
}

# The same, column by column over a haplotype matrix.
.index_minor_af <- function(H) {
  vapply(seq_len(ncol(H)), function(j) .minor_af(H[, j]), numeric(1))
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

# Markers carrying more than two alleles in this group, and how much of the group they cost.
#
# rehh's scan_hh reports one major and one minor allele whatever it is handed: at a marker
# with three, it keeps the two commonest and drops the rest, so the score is a two-allele
# contrast computed on a subset of the haplotypes and the reported FREQ_MIN is the second
# commonest allele's frequency rather than a minor-allele frequency. Nothing in the returned
# table says this happened, which is why it is worth counting here.
.multiallelic_drop <- function(hap, rows) {
  h <- hap$hap[rows, , drop = FALSE]
  k <- apply(h, 2, function(v) length(unique(v[!is.na(v)])))
  hit <- which(k > 2)
  if (!length(hit)) return(NULL)
  drop <- vapply(hit, function(i) {
    u <- h[, i]
    tb <- sort(table(u[!is.na(u)]), decreasing = TRUE)
    sum(tb[-(1:2)]) / sum(tb)
  }, numeric(1))
  ids <- if (!is.null(hap$map$snp_id)) hap$map$snp_id[hit] else names(k)[hit]
  data.frame(snp_id = as.character(ids), n_alleles = as.integer(k[hit]),
             frac_dropped = drop, stringsAsFactors = FALSE)
}

# One warning for the whole scan rather than one per group per marker.
#
# `what` names the statistic that was actually run. The reduction is rehh's and is the same
# for all of them -- it scores each marker on its two commonest alleles -- but the message
# should not tell someone running Rsb that their iHS is affected, and `freq_minor` is a
# column only the iHS scan has.
.warn_multiallelic <- function(ma, what = "iHS") {
  ma <- do.call(rbind, Filter(Negate(is.null), ma))
  if (is.null(ma) || !nrow(ma)) return(invisible(NULL))
  worst <- ma[order(-ma$frac_dropped), , drop = FALSE]
  ids <- unique(worst$snp_id)
  freq_clause <- if (identical(what, "iHS"))
    " and their `freq_minor` is not a minor-allele frequency" else ""
  warning(sprintf(
    paste0("%d marker(s) carry more than two alleles; rehh scores each on its two ",
           "commonest and drops the rest, so their %s is a two-allele contrast on a ",
           "subset of the haplotypes%s. Worst: %s (%.0f%% of haplotypes excluded). ",
           "Recode such a marker to the pairwise contrast you mean -- keep the ",
           "haplotypes carrying either of two alleles and code them 0/1 -- rather than ",
           "reading the score as it stands."),
    length(ids), what, freq_clause, worst$snp_id[1],
    100 * worst$frac_dropped[1]), call. = FALSE)
  invisible(ma)
}

# ---- contrasts -------------------------------------------------------------

# Per-allele iHH at one marker, from rehh's own single-marker integrator.
#
# `calc_ehh()` returns IHH_A, IHH_D1, IHH_D2, ... -- one integral per allele -- where
# `scan_hh()` reports only two. That is the whole difference, and it is a reporting
# difference rather than a limit of the method.
#
# These integrals are **invariant** to which other alleles are present, because EHH for an
# allele class only ever involves that class's haplotypes. Verified against explicit subsets
# to nine decimal places. So every contrast at a marker comes out of one call, and a
# contrast is a choice of which two integrals to divide rather than a separate scan.
.marker_allele_ihh <- function(hap, rows, snp_id, polarized = FALSE, maxgap = NA,
                               scalegap = NA, discard_at_border = NULL) {
  .need_package("rehh", "the haplotype scans")
  j <- match(snp_id, hap$map$snp_id)
  if (is.na(j)) return(NULL)
  objs <- .haplohh_list(hap, rows)
  chr <- as.character(hap$map$chr[j])
  o <- objs[[chr]]
  if (is.null(o)) return(NULL)
  mrk <- match(as.numeric(hap$map$pos[j]), o@positions)
  if (is.na(mrk)) return(NULL)
  e <- try(rehh::calc_ehh(o, mrk = mrk, polarized = polarized, maxgap = maxgap,
                          scalegap = scalegap,
                          discard_integration_at_border =
                            .resolve_border(discard_at_border, maxgap)),
           silent = TRUE)
  if (inherits(e, "try-error") || is.null(e$ihh)) return(NULL)
  # rehh names them A / D1 / D2 ...; the package speaks in allele indices, and the order is
  # positional, so index i is `hap`'s allele i
  list(ihh = stats::setNames(as.numeric(e$ihh), seq_along(e$ihh) - 1L),
       freq = stats::setNames(as.numeric(e$freq), seq_along(e$freq) - 1L))
}

# Which allele pairs a marker is scored on.
#
# `ref` is the default and the one that keeps independent origins apart: each alternate is
# compared against the **reference only**, so a marker carrying three changes gives three
# separate answers rather than one that has quietly pooled or dropped some of them.
# `pairwise` adds alternate-against-alternate, which asks a different question -- how two
# origins compare to each other -- and should look different in the call.
.allele_pairs <- function(alleles, contrast) {
  if (length(alleles) < 2L) return(list())
  if (identical(contrast, "pairwise")) {
    cb <- utils::combn(sort(alleles), 2L, simplify = FALSE)
    return(cb)
  }
  ref <- min(alleles)
  lapply(setdiff(sort(alleles), ref), function(a) c(ref, a))
}

# rehh's own frequency binning, reimplemented so it can be applied to rows that share a
# position -- which contrast rows do, and which every position-keyed lookup in `ihh2ihs`
# assumes cannot happen. `cut()` on `seq(min_maf, 1 - min_maf, freqbin)` then a z-score
# within bin is exactly what `ihh2ihs` does; `freqbin >= 1` means one bin, as there.
.standardise_unihs <- function(unihs, freq, freqbin, min_maf, maf_bands) {
  keep <- is.finite(unihs) & !is.na(freq) & freq >= min_maf & freq <= 1 - min_maf
  z <- rep(NA_real_, length(unihs))
  if (!any(keep)) return(list(ihs = z, logp = z))
  if (!is.null(maf_bands)) {
    raw <- data.frame(CHR = NA, POSITION = seq_along(unihs), UNIHS = unihs)
    out <- .band_standardise(raw[keep, , drop = FALSE], freq[keep], maf_bands)
    z[keep] <- out$IHS
    lp <- rep(NA_real_, length(unihs)); lp[keep] <- out$LOGPVALUE
    return(list(ihs = z, logp = lp))
  }
  fb <- if (freqbin >= 1) (1 - 2 * min_maf) / round(freqbin) else freqbin
  br <- seq(min_maf, 1 - min_maf, fb)
  bins <- if (length(br) < 2L) factor(rep("all", sum(keep)))
          else cut(freq[keep], breaks = br, include.lowest = TRUE)
  u <- unihs[keep]
  m <- tapply(u, bins, mean, na.rm = TRUE)
  sdv <- tapply(u, bins, stats::sd, na.rm = TRUE)
  zz <- (u - m[bins]) / sdv[bins]
  z[keep] <- as.numeric(zz)
  lp <- rep(NA_real_, length(unihs))
  lpk <- -log10(2 * stats::pnorm(-abs(as.numeric(zz))))
  if (any(is.finite(lpk))) lpk[is.infinite(lpk)] <- max(lpk[is.finite(lpk)]) + 1
  lp[keep] <- lpk
  list(ihs = z, logp = lp)
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
  data.frame(CHR = raw$CHR, POSITION = raw$POSITION, UNIHS = as.numeric(u),
             IHS = as.numeric(z), LOGPVALUE = as.numeric(lp), stringsAsFactors = FALSE)
}

# One group's scan -> its standardised iHS, by rehh's own binning or by our frequency bands.
.standardise_ihs <- function(scan, freqbin, min_maf, maf_bands) {
  if (is.null(maf_bands)) {
    std <- rehh::ihh2ihs(scan, freqbin = freqbin, min_maf = min_maf, verbose = FALSE)$ihs
    if (is.null(std) || !nrow(std)) return(std)
    # rehh drops the unstandardised ratio when it standardises, so ask for it separately and
    # match on position: same scan, same min_maf, so the two agree row for row
    raw <- rehh::ihh2ihs(scan, freqbin = 1, min_maf = min_maf, standardize = FALSE,
                         verbose = FALSE)$ihs
    std$UNIHS <- if (is.null(raw)) NA_real_ else
      raw$UNIHS[match(paste(std$CHR, std$POSITION), paste(raw$CHR, raw$POSITION))]
    return(std)
  }
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
# The contrast path: one row per (marker, allele pair), standardised together.
#
# The scan is run once per group as before, for the markers that have one contrast anyway.
# Only the markers carrying more than two alleles cost anything extra, and each costs a
# single `calc_ehh()` -- not a re-scan, because the per-allele integrals do not depend on
# which other alleles are present.
#
# The standardisation is over **all** rows at once. A multiallelic marker contributing two
# rows must have both judged against the same genome-wide distribution, which is why the
# binning is done here rather than by `ihh2ihs()`: its output is keyed by position, and two
# contrasts at one position are two rows with the same key.
# `scan_hh()` reports a two-allele marker as MAJ/MIN -- ordered by *frequency*, with nothing
# in the output saying which allele is which. The contrast rows are labelled by allele
# **index** ("0>1"), so at every marker where allele 0 happens to be the minor one the label
# and the arithmetic disagree and the ratio comes out inverted -- an exact sign flip on
# log(iHH_0 / iHH_1), not an approximation. Measured on a 60-haplotype panel simulated at
# p(alt) = 0.62: 59 of the 60 markers with allele 0 in the minority, |diff| up to 0.66.
#
# The haplotypes do carry the identity, so decide from them: whichever FREQ_ column equals
# allele 0's frequency names allele 0's IHH_ column. At an exact 50/50 tie both columns match
# and there is no way to tell, so return NA and let the caller pay for one `calc_ehh()`, which
# is indexed by allele and cannot be ambiguous. Same for the case where rehh scored a
# different set of haplotypes than we counted -- then neither column matches.
.biallelic_unihs <- function(scan, ihh_cols, frq_cols, i, f0) {
  hit <- abs(c(as.numeric(scan[[frq_cols[1]]][i]),
               as.numeric(scan[[frq_cols[2]]][i])) - f0) < 1e-9
  if (sum(hit, na.rm = TRUE) != 1L) return(NA_real_)
  num <- which(hit); den <- 3L - num
  log(as.numeric(scan[[ihh_cols[num]]][i]) / as.numeric(scan[[ihh_cols[den]]][i]))
}

.run_ihs_contrasts <- function(hap, rows, contrast, polarized, freqbin, min_maf, maf_bands,
                               maxgap, scalegap, discard_at_border, threads, named_group) {
  out <- list()
  na_frac <- numeric(0)
  for (l in names(rows)) {
    r <- rows[[l]]
    scan <- .scan_group(hap, r, polarized, threads, maxgap, scalegap, discard_at_border)
    if (is.null(scan) || !nrow(scan)) next
    na_frac <- c(na_frac, .border_na_frac(scan))

    ihh_cols <- grep("^IHH_", names(scan))
    frq_cols <- grep("^FREQ_", names(scan))
    ids <- paste0(scan$CHR, ":", format(scan$POSITION, scientific = FALSE, trim = TRUE))
    h <- hap$hap[r, , drop = FALSE]
    jj <- match(ids, hap$map$snp_id)
    n_alleles <- vapply(jj, function(j)
      if (is.na(j)) 2L else length(unique(h[!is.na(h[, j]), j])), integer(1))
    # allele 0's frequency among the haplotypes rehh scored, which is what says whether the
    # scan's MAJ column is allele 0 or allele 1
    freq0 <- vapply(jj, function(j)
      if (is.na(j)) NA_real_ else mean(h[!is.na(h[, j]), j] == 0L), numeric(1))

    rec <- list()
    for (i in seq_len(nrow(scan))) {
      if (n_alleles[i] <= 2L) {
        # two alleles, one contrast -- the scan's own numbers, oriented onto the label
        u <- .biallelic_unihs(scan, ihh_cols, frq_cols, i, freq0[i])
        if (!is.na(u)) {
          rec[[length(rec) + 1L]] <- data.frame(
            chr = as.character(scan$CHR[i]), pos = as.numeric(scan$POSITION[i]),
            snp_id = ids[i], contrast = "0>1",
            freq = as.numeric(scan[[frq_cols[2]]][i]),
            unihs = u, stringsAsFactors = FALSE)
          next
        }
        # a 50/50 tie: which allele the scan called MAJ is unknowable from its output, so
        # fall through to the per-allele integrals, which are indexed by allele
      }
      a <- .marker_allele_ihh(hap, r, ids[i], polarized, maxgap, scalegap,
                              discard_at_border)
      if (is.null(a)) next
      for (pr in .allele_pairs(as.integer(names(a$ihh)), contrast)) {
        i1 <- as.character(pr[1]); i2 <- as.character(pr[2])
        rec[[length(rec) + 1L]] <- data.frame(
          chr = as.character(scan$CHR[i]), pos = as.numeric(scan$POSITION[i]),
          snp_id = ids[i], contrast = paste0(pr[1], ">", pr[2]),
          freq = unname(a$freq[i2]),
          unihs = log(unname(a$ihh[i1]) / unname(a$ihh[i2])),
          stringsAsFactors = FALSE)
      }
    }
    if (!length(rec)) next
    df <- do.call(rbind, rec)
    st <- .standardise_unihs(df$unihs, df$freq, freqbin, min_maf, maf_bands)
    df$group <- l
    df$freq_minor <- df$freq
    df$ihs <- st$ihs
    df$neg_log10_p <- st$logp
    out[[length(out) + 1L]] <- df[!is.na(df$ihs),
                                  c("group", "chr", "pos", "snp_id", "contrast",
                                    "freq_minor", "unihs", "ihs", "neg_log10_p")]
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
#' @section Contrasts:
#' `log(iHH_A / iHH_B)` is a ratio, so it needs exactly two terms and there is no k-allele
#' iHS. At a marker carrying three alleles the honest answer is **k-1 contrasts**, each saying
#' which two it compared, and the `contrast` column names them as `"0>1"`, `"0>2"` and so on
#' in allele-index order (`$sites$alt` on the panel says which base each index is).
#'
#' `"ref"` compares each alternate against the reference only. That is what keeps independent
#' origins apart: at a codon where three changes arose separately, pooling them into one
#' "not reference" class merges exactly the distinction the scan exists to draw.
#'
#' What makes this cheap is that **per-allele iHH does not depend on which other alleles are
#' at the marker** -- EHH for an allele class only ever involves that class's haplotypes. So
#' every contrast at a marker comes out of one [rehh::calc_ehh()] call, and the numbers are
#' identical to what an explicitly subset and recoded panel would give (verified to nine
#' decimal places). Only multiallelic markers cost anything extra.
#'
#' A biallelic marker has one contrast, so `"ref"` reproduces `"none"` exactly, value for
#' value. The standardisation is over every row at once, so a marker contributing two rows has
#' both judged against the same genome-wide distribution.
#'
#' Downstream, `contrast` is a grouping key: [ihs_windows()], [ihs_genes()],
#' [selection_peaks()] and [plot_ihs()] all split on it, because two contrasts at one position
#' are two measurements rather than two SNPs at one site.
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
#' @return A tibble with `group`, `chr`, `pos`, `snp_id`, `freq_minor`, `unihs`, `ihs` and
#'   `neg_log10_p`.
#'
#'   `unihs` is the **un**standardised statistic, `log(iHH_major / iHH_minor)` (ancestral
#'   over derived when `polarized = TRUE`), so `exp(unihs)` is the integrated-EHH ratio
#'   itself. `ihs` is that value z-scored within its frequency band, which is what makes
#'   scores comparable along the genome but also throws the scale away -- the band's mean and
#'   sd are not recoverable from `ihs` alone, so keep `unihs` if you ever want the ratio back.
#'   Note that `ihs = 0` does **not** mean a ratio of 1: it means average for that group, and
#'   the average is below 1 wherever minor-allele haplotypes are systematically longer.
#'
#'   `ihs` and `neg_log10_p` are `NA` wherever the integral could not be formed for one of
#'   the two alleles -- see the note on missing scores below.
#' @references
#' Voight, B. F., Kudaravalli, S., Wen, X. & Pritchard, J. K. (2006) A map of recent
#' positive selection in the human genome. \emph{PLoS Biology} 4, e72.
#' \doi{10.1371/journal.pbio.0040072}
#'
#' Gautier, M., Klassmann, A. & Vitalis, R. (2017) rehh 2.0: a reimplementation of the R
#' package rehh to detect positive selection from haplotype structure.
#' \emph{Molecular Ecology Resources} 17, 78-90. \doi{10.1111/1755-0998.12634}
#' @section Why a SNP can appear for one group only, or score `NA`:
#' The scan is per group, so a SNP is tested in a group only where it is polymorphic there
#' and clears `min_maf` there. A variant private to one region therefore has one row, not
#' one per region, and that is a statement about the cohort rather than a fault.
#'
#' A row can be present with `ihs` and `neg_log10_p` both `NA`. That means the SNP passed the
#' frequency filter but \pkg{rehh} could not integrate EHH for at least one of its two
#' alleles, so the log ratio is undefined. The usual causes are too few haplotypes carrying
#' the minor allele for the decay to be estimated, and EHH that never falls below the cutoff
#' before the data runs out -- a chromosome end, or a gap wider than `maxgap`. Both get more
#' common in small groups and at low minor-allele counts, which is the same corner where a
#' score that *is* returned deserves the least trust. Treat `NA` as "not measurable here",
#' not as "no selection".
#' @seealso [ihs_windows()], [ihs_genes()], [plot_ihs()], [run_rsb()], [beta_score()]
#' @examples
#' ps <- example_pop_structure(umap = FALSE)
#' hap <- parasite_haplotypes(ps, maf = 0.05)
#' run_ihs(hap, group = "country")
#' @export
run_ihs <- function(hap, group = NULL, meta = NULL, polarized = FALSE, freqbin = NULL,
                    min_maf = 0.05, maf_bands = NULL, min_samples = 4, maxgap = NA,
                    scalegap = NA, discard_at_border = NULL, threads = 1,
                    contrast = c("ref", "pairwise", "none")) {
  meta <- .normalise_meta(meta)
  .need_package("rehh", "run_ihs()")
  stopifnot(inherits(hap, "parasite_haplotypes"))
  contrast <- match.arg(contrast)
  freqbin <- .resolve_freqbin(freqbin, polarized)
  rows <- .ihs_rows(hap, group, meta, min_samples)
  if (!identical(contrast, "none") && !identical(hap$alleles, "index")) {
    # A dosage panel has one alternate by construction, so there is only ever one contrast to
    # draw and the two paths agree row for row. Fall back rather than making the caller ask.
    contrast <- "none"
  }
  if (!identical(contrast, "none"))
    return(.run_ihs_contrasts(hap, rows, contrast, polarized, freqbin, min_maf, maf_bands,
                              maxgap, scalegap, discard_at_border, threads,
                              named_group = !is.null(group)))

  out <- list()
  na_frac <- numeric(0)
  multi <- list()
  for (l in names(rows)) {
    # per group, because a marker can be biallelic inside one and not inside another
    multi[[l]] <- .multiallelic_drop(hap, rows[[l]])
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
      freq_minor = freq, unihs = res$UNIHS, ihs = res$IHS, neg_log10_p = res$LOGPVALUE,
      stringsAsFactors = FALSE)
  }
  .warn_border_na(na_frac, maxgap, discard_at_border)
  .warn_multiallelic(multi)
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
  # A contrast is a separate measurement that happens to sit at the same position as another,
  # so two of them must not be summarised as two SNPs at one site. Splitting on it keeps each
  # allele's scan its own scan.
  ctr <- if ("contrast" %in% names(df)) as.character(df$contrast) else NULL

  key <- c(if (!is.null(grp)) list(grp), if (!is.null(ctr)) list(ctr), list(df$chr))
  blocks <- split(df, key, drop = TRUE)
  rows <- lapply(blocks, function(d) {
    w <- .window_rows(d, window, step, threshold, metric)
    if (!is.null(ctr)) w <- cbind(contrast = as.character(d$contrast[1]), w,
                                  stringsAsFactors = FALSE)
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
  out <- out[, c(if (!is.null(grp)) "group", if (!is.null(ctr)) "contrast",
                 "chr", "start", "end", "pos", "n_snps",
                 "n_extreme", "frac_extreme", "max_abs")]
  tibble::as_tibble(out[order(out$chr, out$start), , drop = FALSE])
}

# How many alleles each marker carries within one group, aligned to `ids`. Reported rather
# than warned about, so a caller can drop the markers where the two populations disagree.
.alleles_at <- function(hap, rows, ids) {
  j <- match(ids, hap$map$snp_id)
  h <- hap$hap[rows, , drop = FALSE]
  vapply(j, function(k) if (is.na(k)) NA_integer_
         else length(unique(h[!is.na(h[, k]), k])), integer(1))
}

# "run_rsb()" -> "Rsb". Kept beside .cross_pop so the two cannot drift.
.stat_name <- function(label) {
  switch(label, "run_rsb()" = "Rsb", "run_xpehh()" = "XP-EHH", sub("\\(\\)$", "", label))
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
    ids <- paste0(res$CHR, ":", format(res$POSITION, scientific = FALSE, trim = TRUE))
    out[[length(out) + 1L]] <- data.frame(
      pair = paste(a, "vs", b), pop1 = a, pop2 = b,
      chr = as.character(res$CHR), pos = as.numeric(res$POSITION),
      snp_id = ids,
      n_alleles_pop1 = .alleles_at(hap, rows[[a]], ids),
      n_alleles_pop2 = .alleles_at(hap, rows[[b]], ids),
      value = val, neg_log10_p = lp, stringsAsFactors = FALSE)
  }
  # iES is a site-level homozygosity pooled across allele classes, so it falls as the focal is
  # split into more of them -- measured at 36% for three alleles and 44% for four, on
  # identical haplotypes. XP-EHH divides two of those, so a marker carrying different numbers
  # of alleles in the two populations is compared on different footings. iNES normalises by
  # the focal homozygosity and moves by 1% or less, so Rsb is not affected and does not warn.
  if (identical(label, "run_xpehh()") && length(out)) {
    all_out <- do.call(rbind, out)
    hit <- which(all_out$n_alleles_pop1 != all_out$n_alleles_pop2 &
                   pmax(all_out$n_alleles_pop1, all_out$n_alleles_pop2) > 2L)
    if (length(hit)) {
      worst <- hit[which.max(abs(all_out$n_alleles_pop1[hit] - all_out$n_alleles_pop2[hit]))]
      warning(sprintf(paste0("%d marker(s) carry a different allele count in the two ",
                             "populations (worst: %s, %d vs %d). XP-EHH divides their iES, ",
                             "and iES falls with the number of allele classes at the focal ",
                             "-- about 36%% for three and 44%% for four on identical ",
                             "haplotypes -- so those scores partly measure the diversity ",
                             "difference rather than the haplotype-length difference. ",
                             "`n_alleles_pop1`/`n_alleles_pop2` are reported so they can be ",
                             "dropped; Rsb reads iNES instead and is not affected."),
                      length(hit), all_out$snp_id[worst],
                      all_out$n_alleles_pop1[worst], all_out$n_alleles_pop2[worst]),
              call. = FALSE)
    }
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
  # A gene is summarised once per group *and once per contrast*: two contrasts at one
  # position are two measurements, and taking the max over both would report whichever
  # allele happened to score higher as though it were the gene's answer.
  by <- if ("group" %in% names(scan)) "group" else if ("pair" %in% names(scan)) "pair" else
    NULL
  has_contrast <- "contrast" %in% names(scan)
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
    keys <- if (is.null(by)) rep("all", nrow(sub)) else as.character(sub[[by]])
    if (has_contrast) keys <- paste(keys, as.character(sub$contrast), sep = "\r")
    for (l in unique(keys)) {
      s <- sub[keys == l, , drop = FALSE]
      s <- s[is.finite(s$neg_log10_p), , drop = FALSE]
      if (nrow(s) < min_snps) next
      k <- which.max(s$neg_log10_p)
      parts <- strsplit(l, "\r", fixed = TRUE)[[1]]
      rows[[length(rows) + 1L]] <- data.frame(
        by = parts[1], contrast = if (has_contrast) parts[2] else NA_character_,
        gene = as.character(g$name[i]),
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
  if (!has_contrast) out$contrast <- NULL
  if (is.null(by)) {
    out$by <- NULL
    ord <- order(-out$max_neg_log10_p)
  } else {
    names(out)[1] <- by
    out[[by]] <- factor(out[[by]], levels = levels(scan[[by]]) %||% unique(out[[by]]))
    ord <- order(out[[by]], -out$max_neg_log10_p)
  }
  out <- out[ord, , drop = FALSE]
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
#' @section Allele count: not a reduction, and not a problem here:
#' `scan_hh()` reports only two allele frequencies whatever a marker carries, which looks like
#' a reduction and is not one for this statistic. Rsb's input is `iNES`, a **site-level**
#' homozygosity computed over every haplotype: at a four-allele marker `scan_hh()`'s `iNES` is
#' identical to [rehh::calc_ehhs()]'s over the full data. No haplotype is dropped and no allele
#' is ignored; only the reported `FREQ_` columns reduce, and this function does not return them.
#'
#' `iNES` normalises by the focal site's own homozygosity, so splitting the same haplotypes
#' into more allele classes moves it by 1% or less. That is why Rsb carries no multiallelic
#' warning while [run_xpehh()] does -- see there.
#'
#' @section Allele count: iES falls as the focal splits, and XP-EHH divides two of them:
#' No haplotype is dropped and no allele ignored -- `scan_hh()`'s `iES` at a four-allele marker
#' is identical to [rehh::calc_ehhs()]'s over the full data, and only the reported `FREQ_`
#' columns reduce. But `iES` is a site-level homozygosity **pooled across allele classes**, so
#' it falls as the focal is split into more of them. Measured on identical haplotypes:
#'
#' | alleles at the focal | `iES` | `iNES` |
#' |---|---|---|
#' | 2 | 991.8 | 2406.7 |
#' | 3 | 631.5 (-36%) | 2388.6 (-1%) |
#' | 4 | 558.5 (-44%) | 2395.8 (0%) |
#'
#' XP-EHH divides one population's `iES` by the other's, so a marker carrying **different**
#' numbers of alleles in the two populations puts them on different footings, and the score
#' partly measures the diversity difference rather than the haplotype-length difference. It
#' warns when that happens, and reports `n_alleles_pop1` / `n_alleles_pop2` so those markers
#' can be dropped. [run_rsb()] reads `iNES` instead and is not affected.
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
