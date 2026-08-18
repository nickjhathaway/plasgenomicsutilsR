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
                                impute = TRUE, seed = 42, meta = NULL, genotype = NULL) {
  meta <- .normalise_meta(meta)
  het <- match.arg(het)
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

  H <- .haploid_calls(G, "missing")
  rec$n_het_calls <- sum(G == 1, na.rm = TRUE)
  if (identical(het, "sample") && rec$n_het_calls > 0) {
    p <- .snp_freqs(H)$p
    hit <- which(!is.na(G) & G == 1, arr.ind = TRUE)
    pr <- p[hit[, "col"]]
    pr[is.na(pr)] <- 0.5
    H[hit] <- stats::rbinom(nrow(hit), 1, pr)
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
  m <- .snp_freqs(H)
  keep_maf <- !is.na(m$p) & pmin(m$p, 1 - m$p) >= maf
  H <- H[, keep_maf, drop = FALSE]
  rec$n_dropped_maf <- n_snp - ncol(H)
  if (!ncol(H)) stop(sprintf("no SNP passed maf >= %.3f", maf), call. = FALSE)

  rec$n_imputed <- sum(is.na(H))
  if (rec$n_imputed > 0) {
    if (impute) {
      p <- .snp_freqs(H)$p
      hit <- which(is.na(H), arr.ind = TRUE)
      H[hit] <- stats::rbinom(nrow(hit), 1, p[hit[, "col"]])
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
                 meta = meta, filtering = rec, seed = seed, het = het),
            class = "parasite_haplotypes")
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

.scan_group <- function(hap, rows, polarized, threads) {
  objs <- .haplohh_list(hap, rows)
  if (!length(objs)) return(NULL)
  scans <- lapply(objs, function(o)
    rehh::scan_hh(o, polarized = polarized, discard_integration_at_border = FALSE,
                  threads = threads))
  do.call(rbind, scans)
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
#' @param min_samples Skip groups smaller than this.
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
#' @seealso [ihs_genes()], [plot_ihs()], [run_rsb()], [beta_score()]
#' @examples
#' ps <- example_pop_structure(umap = FALSE)
#' hap <- parasite_haplotypes(ps, maf = 0.05)
#' run_ihs(hap, group = "country")
#' @export
run_ihs <- function(hap, group = NULL, meta = NULL, polarized = FALSE, freqbin = NULL,
                    min_maf = 0.05, min_samples = 4, threads = 1) {
  meta <- .normalise_meta(meta)
  .need_package("rehh", "run_ihs()")
  stopifnot(inherits(hap, "parasite_haplotypes"))
  freqbin <- .resolve_freqbin(freqbin, polarized)
  rows <- .ihs_rows(hap, group, meta, min_samples)

  out <- list()
  for (l in names(rows)) {
    scan <- .scan_group(hap, rows[[l]], polarized, threads)
    if (is.null(scan)) next
    res <- rehh::ihh2ihs(scan, freqbin = freqbin, min_maf = min_maf, verbose = FALSE)$ihs
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
  if (!length(out)) {
    warning("no group produced an iHS scan", call. = FALSE)
    return(tibble::tibble())
  }
  df <- do.call(rbind, out)
  df$group <- factor(df$group, levels = names(rows))
  tibble::as_tibble(df)
}

.cross_pop <- function(hap, group, meta, pairs, polarized, min_samples, threads, fn,
                       value_col, label) {
  .need_package("rehh", label)
  stopifnot(inherits(hap, "parasite_haplotypes"))
  rows <- .ihs_rows(hap, group, meta, min_samples)
  levs <- names(rows)
  if (length(levs) < 2) stop("need at least two groups", call. = FALSE)
  if (is.null(pairs)) pairs <- utils::combn(levs, 2, simplify = FALSE)
  scans <- lapply(rows, function(r) .scan_group(hap, r, polarized, threads))

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
                    min_samples = 4, threads = 1) {
  meta <- .normalise_meta(meta)
  .cross_pop(hap, group, meta, pairs, polarized, min_samples, threads,
             rehh::ines2rsb, "^RSB", "run_rsb()")
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
                      min_samples = 4, threads = 1) {
  meta <- .normalise_meta(meta)
  .cross_pop(hap, group, meta, pairs, polarized, min_samples, threads,
             rehh::ies2xpehh, "^XPEHH", "run_xpehh()")
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
