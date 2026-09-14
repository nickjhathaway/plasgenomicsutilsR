# Do two populations carry the *same* haplotype at a locus, or independently derived ones?
#
# The IBD-block family cannot answer this. The IBD caller is least able to reach exactly the
# cross-population comparison the question turns on -- distant pairs are fitted with a low pi
# and their shared segments go uncalled -- so an absence of called IBD between two regions is
# as consistent with a dropped call as with independent origins. This test never touches the
# caller: it works off raw allele agreement (IBS), which is defined for every pair.
#
# The statistic is a difference of differences. Within a window, take the mean pairwise
# identity between two cross-population sets -- carrier(A) x carrier(B), and reference(A) x
# reference(B) -- and subtract. Everything that shifts both strata equally cancels: window SNP
# density, the baseline divergence between the two populations, missingness. Reporting only the
# carrier identity is wrong; the signal is the contrast with the reference stratum, which sits
# below its own genome-wide average because non-carriers at a polymorphic locus carry assorted
# haplotypes.
#
# The null is a set of matched windows tiled off the focal chromosome. Windows are matched on
# SNP *count*, not physical width, because IBS is a mean over the SNPs in a window and its
# sampling variance is set by that count; a physical-width match would make the null
# heteroscedastic. A physical-span bound is applied on top, to keep a 53-SNP window that
# happens to span 400 kb from inflating the null variance.
#
# Ported from the prototype locus_ibs_null.py. One deliberate change: within-group cohesion
# (a_car, b_car) drops self-comparisons, which the prototype left in; the headline statistic is
# cross-group and unaffected.

# Mean pairwise IBS between two index sets over the rows of a window matrix `W` (SNPs x
# samples). A and B are column indices. Self-pairs are dropped; a pair needs `min_sites`
# sites where both samples are called before it counts; `exclude_keys` removes named pairs
# (the "min max" index key) from every comparison.
.ibs_group_mean <- function(W, A, B, min_sites, exclude_keys = NULL) {
  tot <- 0; n <- 0L
  for (bj in B) {
    Aset <- A[A != bj]
    if (!is.null(exclude_keys) && length(Aset)) {
      k <- paste(pmin(Aset, bj), pmax(Aset, bj))
      Aset <- Aset[!(k %in% exclude_keys)]
    }
    if (!length(Aset)) next
    Amat <- W[, Aset, drop = FALSE]
    colb <- W[, bj]
    ok <- !is.na(Amat) & !is.na(colb)              # colb recycles down every column
    eq <- ok & (Amat == colb)
    cnt <- colSums(ok)
    keep <- cnt >= min_sites
    if (any(keep)) { tot <- tot + sum(colSums(eq)[keep] / cnt[keep]); n <- n + sum(keep) }
  }
  if (n) tot / n else NA_real_
}

# Standardise the several genotype inputs to one integer matrix (SNPs x samples) plus aligned
# chr/pos and the sample ids. Missing must survive as NA -- an imputed panel silently fills it.
.as_ibs_panel <- function(x, map = NULL) {
  if (inherits(x, "parasite_haplotypes")) {
    G0 <- x$hap
    mp <- x$map
    if (is.null(mp) || !all(c("chr", "pos") %in% names(mp)))
      stop("this parasite_haplotypes has no $map with chr/pos", call. = FALSE)
    if (!is.null(x$filtering) && isTRUE(x$filtering$n_imputed > 0))
      warning("this parasite_haplotypes was imputed (", x$filtering$n_imputed,
              " calls); IBS will treat imputed calls as observed. Build it with ",
              "impute = FALSE, or pass an allele-index genotype panel, to keep missingness.",
              call. = FALSE)
    chr <- normalise_chr(mp$chr); pos <- as.numeric(mp$pos)
    samples <- rownames(G0)
  } else if (is.list(x) && !is.null(x$genotype)) {
    # An allele_set list is the lossless superset: its matrix holds per-cell *state codes*, not
    # allele indices, so read it through the same derivation the panel machinery uses -- a
    # singleton set is the allele the sample carries, a mixed call has no single identity and
    # is NA (which is what IBS wants: a mixture is not "the same allele" as anything). This
    # makes the index that IBS reads identical to a directly loaded allele_index panel.
    if (identical(x$encoding, "allele_set")) {
      if (is.null(x$state_sets))
        stop("this allele_set panel carries no state-set decomposition, so an allele index ",
             "cannot be derived. Rebuild it with `load_genotypes(encoding = \"allele_set\")`.",
             call. = FALSE)
      G0 <- .index_matrix_from_sets(as.matrix(x$genotype), x$state_sets)
    } else {
      G0 <- x$genotype
    }
    ids <- colnames(G0)
    if (is.null(ids) && !is.null(x$snp.id)) ids <- x$snp.id
    cp <- .snp_id_chr_pos(ids, "the genotype panel's column names")
    chr <- cp$chr; pos <- cp$pos
    samples <- rownames(G0)
    if (identical(x$encoding, "dosage"))
      warning("this genotype panel is dosage-encoded (0/1/2 alt copies), which cannot name ",
              "*which* alternate a call carries; IBS at a multiallelic codon needs the ",
              "allele-index encoding (load_genotypes(encoding = \"allele_index\")).",
              call. = FALSE)
  } else if (is.matrix(x)) {
    G0 <- x
    if (is.null(map) || !all(c("chr", "pos") %in% names(map)) || nrow(map) != ncol(G0))
      stop("with a bare genotype matrix, `map` must be a data frame of chr/pos with one row ",
           "per column of the matrix", call. = FALSE)
    chr <- normalise_chr(map$chr); pos <- as.numeric(map$pos)
    samples <- rownames(G0)
  } else {
    stop("`x` must be a parasite_haplotypes, a load_genotypes() list, or a genotype matrix",
         call. = FALSE)
  }
  if (is.null(samples))
    stop("the genotypes have no sample ids on their rows", call. = FALSE)
  storage.mode(G0) <- "integer"
  Gt <- t(G0)                                       # SNPs x samples
  colnames(Gt) <- .as_id_chr(samples)
  list(Gt = Gt, chr = chr, pos = pos)
}

# Consecutive non-overlapping blocks of `n_snps` SNPs along one chromosome, in position order.
.tile_windows <- function(idx, pos, n_snps) {
  idx <- idx[order(pos[idx])]
  if (length(idx) < n_snps) return(list())
  starts <- seq(1L, length(idx) - n_snps + 1L, by = n_snps)
  lapply(starts, function(s) idx[s:(s + n_snps - 1L)])
}

#' Do two populations carry the same haplotype at a locus? An allele-split IBS test
#'
#' Measures whether the carriers of a variant in two populations are drawn from a common
#' haplotype pool, using pairwise identity (IBS) rather than called IBD, so it never depends
#' on the caller reaching the cross-population pairs the question turns on. The statistic is a
#' difference of differences -- cross-population carrier-vs-carrier identity minus
#' reference-vs-reference identity -- scored against a null of SNP-count-matched windows tiled
#' off the focal chromosome. It answers "one haplotype spreading across the boundary" versus
#' "independent origins", which the IBD block-extension family (the iHS-shaped, within-
#' population question; see [ibd_block_extension_by_allele()]) cannot settle on its own.
#'
#' @section What it does and does not say:
#' A positive, significant result means the carrier chromosomes in the two populations are
#' drawn from a common pool rather than independently derived. It does **not** count how many
#' backgrounds that pool holds; pair it with a tabulation of the other coding variants each
#' carrier also carries. A negative difference means the carriers are *less* alike across the
#' boundary than the reference samples are -- independent backgrounds -- even where each
#' population's carriers cluster tightly on their own (`a_car`, `b_car`).
#'
#' @section Matched windows:
#' Windows are matched on **SNP count** (`n_snps`), not physical width, because IBS is a mean
#' over a window's SNPs and its variance is set by that count. A physical-span bound
#' (`span_tol`) is layered on top to stop a count-matched window that happens to span far more
#' bp from inflating the null variance. The focal chromosome is excluded from the null so the
#' locus cannot seed its own baseline, and a pair needs `min_sites` comparable sites in a
#' window to contribute. With a few hundred windows the smallest achievable two-sided p is a
#' resolution floor, so read a p at the floor as "no matched window reached the observed
#' value" rather than as an exact probability.
#'
#' @param x Genotypes: a [parasite_haplotypes] object, a `load_genotypes()` list (its
#'   `allele-index` encoding, which keeps missingness), or a bare integer matrix
#'   (samples x SNPs) with `map`. Missing calls must be `NA`; an imputed panel is warned
#'   about, since imputation would be read as agreement. For IBS and IBD to be comparable the
#'   calls should be the ones the IBD run used (hmmibd-rs dominant-allele from `FORMAT/AD`);
#'   the allele-index panel is `GT`-derived and agrees with it wherever calls are unmixed.
#' @param locus The focal locus: a one-row data frame or list with `chr`, `start`, `end`
#'   (0-based half-open, used for the `drop_ibd` overlap) and optional `centre` (the focal
#'   window is centred here; defaults to the interval midpoint). A multi-row interval frame is
#'   collapsed to its span.
#' @param allele The variant to split on: a metadata column name or a named `sample -> state`
#'   vector (as [allele_states()] returns).
#' @param carrier,reference State sets (matched with `%in%`), as in
#'   [carrier_genome_wide_relatedness()].
#' @param group_a,group_b The two populations, each a set of `group` values. The contrast is
#'   between them; a sample in neither is unused.
#' @param group Metadata column defining the populations (default `"region"`, or the object's
#'   own group column when it carries one).
#' @param meta Sample metadata (`sample` plus `group` and the `allele` column). Taken from a
#'   `parasite_haplotypes`' `$meta` when not given.
#' @param n_snps Window size in SNPs (default `53`).
#' @param span_tol Keep null windows whose physical span is within this factor of the focal
#'   window's span, either side (default `2`; `NULL` or `0` to not bound).
#' @param exclude_focal_chr Drop the focal chromosome from the null (default `TRUE`).
#' @param min_sites Comparable sites a pair needs in a window to contribute (default `10`).
#' @param drop_ibd Optionally also report the statistic with pairs already IBD across the
#'   locus removed from both strata, which separates "a lineage moved across the border" from
#'   "the carrier populations share a background". Pass an [IbdResults] (its blocks name the
#'   IBD pairs), or `TRUE` to use `x` if it is one. Adds a second result row.
#' @return A tibble, one row per pair set (`"all pairs"`, and `"IBD pairs dropped"` when
#'   `drop_ibd` is given), carrying the per-window null on the `null` attribute:
#'   \describe{
#'     \item{`locus`, `pairs`}{the locus label and which pair set the row is.}
#'     \item{`n_a_car`, `n_b_car`, `n_a_ref`, `n_b_ref`}{samples in each population x stratum.}
#'     \item{`cross_car`, `cross_ref`}{cross-population mean IBS in each stratum.}
#'     \item{`a_car`, `b_car`}{within-population carrier cohesion, each population on its own.}
#'     \item{`d_ind`}{`cross_car - cross_ref`, the statistic.}
#'     \item{`null_n`, `null_mean`, `null_sd`}{the matched-window null.}
#'     \item{`z`}{`(d_ind - null_mean) / null_sd`.}
#'     \item{`p_two`}{two-sided empirical p against the null (a floor when small).}
#'     \item{`focal_span_kb`, `span_tol`, `n_snps`}{the window the null was matched to.}
#'   }
#' @seealso [ibd_block_extension_by_allele()] for the within-population iHS-shaped contrast,
#'   [carrier_genome_wide_relatedness()] for the genome-wide relatedness control,
#'   [parasite_haplotypes()] and `load_genotypes()` for the genotype inputs.
#' @export
locus_ibs_by_allele <- function(x, locus, allele, carrier = NULL, reference = NULL,
                                group_a, group_b, group = NULL, meta = NULL,
                                n_snps = 53L, span_tol = 2, exclude_focal_chr = TRUE,
                                min_sites = 10L, drop_ibd = NULL, map = NULL) {
  if (missing(allele) || is.null(allele)) stop("`allele` is required", call. = FALSE)
  if (missing(group_a) || missing(group_b)) stop("`group_a` and `group_b` are required",
                                                 call. = FALSE)
  n_snps <- as.integer(n_snps)
  if (is.na(n_snps) || n_snps < 2L) stop("`n_snps` must be at least 2", call. = FALSE)

  # ---- genotypes -------------------------------------------------------------------
  panel <- .as_ibs_panel(x, map = map)
  Gt <- panel$Gt; chr <- panel$chr; pos <- panel$pos

  if (is.null(meta)) {
    meta <- if (inherits(x, "parasite_haplotypes")) x$meta
            else if (is.list(x) && !is.null(x$meta)) x$meta else NULL
  }
  meta <- .normalise_meta(meta)
  if (is.null(meta) || !"sample" %in% names(meta))
    stop("`meta` (with a `sample` column) is needed to group and split the samples",
         call. = FALSE)
  if (is.null(group)) group <- if (inherits(x, "IbdResults")) x$get_group_col() else NULL
  if (is.null(group)) group <- if ("region" %in% names(meta)) "region" else
    setdiff(names(meta), "sample")[1]
  if (is.na(group) || !group %in% names(meta))
    stop("meta has no column '", group, "'.\n  columns available: ",
         paste(setdiff(names(meta), "sample"), collapse = ", "), call. = FALSE)

  # ---- locus -----------------------------------------------------------------------
  L <- as.data.frame(locus, stringsAsFactors = FALSE)
  lc <- intersect(c("chr", "chrom"), names(L))[1]
  if (is.na(lc) || !all(c("start", "end") %in% names(L)))
    stop("`locus` needs chr, start and end", call. = FALSE)
  focal_chr <- normalise_chr(L[[lc]][1])
  lstart <- min(as.numeric(L$start)); lend <- max(as.numeric(L$end))
  centre <- if ("centre" %in% names(L)) as.numeric(L$centre[1]) else (lstart + lend) / 2

  # ---- carrier / reference states and the four groups ------------------------------
  m <- meta[!duplicated(meta$sample), , drop = FALSE]
  ss <- .allele_states(allele, m, m$sample)
  st <- ss$state; lv <- ss$levels
  if (is.null(carrier) && is.null(reference)) {
    if (length(lv) != 2)
      stop("`allele` has ", length(lv), " observed states; name `carrier =` and ",
           "`reference =`", call. = FALSE)
    reference <- lv[1]; carrier <- lv[2]
    message("reading `", reference, "` as reference and `", carrier, "` as carrier")
  } else if (is.null(reference) || is.null(carrier)) {
    if (length(lv) > 2) stop("`allele` has ", length(lv), " states; name both `carrier =` ",
                             "and `reference =`", call. = FALSE)
    if (is.null(reference)) reference <- setdiff(lv, carrier) else carrier <- setdiff(lv, reference)
  }
  carrier <- as.character(carrier); reference <- as.character(reference)
  group_a <- as.character(group_a); group_b <- as.character(group_b)
  if (length(intersect(group_a, group_b)))
    stop("`group_a` and `group_b` overlap; the contrast is between two disjoint populations",
         call. = FALSE)

  gp <- stats::setNames(as.character(m[[group]]), m$sample)
  # column index in Gt for each sample carrying a given (group set, state set)
  col_of <- stats::setNames(seq_len(ncol(Gt)), colnames(Gt))
  pick <- function(gset, sset) {
    s <- names(st)[!is.na(st) & st %in% sset & !is.na(gp[names(st)]) & gp[names(st)] %in% gset]
    unname(col_of[intersect(s, colnames(Gt))])
  }
  Ac <- pick(group_a, carrier); Ar <- pick(group_a, reference)
  Bc <- pick(group_b, carrier); Br <- pick(group_b, reference)
  if (!length(Ac) || !length(Bc))
    stop("one of the carrier groups is empty (A: ", length(Ac), ", B: ", length(Bc),
         "); nothing to contrast", call. = FALSE)
  if (length(Ac) < 2 && length(Bc) < 2)
    warning("very few carriers (A: ", length(Ac), ", B: ", length(Bc),
            "); the statistic rests on a handful of pairs", call. = FALSE)

  # ---- drop_ibd: pairs IBD across the locus, as excluded index keys -----------------
  exclude_keys <- NULL
  if (!is.null(drop_ibd) && !isFALSE(drop_ibd)) {
    ibd_obj <- if (isTRUE(drop_ibd)) x else drop_ibd
    if (!inherits(ibd_obj, "IbdResults"))
      stop("`drop_ibd` must be an IbdResults (or TRUE when `x` is one)", call. = FALSE)
    bl <- ibd_obj$get_blocks()
    if (is.null(bl) || !nrow(bl)) stop("`drop_ibd` IbdResults carries no blocks", call. = FALSE)
    bl <- bl[normalise_chr(bl$chr) == focal_chr &
             as.numeric(bl$start) < lend & as.numeric(bl$end) > lstart, , drop = FALSE]
    if (nrow(bl)) {
      i <- col_of[.as_id_chr(bl$sample1)]; j <- col_of[.as_id_chr(bl$sample2)]
      ok <- !is.na(i) & !is.na(j)
      exclude_keys <- unique(paste(pmin(i[ok], j[ok]), pmax(i[ok], j[ok])))
    } else exclude_keys <- character(0)
  }

  # ---- focal window and the matched null -------------------------------------------
  ci <- which(chr == focal_chr)
  if (length(ci) < n_snps)
    stop("the focal chromosome ", focal_chr, " carries only ", length(ci),
         " SNPs, fewer than n_snps = ", n_snps, call. = FALSE)
  focal <- ci[order(abs(pos[ci] - centre))[seq_len(n_snps)]]
  focal <- focal[order(pos[focal])]
  fspan <- max(pos[focal]) - min(pos[focal])

  stat_at <- function(rows, ex) {
    W <- Gt[rows, , drop = FALSE]
    list(cc = .ibs_group_mean(W, Ac, Bc, min_sites, ex),
         cr = .ibs_group_mean(W, Ar, Br, min_sites, ex),
         ac = .ibs_group_mean(W, Ac, Ac, min_sites, ex),
         bc = .ibs_group_mean(W, Bc, Bc, min_sites, ex))
  }

  null_windows <- list()
  for (c in unique(chr)) {
    if (isTRUE(exclude_focal_chr) && c == focal_chr) next
    ww <- .tile_windows(which(chr == c), pos, n_snps)
    for (w in ww) {
      sp <- max(pos[w]) - min(pos[w])
      if (!is.null(span_tol) && span_tol && !(fspan / span_tol <= sp && sp <= fspan * span_tol))
        next
      null_windows[[length(null_windows) + 1L]] <- w
    }
  }
  lab <- L$name[1]; if (is.null(lab) || is.na(lab)) lab <- paste0(focal_chr, ":", round(centre))

  one_pass <- function(ex, pairs_label) {
    obs <- stat_at(focal, ex)
    d_obs <- obs$cc - obs$cr
    nv <- vapply(null_windows, function(w) {
      s <- stat_at(w, ex); s$cc - s$cr
    }, numeric(1))
    v <- nv[is.finite(nv)]
    if (length(v)) {
      pg <- (1 + sum(v >= d_obs)) / (1 + length(v))
      pl <- (1 + sum(v <= d_obs)) / (1 + length(v))
      p_two <- min(1, 2 * min(pg, pl))
      z <- (d_obs - mean(v)) / stats::sd(v)
    } else { p_two <- NA_real_; z <- NA_real_ }
    row <- tibble::tibble(
      locus = lab, pairs = pairs_label,
      n_a_car = length(Ac), n_b_car = length(Bc),
      n_a_ref = length(Ar), n_b_ref = length(Br),
      cross_car = obs$cc, cross_ref = obs$cr, a_car = obs$ac, b_car = obs$bc,
      d_ind = d_obs, null_n = length(v),
      null_mean = if (length(v)) mean(v) else NA_real_,
      null_sd = if (length(v)) stats::sd(v) else NA_real_,
      z = z, p_two = p_two,
      focal_span_kb = fspan / 1000, span_tol = if (is.null(span_tol)) NA_real_ else span_tol,
      n_snps = n_snps)
    attr(row, "null") <- v
    row
  }

  passes <- list(one_pass(NULL, "all pairs"))
  if (!is.null(exclude_keys)) passes[[2]] <- one_pass(exclude_keys, "IBD pairs dropped")
  res <- dplyr::bind_rows(passes)
  attr(res, "null") <- lapply(passes, function(p) attr(p, "null"))
  attr(res, "states") <- list(carrier = carrier, reference = reference)
  attr(res, "groups") <- list(group_a = group_a, group_b = group_b)
  res
}
