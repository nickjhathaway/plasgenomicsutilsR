# Is an IBD segment spanning a locus longer than the pairs sharing it usually carry?
#
# The naive form of that question -- locus median block length over the cohort's genome-wide
# median -- cannot be compared across groups, because a group with less outcrossing carries
# longer blocks *everywhere*, and because the two medians come from different pair sets:
# pairs IBD at a locus are selected for being related enough to share it. The paired form
# makes each pair its own control, which both cancels background relatedness and makes the
# pair the unit of observation, which is what makes a signed-rank test valid.
#
# This is the *length* statistic. XiR,s (Henden et al. 2018, isoRelate) is the *count*
# statistic over the same blocks and is generally better powered; the two are complements
# and the count should usually lead. See `?ibd_block_extension_test` for citations.

# Loci to test: gene names / the object's track go through the gene-track path, while a
# supplied interval frame is merged into one row per feature -- rows sharing a name, a
# chromosome and (where the frame has one) a gene_id become a single interval spanning them,
# so an uncollapsed `aa_intervals()` table for a marker set is one locus. A gene track passed
# as a data frame stays one locus per gene, including where a name repeats across
# chromosomes (paralogues), which is why the chromosome and gene_id are part of the key.
.locus_track_for <- function(x, loci) {
  if (is.null(loci) || is.character(loci)) {
    tr <- .gene_track_for(x, loci)
    tr$.locus <- tr$.label
    return(tr)
  }
  tr <- .as_gene_track(loci)
  if (is.null(tr) || !nrow(tr)) stop("`loci` is empty", call. = FALSE)
  tr$name <- as.character(tr$name)
  gid <- if ("gene_id" %in% names(tr)) as.character(tr$gene_id) else rep(NA_character_, nrow(tr))
  key <- paste(tr$name, tr$chr, ifelse(is.na(gid), "", gid), sep = "\r")
  parts <- split(seq_len(nrow(tr)), factor(key, levels = unique(key)))  # input order
  take <- function(f, mode) unname(vapply(parts, f, mode))
  out <- data.frame(
    name = take(function(ix) tr$name[ix[1]], character(1)),
    chr = take(function(ix) tr$chr[ix[1]], character(1)),
    start = take(function(ix) min(as.numeric(tr$start[ix])), numeric(1)),
    end = take(function(ix) max(as.numeric(tr$end[ix])), numeric(1)),
    gene_id = take(function(ix) gid[ix[1]], character(1)),
    n_parts = take(function(ix) length(ix), integer(1)),
    stringsAsFactors = FALSE)
  out$.locus <- .disambiguate_gene_labels(out$name, out$gene_id)
  if (!identical(out$.locus, out$name)) {
    warning("some locus names repeat across features; disambiguated by gene_id in the ",
            "`locus` column (`name`/`gene_id` kept as columns)", call. = FALSE)
  }
  out
}

# "A\rB" -> group_a / group_b, in place of `group` and in the same column position.
.split_group_pair <- function(df) {
  parts <- if (!nrow(df)) matrix(character(0), 0, 2) else
    do.call(rbind, strsplit(as.character(df$group), "\r", fixed = TRUE))
  i <- match("group", names(df))
  before <- if (i > 1) names(df)[seq_len(i - 1)] else character(0)
  after <- setdiff(names(df), c(before, "group"))
  out <- df
  out$group <- NULL
  out$group_a <- parts[, 1]
  out$group_b <- parts[, 2]
  out[, c(before, "group_a", "group_b", after), drop = FALSE]
}

# Wilcoxon signed-rank against 0, made safe for a genome-wide scan: ties fall back to the
# normal approximation, since thousands of loci must not raise thousands of warnings. A
# stratum of nothing but exact zeros is answered here rather than left to the R version --
# the signed-rank test drops zeros, so recent R reports p = 1 where older R errors, and 1
# (no evidence the locus is longer) is the answer either way.
.signed_rank_p <- function(v, alternative) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  if (all(v == 0)) return(1)
  p <- tryCatch(
    suppressWarnings(stats::wilcox.test(v, mu = 0, alternative = alternative)$p.value),
    error = function(e) NA_real_)
  if (is.null(p)) NA_real_ else as.numeric(p)
}

# One row per grouping of the per-pair log ratios. Shared by the locus test, the window scan
# and the carrier split, so the three cannot drift apart in what `paired_ratio` means.
.summarise_ratios <- function(ratios, by, min_pairs) {
  g <- dplyr::group_by(ratios, dplyr::across(dplyr::all_of(by)))
  dplyr::summarise(
    dplyr::filter(g, dplyr::n() >= min_pairs),
    n_pairs = dplyr::n(),
    locus_median = stats::median(.data$locus_len),
    ref_median = stats::median(.data$ref_median),
    paired_ratio = 2 ^ stats::median(.data$log2_ratio),
    z_paired = stats::median(.data$log2_ratio) /
      (stats::sd(.data$log2_ratio) / sqrt(dplyr::n())),
    p_paired = .signed_rank_p(.data$log2_ratio, "greater"),
    p_two_sided = .signed_rank_p(.data$log2_ratio, "two.sided"),
    .groups = "drop")
}

#' Test whether IBD blocks are longer at a locus than the sharing pairs' own background
#'
#' For each locus and each group of samples, compares the length of the IBD segment a pair
#' shares **across the locus** with the length of that same pair's segments elsewhere in the
#' genome. Long segments at a locus are the footprint of a recent sweep, but raw length is
#' not comparable between groups: a group with less outcrossing carries longer segments
#' genome-wide, so the naive comparison can run backwards. Making each pair its own control
#' removes that, and makes the pair -- not the segment -- the unit of observation, which is
#' what the signed-rank test needs.
#'
#' Two statistics come back, and only one of them is testable:
#'
#' * `naive_ratio` is the locus median segment length over the group's genome-wide median.
#'   Numerator and denominator come from different sets of pairs, because pairs IBD at a
#'   locus are selected for being related enough to share it. There is no error model behind
#'   it, so it gets no p-value. It is reported to be compared against the next one.
#' * `paired_ratio` is `2^median(log2(locus segment / that same pair's median segment on
#'   every other chromosome))`. Background relatedness cancels pair by pair. **This is the
#'   one to quote.** `p_paired` tests it.
#'
#' Where a locus's `ref_median` sits well above `gw_median`, the pairs sharing it are
#' unusually related and `naive_ratio` was counting that background rather than anything
#' local.
#'
#' Each pair's baseline excludes the **whole locus chromosome**, so a sweep's own footprint
#' cannot leak into the reference and flatten the ratio. Segments come from
#' `x$get_blocks()`, which has already dropped non-IBD rows and applied the short-segment
#' floor; taking a filtered numerator against an unfiltered denominator inflates these
#' ratios several-fold, so do not re-read the raw `hmm.txt` for one side of it.
#'
#' Only pairs whose two samples fall in the same group are used, since a cross-group pair
#' has no single group to report under.
#'
#' This measures segment *length*. The count analogue, XiR,s (Henden et al. 2018,
#' \doi{10.1371/journal.pgen.1007279}), is computed by the companion Python pipeline and is
#' generally the better-powered statistic; treat this as its complement rather than its
#' replacement. For background on excess IBD as a selection signal see Albrechtsen et al.
#' (2010, \doi{10.1534/genetics.110.113977}); for the principled coalescent-time version of
#' what this ratio approximates, ASMC (Palamara et al. 2018,
#' \doi{10.1038/s41588-018-0177-x}).
#'
#' @param x An [IbdResults] built with `blocks =` and `meta =`.
#' @param loci Loci to test: gene names from the object's track, `NULL` for every gene in
#'   it, or an interval data frame (`name`, `chr`/`chrom`, `start`, `end`, optionally
#'   `gene_id`). Rows of an interval frame that share a name, a chromosome and a `gene_id`
#'   are merged into one interval spanning them, so an `aa_intervals()` table for a marker
#'   set can be passed as-is once it carries a `name`, while a gene track passed this way
#'   stays one locus per gene even where a name repeats across chromosomes.
#' @param group Metadata column defining the groups. Defaults to the object's declared group
#'   column, then to the first non-`sample` column of `meta`.
#' @param within Pad each locus by this many bp on both sides when deciding whether a
#'   segment counts (default `0`). It applies under either `sharing`, so with padding
#'   `"complete"` asks the segment to cover the padded interval.
#' @param sharing What a pair's segment must do at the locus to be counted, the same choice
#'   [plot_ibd_network()] offers:
#'   \describe{
#'     \item{`"overlap"`}{(default) the segment touches the locus anywhere, which is
#'       [gene_ibd_pairs()]'s rule: `start < end_locus & end > start_locus`. The pair shares
#'       *some* of the locus.}
#'     \item{`"complete"`}{the segment spans the whole locus. Far stricter, and it changes
#'       what the statistic measures -- see the note below on why locus width then matters.}
#'   }
#' @param min_ref_blocks Segments a pair needs off the locus chromosome before its baseline
#'   is usable (default `1`). One is deliberate rather than lax: a pair with a single
#'   segment genome-wide is a distant pair, and that segment is length-biased long because
#'   short ones fell under the floor, so a thin baseline reads long and pushes the ratio
#'   *down*. Raising this is conservative in the wrong direction -- it drops whole groups
#'   while rarely moving a conclusion.
#' @param min_pairs Pairs a group needs at a locus before it is tested (default `5`).
#'   Thinner strata are dropped rather than reported untested.
#' @param adjust Scope of the Benjamini-Hochberg correction, which is a real choice and not
#'   a detail. `"per_locus"` (default) makes one family of the groups tested at each locus,
#'   which is right for a handful of hand-picked loci. `"per_group"` makes one family of the
#'   loci tested in each group, which is what a genome-wide scan wants. `"all"` is a single
#'   family over every row, `"none"` leaves `q_paired` as `NA`.
#' @param pairs Which sample pairs to count. `"within"` (default) uses only pairs whose two
#'   samples share a group and reports a `group` column. `"between"` uses only cross-group
#'   pairs and `"all"` uses both; those two report `group_a` and `group_b` **in place of**
#'   `group`, one row per unordered pair of groups, which is the shape
#'   [plot_pairwise_block_extension()] draws. Nothing about the statistic changes -- a
#'   cross-group pair is still measured against its own segments off the locus chromosome --
#'   so the diagonal of an `"all"` result is exactly the `"within"` result.
#' @param keep_pair_ratios Attach the per-pair log ratios as the `pair_ratios` attribute
#'   (default `TRUE`). Set `FALSE` for a genome-wide scan: on one real cohort that attribute
#'   was 288 MB of a 330 MB result, and nothing downstream of a scan reads it.
#' @param meta Sample metadata; taken from `x` when not given.
#' @return A tibble, one row per locus x group, ordered by locus then descending
#'   `paired_ratio`, carrying the per-pair ratios it was built from on the `pair_ratios`
#'   attribute (`locus`, `group`, `sample1`, `sample2`, `locus_len`, `ref_n`, `ref_median`,
#'   `log2_ratio`):
#'   \describe{
#'     \item{`locus`, `name`, `gene_id`}{locus labels, as [gene_ibd_overlap()] uses them.}
#'     \item{`chr`, `start`, `end`, `span_bp`}{the interval tested, 0-based half-open.
#'       Under `sharing = "complete"`, `span_bp` is a floor on `locus_median`, so compare
#'       loci of like width.}
#'     \item{`group`}{the group both samples of every counted pair belong to.}
#'     \item{`n_pairs`}{pairs IBD across the locus with a usable baseline.}
#'     \item{`gw_median`}{that group's genome-wide median segment length, over all its
#'       within-group pairs.}
#'     \item{`locus_median`}{median segment length across the locus.}
#'     \item{`naive_ratio`}{`locus_median / gw_median`. Descriptive only.}
#'     \item{`ref_median`}{median over the counted pairs of each pair's own baseline. Well
#'       above `gw_median` means the sharing pairs are unusually related.}
#'     \item{`paired_ratio`}{the statistic. `1` means a pair's segment at the locus is no
#'       longer than its segments elsewhere.}
#'     \item{`z_paired`}{`median(log2_ratio)` over the standard error of their mean. A rough
#'       effect size for ranking loci, not the test.}
#'     \item{`p_paired`, `p_two_sided`}{one-sided (`greater`) and two-sided signed-rank
#'       p-values. Read the two-sided column before calling anything shorter: a one-sided p
#'       near `1` means "not longer", never "significantly shorter". Ties use the normal
#'       approximation.}
#'     \item{`q_paired`}{`p_paired` adjusted within the family set by `adjust`.}
#'   }
#' @section Reading it:
#' Do not run the test over segments instead of pairs. One pair contributes several
#' correlated segments and the p-value comes out far too small. Splitting the pairs by
#' carriage of a core variant, which is the IBD analogue of iHS (see [run_ihs()]), turns a
#' statement about a group into one about an allele.
#'
#' A `paired_ratio` of 1 is **not** the null. A locus is a fixed target, and a segment
#' covers it with probability rising in the segment's own length, so the numerator is drawn
#' from a length-biased distribution while each pair's baseline is not. Every locus is
#' inflated by this, whether or not anything happened there. Read a locus against the
#' distribution from scanning every gene in the same group, not against 1.
#'
#' @section Why `sharing` changes what locus width means:
#' Under `"overlap"` a segment counts if it touches the locus, so the selection weight rises
#' with segment length *plus* locus width. Since a gene is small against a typical IBD
#' segment, every gene behaves like a point and takes the same length bias: on one real
#' cohort, locus width explained an R-squared below 0.002 of `log2(paired_ratio)`, and
#' binning the scan by width only added noise.
#'
#' Under `"complete"` the segment must cover the locus, so it cannot be shorter than the
#' locus is wide. Locus width becomes a hard floor under `locus_median`, and the wider the
#' locus the more of the short-segment end of the distribution is cut away. The ratio then
#' rises with width for reasons that have nothing to do with selection, and loci of
#' different widths are no longer on one scale. So with `"complete"`, compare like widths --
#' bin the scan by `span_bp` -- or stay with `"overlap"` when ranking genes against each
#' other.
#' @seealso [gene_ibd_pairs()] for the pairs themselves, [gene_ibd_overlap()] for the
#'   fraction of pairs sharing at a locus (the count side of the same question),
#'   [plasgenomicsutilsR-coordinates] for the interval convention.
#' @examples
#' # five related pairs, each sharing 80 kb across the locus and 20 kb elsewhere
#' pairs <- data.frame(sample1 = paste0("s", seq(1, 9, by = 2)),
#'                     sample2 = paste0("s", seq(2, 10, by = 2)))
#' blocks <- rbind(
#'   data.frame(pairs, chr = "7", start = 3.6e5, end = 4.4e5),   # over the locus
#'   data.frame(pairs, chr = "1", start = 1e5,   end = 1.2e5),   # baseline
#'   data.frame(pairs, chr = "2", start = 1e5,   end = 1.2e5))
#' blocks$different <- 0
#' blocks$Nsnp <- 40
#' meta <- data.frame(sample = paste0("s", 1:10), region = "north")
#' ibd <- ibd_results(blocks = blocks, meta = meta, group_col_in_meta = "region")
#' locus <- data.frame(name = "sweep", chr = "7", start = 4e5, end = 4.1e5)
#' ibd_block_extension_test(ibd, locus)[, c("locus", "group", "n_pairs", "paired_ratio")]
#' @export
ibd_block_extension_test <- function(x, loci = NULL, group = NULL, within = 0,
                                     sharing = c("overlap", "complete"),
                                     pairs = c("within", "all", "between"),
                                     min_ref_blocks = 1L, min_pairs = 5L,
                                     adjust = c("per_locus", "per_group", "all", "none"),
                                     keep_pair_ratios = TRUE, meta = NULL) {
  adjust_given <- !missing(adjust)
  adjust <- match.arg(adjust)
  sharing <- match.arg(sharing)
  pairs <- match.arg(pairs)
  if (!inherits(x, "IbdResults"))
    stop("`x` must be an IbdResults; build one with ibd_results(blocks = , meta = )",
         call. = FALSE)

  blocks <- x$get_blocks()
  if (is.null(blocks) || !nrow(blocks))
    stop("this IbdResults has no IBD blocks; build it with ibd_results(blocks = , meta = )",
         call. = FALSE)
  if (is.null(meta)) meta <- x$get_meta()
  meta <- .normalise_meta(meta)
  if (is.null(meta) || !"sample" %in% names(meta))
    stop("grouping the pairs needs meta with a 'sample' column; pass ibd_results(meta = )",
         call. = FALSE)
  if (is.null(group)) group <- x$get_group_col()
  if (is.null(group)) group <- setdiff(names(meta), "sample")[1]
  if (is.na(group) || !group %in% names(meta))
    stop("meta has no column '", group, "'.\n  columns available: ",
         paste(setdiff(names(meta), "sample"), collapse = ", "), call. = FALSE)

  tr <- .locus_track_for(x, loci)
  if (adjust == "per_locus" && !adjust_given && nrow(tr) > 25) {
    message(nrow(tr), " loci with adjust = \"per_locus\": each locus is its own family, so ",
            "`q_paired` does not correct for scanning them all. A genome-wide scan wants ",
            "adjust = \"per_group\".")
  }

  # ---- within-group pairs only -----------------------------------------------------
  s2g <- stats::setNames(as.character(meta[[group]]), as.character(meta$sample))
  b <- blocks
  b$chr <- normalise_chr(b$chr)
  bs1 <- as.character(b$sample1)
  bs2 <- as.character(b$sample2)
  g1 <- unname(s2g[bs1])
  g2 <- unname(s2g[bs2])
  n_ungrouped <- sum(is.na(g1) | is.na(g2))
  if (n_ungrouped)
    message(n_ungrouped, " of ", nrow(b), " IBD segments touch a sample with no `", group,
            "`, and are left out")
  known <- !is.na(g1) & !is.na(g2)
  same <- known & g1 == g2
  keep <- switch(pairs, within = same, between = known & !same, all = known)
  if (!any(keep))
    stop("no IBD segment joins two samples matching pairs = \"", pairs, "\"", call. = FALSE)
  b <- b[keep, , drop = FALSE]
  # For a cross-group pair the statistic is unchanged -- the baseline is still that pair's own
  # segments off the locus chromosome -- so the only thing that differs is the label it gets
  # reported under: the unordered pair of the two ends' groups instead of a single group.
  blk_group <- if (pairs == "within") g1[keep] else
    paste(pmin(g1[keep], g2[keep]), pmax(g1[keep], g2[keep]), sep = "\r")
  lo1 <- pmin(bs1[keep], bs2[keep])
  hi1 <- pmax(bs1[keep], bs2[keep])
  pair_key <- paste(lo1, hi1, sep = "\r")
  pair_lvl <- unique(pair_key)
  pair_ix <- match(pair_key, pair_lvl)
  first_of_pair <- match(pair_lvl, pair_key)
  pair_group <- blk_group[first_of_pair]
  pair_s1 <- lo1[first_of_pair]
  pair_s2 <- hi1[first_of_pair]
  n_pair_tot <- length(pair_lvl)
  block_len <- as.numeric(b$end) - as.numeric(b$start)

  # ---- each group's genome-wide background (the naive denominator) -------------------
  gw <- dplyr::summarise(
    dplyr::group_by(dplyr::tibble(group = blk_group, len = block_len), .data$group),
    gw_n_blocks = dplyr::n(), gw_median = stats::median(.data$len), .groups = "drop")

  # ---- per-pair baselines, one per excluded chromosome, built once ------------------
  # The exclusion is always the whole locus chromosome, so there are only as many distinct
  # baselines per pair as there are chromosomes carrying loci -- never one per locus. This
  # is what keeps a genome-wide scan from re-deriving the same baseline thousands of times.
  chr_needed <- unique(tr$chr)
  ref_med <- matrix(NA_real_, n_pair_tot, length(chr_needed), dimnames = list(NULL, chr_needed))
  ref_n <- matrix(0L, n_pair_tot, length(chr_needed), dimnames = list(NULL, chr_needed))
  for (j in seq_along(chr_needed)) {
    off <- b$chr != chr_needed[j]
    if (!any(off)) next
    s <- dplyr::summarise(
      dplyr::group_by(dplyr::tibble(p = pair_ix[off], len = block_len[off]), .data$p),
      n = dplyr::n(), m = stats::median(.data$len), .groups = "drop")
    ref_n[s$p, j] <- as.integer(s$n)
    ref_med[s$p, j] <- s$m
  }

  # ---- one pass per locus ----------------------------------------------------------
  by_chr <- split(seq_len(nrow(b)), b$chr)
  ratios <- vector("list", nrow(tr))
  for (i in seq_len(nrow(tr))) {
    idx <- by_chr[[tr$chr[i]]]
    if (is.null(idx) || !length(idx)) next
    ls <- as.numeric(tr$start[i]) - within
    le <- as.numeric(tr$end[i]) + within
    m <- if (sharing == "complete") {
      b$start[idx] <= ls & b$end[idx] >= le           # the segment spans the whole locus
    } else {
      b$start[idx] < le & b$end[idx] > ls             # gene_ibd_pairs()'s overlap rule
    }
    if (!any(m)) next
    ii <- idx[m]
    pk <- pair_ix[ii]
    ln <- block_len[ii]
    # one row per pair: the longest of its segments over the locus, since an interval is
    # occasionally split across two and the longer one is the relevant haplotype
    o <- order(pk, -ln)
    pk <- pk[o]; ln <- ln[o]
    f <- !duplicated(pk)
    up <- pk[f]; ul <- ln[f]
    j <- match(tr$chr[i], chr_needed)
    rn <- ref_n[up, j]
    rm <- ref_med[up, j]
    ok <- is.finite(rm) & rm > 0 & rn >= min_ref_blocks
    if (!any(ok)) next
    up <- up[ok]
    ratios[[i]] <- dplyr::tibble(
      locus = tr$.locus[i], group = pair_group[up],
      sample1 = pair_s1[up], sample2 = pair_s2[up],
      locus_len = ul[ok], ref_n = rn[ok], ref_median = rm[ok],
      log2_ratio = log2(ul[ok] / rm[ok]))
  }
  ratios <- ratios[!vapply(ratios, is.null, logical(1))]
  ratios <- if (length(ratios)) dplyr::bind_rows(ratios) else
    dplyr::tibble(locus = character(), group = character(), sample1 = character(),
                  sample2 = character(), locus_len = numeric(), ref_n = integer(),
                  ref_median = numeric(), log2_ratio = numeric())

  res <- .summarise_ratios(ratios, c("locus", "group"), min_pairs)
  res <- dplyr::left_join(res, gw, by = "group")
  res <- dplyr::mutate(res, naive_ratio = .data$locus_median / .data$gw_median)

  # BH within whichever family `adjust` names -- the scope is the choice, not the method
  res$q_paired <- if (!nrow(res)) numeric(0) else switch(
    adjust,
    per_locus = stats::ave(res$p_paired, res$locus, FUN = function(p) stats::p.adjust(p, "BH")),
    per_group = stats::ave(res$p_paired, res$group, FUN = function(p) stats::p.adjust(p, "BH")),
    all = stats::p.adjust(res$p_paired, "BH"),
    none = NA_real_)

  key <- match(res$locus, tr$.locus)
  res$name <- as.character(tr$name)[key]
  res$gene_id <- if ("gene_id" %in% names(tr)) as.character(tr$gene_id)[key] else NA_character_
  res$chr <- tr$chr[key]
  res$start <- as.numeric(tr$start)[key]
  res$end <- as.numeric(tr$end)[key]
  res$span_bp <- res$end - res$start
  res$locus <- factor(res$locus, levels = tr$.locus[tr$.locus %in% res$locus])
  res <- res[order(res$locus, -res$paired_ratio), , drop = FALSE]

  res <- res[, c("locus", "name", "gene_id", "chr", "start", "end", "span_bp", "group",
                 "n_pairs", "gw_median", "locus_median", "naive_ratio", "ref_median",
                 "paired_ratio", "z_paired", "p_paired", "q_paired", "p_two_sided")]
  # a group-pair label splits back into its two ends, so the result is a matrix to plot
  if (pairs != "within") {
    res <- .split_group_pair(res)
    ratios <- .split_group_pair(ratios)
  }
  if (isTRUE(keep_pair_ratios)) attr(res, "pair_ratios") <- ratios
  attr(res, "adjust") <- adjust
  attr(res, "sharing") <- sharing
  attr(res, "group_col") <- group
  res
}
