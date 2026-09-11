# The IBD analogue of iHS: split the pairs sharing at a locus by their carriage of a core
# variant, and ask whether the carriers' segments are the long ones.
#
# Why it is worth the extra step. Extension at a locus, measured over all pairs, cannot
# separate "this locus is under selection in this group" from "this particular haplotype is".
# Splitting on the core variant does, and if the two strata come out equal in a group, the
# variant is not the target there. iHS is literally the log ratio of integrated EHH between
# the two core alleles; this is the same contrast with IBD segment length in place of EHH,
# so it is named and framed to match rather than to compete (see [run_ihs()]).

# Per-sample allele state, from a metadata column or a vector, as character with NA for
# unknown. Two states is the usual case but nothing here requires exactly two.
# Rows of a (locus, group) frame counted against the result's own key order.
.count_by_key <- function(df, kk) {
  if (!nrow(df)) return(rep(0L, length(kk)))
  tb <- table(paste(df$locus, df$group, sep = "\r"))
  out <- as.integer(tb[kk])
  out[is.na(out)] <- 0L
  out
}

.allele_states <- function(allele, meta, samples) {
  if (is.character(allele) && length(allele) == 1 && allele %in% names(meta)) {
    v <- stats::setNames(as.character(meta[[allele]]), as.character(meta$sample))
    lv <- .levels_of(meta[[allele]])
  } else if (!is.null(names(allele))) {
    v <- stats::setNames(as.character(allele), names(allele))
    lv <- .levels_of(allele)
  } else {
    stop("`allele` must name a column of `meta`, or be a named vector of sample -> state.\n",
         "  columns available: ", paste(setdiff(names(meta), "sample"), collapse = ", "),
         call. = FALSE)
  }
  list(state = v[samples], levels = lv)
}

# Two-sided Mann-Whitney between strata, safe to run thousands of times.
.rank_sum_p <- function(a, b) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (length(a) < 1 || length(b) < 1) return(NA_real_)
  p <- tryCatch(suppressWarnings(stats::wilcox.test(a, b)$p.value), error = function(e) NA_real_)
  if (is.null(p)) NA_real_ else as.numeric(p)
}

#' Split IBD block extension at a locus by carriage of a variant
#'
#' The allele-level form of [ibd_block_extension_test()]. Within each group, the pairs are
#' split by their state at a core variant into carrier/carrier, reference/reference and
#' discordant, the same per-pair extension statistic is computed in each stratum, and the two
#' matching strata are compared. Extension measured over all pairs says a locus is unusual in
#' a group; this says whether it is the *haplotype* that is unusual, which is a different and
#' stronger claim.
#'
#' Two statistics come back per locus and group, and they answer different questions:
#'
#' * **length**, `ratio_carrier` against `ratio_reference`, tested by `p_length`. Among pairs
#'   that share the locus at all, are the carriers' segments the longer ones?
#' * **fraction**, `frac_carrier` against `frac_reference`, tested by `p_fraction`. Of all
#'   the pairs that *could* share, how many do? This is usually the better powered of the
#'   two, because it uses every pair rather than only the sharing ones, so read it first.
#'
#' @section Which states are contrasted:
#' The non-carrier class is the **named reference state only**: a sample carrying some third
#' state is excluded from both strata, not pooled into the reference. That is deliberate and
#' is the whole point at a multiallelic site. At *pfpx1* codon 384, D384A, D384G and D384Y
#' arose independently, so pooling the alternates into one "not reference" class would merge
#' origins the analysis exists to separate. To contrast two alternates directly, name them as
#' `carrier` and `reference`; that is a different question and it should look different in
#' the call. The pairs left out either way are counted in `n_excluded_other_allele`.
#'
#' @section The discordant stratum is a check, not a result:
#' A pair sharing an interval by descent shares whatever allele sits in it, so pairs that are
#' IBD across the locus and discordant at the variant should be rare. A large
#' `n_discordant` points at genotyping error, at a recombination inside the interval, or at
#' the interval being wider than the haplotype -- look there before reading anything else.
#'
#' @param x An [IbdResults] built with `blocks =` and `meta =`.
#' @param loci Loci to test, as [ibd_block_extension_test()] takes them.
#' @param allele The variant to split on: the name of a metadata column, or a named vector
#'   of `sample -> state`. Samples with `NA` are left out of every stratum.
#' @param carrier,reference Which states count as carrying and as reference. With exactly
#'   two states either may be left out and the other is implied; leaving both out reads the
#'   column's first two levels, `reference` first, with a message saying which way round it
#'   was read. With **three or more states both must be named** -- there is nothing to imply
#'   the second one from, and guessing would contrast one alternate against another. Name
#'   them explicitly when the labels are not self-evident.
#' @param min_pairs Pairs a stratum needs at a locus before its length statistic is computed
#'   (default `5`). The fraction statistic is reported whatever the count, since it has a
#'   denominator either way.
#' @param adjust Benjamini-Hochberg scope for `p_length` and `p_fraction`, as in
#'   [ibd_block_extension_test()].
#' @param group,within,sharing,min_ref_blocks,meta Passed to [ibd_block_extension_test()].
#' @return A tibble, one row per locus x group, ordered by locus then by `p_fraction`:
#'   \describe{
#'     \item{`locus`, `name`, `gene_id`, `chr`, `start`, `end`, `span_bp`, `group`}{as in
#'       [ibd_block_extension_test()].}
#'     \item{`n_carrier`, `n_reference`, `n_discordant`}{pairs IBD across the locus in each
#'       stratum. Report all three; see the note on the discordant one.}
#'     \item{`n_excluded_other_allele`}{pairs IBD across the locus that were left out
#'       because an end carries a state that is neither `carrier` nor `reference`. Zero at a
#'       two-state locus. These are not discordant and not unknown, so they would otherwise
#'       be invisible -- and at a multiallelic locus they can be the pairs that would show
#'       the interval is wider than the haplotype.}
#'     \item{`ratio_carrier`, `ratio_reference`}{`paired_ratio` within each stratum.}
#'     \item{`ratio_contrast`}{`ratio_carrier / ratio_reference` on the per-pair scale, so
#'       `1` means carriage makes no difference to segment length.}
#'     \item{`p_length`}{two-sided rank-sum test of the two strata's per-pair log ratios.}
#'     \item{`n_carrier_possible`, `n_reference_possible`}{pairs that *could* share, from the
#'       analyzed samples of each stratum.}
#'     \item{`frac_carrier`, `frac_reference`}{sharing pairs over possible pairs.}
#'     \item{`odds_ratio`, `p_fraction`}{Fisher's exact test of those two fractions. The
#'       better-powered comparison.}
#'     \item{`q_length`, `q_fraction`}{adjusted within the family set by `adjust`.}
#'   }
#'   The per-stratum full statistics are on the `strata` attribute, and the per-pair ratios,
#'   labelled by stratum, on `pair_ratios`.
#' @seealso [ibd_block_extension_test()], [ibd_block_extension_scan()] for the genome-wide
#'   background, [run_ihs()] for the haplotype-length statistic this mirrors.
#' @examples
#' \dontrun{
#' ibd <- ibd_results(blocks = "hmm.txt", meta = meta, group_col_in_meta = "region")
#' ibd_block_extension_by_allele(ibd, loci = "pfpx1", allele = "PIN_Haplotype",
#'                               carrier = "Present", reference = "Absent")
#' }
#' @export
ibd_block_extension_by_allele <- function(x, loci, allele, carrier = NULL, reference = NULL,
                                          group = NULL, within = 0,
                                          sharing = c("overlap", "complete"),
                                          min_ref_blocks = 1L, min_pairs = 5L,
                                          adjust = c("per_locus", "per_group", "all", "none"),
                                          meta = NULL) {
  sharing <- match.arg(sharing)
  adjust <- match.arg(adjust)
  if (!inherits(x, "IbdResults"))
    stop("`x` must be an IbdResults; build one with ibd_results(blocks = , meta = )",
         call. = FALSE)
  if (missing(allele) || is.null(allele)) stop("`allele` is required", call. = FALSE)
  if (is.null(meta)) meta <- x$get_meta()
  meta <- .normalise_meta(meta)
  if (is.null(meta) || !"sample" %in% names(meta))
    stop("splitting the pairs needs meta with a 'sample' column", call. = FALSE)

  # the pair-level ratios, computed exactly as the locus test computes them -- min_pairs is
  # applied per stratum below, so nothing is filtered out here
  base <- ibd_block_extension_test(
    x, loci = loci, group = group, within = within, sharing = sharing,
    min_ref_blocks = min_ref_blocks, min_pairs = 1L, adjust = "none", meta = meta)
  pr <- attr(base, "pair_ratios")
  if (!nrow(pr)) stop("no pair shares any of these loci", call. = FALSE)

  st <- .allele_states(allele, meta, unique(c(pr$sample1, pr$sample2, x$get_analyzed_samples())))
  lv <- st$levels
  if (is.null(reference) && is.null(carrier)) {
    if (length(lv) < 2)
      stop("`allele` has fewer than two observed states, so there is nothing to contrast",
           call. = FALSE)
    if (length(lv) > 2)
      stop("`allele` has ", length(lv), " states (", paste(lv, collapse = ", "),
           "); name `carrier =` and `reference =`", call. = FALSE)
    reference <- lv[1]; carrier <- lv[2]
    message("reading `", reference, "` as reference and `", carrier, "` as carrier; ",
            "pass carrier=/reference= to swap them")
  } else if (is.null(reference) || is.null(carrier)) {
    # Only one side was named. With two states the other is implied; with three or more
    # there is nothing to imply it, and taking the first remaining level in sort order
    # would silently contrast one alternate against another -- at a codon carrying D384A,
    # D384G and D384Y, `carrier = "D384A"` would be measured against D384G, two independent
    # origins, with no message. The no-argument path above already refuses this; so does
    # this one.
    if (length(lv) > 2) {
      named <- if (is.null(reference)) "carrier" else "reference"
      missing <- if (is.null(reference)) "reference" else "carrier"
      stop("`allele` has ", length(lv), " states (", paste(lv, collapse = ", "),
           "); with `", named, " =` given, `", missing,
           " =` cannot be inferred -- name it too. Contrasting one alternate against ",
           "another is a different question from contrasting it against the reference.",
           call. = FALSE)
    }
    if (is.null(reference)) reference <- setdiff(lv, carrier)[1]
    else carrier <- setdiff(lv, reference)[1]
  }
  if (is.na(carrier) || is.na(reference) || identical(carrier, reference))
    stop("`carrier` and `reference` must be two different states of `allele`", call. = FALSE)

  code <- function(s) {
    a <- unname(st$state[s])
    ifelse(is.na(a), NA_character_,
           ifelse(a == carrier, "carrier", ifelse(a == reference, "reference", NA_character_)))
  }
  pr$a1 <- code(pr$sample1)
  pr$a2 <- code(pr$sample2)
  pr$stratum <- ifelse(is.na(pr$a1) | is.na(pr$a2), NA_character_,
                       ifelse(pr$a1 == pr$a2, pr$a1, "discordant"))

  # Pairs dropped because an end carries a *third* state -- both states known, neither the
  # carrier nor the reference. They are not discordant (that column is carrier-vs-reference)
  # and they are not unknown, so without a count of their own they vanish from a result the
  # documentation tells the reader to judge by `n_discordant`. At a multiallelic locus these
  # can be exactly the pairs that would show the interval is wider than the haplotype.
  known <- !is.na(unname(st$state[pr$sample1])) & !is.na(unname(st$state[pr$sample2]))
  excl <- pr[is.na(pr$stratum) & known, c("locus", "group"), drop = FALSE]

  pr <- pr[!is.na(pr$stratum), , drop = FALSE]
  if (!nrow(pr))
    stop("no sharing pair has a known `allele` state at both ends", call. = FALSE)

  strata <- .summarise_ratios(pr, c("locus", "group", "stratum"), min_pairs)

  # Denominators: the pairs that *could* share. Taken from the analyzed samples -- the set
  # before the IBD filter -- so a stratum whose pairs never share still has a denominator and
  # the fraction does not quietly condition on sharing.
  gcol <- attr(base, "group_col")
  gmap <- stats::setNames(as.character(meta[[gcol]]), as.character(meta$sample))
  an <- x$get_analyzed_samples()
  cnt <- table(group = gmap[an], stratum = code(an), useNA = "no")
  possible <- function(g, s) {
    if (!(g %in% rownames(cnt)) || !(s %in% colnames(cnt))) return(0)
    n <- as.numeric(cnt[g, s])
    n * (n - 1) / 2
  }

  key <- unique(pr[, c("locus", "group")])
  kk <- paste(key$locus, key$group, sep = "\r")
  n_by <- stats::aggregate(list(n = pr$log2_ratio), pr[, c("locus", "group", "stratum")], length)

  pick <- function(nm) {
    sub <- n_by[n_by$stratum == nm, , drop = FALSE]
    out <- sub$n[match(kk, paste(sub$locus, sub$group, sep = "\r"))]
    out[is.na(out)] <- 0
    out
  }
  res <- dplyr::tibble(
    locus = key$locus, group = key$group,
    n_carrier = as.integer(pick("carrier")),
    n_reference = as.integer(pick("reference")),
    n_discordant = as.integer(pick("discordant")),
    n_excluded_other_allele = as.integer(.count_by_key(excl, kk)))

  from_strata <- function(nm, col) {
    sub <- strata[strata$stratum == nm, , drop = FALSE]
    sub[[col]][match(kk, paste(sub$locus, sub$group, sep = "\r"))]
  }
  res$ratio_carrier <- from_strata("carrier", "paired_ratio")
  res$ratio_reference <- from_strata("reference", "paired_ratio")
  res$ratio_contrast <- res$ratio_carrier / res$ratio_reference

  res$p_length <- vapply(seq_len(nrow(res)), function(i) {
    sel <- pr$locus == res$locus[i] & pr$group == res$group[i]
    a <- pr$log2_ratio[sel & pr$stratum == "carrier"]
    b <- pr$log2_ratio[sel & pr$stratum == "reference"]
    if (length(a) < min_pairs || length(b) < min_pairs) return(NA_real_)
    .rank_sum_p(a, b)
  }, numeric(1))

  # unname: vapply over a character vector carries its values through as names, which would
  # otherwise ride along on every fraction column
  res$n_carrier_possible <- unname(vapply(as.character(res$group),
                                          function(g) possible(g, "carrier"), numeric(1)))
  res$n_reference_possible <- unname(vapply(as.character(res$group),
                                            function(g) possible(g, "reference"), numeric(1)))
  res$frac_carrier <- ifelse(res$n_carrier_possible > 0,
                             res$n_carrier / res$n_carrier_possible, NA_real_)
  res$frac_reference <- ifelse(res$n_reference_possible > 0,
                               res$n_reference / res$n_reference_possible, NA_real_)
  ft <- lapply(seq_len(nrow(res)), function(i) {
    m <- matrix(c(res$n_carrier[i], max(0, res$n_carrier_possible[i] - res$n_carrier[i]),
                  res$n_reference[i], max(0, res$n_reference_possible[i] - res$n_reference[i])),
                nrow = 2)
    if (any(is.na(m)) || any(m < 0) || sum(m) == 0) return(c(NA_real_, NA_real_))
    r <- tryCatch(stats::fisher.test(m), error = function(e) NULL)
    if (is.null(r)) c(NA_real_, NA_real_) else c(unname(r$estimate), r$p.value)
  })
  res$odds_ratio <- vapply(ft, `[`, numeric(1), 1)
  res$p_fraction <- vapply(ft, `[`, numeric(1), 2)

  fam <- switch(adjust, per_locus = res$locus, per_group = res$group,
                all = rep("all", nrow(res)), none = NULL)
  adj <- function(p) if (is.null(fam)) NA_real_ else
    stats::ave(p, fam, FUN = function(z) stats::p.adjust(z, "BH"))
  res$q_length <- adj(res$p_length)
  res$q_fraction <- adj(res$p_fraction)

  k <- match(res$locus, as.character(base$locus))
  for (cc in c("name", "gene_id", "chr", "start", "end", "span_bp")) res[[cc]] <- base[[cc]][k]
  res$locus <- factor(res$locus, levels = levels(base$locus))
  res <- res[order(res$locus, res$p_fraction), , drop = FALSE]
  res <- res[, c("locus", "name", "gene_id", "chr", "start", "end", "span_bp", "group",
                 "n_carrier", "n_reference", "n_discordant", "n_excluded_other_allele",
                 "ratio_carrier", "ratio_reference", "ratio_contrast", "p_length", "q_length",
                 "n_carrier_possible", "n_reference_possible", "frac_carrier",
                 "frac_reference", "odds_ratio", "p_fraction", "q_fraction")]
  attr(res, "strata") <- strata
  attr(res, "pair_ratios") <- pr
  attr(res, "states") <- c(carrier = carrier, reference = reference)
  attr(res, "sharing") <- sharing
  res
}
