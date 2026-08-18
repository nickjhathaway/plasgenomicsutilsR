# Group-level summaries of the genome-wide per-pair IBD fraction: one row per pair of
# metadata groups, over every sample pair spanning them.

#' Summarise genome-wide IBD sharing between metadata groups
#'
#' Reduces the per-pair IBD fraction table to one row per pair of `group` levels, over
#' every sample pair that spans them: two groups' row covers all `n_a x n_b` cross-group
#' pairs, and a group's row against itself covers its `choose(n, 2)` within-group pairs.
#' Within-group rows are the ones that say whether a group is internally related at all,
#' so they are included by default and the comparison to the between-group rows is the
#' point.
#'
#' The summary is over *pairs*, not samples, so a highly related cluster inside one group
#' pulls its mean up through every pair it takes part in. The median and the quartiles are
#' there because relatedness is skewed: most pairs share almost nothing and a few share a
#' great deal, so a mean on its own reads as a much more related cohort than any pair
#' actually is.
#'
#' @param x An [IbdResults] holding a pair table (`ibd_results(pair_fraction = )` or
#'   `$set_pair_fraction()`), or the pair table itself as a data frame or path.
#' @param group Metadata column defining the groups. Defaults to the object's declared
#'   group column, then to the first non-`sample` column of `meta`.
#' @param meta Sample metadata (`sample` plus `group`). Taken from `x` when it is an
#'   [IbdResults].
#' @param value Column holding the sharing measure (default `"ibd_fraction_accessible"`,
#'   what `plasgenomicsutils ibd_fraction_and_snp_density` writes). Pass
#'   `"ibd_fraction_full_genome"` to summarise against the full genome instead.
#' @param probs Quantiles to report, as proportions (default the quartiles,
#'   `c(0.25, 0.75)`). Each becomes a `q<pct>` column.
#' @param within Include each group's row against itself (default `TRUE`).
#' @param min_pairs Groups contributing fewer than this many pairs are dropped, with a
#'   note (default `1`, so only empty combinations go).
#' @return A tibble, one row per group pair, ordered by the object's group order (or the
#'   column's factor levels):
#'   \describe{
#'     \item{`group_a`, `group_b`}{the two groups, `group_a == group_b` on a within-group
#'       row.}
#'     \item{`n_samples_a`, `n_samples_b`}{samples of each group present in the pair
#'       table.}
#'     \item{`n_pairs`}{sample pairs the summary is over. On a complete table this is
#'       `n_samples_a * n_samples_b` between groups and `choose(n_samples_a, 2)` within
#'       one; short of that, the table is missing pairs.}
#'     \item{`mean`, `median`, `sd`, `min`, `max`}{over those pairs.}
#'     \item{`q25`, `q75`}{the requested quantiles (type 7, `stats::quantile()`'s
#'       default).}
#'   }
#' @seealso [plot_ibd_pair_network()], which draws the same table pair by pair.
#' @examples
#' pairs <- data.frame(
#'   sample1 = c("s1", "s1", "s1", "s2", "s2", "s3"),
#'   sample2 = c("s2", "s3", "s4", "s3", "s4", "s4"),
#'   ibd_fraction_accessible = c(0.02, 0.31, 0.04, 0.28, 0.03, 0.05))
#' meta <- data.frame(sample = paste0("s", 1:4),
#'                    region = c("north", "north", "south", "south"))
#' pair_fraction_summary(pairs, group = "region", meta = meta)
#' @export
pair_fraction_summary <- function(x, group = NULL, meta = NULL,
                                  value = .PAIR_FRACTION_COL, probs = c(0.25, 0.75),
                                  within = TRUE, min_pairs = 1L) {
  meta <- .normalise_meta(meta)
  if (inherits(x, "IbdResults")) {
    if (is.null(meta)) meta <- x$get_meta()
    if (is.null(group)) group <- x$get_group_col()
    # the stored order belongs to the declared group column; summarising by any other one
    # takes that column's own levels, or the order is a list of labels this grouping has
    # never heard of
    ord <- if (identical(group, x$get_group_col())) x$get_group_order() else NULL
    pf <- x$get_pair_fraction()
    if (is.null(pf))
      stop("this IbdResults has no pair table; build it with ",
           "ibd_results(pair_fraction = ) or call $set_pair_fraction()", call. = FALSE)
    x <- pf
  } else {
    ord <- NULL
    if (is.character(x) && length(x) == 1) x <- .read_maybe(x, "pair table")
  }
  df <- .pair_endpoints(as.data.frame(x, stringsAsFactors = FALSE))
  if (!nrow(df)) stop("the pair table is empty", call. = FALSE)
  if (!value %in% names(df))
    stop("the pair table has no column '", value, "'. Available: ",
         paste(names(df), collapse = ", "), call. = FALSE)

  if (is.null(meta)) stop("`meta` is needed to group the pairs", call. = FALSE)
  meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  if (!"sample" %in% names(meta)) stop("`meta` needs a `sample` column", call. = FALSE)
  if (is.null(group)) group <- setdiff(names(meta), "sample")[1]
  if (is.na(group) || !group %in% names(meta))
    stop("no metadata column `", group, "` to group by.\n  columns available: ",
         paste(setdiff(names(meta), "sample"), collapse = ", "), call. = FALSE)
  if (!is.numeric(probs) || any(probs < 0 | probs > 1))
    stop("`probs` are proportions in [0, 1]", call. = FALSE)

  # a declared order sets the order, but never the membership: a group the order does not
  # mention still gets a row, after the ones it does
  observed <- .levels_of(meta[[group]])
  levs <- if (is.null(ord)) observed
          else c(intersect(ord, observed), setdiff(observed, ord))
  g <- stats::setNames(as.character(meta[[group]]), meta$sample)
  ga <- g[df$sample1]
  gb <- g[df$sample2]
  v <- suppressWarnings(as.numeric(df[[value]]))

  # A pair only belongs to a group pair if both of its samples are grouped, and only
  # contributes if it has a value; both losses are worth naming, since either silently
  # shrinks the denominator every mean below is divided by.
  ungrouped <- sum(is.na(ga) | is.na(gb))
  if (ungrouped)
    message(ungrouped, " of ", nrow(df), " pairs touch a sample with no `", group,
            "`, and are left out")
  ok <- !is.na(ga) & !is.na(gb) & !is.na(v)
  dropped_na <- sum(!is.na(ga) & !is.na(gb) & is.na(v))
  if (dropped_na) message(dropped_na, " pairs have no `", value, "` value, and are left out")
  ga <- ga[ok]; gb <- gb[ok]; v <- v[ok]
  if (!length(v)) stop("no pair has both samples grouped and a value", call. = FALSE)

  # order each pair's groups by the group order, so A-vs-B and B-vs-A land on one row
  ia <- match(ga, levs); ib <- match(gb, levs)
  lo <- ifelse(ia <= ib, ga, gb)
  hi <- ifelse(ia <= ib, gb, ga)
  key <- paste(lo, hi, sep = "\r")

  present <- intersect(levs, unique(c(ga, gb)))
  # counted over the whole pair table, so `n_pairs` short of n_a * n_b means the table is
  # missing pairs rather than the groups being smaller than they look
  in_table <- unique(c(df$sample1, df$sample2))
  n_by <- vapply(present, function(l) sum(!is.na(g[in_table]) & g[in_table] == l),
                 integer(1))

  # one group is a legitimate cohort -- it just has no between-group row
  combos <- if (length(present) >= 2) utils::combn(present, 2)
            else matrix(character(0), nrow = 2)
  wanted <- rbind(data.frame(group_a = combos[1, ], group_b = combos[2, ],
                             stringsAsFactors = FALSE),
                  if (within) data.frame(group_a = present, group_b = present,
                                         stringsAsFactors = FALSE))
  if (!nrow(wanted))
    stop("every pair is within the one `", group, "` group (",
         paste(present, collapse = ", "),
         "); there is no between-group comparison to make (use `within = TRUE`)",
         call. = FALSE)
  wanted <- wanted[order(match(wanted$group_a, levs), match(wanted$group_b, levs)), ,
                   drop = FALSE]

  qn <- paste0("q", .pct_label(probs))
  out <- lapply(seq_len(nrow(wanted)), function(i) {
    a <- wanted$group_a[i]; b <- wanted$group_b[i]
    vals <- v[key == paste(a, b, sep = "\r")]
    q <- if (length(vals)) stats::quantile(vals, probs, names = FALSE)
         else rep(NA_real_, length(probs))
    row <- data.frame(group_a = a, group_b = b,
                      n_samples_a = unname(n_by[a]), n_samples_b = unname(n_by[b]),
                      n_pairs = length(vals),
                      mean = if (length(vals)) mean(vals) else NA_real_,
                      median = if (length(vals)) stats::median(vals) else NA_real_,
                      sd = if (length(vals) > 1) stats::sd(vals) else NA_real_,
                      min = if (length(vals)) min(vals) else NA_real_,
                      max = if (length(vals)) max(vals) else NA_real_,
                      stringsAsFactors = FALSE)
    row[qn] <- as.list(q)
    row
  })
  out <- do.call(rbind, out)

  thin <- out$n_pairs < min_pairs
  if (any(thin)) {
    which_ones <- paste(paste(out$group_a[thin], out$group_b[thin], sep = " vs "),
                        collapse = ", ")
    message(sum(thin), " group pair(s) ",
            if (min_pairs <= 1L) "have no pairs in the table"
            else paste("have fewer than", min_pairs, "pairs"),
            " and are left out: ", which_ones)
    out <- out[!thin, , drop = FALSE]
  }
  out$group_a <- factor(out$group_a, levels = levs[levs %in% out$group_a])
  out$group_b <- factor(out$group_b, levels = levs[levs %in% out$group_b])
  rownames(out) <- NULL
  attr(out, "value") <- value
  attr(out, "group_col") <- group
  tibble::as_tibble(out)
}

# quantile column labels: 0.25 -> "25", 0.5 -> "50", 0.025 -> "2.5". Only a fractional part
# is trimmed -- stripping trailing zeros unconditionally turns 50% into "5", which then
# collides with a genuine 5% column.
.pct_label <- function(probs) {
  lab <- format(probs * 100, trim = TRUE, scientific = FALSE)
  ifelse(grepl(".", lab, fixed = TRUE), sub("\\.?0+$", "", lab), lab)
}
