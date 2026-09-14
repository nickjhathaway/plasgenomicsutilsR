# The complement of ibd_pair_clusters(): that one reports who is related to whom, this one
# reports the largest group in which nobody is related to anybody. Same graph, same edges
# (.pair_edges), read the other way -- so the two can never disagree about which pairs are
# linked, and neither can disagree with plot_ibd_pair_network().

# Greedy minimum-degree independent set. Take the surviving node with the fewest surviving
# neighbours, keep it, delete it and its neighbours. Minimum-degree is the right rule here
# because it spends the fewest exclusions per node kept. Degrees are maintained
# incrementally -- recomputing them per pick is O(n^2) a pick and dominates on big cohorts.
.mis_greedy <- function(adj, n, random = FALSE) {
  alive <- rep(TRUE, n)
  deg <- lengths(adj)
  chosen <- logical(n)
  remaining <- n
  while (remaining > 0L) {
    live <- which(alive)
    d <- deg[live]
    cand <- live[d == min(d)]
    v <- if (length(cand) == 1L) cand
         else if (random) cand[sample.int(length(cand), 1L)]
         else cand[1L]
    chosen[v] <- TRUE
    nb <- adj[[v]]
    drop <- c(v, nb[alive[nb]])
    alive[drop] <- FALSE
    remaining <- remaining - length(drop)
    # only the surviving neighbours of what we just removed change degree
    for (u in drop) {
      un <- adj[[u]]
      un <- un[alive[un]]
      if (length(un)) deg[un] <- deg[un] - 1L
    }
  }
  chosen
}

# Add back any node with no chosen neighbour, to a fixpoint. The greedy pass can strand a
# node that becomes addable once its blockers are removed; this guarantees the set is
# maximal, so "unselected" always means "genuinely conflicts with something selected".
.mis_maximalise <- function(chosen, adj, n) {
  repeat {
    added <- FALSE
    for (v in which(!chosen)) {
      if (!any(chosen[adj[[v]]])) {
        chosen[v] <- TRUE
        added <- TRUE
      }
    }
    if (!added) break
  }
  chosen
}

#' The largest mutually unrelated set of samples
#'
#' The complement of [ibd_pair_clusters()]: instead of who is related to whom, the largest
#' group found in which **no** pair shares more than `max_ibd` of the genome. Use it to pick
#' founder haplotypes for simulation, to de-duplicate a cohort before statistics that assume
#' independent samples (allele frequencies, \eqn{F_{ST}}, PCA, admixture), or to report how
#' much independent signal a cohort really carries.
#'
#' @details
#' Keeping one sample per single-linkage cluster is **not** the same thing, and is usually
#' much worse. The two coincide only when clusters are cliques; at a low cutoff the
#' relatedness graph is typically one sprawling component, so one-per-cluster returns a
#' handful of samples where this returns many.
#'
#' Finding the true maximum is NP-hard, so this is a heuristic: greedy minimum-degree with
#' randomised restarts, then a pass that adds back any sample that turns out to fit. The
#' result is always **valid** (no pair inside it exceeds `max_ibd`) and **maximal** (no
#' sample could be added), but is not guaranteed to be the largest such set that exists.
#' More `restarts` searches harder. Samples sharing with nobody are always selected.
#'
#' `weight` should hold a fraction computed from **length-filtered** IBD segments. A raw
#' posterior site fraction carries a background floor -- routinely a median near 0.015 across
#' unrelated pairs -- so a cutoff near 1% on one of those counts almost every pair as related
#' and collapses the set to nearly nothing.
#'
#' Set size falls away steeply as `max_ibd` tightens, so scan a few cutoffs rather than
#' trusting one. The cutoff is an inclusive ceiling: a pair sharing exactly `max_ibd` may
#' stay, matching the strictly-greater comparison [plot_ibd_pair_network()] draws with.
#'
#' @inheritParams ibd_pair_clusters
#' @param max_ibd Allow no pair in the set to share more than this fraction (default `0.01`).
#' @param restarts Randomised restarts to try (default `200L`); the first is deterministic.
#' @param seed Seed for tie-breaking, so a run reproduces. The caller's random state is
#'   restored afterwards.
#' @return A tibble with one row per analysed sample, selected first then by name:
#'   \describe{
#'     \item{`sample`}{the sample.}
#'     \item{`selected`}{whether it is in the unrelated set.}
#'     \item{`n_links`}{how many samples it shares more than `max_ibd` with.}
#'     \item{`max_ibd`}{the most it shares with any of them, `NA` when none.}
#'     \item{`blocked_by`}{for an unselected sample, a selected sample it conflicts with, so
#'       the result can be audited; `NA` when selected.}
#'   }
#'   plus one column for each of `add_meta_cols`. Carries a `set_size` attribute and a
#'   `max_ibd_within_set` attribute -- the most any two selected samples actually share, which
#'   is at most `max_ibd` and says how unrelated the set really is.
#' @seealso [ibd_pair_clusters()] for the related-to-whom view, [ibd_pair_links()] for the
#'   edges themselves, [plot_ibd_pair_network()] for the picture.
#' @examples
#' \dontrun{
#' u <- ibd_unrelated_set(ibd_all, max_ibd = 0.01)
#' sum(u$selected)
#' attr(u, "max_ibd_within_set")   # <= 0.01 by construction
#'
#' # how the set shrinks as the cutoff tightens
#' sapply(c(0.005, 0.01, 0.02, 0.05),
#'        function(t) sum(ibd_unrelated_set(ibd_all, max_ibd = t)$selected))
#'
#' # founder haplotypes for a simulation
#' writeLines(u$sample[u$selected], "founders.txt")
#'
#' # does the unrelated set still span every site?
#' u <- ibd_all$ibd_unrelated_set(max_ibd = 0.01, add_meta_cols = "region")
#' table(subset(u, selected)$region)
#' }
#' @export
ibd_unrelated_set <- function(pairs, weight = NULL, max_ibd = 0.01, samples = NULL,
                              restarts = 200L, seed = 1L,
                              add_meta_cols = NULL, meta = NULL) {
  if (!is.null(add_meta_cols) && is.null(meta) && inherits(pairs, "IbdResults"))
    meta <- pairs$get_meta()
  meta <- .normalise_meta(meta)

  # read the table once: .pair_edges() keeps only the pairs above the cutoff, but reporting
  # how related the chosen set actually is needs the pairs below it too
  df_full <- as.data.frame(.pair_fraction_of(pairs), stringsAsFactors = FALSE)

  # edges above the cutoff are exactly the conflicts; `analyzed` is every sample in the table
  pe <- .pair_edges(df_full, weight, min_ibd = max_ibd, samples)
  analyzed <- sort(pe$analyzed)
  n <- length(analyzed)
  e <- pe$edges
  from <- match(e$from, analyzed)
  to <- match(e$to, analyzed)

  adj <- replicate(n, integer(0), simplify = FALSE)
  if (nrow(e)) {
    ends <- c(from, to)
    other <- c(to, from)
    adj <- lapply(split(other, factor(ends, levels = seq_len(n))),
                  function(z) unique(as.integer(z)))
  }

  best <- logical(n)
  if (n) {
    old <- if (exists(".Random.seed", envir = globalenv()))
      get(".Random.seed", envir = globalenv()) else NULL
    on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()), add = TRUE)
    set.seed(seed)
    for (r in seq_len(max(1L, as.integer(restarts)))) {
      ch <- .mis_maximalise(.mis_greedy(adj, n, random = r > 1L), adj, n)
      if (sum(ch) > sum(best)) best <- ch
    }
  }

  out <- data.frame(sample = analyzed, selected = best,
                    n_links = lengths(adj), max_ibd = NA_real_,
                    blocked_by = NA_character_, stringsAsFactors = FALSE)
  if (nrow(e)) {
    ends <- c(from, to)
    w <- c(e$weight, e$weight)
    out$max_ibd <- as.numeric(tapply(w, factor(ends, levels = seq_len(n)), max))
    for (v in which(!best)) {
      hit <- adj[[v]][best[adj[[v]]]]
      if (length(hit)) out$blocked_by[v] <- analyzed[hit[1L]]
    }
  }
  out$n_links <- as.integer(out$n_links)
  out <- out[order(!out$selected, out$sample), , drop = FALSE]
  out <- .add_sample_meta(out, add_meta_cols, meta)
  rownames(out) <- NULL
  out <- tibble::as_tibble(out)

  # the most any two selected samples actually share -- <= max_ibd by construction, and worth
  # reporting: it says how unrelated the set really is, not just that the cutoff held
  sel <- analyzed[best]
  dfp <- .pair_endpoints(df_full)
  keep <- dfp$sample1 %in% sel & dfp$sample2 %in% sel
  within <- suppressWarnings(as.numeric(dfp[[pe$weight_col]][keep]))
  within <- within[!is.na(within)]
  attr(out, "set_size") <- length(sel)
  attr(out, "max_ibd_within_set") <- if (length(within)) max(within) else NA_real_
  out
}
