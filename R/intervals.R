# A small, general BED-style interval-overlap helper: intersect two sets of genomic
# intervals and partition them (overlapping pairs, and each set's non-overlapping rows).

#' Intersect two sets of genomic intervals
#'
#' A lightweight, dependency-free "bedtools intersect"-style overlap between two interval
#' tables (e.g. genes vs. core regions, SNPs vs. paralog masks). Two intervals overlap when
#' they are on the same chromosome and share at least one base
#' (`start1 <= end2 & end1 >= start2`); chromosome names are matched via [normalise_chr()]
#' so `"Pf3D7_07_v3"` and `"7"` agree.
#'
#' @param locs1,locs2 Interval tables (data frames / tibbles).
#' @param chrom1,start1,end1 Column names for the chromosome / start / end in `locs1`
#'   (defaults `"chr"`, `"start"`, `"end"`; `"chrom"` is accepted as a chromosome alias).
#' @param chrom2,start2,end2 As above for `locs2`.
#' @return A list of three tibbles:
#'   \describe{
#'     \item{`overlap`}{one row per overlapping `locs1`x`locs2` pair -- every `locs1`
#'       column, then every `locs2` column suffixed `.2`, then `overlap_start`,
#'       `overlap_end`, `overlap_bp`.}
#'     \item{`only1`}{`locs1` rows overlapping nothing in `locs2`.}
#'     \item{`only2`}{`locs2` rows overlapping nothing in `locs1`.}
#'   }
#' @examples
#' genes <- data.frame(name = c("a", "b"), chr = c("7", "7"),
#'                     start = c(100, 5000), end = c(200, 5200))
#' core  <- data.frame(chr = "7", start = 50, end = 1000)
#' hit <- bed_intersect(genes, core)
#' hit$overlap$name   # "a"  (inside core)
#' hit$only1$name     # "b"  (subtelomeric)
#' @export
bed_intersect <- function(locs1, locs2,
                          chrom1 = "chr", start1 = "start", end1 = "end",
                          chrom2 = "chr", start2 = "start", end2 = "end") {
  locs1 <- as.data.frame(locs1, stringsAsFactors = FALSE)
  locs2 <- as.data.frame(locs2, stringsAsFactors = FALSE)
  pick <- function(df, col, what) {
    if (col %in% names(df)) return(col)
    if (col == "chr" && "chrom" %in% names(df)) return("chrom")   # accept chrom alias
    stop(sprintf("%s has no '%s' column", what, col), call. = FALSE)
  }
  chrom1 <- pick(locs1, chrom1, "locs1"); chrom2 <- pick(locs2, chrom2, "locs2")
  for (cc in list(list(locs1, start1, end1, "locs1"), list(locs2, start2, end2, "locs2"))) {
    miss <- setdiff(c(cc[[2]], cc[[3]]), names(cc[[1]]))
    if (length(miss)) stop(sprintf("%s has no column(s): %s", cc[[4]],
                                    paste(miss, collapse = ", ")), call. = FALSE)
  }
  c1 <- normalise_chr(locs1[[chrom1]]); s1 <- as.numeric(locs1[[start1]]); e1 <- as.numeric(locs1[[end1]])
  c2 <- normalise_chr(locs2[[chrom2]]); s2 <- as.numeric(locs2[[start2]]); e2 <- as.numeric(locs2[[end2]])

  hit1 <- logical(nrow(locs1)); hit2 <- logical(nrow(locs2))
  parts <- list()
  for (ch in intersect(unique(c1), unique(c2))) {
    i1 <- which(c1 == ch); i2 <- which(c2 == ch)
    for (i in i1) {
      ov <- i2[s2[i2] <= e1[i] & e2[i2] >= s1[i]]
      if (length(ov)) {
        hit1[i] <- TRUE; hit2[ov] <- TRUE
        parts[[length(parts) + 1L]] <- data.frame(
          .i1 = i, .i2 = ov,
          overlap_start = pmax(s1[i], s2[ov]), overlap_end = pmin(e1[i], e2[ov]),
          stringsAsFactors = FALSE)
      }
    }
  }
  if (length(parts)) {
    pp <- do.call(rbind, parts)
    o2 <- locs2[pp$.i2, , drop = FALSE]; names(o2) <- paste0(names(o2), ".2")
    overlap <- cbind(locs1[pp$.i1, , drop = FALSE], o2,
                     overlap_start = pp$overlap_start, overlap_end = pp$overlap_end,
                     overlap_bp = pmax(0, pp$overlap_end - pp$overlap_start))
    rownames(overlap) <- NULL
  } else {
    overlap <- cbind(locs1[0, , drop = FALSE],
                     stats::setNames(locs2[0, , drop = FALSE], paste0(names(locs2), ".2")),
                     overlap_start = numeric(0), overlap_end = numeric(0), overlap_bp = numeric(0))
  }
  list(overlap = tibble::as_tibble(overlap),
       only1 = tibble::as_tibble(locs1[!hit1, , drop = FALSE]),
       only2 = tibble::as_tibble(locs2[!hit2, , drop = FALSE]))
}
