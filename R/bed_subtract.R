# The complement of bed_intersect(): what is left of one interval set once another is
# removed from it, at base resolution.

# Anything that can name a set of covered positions -> a 0-based half-open interval table.
.as_interval_table <- function(x, chrom = "chr", start = "start", end = "end",
                               what = "locs2") {
  if (inherits(x, "PopStructure")) x <- x$genotype()
  if (is.matrix(x)) x <- colnames(x)
  if (is.list(x) && !is.data.frame(x) && !is.null(x$genotype))
    x <- colnames(as.matrix(x$genotype))
  if (is.character(x)) {
    if (!length(x)) return(data.frame(chr = character(), start = numeric(),
                                      end = numeric(), stringsAsFactors = FALSE))
    loci <- .parse_snp_ids(x)                       # `chr:pos`, pos 0-based
    return(data.frame(chr = loci$chr, start = loci$pos, end = loci$pos + 1,
                      stringsAsFactors = FALSE))
  }
  df <- as.data.frame(x, stringsAsFactors = FALSE)
  ch <- if (chrom %in% names(df)) chrom else if ("chrom" %in% names(df)) "chrom" else
    stop(sprintf("%s has no '%s' column", what, chrom), call. = FALSE)
  miss <- setdiff(c(start, end), names(df))
  if (length(miss))
    stop(sprintf("%s has no column(s): %s", what, paste(miss, collapse = ", ")),
         call. = FALSE)
  data.frame(chr = df[[ch]], start = as.numeric(df[[start]]),
             end = as.numeric(df[[end]]), stringsAsFactors = FALSE)
}

# merge overlapping/abutting intervals of one chromosome, given sorted starts
.merge_spans <- function(s, e) {
  if (!length(s)) return(list(s = numeric(), e = numeric()))
  o <- order(s, e); s <- s[o]; e <- e[o]
  ms <- s[1]; me <- e[1]; out_s <- numeric(); out_e <- numeric()
  for (i in seq_along(s)[-1]) {
    if (s[i] <= me) {
      me <- max(me, e[i])
    } else {
      out_s <- c(out_s, ms); out_e <- c(out_e, me); ms <- s[i]; me <- e[i]
    }
  }
  list(s = c(out_s, ms), e = c(out_e, me))
}

#' Subtract one set of genomic intervals from another
#'
#' The complement of [bed_intersect()], and a dependency-free `bedtools subtract`: what is
#' left of `locs1` once everything in `locs2` is removed, cut at **base** resolution. An
#' interval partly covered comes back as the pieces that are not, one row each, carrying its
#' original columns; one covered end to end disappears.
#'
#' The use it was written for is filling gaps in a callset. [aa_intervals()] gives the
#' genomic span of each codon you care about; subtracting the SNPs a panel already has leaves
#' the bases that were never called, which is what to hand `bcftools mpileup -R`. Because the
#' cut is per base, a codon with one of its three bases in the panel still returns the other
#' two -- which is the point, since a residue cannot be read from one base.
#'
#' @section Coordinates: Intervals are **0-based half-open** `[start, end)`, as in BED and
#'   throughout this package, so an interval abutting another (`end1 == start2`) loses
#'   nothing. See [plasgenomicsutilsR-coordinates].
#'
#' @param locs1 Interval table to subtract from (a data frame / tibble).
#' @param locs2 What to remove. An interval table, a character vector of `chr:pos` SNP ids,
#'   a genotype matrix or [load_genotypes()] list (its column names are the ids), or a
#'   [PopStructure] (its genotype panel).
#' @param chrom1,start1,end1 Column names in `locs1` (defaults `"chr"`, `"start"`, `"end"`;
#'   `"chrom"` is accepted as a chromosome alias).
#' @param chrom2,start2,end2 As above for `locs2` when it is a table.
#' @param min_width Drop leftover pieces narrower than this (default `1`, i.e. keep every
#'   base). Raise it to ignore slivers.
#' @return A tibble of the uncovered pieces: every `locs1` column, with `start`/`end`
#'   replaced by the piece's own bounds, plus `piece` (which piece of that row this is) and
#'   `width`. Rows of `locs1` covered completely are absent. Input row order is kept.
#' @seealso [bed_intersect()] for the overlap, [write_bed()] to write the result out.
#' @examples
#' cds <- read_gff_cds(system.file("extdata", "pf3d7_drug_gene_cds.gff",
#'                                 package = "plasgenomicsutilsR"))
#' want <- aa_intervals(data.frame(transcript_id = c("pfcrt", "pfcrt", "pfkelch13"),
#'                                 aa_position = c(72, 76, 580)), cds)
#'
#' # a panel that happens to carry only the middle base of the K76T codon
#' have <- paste0("Pf3D7_07_v3:", want$start[want$aa_position == 76] + 1)
#' gaps <- bed_subtract(want, have)
#' gaps[, c("name", "aa_position", "start", "end", "width")]
#' @export
bed_subtract <- function(locs1, locs2,
                         chrom1 = "chr", start1 = "start", end1 = "end",
                         chrom2 = "chr", start2 = "start", end2 = "end",
                         min_width = 1) {
  a <- as.data.frame(locs1, stringsAsFactors = FALSE)
  if (!nrow(a)) return(tibble::as_tibble(a))
  ch1 <- if (chrom1 %in% names(a)) chrom1 else if ("chrom" %in% names(a)) "chrom" else
    stop(sprintf("locs1 has no '%s' column", chrom1), call. = FALSE)
  miss <- setdiff(c(start1, end1), names(a))
  if (length(miss))
    stop(sprintf("locs1 has no column(s): %s", paste(miss, collapse = ", ")), call. = FALSE)

  b <- .as_interval_table(locs2, chrom2, start2, end2)
  ac <- normalise_chr(a[[ch1]])
  as_ <- as.numeric(a[[start1]]); ae <- as.numeric(a[[end1]])
  bc <- normalise_chr(b$chr)

  # per chromosome, merge what is being removed once, then cut each row against it
  merged <- lapply(split(seq_len(nrow(b)), bc), function(i) .merge_spans(b$start[i], b$end[i]))

  keep_row <- integer(); piece_no <- integer()
  new_s <- numeric(); new_e <- numeric()
  for (i in seq_len(nrow(a))) {
    m <- merged[[ac[i]]]
    s <- as_[i]; e <- ae[i]
    if (is.null(m)) {
      cuts_s <- numeric(); cuts_e <- numeric()
    } else {
      hit <- which(m$s < e & m$e > s)
      cuts_s <- m$s[hit]; cuts_e <- m$e[hit]
    }
    left <- s
    starts <- numeric(); ends <- numeric()
    for (k in seq_along(cuts_s)) {
      if (cuts_s[k] > left) { starts <- c(starts, left); ends <- c(ends, cuts_s[k]) }
      left <- max(left, cuts_e[k])
    }
    if (left < e) { starts <- c(starts, left); ends <- c(ends, e) }
    ok <- (ends - starts) >= min_width
    starts <- starts[ok]; ends <- ends[ok]
    if (!length(starts)) next
    keep_row <- c(keep_row, rep(i, length(starts)))
    piece_no <- c(piece_no, seq_along(starts))
    new_s <- c(new_s, starts); new_e <- c(new_e, ends)
  }
  if (!length(keep_row)) {
    out <- a[0, , drop = FALSE]
    out$piece <- integer(); out$width <- numeric()
    return(tibble::as_tibble(out))
  }
  out <- a[keep_row, , drop = FALSE]
  out[[start1]] <- new_s
  out[[end1]] <- new_e
  out$piece <- piece_no
  out$width <- new_e - new_s
  rownames(out) <- NULL
  tibble::as_tibble(out)
}

#' Write an interval table as a BED file
#'
#' Three columns, tab separated, no header, `start` 0-based half-open -- what `bedtools`
#' and `bcftools mpileup -R` expect. A fourth `name` column is written when the table has
#' one, since a BED that says what each interval is survives being looked at later.
#'
#' **`chrom` is preferred over `chr`.** Tables in this package carry both: `chr` normalised
#' for matching (`"7"`), and `chrom` as the source file spells it (`"Pf3D7_07_v3"`). A BED is
#' read by other tools against a real reference, so it has to carry the name the FASTA and
#' the BAMs use -- writing the normalised one produces a file that matches nothing, silently.
#'
#' @param x An interval table (`chrom`/`chr`, `start`, `end`), e.g. from [bed_subtract()].
#' @param file Path to write.
#' @param chrom Column holding the chromosome name to write. Defaults to `"chrom"` when the
#'   table has it, else `"chr"`.
#' @param name Column to use as the BED name field, or `NULL` for none. Defaults to `"name"`
#'   when the table has it.
#' @param sort Sort by chromosome and start (default `TRUE`), which is what the tools want.
#' @return `file`, invisibly.
#' @seealso [bed_subtract()]
#' @examples
#' iv <- data.frame(chr = c("Pf3D7_07_v3", "Pf3D7_07_v3"), start = c(403623, 403700),
#'                  end = c(403626, 403703), name = c("pfcrt-76", "pfcrt-102"))
#' write_bed(iv, file.path(tempdir(), "targets.bed"))
#' @export
write_bed <- function(x, file, name = NULL, sort = TRUE, chrom = NULL) {
  df <- as.data.frame(x, stringsAsFactors = FALSE)
  ch <- chrom %||% if ("chrom" %in% names(df)) "chrom" else
    if ("chr" %in% names(df)) "chr" else
      stop("`x` needs a `chrom` (or `chr`) column", call. = FALSE)
  if (!ch %in% names(df))
    stop("no `", ch, "` column to use as the chromosome", call. = FALSE)
  miss <- setdiff(c("start", "end"), names(df))
  if (length(miss))
    stop("`x` needs column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  if (missing(name) && "name" %in% names(df)) name <- "name"
  if (!is.null(name) && !name %in% names(df))
    stop("no `", name, "` column to use as the BED name field", call. = FALSE)

  out <- data.frame(chrom = as.character(df[[ch]]),
                    start = format(as.numeric(df$start), scientific = FALSE, trim = TRUE),
                    end = format(as.numeric(df$end), scientific = FALSE, trim = TRUE),
                    stringsAsFactors = FALSE)
  if (!is.null(name)) out$name <- as.character(df[[name]])
  if (sort) out <- out[order(out$chrom, as.numeric(out$start)), , drop = FALSE]
  utils::write.table(out, file, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = FALSE)
  invisible(file)
}
