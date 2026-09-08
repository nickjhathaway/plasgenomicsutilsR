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

# merge overlapping/abutting intervals of one chromosome (or ones within `gap` bases).
# Vectorised: a running maximum of the end handles a short interval nested in a long one,
# and a mask of tens of thousands of blocks (a genome's tandem repeats) merges in
# milliseconds where an append loop took minutes.
.merge_spans <- function(s, e, gap = 0) {
  if (!length(s)) return(list(s = numeric(), e = numeric()))
  o <- order(s, e); s <- s[o]; e <- e[o]
  cm <- cummax(e)
  n <- length(s)
  new_block <- c(TRUE, s[-1] > cm[-n] + gap)
  starts <- which(new_block)
  list(s = s[starts], e = cm[c(starts[-1] - 1L, n)])
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
#' @param pad Widen every `locs2` interval by this many bases on each side before
#'   subtracting (default `0`). Masking tandem repeats for a primer design, say, wants the
#'   repeat plus a few bases of clearance, so nothing ends right at a repeat's edge.
#'   Clamped at the start of the chromosome.
#' @param min_width Drop leftover pieces narrower than this (default `1`, i.e. keep every
#'   base). Raise it to ignore slivers.
#' @return A tibble of the uncovered pieces: every `locs1` column, with `start`/`end`
#'   replaced by the piece's own bounds, plus `piece` (which piece of that row this is) and
#'   `width`. Rows of `locs1` covered completely are absent. Input row order is kept.
#' @seealso [bed_intersect()] for the overlap, [bed_merge()] for joining intervals,
#'   [tandem_repeats_to_avoid()] for a mask worth subtracting, [write_bed()] to write the
#'   result out.
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
#'
#' # a gene minus the tandem repeats a primer should avoid, with 10 bp of clearance
#' crt <- PF3D7_GENES[PF3D7_GENES$name == "pfcrt", ]
#' mask <- tandem_repeats_to_avoid(pf3d7_tandem_repeats())
#' bed_subtract(crt, mask, pad = 10)[, c("start", "end", "piece", "width")]
#' @export
bed_subtract <- function(locs1, locs2,
                         chrom1 = "chr", start1 = "start", end1 = "end",
                         chrom2 = "chr", start2 = "start", end2 = "end",
                         pad = 0, min_width = 1) {
  a <- as.data.frame(locs1, stringsAsFactors = FALSE)
  if (!nrow(a)) return(tibble::as_tibble(a))
  ch1 <- if (chrom1 %in% names(a)) chrom1 else if ("chrom" %in% names(a)) "chrom" else
    stop(sprintf("locs1 has no '%s' column", chrom1), call. = FALSE)
  miss <- setdiff(c(start1, end1), names(a))
  if (length(miss))
    stop(sprintf("locs1 has no column(s): %s", paste(miss, collapse = ", ")), call. = FALSE)

  b <- .as_interval_table(locs2, chrom2, start2, end2)
  if (!is.numeric(pad) || length(pad) != 1L || is.na(pad) || pad < 0)
    stop("`pad` must be one number >= 0", call. = FALSE)
  if (pad > 0) {
    b$start <- pmax(0, b$start - pad)
    b$end <- b$end + pad
  }
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

# the column carrying the reference's own chromosome spelling, by the package's conventions
.bed_chrom_col <- function(nms) {
  ref <- grep("_chrom$", nms, value = TRUE)
  if (length(ref)) return(ref[1])
  if ("chrom" %in% nms) return("chrom")
  if ("chr" %in% nms) return("chr")
  stop("`x` needs a `chrom` (or `chr`) column", call. = FALSE)
}

# `pfcrt`, `pfcrt`, `pfcrt` -> `pfcrt_1`, `pfcrt_2`, `pfcrt_3`; twelve of one name get
# `_01` .. `_12`: zero-padded to the width of that name's own count, in the order given
.make_names_unique <- function(nm) {
  key <- ifelse(is.na(nm), "\001NA", nm)
  dup <- key %in% key[duplicated(key)]
  if (!any(dup)) return(nm)
  idx <- stats::ave(seq_along(key), key, FUN = seq_along)
  cnt <- stats::ave(seq_along(key), key, FUN = length)
  nm[dup] <- paste0(nm[dup], "_", sprintf(paste0("%0", nchar(cnt[dup]), "d"), idx[dup]))
  nm
}

# the BED name field: a column, or `chrom-start-end` when asked for (or when uniqueness is
# wanted and there is no column to make unique)
.bed_names <- function(df, name_col, chrom_v, start_v, end_v, name_is_coords,
                       make_names_unique, fmt) {
  if (isTRUE(name_is_coords) || (isTRUE(make_names_unique) && is.null(name_col)) ||
      (is.null(name_col) && !is.null(fmt$always) && fmt$always)) {
    nm <- paste0(chrom_v, "-", fmt$num(start_v), "-", fmt$num(end_v))
  } else if (!is.null(name_col)) {
    nm <- as.character(df[[name_col]])
  } else {
    return(NULL)
  }
  nm
}

#' Write an interval table as a BED file
#'
#' Three columns, tab separated, no header, `start` 0-based half-open -- what `bedtools`
#' and `bcftools mpileup -R` expect. A fourth `name` column is written when the table has
#' one, since a BED that says what each interval is survives being looked at later.
#'
#' **The reference's own spelling is preferred.** A BED is read by other tools against a
#' real reference, so it has to carry the name the FASTA and the BAMs use -- writing a
#' normalised one produces a file that matches nothing, silently. Tables in this package
#' carry two kinds of chromosome column, and the default picks the one with the full name:
#' a `<assembly>_chrom` column first (`Pf3D7_chrom` in [PF3D7_GENES] and the other bundled
#' datasets, where `chrom` is the short `"7"`), then `chrom` (which [aa_intervals()],
#' [tandem_repeats()] and friends fill with the source spelling, keeping `chr` for the
#' normalised `"7"`), then `chr`. Pieces from [bed_subtract()] keep their `locs1` columns,
#' so a gene table cut by a mask still writes the right names.
#'
#' @param x An interval table (`chrom`/`chr`, `start`, `end`), e.g. from [bed_subtract()].
#' @param file Path to write.
#' @param chrom Column holding the chromosome name to write. Defaults to the first of a
#'   `<assembly>_chrom` column, `"chrom"`, and `"chr"` that the table has (see above).
#' @param name Column to use as the BED name field, or `NULL` for none. Defaults to `"name"`
#'   when the table has it.
#' @param name_is_coords Write `chrom-start-end` (`Pf3D7_07_v3-403221-403363`) as the name
#'   field instead of a column (default `FALSE`).
#' @param make_names_unique Make the name field unique within the file (default `FALSE`):
#'   a name that occurs more than once gets `_1`, `_2`, ... in file order, zero-padded to the
#'   width of that name's count (`_01` to `_12` for twelve). Combines with `name_is_coords`,
#'   so identical coordinates written twice get distinct names. With no name column and
#'   `name_is_coords = FALSE`, the coordinates are used as the names to make unique.
#' @param sort Sort by chromosome and start (default `TRUE`), which is what the tools want.
#'   Numbering for uniqueness follows the written order.
#' @return `file`, invisibly.
#' @seealso [bed_subtract()]
#' @examples
#' iv <- data.frame(chr = c("Pf3D7_07_v3", "Pf3D7_07_v3"), start = c(403623, 403700),
#'                  end = c(403626, 403703), name = c("pfcrt-76", "pfcrt-102"))
#' write_bed(iv, file.path(tempdir(), "targets.bed"))
#' @export
write_bed <- function(x, file, name = NULL, sort = TRUE, chrom = NULL,
                      name_is_coords = FALSE, make_names_unique = FALSE) {
  df <- as.data.frame(x, stringsAsFactors = FALSE)
  ch <- chrom %||% .bed_chrom_col(names(df))
  if (!ch %in% names(df))
    stop("no `", ch, "` column to use as the chromosome", call. = FALSE)
  miss <- setdiff(c("start", "end"), names(df))
  if (length(miss))
    stop("`x` needs column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  if (missing(name) && "name" %in% names(df)) name <- "name"
  if (!is.null(name) && !name %in% names(df))
    stop("no `", name, "` column to use as the BED name field", call. = FALSE)

  num <- function(v) format(as.numeric(v), scientific = FALSE, trim = TRUE)
  chrom_v <- as.character(df[[ch]])
  start_v <- as.numeric(df$start); end_v <- as.numeric(df$end)
  out <- data.frame(chrom = chrom_v, start = num(start_v), end = num(end_v),
                    stringsAsFactors = FALSE)
  nm <- .bed_names(df, name, chrom_v, start_v, end_v, name_is_coords, make_names_unique,
                   list(num = num))
  if (!is.null(nm)) out$name <- nm
  if (sort) out <- out[order(out$chrom, as.numeric(out$start)), , drop = FALSE]
  if (isTRUE(make_names_unique)) out$name <- .make_names_unique(out$name)
  utils::write.table(out, file, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = FALSE)
  invisible(file)
}

#' Write an interval table as a six-column BED, with optional metadata
#'
#' The BED6 layout: `chrom`, `start`, `end`, `name`, `score`, `strand`, tab separated, no
#' header, `start` 0-based half-open. An optional seventh column carries any other columns
#' of the table as `[field=value;field=value;]`, the form `elucidator` reads metadata in
#' (`meta_cols`; the name `meta` is kept for sample-metadata tables throughout the package).
#'
#' Defaults follow the table: `name` is its `name` column, else `chrom-start-end`; `score`
#' is its `score` column, else the interval's width in bp; `strand` is its `strand` column,
#' else `+`. The chromosome column is chosen as [write_bed()] does, preferring the
#' reference's own spelling. Missing values in `name`, `score` and `strand` are written as
#' `.`.
#'
#' Metadata fields and values are written as-is, whitespace included, except that `;` --
#' the separator -- is replaced by `semicolon` (default `:`) so a value can never split a
#' field. Numbers are written in full, never in scientific notation.
#'
#' @param x,file,chrom,sort As for [write_bed()].
#' @param name,score,strand Columns to write in those fields, or `NULL` for the defaults
#'   above.
#' @inheritParams write_bed
#' @param meta_cols Columns to carry in the seventh column, as a character vector of names,
#'   or `TRUE` for every column not already written. `NULL` (the default) writes six columns.
#' @param semicolon What replaces a `;` inside a metadata field or value.
#' @return `file`, invisibly.
#' @seealso [write_bed()] for the three- or four-column form.
#' @examples
#' targets <- bed_subtract(PF3D7_GENES[PF3D7_GENES$name == "pfcrt", ],
#'                         tandem_repeats_to_avoid(pf3d7_tandem_repeats()), pad = 10)
#' f <- file.path(tempdir(), "pfcrt_pieces.bed")
#' write_bed6(targets, f, meta_cols = c("gene_id", "piece"))
#' readLines(f)[1:2]
#'
#' # every piece is called "pfcrt": number them, or name them by their coordinates
#' write_bed6(targets, f, make_names_unique = TRUE)
#' readLines(f)[1:2]
#' write_bed6(targets, f, name_is_coords = TRUE)
#' readLines(f)[1:2]
#' @export
write_bed6 <- function(x, file, name = NULL, score = NULL, strand = NULL, meta_cols = NULL,
                       semicolon = ":", chrom = NULL, sort = TRUE,
                       name_is_coords = FALSE, make_names_unique = FALSE) {
  df <- as.data.frame(x, stringsAsFactors = FALSE)
  ch <- chrom %||% .bed_chrom_col(names(df))
  if (!ch %in% names(df))
    stop("no `", ch, "` column to use as the chromosome", call. = FALSE)
  miss <- setdiff(c("start", "end"), names(df))
  if (length(miss))
    stop("`x` needs column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  pick <- function(arg, default_col, what) {
    if (is.null(arg)) return(if (default_col %in% names(df)) default_col else NULL)
    if (!is.character(arg) || length(arg) != 1L || !arg %in% names(df))
      stop("no `", arg, "` column to use as the BED ", what, " field", call. = FALSE)
    arg
  }
  name_col <- pick(name, "name", "name")
  score_col <- pick(score, "score", "score")
  strand_col <- pick(strand, "strand", "strand")

  num <- function(v) format(as.numeric(v), scientific = FALSE, trim = TRUE)
  chrom_v <- as.character(df[[ch]])
  start_v <- as.numeric(df$start); end_v <- as.numeric(df$end)
  out <- data.frame(
    chrom = chrom_v, start = num(start_v), end = num(end_v),
    name = .bed_names(df, name_col, chrom_v, start_v, end_v, name_is_coords,
                      make_names_unique, list(num = num, always = TRUE)),
    score = if (is.null(score_col)) num(end_v - start_v) else {
      v <- df[[score_col]]; if (is.numeric(v)) num(v) else as.character(v) },
    strand = if (is.null(strand_col)) rep("+", nrow(df)) else as.character(df[[strand_col]]),
    stringsAsFactors = FALSE)
  for (col in c("name", "score", "strand")) out[[col]][is.na(out[[col]])] <- "."

  if (!is.null(meta_cols) && !isFALSE(meta_cols)) {
    used <- c(ch, "start", "end", name_col, score_col, strand_col)
    cols <- if (isTRUE(meta_cols)) setdiff(names(df), used) else as.character(meta_cols)
    bad <- setdiff(cols, names(df))
    if (length(bad))
      stop("no column(s) to carry as metadata: ", paste(bad, collapse = ", "), call. = FALSE)
    if (!is.character(semicolon) || length(semicolon) != 1L)
      stop("`semicolon` must be one string", call. = FALSE)
    clean <- function(v) gsub(";", semicolon, v, fixed = TRUE)
    if (length(cols)) {
      vals <- lapply(cols, function(col) {
        v <- df[[col]]
        v <- if (is.numeric(v)) num(v) else as.character(v)
        v[is.na(v)] <- "NA"
        paste0(clean(col), "=", clean(v), ";")
      })
      out$meta <- paste0("[", do.call(paste0, vals), "]")
    } else {
      out$meta <- "[]"
    }
  }
  if (sort) out <- out[order(out$chrom, as.numeric(out$start)), , drop = FALSE]
  if (isTRUE(make_names_unique)) out$name <- .make_names_unique(out$name)
  utils::write.table(out, file, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = FALSE)
  invisible(file)
}
