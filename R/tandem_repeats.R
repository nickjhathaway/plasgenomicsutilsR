# Short tandem repeats: read a repeat-finder BED, class each repeat by the period of its
# unit, flag the ones a polymerase is likely to slip on, and merge runs of repeats that flow
# into one another -- so a target design can subtract them (bed_subtract(..., pad = )).
#
# The rules and the bundled Pf3D7 table come from the HEOME redesign work
# (Initial_Windows.qmd): repeats from `trf` plus a custom simple-repeat finder, then
# per-record thresholds by unit size, then everything merged and anything long kept.

.tr_cache <- new.env(parent = emptyenv())

# --- the period of a repeat unit -----------------------------------------------------------
# The smallest p such that the unit is its first p bases repeated: "ATATAT" -> 2, "AAT" -> 3,
# "AAAA" -> 1, "TATT" -> 4. Repeat finders spell the same run as "AT", "ATAT" or "ATATAT"
# depending on the alignment, and slippage follows the true period, not the spelling.
.repeat_period <- function(unit) {
  out <- rep(NA_integer_, length(unit))
  ok <- !is.na(unit) & nzchar(unit)
  if (!any(ok)) return(out)
  u <- unique(unit[ok])
  n <- nchar(u)
  p <- n                                   # a unit with no shorter period is its own
  found <- n <= 1L
  for (k in seq_len(max(n) %/% 2)) {
    i <- which(!found & n %% k == 0L & n > k)
    if (!length(i)) next
    hit <- strrep(substr(u[i], 1L, k), n[i] %/% k) == u[i]
    p[i[hit]] <- k
    found[i[hit]] <- TRUE
  }
  out[ok] <- p[match(unit[ok], u)]
  out
}

# `UNIT_xCOPIES`, with or without a `chrom-start-end__` location prefix (the prefix only
# restates the coordinates) -> unit and copy number; NA where the name is not of that form
.parse_repeat_name <- function(name) {
  code <- sub("^.*__", "", as.character(name))
  ok <- !is.na(code) & grepl("^[A-Za-z]+_x[0-9]+(\\.[0-9]+)?$", code)
  unit <- ifelse(ok, toupper(sub("_x.*$", "", code)), NA_character_)
  copies <- ifelse(ok, suppressWarnings(as.numeric(sub("^.*_x", "", code))), NA_real_)
  list(unit = unit, copies = copies)
}

# a table that already went through tandem_repeats() is passed straight through
.is_tandem_table <- function(x) {
  is.data.frame(x) && all(c("chr", "start", "end", "width", "repeat_unit", "period") %in% names(x))
}

.as_tandem_repeats <- function(x) {
  if (.is_tandem_table(x)) tibble::as_tibble(x) else tandem_repeats(x)
}

# the first four BED fields out of a file or a data frame, whatever they are called
.bed4 <- function(x) {
  if (is.character(x)) {
    if (length(x) != 1L)
      stop("`x` must be one file path, or a data frame", call. = FALSE)
    if (!file.exists(x)) stop("no such file: ", x, call. = FALSE)
    df <- utils::read.delim(x, header = FALSE, comment.char = "#", quote = "",
                            colClasses = "character", stringsAsFactors = FALSE)
    df <- df[!grepl("^(track|browser)\\b", df[[1]]), , drop = FALSE]
    if (ncol(df) < 4L)
      stop("`", basename(x), "` has ", ncol(df), " columns; a tandem-repeat BED needs ",
           "chrom, start, end and a name carrying the repeat unit (`UNIT_xCOPIES`)",
           call. = FALSE)
    return(data.frame(chrom = df[[1]], start = as.numeric(df[[2]]), end = as.numeric(df[[3]]),
                      name = df[[4]], stringsAsFactors = FALSE))
  }
  if (!is.data.frame(x))
    stop("`x` must be a BED file path or a data frame", call. = FALSE)
  df <- as.data.frame(x, stringsAsFactors = FALSE)
  nms <- names(df)
  by_name <- all(c("start", "end", "name") %in% nms) && any(c("chrom", "chr") %in% nms)
  if (by_name) {
    ch <- if ("chrom" %in% nms) "chrom" else "chr"
    return(data.frame(chrom = as.character(df[[ch]]), start = as.numeric(df$start),
                      end = as.numeric(df$end), name = as.character(df$name),
                      stringsAsFactors = FALSE))
  }
  if (ncol(df) < 4L)
    stop("`x` needs `chrom`/`chr`, `start`, `end` and `name` columns, or at least four ",
         "unnamed BED columns", call. = FALSE)
  data.frame(chrom = as.character(df[[1]]), start = as.numeric(df[[2]]),
             end = as.numeric(df[[3]]), name = as.character(df[[4]]),
             stringsAsFactors = FALSE)
}

#' Read a table of short tandem repeats
#'
#' Turns the BED a repeat finder writes into a table that says, for each repeat, what its
#' unit is, how long that unit really is (its *period*), and how many bases the run covers
#' -- the three things slippage depends on. [flag_tandem_repeats()] applies the thresholds,
#' [merge_tandem_repeats()] joins runs that flow into one another, and
#' [tandem_repeats_to_avoid()] does both.
#'
#' Only the first four BED fields are read. The name field must be `UNIT_xCOPIES`
#' (`AACCCTA_x50.5714`), optionally preceded by a `chrom-start-end__` location, which is
#' what `elucidator TandemRepeatFinderOutputToBed` and `findSimpleTandemRepeatLocations`
#' write (`Pf3D7_01_v3-2-356__AACCCTA_x50.5714`); the location only restates the
#' coordinates and is ignored. A name of any other form gives `NA` for the unit, period and
#' copies, and such a row can still be flagged on its total width alone.
#'
#' **Period, not unit size.** A finder spells the same run as `AT`, `ATAT` or `ATATAT`
#' depending on where its alignment started. `period` is the smallest number of bases the
#' unit is a repetition of (2 for all three spellings; 1 for `AAAA`; 4 for `TATT`), so a
#' dinucleotide threshold reaches every dinucleotide run however it was reported.
#' `unit_size` keeps the reported spelling's length. The same span is sometimes reported
#' twice with different units (the two finders disagree on the alignment); both rows are
#' kept, and merging collapses them.
#'
#' @section Coordinates: 0-based half-open `[start, end)`, as in BED and throughout the
#'   package; see [plasgenomicsutilsR-coordinates]. `chr` is the chromosome normalised for
#'   matching ([normalise_chr()]), `chrom` as the file spells it.
#'
#' @param x Path to a BED file (plain or gzipped), or a data frame holding the first four
#'   BED fields -- by name (`chrom`/`chr`, `start`, `end`, `name`) or, failing that, by
#'   position.
#' @return A tibble with one row per repeat: `chrom`, `chr`, `start`, `end`, `width`,
#'   `name`, `repeat_unit`, `unit_size`, `period`, `copies`. Input order is kept.
#' @seealso [pf3d7_tandem_repeats()] for the bundled Pf3D7 table,
#'   [tandem_repeats_to_avoid()] for the whole recipe.
#' @examples
#' bed <- data.frame(chrom = "Pf3D7_07_v3", start = c(100, 130, 500),
#'                   end = c(120, 160, 560),
#'                   name = c("Pf3D7_07_v3-100-120__A_x20",
#'                            "Pf3D7_07_v3-130-160__ATAT_x7.5",
#'                            "Pf3D7_07_v3-500-560__TATT_x15"))
#' tandem_repeats(bed)[, c("start", "end", "width", "repeat_unit", "period")]
#' @export
tandem_repeats <- function(x) {
  df <- .bed4(x)
  bad <- is.na(df$start) | is.na(df$end) | df$end < df$start
  if (any(bad))
    stop(sum(bad), " row(s) have a missing or inverted start/end (first at row ",
         which(bad)[1], ")", call. = FALSE)
  parsed <- .parse_repeat_name(df$name)
  out <- tibble::tibble(
    chrom = df$chrom,
    chr = normalise_chr(df$chrom),
    start = df$start,
    end = df$end,
    width = df$end - df$start,
    name = df$name,
    repeat_unit = parsed$unit,
    unit_size = nchar(parsed$unit),
    period = .repeat_period(parsed$unit),
    copies = parsed$copies
  )
  n_na <- sum(is.na(out$repeat_unit))
  if (n_na)
    message(n_na, " of ", nrow(out), " repeat names are not `UNIT_xCOPIES` (with or ",
            "without a `chrom-start-end__` prefix); their unit and period are NA ",
            "(flagged on width alone)")
  out
}

# The stored table keeps the finder's copy number only where it is not width / unit_size
# (to the six significant digits the finder printed); the rest is recomputed here. The
# build script checks the names come back byte-identical to the source BED.
.rebuild_repeat_names <- function(z) {
  chrom <- as.character(z$chrom)
  unit <- as.character(z$repeat_unit)
  copies <- z$copies
  na <- is.na(copies)
  copies[na] <- signif((z$end[na] - z$start[na]) / nchar(unit[na]), 6)
  paste0(chrom, "-", z$start, "-", z$end, "__", unit, "_x", as.character(copies))
}

#' Tandem repeats of the Pf3D7 reference genome
#'
#' Every short tandem repeat found in the *P. falciparum* 3D7 reference, as a
#' [tandem_repeats()] table, so nobody has to run the finders again. Feed it to
#' [tandem_repeats_to_avoid()] and subtract the result from a set of targets with
#' [bed_subtract()].
#'
#' @details
#' Built once (April 2021) from `Pf3D7.fasta` (PlasmoDB, the assembly [PF3D7_GENES] is
#' annotated on) by the union of two finders:
#'
#' * Tandem Repeats Finder (Benson 1999), `trf Pf3D7.fasta 2 7 7 80 10 50 1000 -f -d -m -h`
#'   -- match 2, mismatch 7, indel 7, match / indel probabilities 80 / 10, minimum score 50,
#'   maximum period 1000 -- converted to BED with
#'   `elucidator TandemRepeatFinderOutputToBed`;
#' * `elucidator findSimpleTandemRepeatLocations --maxRepeatUnitSize 10`, an exhaustive
#'   search for perfect repeats of units up to 10 bp, which `trf` is inconsistent about
#'   around short simple repeats.
#'
#' `copies` is the finder's own copy number, which for `trf` reflects an alignment with
#' indels and so is not always `width / unit_size`. All 16 sequences are covered, the
#' apicoplast and mitochondrion included.
#'
#' The table is stored compactly under `inst/extdata` (about 1.1 MB) and expanded on first
#' use; later calls return the cached copy. `data-raw/PF3D7_TANDEM_REPEATS.R` rebuilds it.
#'
#' @return The [tandem_repeats()] table for Pf3D7: about 272,000 rows.
#' @references Benson G (1999). Tandem repeats finder: a program to analyze DNA sequences.
#'   *Nucleic Acids Research* 27(2):573-580.
#' @seealso [tandem_repeats_to_avoid()], [bed_subtract()]
#' @examples
#' tr <- pf3d7_tandem_repeats()
#' table(pmin(tr$period, 5), useNA = "ifany")     # 1..4, and "5" for anything longer
#' @export
pf3d7_tandem_repeats <- function() {
  if (!is.null(.tr_cache$pf3d7)) return(.tr_cache$pf3d7)
  f <- system.file("extdata", "pf3d7_tandem_repeats.rds", package = "plasgenomicsutilsR")
  z <- readRDS(f)
  df <- data.frame(chrom = as.character(z$chrom), start = as.numeric(z$start),
                   end = as.numeric(z$end), name = .rebuild_repeat_names(z),
                   stringsAsFactors = FALSE)
  out <- tandem_repeats(df)
  attr(out, "source_version") <- attr(z, "source_version")
  .tr_cache$pf3d7 <- out
  out
}

# `c(11, 12, 21)` or `c("1" = 11, "2" = 12, "3" = 21)` -> named by period
.period_thresholds <- function(x) {
  if (is.null(x) || !length(x)) return(stats::setNames(numeric(), character()))
  x <- unlist(x)
  if (!is.numeric(x) || any(is.na(x)))
    stop("`min_width_by_period` must be numeric: one minimum width per period", call. = FALSE)
  nms <- names(x)
  if (is.null(nms) || any(!nzchar(nms))) nms <- as.character(seq_along(x))
  p <- suppressWarnings(as.integer(nms))
  if (any(is.na(p)) || any(p < 1L) || any(p != as.numeric(nms)))
    stop("names of `min_width_by_period` must be periods (positive whole numbers), got: ",
         paste(nms, collapse = ", "), call. = FALSE)
  if (anyDuplicated(p))
    stop("`min_width_by_period` names a period twice", call. = FALSE)
  stats::setNames(as.numeric(x), as.character(p))
}

#' Flag the tandem repeats a polymerase is likely to slip on
#'
#' Adds an `avoid` column to a [tandem_repeats()] table. A repeat is flagged when it is at
#' least `min_width` bases long whatever its unit, or when its **period** has its own,
#' shorter, threshold and the run reaches that. Short units slip at shorter lengths, so
#' the defaults descend with the period:
#'
#' | period | flagged from | i.e. |
#' | --- | --- | --- |
#' | 1 (homopolymer) | 11 bp | `AAAAAAAAAAA` |
#' | 2 (dinucleotide) | 12 bp | six copies of `AT` |
#' | 3 (trinucleotide) | 21 bp | seven copies of `AAT` |
#' | any | 50 bp | |
#'
#' These are the thresholds the HEOME Pf3D7 design used. Tighten or loosen them per
#' period, add periods (`c("1" = 11, "2" = 12, "3" = 21, "4" = 30)`), or drop the period
#' rules altogether (`min_width_by_period = NULL`) to flag on total length only.
#'
#' **Records or blocks.** Given a [tandem_repeats()] table, each repeat is judged on its own
#' width. Given the blocks from [merge_tandem_repeats()], each block is judged on its
#' *combined* width, held to the lowest threshold among the periods it contains: a block
#' holding a homopolymer is flagged from 11 bp however much of it is something else, so an
#' 8 bp `A` run flowing into a 9 bp `ATA` run is avoided as one stretch although neither
#' would be on its own. [tandem_repeats_to_avoid()] offers both as `rule`.
#'
#' @param x A [tandem_repeats()] table, or anything [tandem_repeats()] accepts; or the block
#'   table from [merge_tandem_repeats()].
#' @param min_width_by_period Minimum total width (bp) to flag, per period. A named
#'   numeric vector whose names are periods, or an unnamed one read as periods 1, 2, 3, ...
#'   in order. `NULL` for no per-period rule.
#' @param min_width Minimum total width (bp) to flag a repeat of any period. Also the length
#'   from which a *merged* run is avoided in [tandem_repeats_to_avoid()].
#' @return `x` as a tibble with a logical `avoid` column.
#' @seealso [tandem_repeats_to_avoid()], [merge_tandem_repeats()]
#' @examples
#' bed <- data.frame(chrom = "Pf3D7_07_v3", start = c(100, 200, 300, 400),
#'                   end = c(110, 212, 320, 460),
#'                   name = c("Pf3D7_07_v3-100-110__A_x10",      # 10 bp homopolymer: kept
#'                            "Pf3D7_07_v3-200-212__AT_x6",      # 12 bp dinucleotide: flagged
#'                            "Pf3D7_07_v3-300-320__AAT_x6.67",  # 20 bp trinucleotide: kept
#'                            "Pf3D7_07_v3-400-460__TATTG_x12")) # 60 bp: flagged
#' flag_tandem_repeats(bed)[, c("repeat_unit", "period", "width", "avoid")]
#'
#' # merged first: an 8 bp A run and a 9 bp ATA run 3 bp apart make one 20 bp block, which
#' # the homopolymer rule (11 bp) catches although neither repeat is flagged alone
#' close <- data.frame(chrom = "Pf3D7_07_v3", start = c(100, 111), end = c(108, 120),
#'                     name = c("A_x8", "ATA_x3"))
#' flag_tandem_repeats(close)$avoid
#' flag_tandem_repeats(merge_tandem_repeats(close, gap = 10))[, c("width", "periods", "avoid")]
#' @export
flag_tandem_repeats <- function(x, min_width_by_period = c("1" = 11, "2" = 12, "3" = 21),
                                min_width = 50) {
  is_block <- is.data.frame(x) && all(c("periods", "width") %in% names(x)) &&
    !"period" %in% names(x)
  if (!is_block) x <- .as_tandem_repeats(x)
  thr <- .period_thresholds(min_width_by_period)
  if (!is.numeric(min_width) || length(min_width) != 1L || is.na(min_width))
    stop("`min_width` must be one number", call. = FALSE)
  w <- x$width
  avoid <- !is.na(w) & w >= min_width
  if (length(thr)) {
    t <- if (is_block) .block_threshold(x$periods, thr) else unname(thr[as.character(x$period)])
    avoid <- avoid | (!is.na(t) & w >= t)
  }
  x$avoid <- avoid
  x
}

# a merged block is judged by the lowest threshold among the periods it contains: a run
# holding a homopolymer is held to the homopolymer's length, whatever else is in it
.block_threshold <- function(periods, thr) {
  vapply(strsplit(as.character(periods), ",", fixed = TRUE), function(p) {
    t <- thr[p]
    if (!length(t) || all(is.na(t))) NA_real_ else min(t, na.rm = TRUE)
  }, numeric(1))
}

# --- merging --------------------------------------------------------------------------------
# Group intervals that overlap or lie within `gap` bases of each other, per chromosome.
# Returns an integer group id per input row, numbered in (chromosome, start) order. A
# running maximum of the end handles a short interval nested inside a long one: the next
# interval is compared to the furthest end seen so far, not the previous row's.
.merge_groups <- function(chr, start, end, gap = 0) {
  n <- length(start)
  if (!n) return(integer())
  if (!is.numeric(gap) || length(gap) != 1L || is.na(gap) || gap < 0)
    stop("`gap` must be one number >= 0", call. = FALSE)
  chr_f <- factor(chr, levels = unique(chr))          # keep the input's chromosome order
  o <- order(chr_f, start, end)
  cs <- chr_f[o]; ss <- start[o]; es <- end[o]
  cm <- stats::ave(es, cs, FUN = cummax)
  new_block <- c(TRUE, cs[-1] != cs[-n] | ss[-1] > cm[-n] + gap)
  grp <- integer(n)
  grp[o] <- cumsum(new_block)
  grp
}

#' Merge overlapping or nearby intervals
#'
#' A dependency-free `bedtools merge`: intervals on one chromosome that overlap, abut, or
#' lie within `gap` bases of each other become one. The complement of what [bed_subtract()]
#' does to a mask internally, exposed for when the merged blocks are what you want to look
#' at or write out.
#'
#' @section Coordinates: 0-based half-open `[start, end)`; abutting intervals
#'   (`end1 == start2`) are merged, as `bedtools merge` does by default. See
#'   [plasgenomicsutilsR-coordinates].
#'
#' @param x An interval table.
#' @param gap Intervals separated by at most this many bases are merged too (default `0`:
#'   only overlapping or abutting ones).
#' @param chrom,start,end Column names in `x` (`"chrom"` is accepted as an alias of the
#'   default `"chr"`).
#' @return A tibble of merged blocks in chromosome / start order: `chrom` (as `x` spells
#'   it, taken from the block's first interval), `chr` (normalised), `start`, `end`,
#'   `width`, and `n`, how many input intervals the block absorbed.
#' @seealso [merge_tandem_repeats()], which also summarises what was merged;
#'   [bed_subtract()], [bed_intersect()]
#' @examples
#' iv <- data.frame(chr = "7", start = c(100, 150, 300, 320), end = c(200, 250, 310, 330))
#' bed_merge(iv)              # 100-250 and two singletons
#' bed_merge(iv, gap = 10)    # 300-330 as well
#' @export
bed_merge <- function(x, gap = 0, chrom = "chr", start = "start", end = "end") {
  b <- .as_interval_table(x, chrom, start, end, what = "x")
  if (!nrow(b))
    return(tibble::tibble(chrom = character(), chr = character(), start = numeric(),
                          end = numeric(), width = numeric(), n = integer()))
  chr <- normalise_chr(b$chr)
  grp <- .merge_groups(chr, b$start, b$end, gap)
  first <- match(seq_len(max(grp)), grp)
  out <- tibble::tibble(
    chrom = b$chr[first], chr = chr[first],
    start = as.numeric(tapply(b$start, grp, min)),
    end = as.numeric(tapply(b$end, grp, max)),
    n = tabulate(grp))
  out$width <- out$end - out$start
  out[, c("chrom", "chr", "start", "end", "width", "n")]
}

# per-block concatenation of the distinct values, only doing work for blocks with > 1 member
.collapse_by <- function(v, grp, n_grp, sort_num = FALSE) {
  out <- rep(NA_character_, n_grp)
  size <- tabulate(grp, n_grp)
  single <- which(size == 1L)
  out[single] <- as.character(v[match(single, grp)])
  multi <- which(size > 1L)
  if (length(multi)) {
    keep <- grp %in% multi
    parts <- split(v[keep], grp[keep])
    out[multi] <- vapply(parts[as.character(multi)], function(z) {
      z <- unique(z[!is.na(z)])
      if (!length(z)) return(NA_character_)
      if (sort_num) z <- sort(z)
      paste(z, collapse = ",")
    }, character(1))
  }
  out
}

# The unit reduced to its period, then to the smallest of its rotations and of its reverse
# complement's rotations, so every spelling of one run compares equal: "ATAT" and "TA" ->
# "AT"; "ATA", "TAA" and the reverse complement "TAT" -> "AAT".
.canonical_unit <- function(unit, period) {
  core <- substr(unit, 1L, period)
  uc <- unique(core[!is.na(core)])
  if (!length(uc)) return(rep(NA_character_, length(unit)))
  rotations <- function(u) {
    n <- nchar(u)
    vapply(seq_len(n), function(i) paste0(substr(u, i, n), substr(u, 1L, i - 1L)), "")
  }
  canon <- vapply(uc, function(u) {
    rc <- .complement(paste(rev(strsplit(u, "", fixed = TRUE)[[1]]), collapse = ""))
    min(c(rotations(u), rotations(rc)))
  }, "")
  unname(canon[core])
}

# TRUE for a repeat lying entirely inside another on the same chromosome with the same
# canonical unit (an identical span counts, keeping the first). Such a record cannot change
# a flag or a block: the container has the same period and at least the width.
.contained_in_same_unit <- function(x) {
  canon <- .canonical_unit(x$repeat_unit, x$period)
  canon[is.na(canon)] <- ""                       # unparsed names only match each other
  o <- order(canon, x$chr, x$start, -x$end)
  k <- canon[o]; ch <- x$chr[o]; e <- x$end[o]
  n <- length(o)
  same <- c(FALSE, k[-1] == k[-n] & ch[-1] == ch[-n])
  cm <- stats::ave(e, cumsum(!same), FUN = cummax)
  prev_max <- c(-Inf, cm[-n]); prev_max[!same] <- -Inf
  drop <- logical(n); drop[o] <- e <= prev_max
  drop
}

#' Merge tandem repeats that run into one another
#'
#' Repeat finders report a stretch like `AAAAAAAAAAATATATATAT` as two records, an `A`
#' homopolymer and an `AT` run, that overlap or abut; one trinucleotide often flows into
#' another the same way. For avoiding slippage the whole stretch is one region, so this
#' joins records that overlap, abut, or lie within `gap` bases, and keeps a record of what
#' went into each block.
#'
#' @param x A [tandem_repeats()] table, or anything [tandem_repeats()] accepts. An `avoid`
#'   column ([flag_tandem_repeats()]) is carried through as *any member flagged*.
#' @param gap Records separated by at most this many bases are merged too (default `0`:
#'   only overlapping or abutting ones).
#' @param dedupe Drop a repeat that lies entirely inside another on the same chromosome
#'   whose unit is the same up to rotation or reverse complement -- `TAT` within a longer
#'   `ATA` run, `AT` within `ATAT` -- before merging (default `FALSE`). The two finders often report one run twice
#'   that way. Such a record adds nothing: the container has the same period and at least
#'   the width, so the blocks and every flag come out the same; only `n_repeats` and the
#'   `repeat_units` list, which no longer names both spellings, change.
#' @return A tibble of blocks in chromosome / start order: `chrom`, `chr`, `start`, `end`,
#'   `width`, `name` (`chrom-start-end`), `n_repeats` (records merged), `repeat_units`
#'   (the distinct units, comma-separated, left to right along the block), `periods` (the distinct
#'   periods, ascending), `min_period`, and `avoid` when `x` had it.
#' @seealso [tandem_repeats_to_avoid()], [bed_merge()] for plain intervals
#' @examples
#' bed <- data.frame(chrom = "Pf3D7_07_v3", start = c(100, 115, 300),
#'                   end = c(116, 130, 330),
#'                   name = c("Pf3D7_07_v3-100-116__A_x16",
#'                            "Pf3D7_07_v3-115-130__AT_x7.5",
#'                            "Pf3D7_07_v3-300-330__AAT_x10"))
#' merge_tandem_repeats(bed)[, c("start", "end", "n_repeats", "repeat_units", "periods")]
#'
#' # the same run reported as `A` and `AA`: dedupe keeps one spelling
#' twice <- rbind(bed, data.frame(chrom = "Pf3D7_07_v3", start = 100, end = 116,
#'                                name = "Pf3D7_07_v3-100-116__AA_x8"))
#' merge_tandem_repeats(twice)$repeat_units[1]
#' merge_tandem_repeats(twice, dedupe = TRUE)$repeat_units[1]
#' @export
merge_tandem_repeats <- function(x, gap = 0, dedupe = FALSE) {
  x <- .as_tandem_repeats(x)
  if (isTRUE(dedupe) && nrow(x)) x <- x[!.contained_in_same_unit(x), , drop = FALSE]
  if (!nrow(x)) {
    out <- tibble::tibble(chrom = character(), chr = character(), start = numeric(),
                          end = numeric(), width = numeric(), name = character(),
                          n_repeats = integer(), repeat_units = character(),
                          periods = character(), min_period = integer())
    if ("avoid" %in% names(x)) out$avoid <- logical()
    return(out)
  }
  grp <- .merge_groups(x$chr, x$start, x$end, gap)
  o <- order(grp, x$start, x$end)          # members left to right, so the unit list reads so
  x <- x[o, , drop = FALSE]; grp <- grp[o]
  n_grp <- max(grp)
  first <- match(seq_len(n_grp), grp)
  start <- as.numeric(tapply(x$start, grp, min))
  end <- as.numeric(tapply(x$end, grp, max))
  min_p <- suppressWarnings(as.integer(tapply(x$period, grp, min, na.rm = TRUE)))
  min_p[!is.finite(min_p)] <- NA_integer_
  out <- tibble::tibble(
    chrom = x$chrom[first], chr = x$chr[first], start = start, end = end,
    width = end - start,
    name = paste0(x$chrom[first], "-", format(start, scientific = FALSE, trim = TRUE), "-",
                  format(end, scientific = FALSE, trim = TRUE)),
    n_repeats = tabulate(grp, n_grp),
    repeat_units = .collapse_by(x$repeat_unit, grp, n_grp),
    periods = .collapse_by(x$period, grp, n_grp, sort_num = TRUE),
    min_period = min_p)
  if ("avoid" %in% names(x))
    out$avoid <- as.logical(tapply(x$avoid, grp, any))
  out
}

#' The tandem repeats a target design should avoid
#'
#' The whole recipe in one call: merge repeats that run into one another
#' ([merge_tandem_repeats()]), apply the period and length thresholds
#' ([flag_tandem_repeats()]), keep what is flagged. Subtract the result from your targets
#' with [bed_subtract()], padding by a few bases so a primer cannot end right at a repeat's
#' edge.
#'
#' `rule` says what the thresholds are applied to:
#'
#' * `"record"` (the default): each repeat is judged on its own width, then merged; a
#'   merged run is kept when any repeat in it was flagged, or when the run as a whole
#'   reaches `min_width`, since a chain of individually harmless short repeats adds up to a
#'   stretch that slips like a long one. This reproduces the mask the HEOME Pf3D7 design
#'   used (per-record thresholds; every repeat merged and runs of 50 bp or more kept; the two
#'   sets merged), except that a short unflagged repeat abutting a flagged one is absorbed
#'   into its block, extending it by at most those few bases.
#' * `"block"`: repeats are merged first (with `gap`) and each block is judged on its
#'   *combined* width, held to the lowest threshold among the periods it contains. An 8 bp
#'   `A` run within 10 bp of a 9 bp `ATA` run is neither flagged alone, but together with
#'   the bases between them they are a 20 bp stretch holding a homopolymer, and 20 is past
#'   the homopolymer's 11. The combined width is the block's span, gap included.
#'
#' @inheritParams flag_tandem_repeats
#' @inheritParams merge_tandem_repeats
#' @param rule `"record"` or `"block"`; see above.
#' @return The [merge_tandem_repeats()] table restricted to the blocks to avoid, without the
#'   `avoid` column.
#' @seealso [pf3d7_tandem_repeats()] for the bundled input, [bed_subtract()] for what to do
#'   with the output, [write_bed()] to write it.
#' @examples
#' avoid <- tandem_repeats_to_avoid(pf3d7_tandem_repeats())
#' nrow(avoid); sum(avoid$width)                # blocks, and bases masked
#'
#' # stricter: merge repeats within 10 bp and judge each block on its combined width,
#' # with homopolymers and dinucleotides held to the trinucleotide length
#' strict <- tandem_repeats_to_avoid(pf3d7_tandem_repeats(), rule = "block", gap = 10,
#'                                   min_width_by_period = c("1" = 21, "2" = 21, "3" = 21))
#' nrow(strict); sum(strict$width)
#'
#' # the targets: a few genes minus the repeats, with 10 bp of clearance around each
#' genes <- PF3D7_GENES[PF3D7_GENES$name %in% c("pfcrt", "pfdhfr", "pfkelch13"), ]
#' targets <- bed_subtract(genes, avoid, pad = 10)
#' targets[, c("name", "start", "end", "piece", "width")]
#' \dontrun{
#' write_bed(targets, "targets_no_repeats.bed")
#' }
#' @export
tandem_repeats_to_avoid <- function(x, min_width_by_period = c("1" = 11, "2" = 12, "3" = 21),
                                    min_width = 50, gap = 0, dedupe = FALSE,
                                    rule = c("record", "block")) {
  rule <- match.arg(rule)
  if (rule == "record") {
    flagged <- flag_tandem_repeats(x, min_width_by_period = min_width_by_period,
                                   min_width = min_width)
    blocks <- merge_tandem_repeats(flagged, gap = gap, dedupe = dedupe)
    keep <- blocks$avoid | blocks$width >= min_width
  } else {
    blocks <- merge_tandem_repeats(x, gap = gap, dedupe = dedupe)
    blocks <- flag_tandem_repeats(blocks, min_width_by_period = min_width_by_period,
                                  min_width = min_width)
    keep <- blocks$avoid
  }
  out <- blocks[keep, , drop = FALSE]
  out$avoid <- NULL
  out
}
