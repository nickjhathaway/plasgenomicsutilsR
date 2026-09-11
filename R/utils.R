#' Require an optional (Suggests) package
#'
#' Stops with an actionable message if a package listed in Suggests is needed at
#' runtime but not installed.
#'
#' @param pkg Package name.
#' @param what Short description of what needs it (used in the error message).
#' @noRd
.need_package <- function(pkg, what = "this function") {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required for %s but is not installed.", pkg, what),
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Accept either spelling of a colour argument
#'
#' The package's plotting arguments were a mix of `colour`/`colours` and `color`/`colors`,
#' sometimes both inside one function. Every such argument now has a formal for each
#' spelling: the canonical one carries the default, the alias defaults to `NULL`, and this
#' resolves them by looking at how the caller actually wrote the call -- so an argument
#' whose default is a real colour (`border_colour = "black"`) still works, which comparing
#' the two values could not tell apart from the caller passing that colour by hand.
#'
#' Reading the call is what distinguishes "the caller wrote the default" from "the caller
#' wrote nothing", but a call reached through `...` (an `R6` `$new()`, say) cannot be
#' matched; there the non-`NULL` alias simply wins, which is the same answer in every case
#' but the one where both spellings were named at once.
#'
#' @param canonical,alias The two argument names, as strings.
#' @param env The calling function's frame.
#' @return The value of whichever spelling was used, or the canonical default.
#' @noRd
.alias_arg <- function(canonical, alias, env = parent.frame()) {
  used <- tryCatch(names(match.call(sys.function(sys.parent()), sys.call(sys.parent()))),
                   error = function(e) NULL)
  if (!is.null(used)) {
    if (canonical %in% used && alias %in% used)
      stop("give either `", canonical, "` or `", alias, "`, not both", call. = FALSE)
    return(get(if (alias %in% used) alias else canonical, envir = env))
  }
  a <- get(alias, envir = env)
  if (is.null(a)) get(canonical, envir = env) else a
}

#' Accept any capitalisation of a metadata table's `sample` column
#'
#' Metadata arrives from whoever assembled it, and `Sample`, `sample` and `SAMPLE` are the
#' same column to everyone except a string comparison. Renames whichever case-variant a
#' table uses to the canonical name, so nothing downstream has to ask twice.
#'
#' Two columns differing only in case is an error rather than a coin toss: which one holds
#' the ids is not ours to guess. A table with neither is left alone -- the function that
#' needs the column raises its own error, which knows what it is for.
#'
#' @param meta A data frame, or `NULL`.
#' @param want The canonical column name (default `"sample"`).
#' @return `meta`, with the column renamed if it needed it.
#' @noRd
# A sample id is a label, not a number. An all-numeric cohort -- micronix ids, barcodes,
# anything without a letter -- reads back from a file as numeric, and R then does something
# worse than failing: `named_vector[c(3, 1)]` indexes by *position*, so a lookup meant to be
# by name silently returns whoever happens to sit in those rows. Ids are made character at
# every point one enters the package, so that can never arise.
.as_id_chr <- function(x) {
  if (is.character(x)) return(x)
  if (is.factor(x)) return(as.character(x))
  if (!is.numeric(x)) return(as.character(x))
  # format(), not as.character(): a large id would otherwise come back as "4.06e+09"
  out <- format(x, scientific = FALSE, trim = TRUE)
  out[is.na(x)] <- NA_character_
  out
}

.normalise_meta <- function(meta, want = "sample") {
  if (is.null(meta) || !is.data.frame(meta)) return(meta)
  nms <- names(meta)
  if (want %in% nms) {
    meta[[want]] <- .as_id_chr(meta[[want]])
    return(meta)
  }
  hit <- which(tolower(nms) == tolower(want))
  if (!length(hit)) return(meta)
  if (length(hit) > 1)
    stop("`meta` has ", length(hit), " columns differing only in case for `", want, "`: ",
         paste(nms[hit], collapse = ", "), ". Rename all but one.", call. = FALSE)
  message("reading metadata column `", nms[hit], "` as `", want, "`")
  names(meta)[hit] <- want
  meta[[want]] <- .as_id_chr(meta[[want]])
  meta
}


# --- carrying sample metadata onto a result table ------------------------------------------
# Shared by ibd_pair_clusters(), ibd_pair_links() and gene_ibd_pairs(), so a bad column name
# is reported the same way wherever `add_meta_cols` is offered. `meta` arrives normalised --
# each exported function does that itself.
.meta_cols_for <- function(cols, meta, new, out) {
  cols <- as.character(cols)
  if (is.null(meta) || !is.data.frame(meta))
    stop("`add_meta_cols` needs meta; build with ibd_results(meta = ), or pass meta = ",
         call. = FALSE)
  if (!"sample" %in% names(meta))
    stop("`meta` has no `sample` column to match the rows on", call. = FALSE)
  miss <- setdiff(cols, names(meta))
  if (length(miss))
    stop("meta has no column ", paste0("'", miss, "'", collapse = ", "), ". Available: ",
         paste(setdiff(names(meta), "sample"), collapse = ", "), call. = FALSE)
  clash <- intersect(new(cols), names(out))
  if (length(clash))
    stop("`add_meta_cols` would overwrite ", paste0("'", clash, "'", collapse = ", "),
         ", which this table already reports; rename the column in `meta` first.",
         call. = FALSE)
  cols
}

# One column per metadata column, keyed on the table's `sample`. Indexing by match() rather
# than merging keeps the row order and a factor's level order, and gives a sample that is not
# in the metadata an NA instead of dropping its row.
.add_sample_meta <- function(out, cols, meta) {
  if (is.null(cols)) return(out)
  cols <- .meta_cols_for(cols, meta, identity, out)
  i <- match(out$sample, meta$sample)
  for (col in cols) out[[col]] <- meta[[col]][i]
  out
}

# The same, but for a table whose rows are pairs: each column lands on both endpoints.
.add_endpoint_meta <- function(out, cols, meta) {
  if (is.null(cols)) return(out)
  cols <- .meta_cols_for(cols, meta,
                         function(z) c(paste0("sample1_", z), paste0("sample2_", z)), out)
  i1 <- match(out$sample1, meta$sample)
  i2 <- match(out$sample2, meta$sample)
  for (col in cols) {
    out[[paste0("sample1_", col)]] <- meta[[col]][i1]
    out[[paste0("sample2_", col)]] <- meta[[col]][i2]
  }
  out
}

# ---- genotype encodings ----------------------------------------------------

#' @keywords internal
#' @noRd
# A dosage matrix (0/1/2 alt copies) and an allele-index matrix (0..k-1, naming *which*
# allele) are both integer matrices, and at a triallelic marker coded 0/1/2 they are
# indistinguishable by inspection. So the object carries its own `encoding` and this reads
# it, falling back to a value check for a bare matrix -- which catches the unambiguous half,
# an index above 2.
#
# The point is to refuse rather than to mangle. `.haploid_calls()` maps 0 -> 0 and 2 -> 1 and
# leaves the rest NA, so an index panel comes back with allele 2 renamed allele 1, allele 1
# turned to missing and allele 3 turned to missing: a matrix of the right shape and the wrong
# contents, which nothing downstream can detect.
.require_dosage <- function(x, what) {
  enc <- if (is.list(x) && !is.null(x$encoding)) x$encoding else NULL
  mat <- if (is.list(x) && !is.null(x$genotype)) x$genotype else x
  if (!is.null(enc) && !identical(enc, "dosage"))
    stop(what, " needs alt-allele dosages (0/1/2), and this panel is encoded \"", enc,
         "\". A dosage cannot say which of several alternates a call carries, so there is ",
         "no faithful conversion -- contrast one allele at a time instead.", call. = FALSE)
  if (is.null(enc) && is.numeric(mat) && length(mat)) {
    mx <- suppressWarnings(max(mat, na.rm = TRUE))
    if (is.finite(mx) && mx > 2)
      stop(what, " needs alt-allele dosages (0/1/2); this matrix holds ", mx,
           ", which is an allele index rather than a dosage. Pass a dosage panel, or ",
           "record the panel's `encoding` if it really is one.", call. = FALSE)
  }
  invisible(TRUE)
}

# `snp_id` -> chr and pos, for a key with any number of trailing fields.
#
# The one-hot expansion names a split multiallelic column `chr:pos:allele`, and the Python
# package writes `chr:pos:ref:alt` wherever an allele has to be named. Parsing from the RIGHT
# ("everything before the last colon is the chromosome") reads `Pf3D7_01_v3:100:A:G` as
# chromosome "Pf3D7_01_v3:100:A" at position NA, so the whole scan is refused. Parse from the
# LEFT instead: field 1 is the chromosome, field 2 is the position, and anything after names
# an allele. No Pf contig name contains a colon, so nothing is lost by fixing the split.
.snp_id_chr_pos <- function(id, what = "`snp_id`") {
  parts <- strsplit(as.character(id), ":", fixed = TRUE)
  if (any(lengths(parts) < 2L))
    stop("could not read a position out of ", what,
         "; expected \"chr:pos\" (optionally followed by allele fields)", call. = FALSE)
  pos <- suppressWarnings(as.numeric(vapply(parts, `[`, character(1), 2L)))
  if (anyNA(pos))
    stop("could not read a position out of ", what,
         "; expected \"chr:pos\" (optionally followed by allele fields)", call. = FALSE)
  list(chr = normalise_chr(vapply(parts, `[`, character(1), 1L)), pos = pos)
}
