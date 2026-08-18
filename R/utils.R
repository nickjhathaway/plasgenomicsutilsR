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
.normalise_meta <- function(meta, want = "sample") {
  if (is.null(meta) || !is.data.frame(meta)) return(meta)
  nms <- names(meta)
  if (want %in% nms) return(meta)
  hit <- which(tolower(nms) == tolower(want))
  if (!length(hit)) return(meta)
  if (length(hit) > 1)
    stop("`meta` has ", length(hit), " columns differing only in case for `", want, "`: ",
         paste(nms[hit], collapse = ", "), ". Rename all but one.", call. = FALSE)
  message("reading metadata column `", nms[hit], "` as `", want, "`")
  names(meta)[hit] <- want
  meta
}
