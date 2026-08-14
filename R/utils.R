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
