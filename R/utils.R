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
