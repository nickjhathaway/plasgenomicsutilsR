#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom dplyr arrange mutate
#' @importFrom rlang .data
#' @importFrom tibble tibble
## usethis namespace: end
NULL

# bare column names used in a formula (aggregate) and a facet spec
utils::globalVariables(c("frac_pairs_ibd", "region_a", "region_b", "region"))
