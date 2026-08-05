#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom dplyr arrange mutate
#' @importFrom rlang .data
#' @importFrom tibble tibble
## usethis namespace: end
NULL

# bare column names used in a formula (aggregate) / facet spec, and lazy-loaded datasets
utils::globalVariables(c("frac_pairs_ibd", "group_a", "group_b", "group",
                         "PF_EXAMPLE_DRUG_GENES", "PF3D7_GENES",
                         "PF3D7_CORE_REGIONS", "PF3D7_PARALOG_GENES"))
