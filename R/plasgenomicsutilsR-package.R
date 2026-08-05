#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom dplyr arrange mutate
#' @importFrom rlang .data
#' @importFrom tibble tibble
## usethis namespace: end
NULL

#' Genomic coordinate conventions
#'
#' @description
#' **Everything is 0-based.** There is one rule, so nothing has to be converted per call
#' and there is no part of the package you have to remember an exception for:
#'
#' \itemize{
#'   \item **Intervals are half-open `[start, end)`** -- the BED convention. This covers the
#'     bundled datasets ([PF3D7_GENES], [PF_EXAMPLE_DRUG_GENES], [PF3D7_CORE_REGIONS],
#'     [PF3D7_PARALOG_GENES]), any `genes` track you supply, the IBD `blocks`, `locus =`
#'     arguments, and [bed_intersect()]. `end - start` is the width in bp, and intervals
#'     that merely touch (`end1 == start2`) do not overlap.
#'   \item **Variant positions are 0-based too**, including the `pos` column and the
#'     `chr:pos` `snp_id` of the per-SNP tables, and the `snps =` arguments that select
#'     them. A variant sits inside a gene when `pos >= start & pos < end`.
#' }
#'
#' Formats that number differently are converted once, at the boundary, and never leak
#' inward: the PlasmoDB GFF (1-based inclusive) is shifted where the gene datasets are
#' built in `data-raw/`; `hmmibd-rs` blocks (0-based, both endpoints inclusive) get
#' `end + 1` when an [IbdResults] reads them; and VCF `POS` becomes `POS - 1` in the
#' companion Python package before any table reaches R.
#'
#' A `snp_id` is therefore `chr:pos0`, one less than the position the VCF or a genome
#' browser shows. Pass `--with-pos-vcf` to the Python tools to carry the 1-based position
#' alongside as `pos_vcf` when you want to look variants up by eye. Ids already present in
#' an input file are never trusted as keys, since `bcftools annotate --set-id` may have
#' written either `%POS` or `%POS0` and the file does not record which.
#'
#' @name plasgenomicsutilsR-coordinates
#' @keywords internal
NULL

# bare column names used in a formula (aggregate) / facet spec, and lazy-loaded datasets
utils::globalVariables(c("frac_pairs_ibd", "group_a", "group_b", "group",
                         "PF_EXAMPLE_DRUG_GENES", "PF3D7_GENES",
                         "PF3D7_CORE_REGIONS", "PF3D7_PARALOG_GENES"))
