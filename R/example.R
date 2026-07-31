# Bundled example IBD dataset (see inst/extdata) and a small drug-gene track, for
# documentation, tests, and quick plotting.

#' Pf3D7 drug-resistance gene coordinates
#'
#' A small `data.frame` of well-known *P. falciparum* drug-resistance genes
#' (`name`, `chr`, `start`, `end`; 1-based, Pf3D7 assembly) for use as the `genes`
#' track of an [IbdResults] object (gene reference lines on Manhattan plots, and
#' the drug-gene triangles). These are public reference-genome coordinates.
#'
#' @format A data frame with columns `name`, `chr`, `start`, `end`.
#' @export
EXAMPLE_DRUG_GENES <- data.frame(
  name  = c("crt", "dhfr", "mdr1", "dhps", "kelch13"),
  chr   = c("7", "4", "5", "8", "13"),
  start = c(403222L, 748088L, 957890L, 548200L, 1724817L),
  end   = c(406317L, 749914L, 962149L, 550616L, 1726997L),
  stringsAsFactors = FALSE
)

#' Load the bundled example IBD results
#'
#' Returns an [IbdResults] built from the small public example dataset shipped in
#' `inst/extdata`: per-SNP-per-region IBD, pairwise-region IBD, and the per-region
#' selection statistic (with thresholds) for five African countries, plus the
#' [EXAMPLE_DRUG_GENES] track. Use it to try the `plot_*()` functions without your
#' own data.
#'
#' The data are derived from publicly available *P. falciparum* whole-genome samples
#' (five countries), run through the `plasgenomicsutils` IBD tools and downsampled
#' to a small SNP panel for a compact fixture.
#'
#' @return An [IbdResults] object.
#' @examples
#' ibd <- example_ibd_results()
#' ibd
#' @export
example_ibd_results <- function() {
  ex <- function(f) {
    p <- system.file("extdata", f, package = "plasgenomicsutilsR")
    if (!nzchar(p)) stop("example data file not found: ", f, call. = FALSE)
    p
  }
  thr <- utils::read.delim(ex("ibd_selection_per_region_threshold.tsv"),
                           stringsAsFactors = FALSE)
  ibd_results(
    per_snp_region  = ex("ibd_per_snp_region.tsv.gz"),
    pairwise_region = ex("ibd_pairwise_region.tsv.gz"),
    selection       = ex("ibd_selection_per_region.tsv.gz"),
    threshold       = thr,
    genes           = EXAMPLE_DRUG_GENES,
    reference       = "pf3d7"
  )
}
