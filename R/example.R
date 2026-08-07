# Bundled example IBD dataset (see inst/extdata), for documentation, tests, and quick
# plotting. The gene tracks (PF3D7_GENES, PF_EXAMPLE_DRUG_GENES) are lazy-loaded
# datasets documented in data.R and built by data-raw/PF3D7_GENES.R.

#' Load the bundled example IBD results
#'
#' Returns an [IbdResults] built from the small public example dataset shipped in
#' `inst/extdata`: per-SNP-per-group IBD, pairwise-group IBD, and the per-group
#' selection statistic (with thresholds) for five African countries, plus the
#' [PF_EXAMPLE_DRUG_GENES] track. Use it to try the `plot_*()` functions without your
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
  thr <- utils::read.delim(ex("ibd_selection_per_group_threshold.tsv"),
                           stringsAsFactors = FALSE)
  ibd_results(
    per_snp_group  = ex("ibd_per_snp_group.tsv.gz"),
    pairwise_group = ex("ibd_pairwise_group.tsv.gz"),
    selection      = ex("ibd_selection_per_group.tsv.gz"),
    threshold      = thr,
    genes          = PF_EXAMPLE_DRUG_GENES,
    reference      = "pf3d7"
  )
}
