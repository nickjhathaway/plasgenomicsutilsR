# Per-sample allele states read straight from a callset.

#' What each sample carries at one variant
#'
#' The set of alleles each sample's genotype names at one position, as a named character
#' vector -- exactly the `sample -> state` shape [ibd_block_extension_by_allele()] takes for
#' its `allele =` argument.
#'
#' This closes the gap that made a carrier contrast awkward to reach. `allele =` has always
#' accepted a named vector, but nothing in the package produced one: the only code that read a
#' per-sample allele set was private and wired into [plot_region_haplotypes()], so "the D384A
#' carriers" was reachable from a hand-built metadata column and from nowhere else.
#'
#' @section Why a set and not a dosage:
#' A dosage says how many alternate copies a sample has, not **which** alternate. At a codon
#' carrying three independently arisen changes that is the whole question, so the state is the
#' allele set: `"C"`, `"G"`, or `"C + G"` for a mixed infection carrying both. A sample with no
#' call is `NA` and is left out of every stratum downstream.
#'
#' @param vcf Path to a VCF/BCF. Needs `bcftools` on `PATH`.
#' @param position `"chr:pos"`. **0-based**, like every position in this package
#'   (`?"plasgenomicsutilsR-coordinates"`); pass `one_based = TRUE` to give a VCF `POS`
#'   instead. The chromosome spelling is normalised, so `13` and `Pf3D7_13_v3` both work.
#' @param names How to label a state. `"base"` (default) uses the alleles themselves, which is
#'   what a carrier contrast wants -- `carrier = "C"` names something the reader can check
#'   against the callset. `"index"` uses the positional wording a legend wants
#'   (`"reference"`, `"alternate 1"`), which is what [plot_region_haplotypes()] shows.
#' @param one_based Read `position` as a 1-based VCF `POS`.
#' @return A named character vector, one entry per sample in the callset, `NA` where the
#'   genotype is missing.
#' @seealso [ibd_block_extension_by_allele()], [plot_region_haplotypes()]
#' @examples
#' \dontrun{
#' st <- allele_states("px1.bcf", "Pf3D7_13_v3:1725591")
#' table(st)
#' ibd_block_extension_by_allele(ibd, loci, allele = st, carrier = "C", reference = "A")
#' }
#' @export
allele_states <- function(vcf, position, names = c("base", "index"),
                          one_based = FALSE) {
  names <- match.arg(names)
  loc <- .parse_snp_ids(position)
  if (nrow(loc) != 1L)
    stop("`position` must be a single \"chr:pos\" id", call. = FALSE)
  pos1 <- as.integer(loc$pos) + if (isTRUE(one_based)) 0L else 1L

  sets <- .read_genotype_sets(vcf, names = names)
  want <- paste0(loc$chr, ":", pos1 - 1L)
  # the callset's own chromosome spelling need not match the caller's, so match on the
  # normalised name rather than the literal one
  have <- .parse_snp_ids(colnames(sets$codes))
  hit <- which(normalise_chr(have$chr) == normalise_chr(loc$chr) &
                 as.integer(have$pos) == pos1 - 1L)
  if (!length(hit))
    stop("no record at ", want, " in ", basename(vcf),
         ". Positions in this package are 0-based; pass one_based = TRUE for a VCF POS.",
         call. = FALSE)

  lv <- sets$levels[[hit[1]]]
  code <- sets$codes[, hit[1]]
  out <- ifelse(is.na(code), NA_character_, lv[code + 1L])
  stats::setNames(out, rownames(sets$codes))
}
