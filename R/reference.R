# Reference-genome registry (mirrors the Python plasgenomicsutils lib/reference.py).
# The only species/assembly-specific facts live here, keyed by a short id, so the
# plotting/analysis code stays species-agnostic.

#' Pf3D7 core chromosome lengths (bp)
#'
#' Named numeric vector of the 14 nuclear chromosome lengths for the Pf3D7
#' reference assembly, chromosome ids normalised to "1".."14".
#'
#' @export
PF3D7_CORE_CHROM_LENGTHS_BP <- c(
  "1" = 640851, "2" = 947102, "3" = 1067971, "4" = 1200490, "5" = 1343557,
  "6" = 1418242, "7" = 1445207, "8" = 1472805, "9" = 1541735, "10" = 1687656,
  "11" = 2038340, "12" = 2271494, "13" = 2925236, "14" = 3291936
)

#' Pf3D7 constant genetic-map rate (bp/cM)
#' @export
PF3D7_BP_PER_CM <- 15000

.REFERENCES <- list(
  pf3d7 = list(
    ref_id = "pf3d7",
    species = "Plasmodium falciparum",
    assembly = "Pf3D7",
    core_chrom_lengths_bp = PF3D7_CORE_CHROM_LENGTHS_BP,
    bp_per_cm = PF3D7_BP_PER_CM
  )
)

#' Default reference id
#' @export
DEFAULT_REFERENCE <- "pf3d7"

#' List available reference ids
#' @return Character vector of registered reference ids.
#' @examples
#' available_references()
#' @export
available_references <- function() {
  sort(names(.REFERENCES))
}

#' Look up a reference genome's facts by id
#'
#' @param ref_id Reference id (case-insensitive), e.g. "pf3d7".
#' @return A list with `ref_id`, `species`, `assembly`, `core_chrom_lengths_bp`,
#'   and `bp_per_cm`.
#' @examples
#' ref <- get_reference("pf3d7")
#' ref$bp_per_cm
#' head(ref$core_chrom_lengths_bp)
#' @export
get_reference <- function(ref_id = DEFAULT_REFERENCE) {
  key <- tolower(ref_id)
  if (!key %in% names(.REFERENCES)) {
    stop(sprintf("Unknown reference '%s'. Available: %s",
                 ref_id, paste(available_references(), collapse = ", ")),
         call. = FALSE)
  }
  .REFERENCES[[key]]
}

#' Normalise a chromosome name to a bare number string
#'
#' `"Pf3D7_07_v3"`, `"chr7"`, `"07"`, `7` all become `"7"`.
#'
#' @param c A chromosome name (character or numeric), scalar or vector.
#' @return Character vector of normalised chromosome ids.
#' @examples
#' # every spelling of a chromosome collapses to the same key, so tables from
#' # different tools join
#' normalise_chr(c("Pf3D7_07_v3", "chr7", "7", "Pf3D7_API_v3"))
#' @export
normalise_chr <- function(c) {
  s <- as.character(c)
  s <- sub("^(Pf3D7_|PvP01_|Pf_)0*([0-9]+).*$", "\\2", s)
  s <- sub("^chr", "", s)
  s <- sub("^0+([0-9])", "\\1", s)
  s[s == ""] <- "0"
  s
}
