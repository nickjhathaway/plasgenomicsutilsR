# Build inst/extdata/pf3d7_tandem_repeats.rds, the compact form behind
# pf3d7_tandem_repeats(): every short tandem repeat of the Pf3D7 reference.
#
# Source: the `combined.bed` of the HEOME redesign
# (Pf_Determine_Regions_Of_Interest/Initial_Windows.qmd, April 2021) -- the union of
#
#   trf Pf3D7.fasta 2 7 7 80 10 50 1000 -f -d -m -h
#     | elucidator TandemRepeatFinderOutputToBed
#   elucidator findSimpleTandemRepeatLocations --fasta Pf3D7.fasta --maxRepeatUnitSize 10
#
# sorted and de-duplicated. Only the first four BED fields carry information: the name is
# `chrom-start-end__UNIT_xCOPIES`, so the file is stored as chrom / start / end / unit /
# copies (chrom and unit as factors, coordinates as integers) and the name is rebuilt on
# load. `copies` is the finder's own copy number; for 96% of records it is width /
# unit_size to the six significant digits the finder printed, and those are stored as NA
# and recomputed on load. The rest (trf records whose alignment had indels) are kept as
# written. xz then brings 18 MB of BED down to about 1.1 MB.

bed <- file.path("~/Dropbox/ownCloud/documents/plasmodium/falciparum/newRedesignHeome_2021_04",
                 "Pf_Determine_Regions_Of_Interest/tandemRepeatRegions/combined_1-4.bed")
PF3D7_TANDEM_REPEATS_VERSION <- "2021-04-06"   # date combined.bed was generated

x <- utils::read.delim(bed, header = FALSE, quote = "", comment.char = "",
                       colClasses = "character",
                       col.names = c("chrom", "start", "end", "name"))
stopifnot(all(grepl("__", x$name, fixed = TRUE)))
code <- sub("^.*__", "", x$name)
unit <- toupper(sub("_x.*$", "", code))
copies <- as.numeric(sub("^.*_x", "", code))
stopifnot(!anyNA(copies), all(nzchar(unit)))

chrom_levels <- unique(x$chrom)                      # file order: 01..14, API, MIT
start <- as.integer(x$start); end <- as.integer(x$end)
implied <- signif((end - start) / nchar(unit), 6)
derivable <- as.character(implied) == as.character(copies)
cat(sum(derivable), "of", length(copies), "copy numbers are width / unit_size; storing the rest\n")
z <- data.frame(
  chrom = factor(x$chrom, levels = chrom_levels),
  start = start,
  end = end,
  repeat_unit = factor(unit),
  copies = ifelse(derivable, NA_real_, copies),
  stringsAsFactors = FALSE
)
stopifnot(!anyNA(z$start), !anyNA(z$end), all(z$end > z$start))
attr(z, "source_version") <- PF3D7_TANDEM_REPEATS_VERSION

out <- file.path("inst", "extdata", "pf3d7_tandem_repeats.rds")
saveRDS(z, out, compress = "xz")
cat(nrow(z), "repeats ->", out, sprintf("(%.2f MB)\n", file.size(out) / 1e6))

# round trip: the rebuilt names must equal the file's, through the same code the loader uses
devtools::load_all(quiet = TRUE)
z2 <- readRDS(out)
name2 <- plasgenomicsutilsR:::.rebuild_repeat_names(z2)
stopifnot(identical(name2, x$name))
cat("names round-trip exactly\n")
