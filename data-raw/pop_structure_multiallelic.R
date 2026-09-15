# Builds `inst/extdata/pop_structure_multiallelic.rds`, the example dataset that actually has
# multiallelic SNPs in it.
#
# The other two fixtures predate the SeqArray backend and were made ad hoc, with no script;
# they are biallelic dosage matrices and cannot demonstrate anything this work is about. This
# one is built from the same public Pf7 Ghana/Cambodia callset the companion Python package
# ships, through `load_genotypes()` itself, so it carries a `sites` table and an allele-index
# panel and stays reproducible.
#
# Run from the package root with the sibling Python repo checked out beside it:
#   Rscript data-raw/pop_structure_multiallelic.R

suppressMessages({
  library(SeqArray)
  devtools::load_all(".", quiet = TRUE)
})

bcf <- normalizePath(file.path("..", "plasgenomicsutils", "tests", "data",
                               "ghana_cambodia.pf7.tiny.bcf"), mustWork = TRUE)
stopifnot(nzchar(Sys.which("bcftools")))

tmp <- tempfile(fileext = ".vcf.gz")
# Trim alleles no genotype carries, then keep SNP records only -- the same reading the
# companion package's `biallelic_snp_filter --snps-only --no-biallelic` uses, so the fixture
# is what that pipeline would actually hand to R. A `*` disqualifies a record, so what
# survives is SNPs whose every allele is a real base.
system2("bash", c("-c", shQuote(paste(
  "bcftools view --trim-alt-alleles", shQuote(bcf), "-Ou |",
  "bcftools view -V indels,mnps,other,ref,bnd -m2 -e 'ALT=\"*\"' -Oz -o", shQuote(tmp)))))

gds <- tempfile(fileext = ".gds")
idx <- load_genotypes(tmp, gds = gds, prune = FALSE, variants = "all",
                      encoding = "allele_index")
dos <- load_genotypes(tmp, gds = gds, prune = FALSE)          # the biallelic dosage view

meta <- data.frame(sample = idx$sample.id,
                   country = ifelse(seq_along(idx$sample.id) <= 30, "Ghana", "Cambodia"),
                   stringsAsFactors = FALSE)

message(sprintf("allele_index panel: %d SNPs, %d multiallelic",
                ncol(idx$genotype), sum(idx$sites$n_alt_real > 1)))
message(sprintf("dosage panel:       %d SNPs", ncol(dos$genotype)))
stopifnot(sum(idx$sites$n_alt_real > 1) >= 20)

saveRDS(list(index = idx, dosage = dos, meta = meta),
        file.path("inst", "extdata", "pop_structure_multiallelic.rds"),
        compress = "xz")
message("wrote inst/extdata/pop_structure_multiallelic.rds")
