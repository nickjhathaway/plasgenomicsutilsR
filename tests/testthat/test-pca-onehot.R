# PCA and UMAP need a matrix whose columns each count copies of ONE allele. An allele-index
# panel is not that -- allele 2 is a different base, not two copies -- so something has to be
# derived, and the object used to derive the *biallelic dosage* view: every multiallelic site
# dropped. On the loci this work exists for (pfpx1 384, where two alternates of independent
# origin sit at one position) that drops exactly the site the analysis is about, and a panel
# that is entirely multiallelic could not be built at all.
#
# One-hot is the derivation that keeps them: one indicator column per ALT, reference column
# dropped. On a biallelic site that indicator column IS the alt-dosage column, so a wholly
# biallelic panel's PCA does not move -- which is what makes this safe as the default.

.mk_vcf <- function(n_bi = 6, n_multi = 2, n = 10) {
  vcf <- tempfile(fileext = ".vcf")
  rows <- character(0); p <- 1000
  for (i in seq_len(n_bi)) {
    rows <- c(rows, paste0("Pf3D7_01_v3\t", p, "\t.\tA\tG\t.\tPASS\t.\tGT\t",
                           paste(rep(c("0/0", "1/1"), length.out = n), collapse = "\t")))
    p <- p + 500
  }
  for (i in seq_len(n_multi)) {
    rows <- c(rows, paste0("Pf3D7_01_v3\t", p, "\t.\tA\tC,G\t.\tPASS\t.\tGT\t",
                           paste(rep(c("0/0", "1/1", "2/2"), length.out = n), collapse = "\t")))
    p <- p + 500
  }
  writeLines(c("##fileformat=VCFv4.2", "##contig=<ID=Pf3D7_01_v3,length=640851>",
               '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">',
               paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",
                      paste(sprintf("s%d", seq_len(n)), collapse = "\t")), rows), vcf)
  vcf
}

.idx_geno <- function(...) suppressMessages(
  load_genotypes(.mk_vcf(...), gds = tempfile(fileext = ".gds"), prune = FALSE,
                 variants = "all", encoding = "allele_index"))

test_that("a panel with nothing but multiallelic sites can be built at all", {
  skip_if_not_installed("SeqArray")
  # before: "panel \"full\" has no biallelic site in it, so no dosage view exists"
  ps <- suppressMessages(PopStructure$new(.idx_geno(n_bi = 0, n_multi = 6), meta = NULL))
  expect_equal(nrow(ps$pca_scores()), 10L)
})

test_that("a multiallelic site reaches PCA instead of being dropped", {
  skip_if_not_installed("SeqArray")
  ps <- suppressMessages(PopStructure$new(.idx_geno(), meta = NULL))
  auto <- ps$pca_scores()
  drop <- suppressMessages(PopStructure$new(.idx_geno(), meta = NULL,
                                            pca_panel = "dosage"))$pca_scores()
  # the two multiallelic sites split into four columns rather than vanishing, so the
  # ordinations are genuinely different -- this is not a relabelling
  expect_false(isTRUE(all.equal(unname(auto[, 1]), unname(drop[, 1]))))
})

test_that("a wholly biallelic panel's PCA does not move", {
  skip_if_not_installed("SeqArray")
  # the safety property behind making one-hot the default: on a biallelic site the indicator
  # column is the alt-dosage column, so there is nothing for one-hot to change
  g <- .idx_geno(n_bi = 8, n_multi = 0)
  a <- suppressMessages(PopStructure$new(g, meta = NULL, pca_panel = "onehot"))$pca_scores()
  b <- suppressMessages(PopStructure$new(g, meta = NULL, pca_panel = "dosage"))$pca_scores()
  expect_equal(unname(a), unname(b), tolerance = 1e-10)
})

test_that("UMAP runs on the same matrix as the PCA, not the raw index matrix", {
  skip_if_not_installed("SeqArray")
  skip_if_not_installed("uwot")
  ps <- suppressMessages(PopStructure$new(.idx_geno(), meta = NULL, n_pcs = 4))
  expect_no_error(suppressWarnings(ps$run_umap(pca_components = 2, n_neighbors = 3)))
  u <- ps$umap_df()
  expect_equal(nrow(u), 10L)
  expect_false(anyNA(u$UMAP1))
})

test_that("pca_panel = \"dosage\" is the way back to the old behaviour", {
  skip_if_not_installed("SeqArray")
  ps <- suppressMessages(PopStructure$new(.idx_geno(), meta = NULL, pca_panel = "dosage"))
  expect_equal(ncol(ps$genotype(needs = "dosage")), 6L)   # the biallelic sites only
})
