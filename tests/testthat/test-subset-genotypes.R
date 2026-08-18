# Reusing the cohort a haplotype set kept. The scans run on monoclonals; a diversity or
# differentiation run that should describe the same samples needs the same set, and the
# full SNP panel rather than the scans' MAF-filtered one.

.hap_fixture <- function() {
  ps <- example_pop_structure("africa", umap = FALSE)
  G <- ps$genotype()
  set.seed(1)
  fws <- data.frame(sample = rownames(G),
                    fws = c(stats::runif(200, 0.93, 1), stats::runif(58, 0.4, 0.9)),
                    stringsAsFactors = FALSE)
  list(ps = ps, G = G, meta = ps$get_meta(),
       hap = parasite_haplotypes(G, meta = ps$get_meta(), fws = fws,
                                 min_fws = 0.92, maf = 0.03))
}

test_that("haplotype_samples returns what survived every filter, not just the Fws gate", {
  f <- .hap_fixture()
  ids <- haplotype_samples(f$hap)
  expect_identical(ids, rownames(f$hap$hap))
  expect_length(ids, nrow(f$G) - f$hap$filtering$n_dropped_polyclonal
                     - f$hap$filtering$n_dropped_sample_missing)
  expect_true(all(ids %in% rownames(f$G)))
  expect_error(haplotype_samples(f$G), "not a parasite_haplotypes")
})

test_that("subset_genotypes narrows samples and keeps every SNP", {
  f <- .hap_fixture()
  m <- suppressMessages(subset_genotypes(f$G, f$hap))
  expect_equal(nrow(m), length(haplotype_samples(f$hap)))
  expect_equal(ncol(m), ncol(f$G))                 # SNPs untouched -- this is not hap$hap
  expect_setequal(rownames(m), haplotype_samples(f$hap))
  # and the values are the original genotypes, not the imputed 0/1 haplotypes
  expect_equal(m[rownames(m)[1], ], f$G[rownames(m)[1], ])
  expect_gt(ncol(f$G), ncol(f$hap$hap))            # the panels really do differ
})

test_that("a load_genotypes list comes back a list, with its slots intact", {
  f <- .hap_fixture()
  geno <- list(genotype = f$G, sample.id = rownames(f$G), snp.id = colnames(f$G),
               allele = "alt", pruned = TRUE, positions = "0-based",
               variants = "biallelic_snvs")
  out <- suppressMessages(subset_genotypes(geno, f$hap))

  expect_type(out, "list")
  expect_equal(nrow(out$genotype), length(haplotype_samples(f$hap)))
  # sample.id has to move with the matrix or every later join is off
  expect_identical(out$sample.id, rownames(out$genotype))
  for (k in c("snp.id", "allele", "pruned", "positions", "variants"))
    expect_identical(out[[k]], geno[[k]])
})

test_that("row order follows the panel, so metadata joined to it stays in step", {
  f <- .hap_fixture()
  asked <- rev(haplotype_samples(f$hap))            # ask in a different order
  m <- suppressMessages(subset_genotypes(f$G, asked))
  expect_identical(rownames(m), rownames(f$G)[rownames(f$G) %in% asked])
})

test_that("samples can be named by any of the things that know them", {
  f <- .hap_fixture()
  ids <- haplotype_samples(f$hap)
  by_hap <- suppressMessages(subset_genotypes(f$G, f$hap))
  for (spec in list(ids,
                    f$ps$subset(samples = ids),
                    list(genotype = f$G[ids, , drop = FALSE], sample.id = ids),
                    data.frame(sample = ids, stringsAsFactors = FALSE))) {
    expect_equal(rownames(suppressMessages(subset_genotypes(f$G, spec))), rownames(by_hap))
  }
  expect_error(subset_genotypes(f$G, 1:3), "cannot read sample ids")
})

test_that("a PopStructure is narrowed as a PopStructure", {
  f <- .hap_fixture()
  out <- subset_genotypes(f$ps, f$hap)
  expect_s3_class(out, "PopStructure")
  expect_setequal(out$get_samples(), haplotype_samples(f$hap))
  expect_length(f$ps$get_samples(), nrow(f$G))     # the original is untouched
})

test_that("unknown samples warn, or error under strict", {
  f <- .hap_fixture()
  one <- rownames(f$G)[1]
  expect_warning(subset_genotypes(f$G, c(one, "ghost")), "not in the panel")
  expect_error(subset_genotypes(f$G, c(one, "ghost"), strict = TRUE), "not in the panel")
  expect_error(suppressWarnings(subset_genotypes(f$G, "ghost")), "none of those samples")
  expect_error(subset_genotypes(f$G, character(0)), "names no samples")
})

test_that("the subset flows into pop_diversity and changes the cohort it describes", {
  f <- .hap_fixture()
  mono <- suppressMessages(subset_genotypes(f$G, f$hap))
  a <- pop_diversity(mono, group = "region", meta = f$meta, by = "genome")
  b <- pop_diversity(f$G, group = "region", meta = f$meta, by = "genome")
  expect_true(all(a$n_samples <= b$n_samples))
  expect_true(any(a$n_samples < b$n_samples))
  # metadata needs no subsetting: it is matched by sample
  c2 <- pop_diversity(mono, group = "region",
                      meta = f$meta[f$meta$sample %in% rownames(mono), ], by = "genome")
  expect_equal(a, c2)
})
