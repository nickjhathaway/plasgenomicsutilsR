# Rsb and XP-EHH at a marker with more than two alleles.
#
# The Phase 0 warning on these two said their score was "a two-allele contrast on a subset of
# the haplotypes", copied from the iHS wording. **That was wrong**, and these tests are what
# establishes it: `scan_hh()`'s `IES` and `INES` at a four-allele marker are identical to
# `calc_ehhs()` computed over every haplotype. The two-allele reduction in `scan_hh` touches
# only the reported FREQ_ columns, which `.cross_pop()` never returns. No haplotype is
# dropped and no allele is ignored.
#
# There *is* a real multiallelic problem here, and it is a different one. `iES` is a
# site-level homozygosity pooled across allele classes, so it falls as the focal is split into
# more of them -- measured at 36% for three alleles and 44% for four, on identical haplotypes.
# XP-EHH reads `iES`, so comparing a marker between two populations that carry different
# numbers of alleles there puts the two integrals on different footings. Rsb reads `iNES`,
# which normalises by the focal homozygosity and moves by 1% or less.

.xp_hap <- function(counts, m = 21, seed = 21, focal_col = NULL) {
  set.seed(seed)
  n <- if (is.null(focal_col)) sum(counts) else length(focal_col)
  G <- matrix(stats::rbinom(n * m, 1, 0.4), n, m)
  G[, (m + 1) %/% 2] <- if (is.null(focal_col))
    rep(seq_along(counts) - 1L, times = counts) else focal_col
  rownames(G) <- sprintf("s%02d", seq_len(n))
  colnames(G) <- paste0("Pf3D7_01_v3:", seq(1000, by = 1000, length.out = m))
  h <- parasite_haplotypes(G, maf = 0.02, alleles = "index", impute = FALSE)
  # `parasite_haplotypes()` leaves `meta` NULL when none is supplied, so the two populations
  # the cross-population scans need are attached here
  h$meta <- data.frame(sample = rownames(h$hap),
                       grp = rep(c("a", "b"), each = nrow(h$hap) / 2),
                       stringsAsFactors = FALSE)
  h
}

# Group a carries two alleles at the focal, group b carries four -- which is the shape the
# XP-EHH confound is about, and which alternating group labels would hide because both halves
# would then span every allele.
.xp_split <- function(seed = 21) {
  .xp_hap(NULL, seed = seed,
          focal_col = c(rep(0:1, each = 20), rep(0:3, each = 10)))
}

test_that("scan_hh's iES and iNES use every haplotype, whatever the allele count", {
  skip_if_not_installed("rehh")
  # the fact the corrected warning rests on
  hap <- .xp_hap(c(40, 20, 10, 10))
  rows <- seq_len(nrow(hap$hap))
  objs <- plasgenomicsutilsR:::.haplohh_list(hap, rows)
  o <- objs[["Pf3D7_01_v3"]]
  s <- rehh::scan_hh(o, discard_integration_at_border = FALSE)
  mrk <- match(11000, o@positions)
  e <- rehh::calc_ehhs(o, mrk = mrk, discard_integration_at_border = FALSE)
  r <- which(s$POSITION == 11000)
  expect_equal(s$IES[r], e$IES, tolerance = 1e-9)
  expect_equal(s$INES[r], e$INES, tolerance = 1e-9)
  # and the reported frequencies really do sum to less than one, which is what misled us
  fq <- as.numeric(s[r, grep("^FREQ_", names(s))])
  expect_lt(sum(fq), 0.95)
})

test_that("iES falls sharply with the allele count and iNES does not", {
  skip_if_not_installed("rehh")
  # identical haplotypes throughout; the only thing that changes is how many classes the
  # focal marker is split into
  ies <- ines <- numeric(0)
  for (counts in list(c(40, 40), c(40, 20, 20), c(40, 20, 10, 10))) {
    hap <- .xp_hap(counts)
    objs <- plasgenomicsutilsR:::.haplohh_list(hap, seq_len(nrow(hap$hap)))
    s <- rehh::scan_hh(objs[["Pf3D7_01_v3"]], discard_integration_at_border = FALSE)
    r <- which(s$POSITION == 11000)
    ies <- c(ies, s$IES[r]); ines <- c(ines, s$INES[r])
  }
  expect_lt(ies[2] / ies[1], 0.8)          # three alleles: a large drop
  expect_lt(ies[3] / ies[1], 0.7)          # four: larger still
  expect_equal(ines[2] / ines[1], 1, tolerance = 0.05)
  expect_equal(ines[3] / ines[1], 1, tolerance = 0.05)
})

test_that("the warning says what is actually wrong, and only for XP-EHH", {
  skip_if_not_installed("rehh")
  hap <- .xp_split()
  w <- tryCatch(run_xpehh(hap, group = "grp"), warning = function(e) conditionMessage(e))
  expect_true(is.character(w))
  expect_match(w, "allele count")
  expect_false(grepl("subset of the haplotypes", w))
  expect_match(w, "XP-EHH")

  # Rsb reads iNES, which the allele count barely moves, so it must not carry the warning
  expect_no_warning(run_rsb(hap, group = "grp"))
})

test_that("the scan reports each population's allele count so the confound is filterable", {
  skip_if_not_installed("rehh")
  hap <- .xp_split()
  r <- suppressWarnings(run_xpehh(hap, group = "grp"))
  expect_true(all(c("n_alleles_pop1", "n_alleles_pop2") %in% names(r)))
  focal <- r[r$snp_id == "Pf3D7_01_v3:11000", ]
  expect_equal(focal$n_alleles_pop1, 2L)
  expect_equal(focal$n_alleles_pop2, 4L)
  # a caller can drop the markers where the two populations disagree
  expect_true(is.integer(r$n_alleles_pop1))
})

test_that("a biallelic panel is silent and unchanged", {
  skip_if_not_installed("rehh")
  hap <- .xp_hap(c(40, 40))
  expect_no_warning(run_xpehh(hap, group = "grp"))
  expect_no_warning(run_rsb(hap, group = "grp"))
  r <- run_xpehh(hap, group = "grp")
  expect_true(all(r$n_alleles_pop1 <= 2L, na.rm = TRUE))
})
