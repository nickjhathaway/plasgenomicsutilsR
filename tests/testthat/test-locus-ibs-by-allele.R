# A controllable two-population panel. Focal chromosome "7"; background chromosomes 1-3 carry
# the null windows. Carriers of both populations can be made to share one haplotype (positive
# difference of differences) or to carry population-specific haplotypes (negative).

AC <- sprintf("ac%02d", 1:6); AR <- sprintf("ar%02d", 1:6)
BC <- sprintf("bc%02d", 1:6); BR <- sprintf("br%02d", 1:6)
SAMPLES <- c(AC, AR, BC, BR)

# 4 chromosomes x 40 SNPs, positions uniform so every count-matched window has the same span
CHR <- rep(c("7", "1", "2", "3"), each = 40)
POS <- rep(seq(1000, 40000, by = 1000), 4)
SNPIDS <- paste0(CHR, ":", POS)
FOCAL_COLS <- which(CHR == "7")

mk_panel <- function(carrier_focal, seed = 1) {
  set.seed(seed)
  G <- matrix(sample(0:1, length(SAMPLES) * length(CHR), replace = TRUE),
              nrow = length(SAMPLES), dimnames = list(SAMPLES, SNPIDS))
  # overwrite the focal chromosome for carriers per the scenario; references stay random
  G[AC, FOCAL_COLS] <- carrier_focal$A
  G[BC, FOCAL_COLS] <- carrier_focal$B
  G
}

META <- data.frame(sample = SAMPLES,
                   region = rep(c("A", "B"), each = 12),
                   pin = rep(c("Present", "Absent", "Present", "Absent"), each = 6),
                   stringsAsFactors = FALSE)
MAP <- data.frame(chr = CHR, pos = POS, stringsAsFactors = FALSE)
LOCUS <- data.frame(name = "sweep", chr = "7", start = 20000, end = 21000,
                    stringsAsFactors = FALSE)

run_ibs <- function(G, ...) {
  locus_ibs_by_allele(G, LOCUS, allele = "pin", carrier = "Present", reference = "Absent",
                      group_a = "A", group_b = "B", meta = META, map = MAP,
                      n_snps = 8L, min_sites = 4L, ...)
}

test_that("a haplotype shared across both populations gives a large positive difference", {
  # both populations' carriers all-0 across the focal chromosome: cross-carrier identity = 1
  G <- mk_panel(list(A = 0L, B = 0L))
  r <- run_ibs(G)
  expect_equal(nrow(r), 1L)
  expect_equal(r$pairs, "all pairs")
  expect_equal(r$cross_car, 1, tolerance = 1e-9)   # identical carrier haplotype across A and B
  expect_lt(r$cross_ref, 0.9)                       # references are assorted
  expect_gt(r$d_ind, 0.3)
  expect_gt(r$z, 3)
  # the observed difference beats every null window, so p sits at the two-sided floor
  # 2/(null_n + 1) -- a resolution limit, not an exact probability
  expect_equal(r$p_two, 2 / (r$null_n + 1))
  expect_equal(r$a_car, 1, tolerance = 1e-9)        # each population's carriers cohere internally
  expect_equal(r$b_car, 1, tolerance = 1e-9)
  expect_gt(r$null_n, 5)
  expect_equal(length(attr(r, "null")[[1]]), r$null_n)
})

test_that("population-specific carrier haplotypes give a negative difference", {
  # A carriers all-0, B carriers all-1: internally tight, but nothing alike across the border
  G <- mk_panel(list(A = 0L, B = 1L))
  r <- run_ibs(G)
  expect_equal(r$cross_car, 0, tolerance = 1e-9)    # 0 vs 1 never agree
  expect_equal(r$a_car, 1, tolerance = 1e-9)
  expect_equal(r$b_car, 1, tolerance = 1e-9)
  expect_lt(r$d_ind, 0)
  expect_lt(r$z, 0)
})

test_that("the difference of differences cancels a shared population baseline", {
  # carriers and references drawn from the same pool (no focal signal): d_ind near 0
  set.seed(7)
  G <- matrix(sample(0:1, length(SAMPLES) * length(CHR), replace = TRUE),
              nrow = length(SAMPLES), dimnames = list(SAMPLES, SNPIDS))
  r <- run_ibs(G)
  # no focal signal, so the carrier and reference strata shift together and the difference
  # of differences cancels to near zero (z on 15 random windows is too noisy to bound)
  expect_lt(abs(r$d_ind), 0.25)
})

test_that("a load_genotypes-style list is accepted and matches the matrix path", {
  G <- mk_panel(list(A = 0L, B = 0L))
  panel <- list(genotype = G, encoding = "allele_index")   # colnames already chr:pos
  r_list <- locus_ibs_by_allele(panel, LOCUS, allele = "pin", carrier = "Present",
                                reference = "Absent", group_a = "A", group_b = "B",
                                meta = META, n_snps = 8L, min_sites = 4L)
  r_mat <- run_ibs(G)
  expect_equal(r_list$d_ind, r_mat$d_ind, tolerance = 1e-9)
  expect_equal(r_list$cross_car, r_mat$cross_car, tolerance = 1e-9)
})

test_that("drop_ibd adds a second pair set and removes the named pairs", {
  G <- mk_panel(list(A = 0L, B = 0L))
  # mark every A-carrier x B-carrier pair as IBD across the locus
  grid <- expand.grid(sample1 = AC, sample2 = BC, stringsAsFactors = FALSE)
  # >= 15 kb so the block survives the IbdResults short-segment filter, and overlapping the locus
  blocks <- data.frame(grid, chr = "7", start = 10000, end = 30000,
                       different = 0, Nsnp = 40, stringsAsFactors = FALSE)
  ibd <- ibd_results(blocks = blocks, meta = META, group_col_in_meta = "region")
  r <- run_ibs(G, drop_ibd = ibd)
  expect_equal(nrow(r), 2L)
  expect_setequal(r$pairs, c("all pairs", "IBD pairs dropped"))
  # with every cross-carrier pair excluded there is no carrier contrast left to compute
  dropped <- r[r$pairs == "IBD pairs dropped", ]
  expect_true(is.na(dropped$cross_car))
})

test_that("it refuses to guess a second state at a multiallelic locus", {
  G <- mk_panel(list(A = 0L, B = 0L))
  st <- stats::setNames(rep(c("X", "Y", "Z"), length.out = length(SAMPLES)), SAMPLES)
  expect_error(
    locus_ibs_by_allele(G, LOCUS, allele = st, carrier = "X", group_a = "A", group_b = "B",
                        meta = META, map = MAP, n_snps = 8L),
    "name both")
})
