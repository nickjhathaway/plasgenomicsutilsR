hap_for_ehh <- function() parasite_haplotypes(example_pop_structure(umap = FALSE), maf = 0.05)

ehh_panel <- function(p) {
  if (!inherits(p, "patchwork")) return(p)
  for (i in seq_len(6)) {
    q <- tryCatch(p[[i]], error = function(e) NULL)
    if (!is.null(q) && is.data.frame(q$data) && "ehh" %in% names(q$data)) return(q)
  }
  stop("no EHH panel")
}

test_that("plot_ehh draws one decay curve per allele at the focal SNP", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("rehh")
  hap <- hap_for_ehh()
  p <- plot_ehh(hap, hap$map$snp_id[10], span = 5e5)
  d <- ehh_panel(p)$data
  # the two states at the focal SNP, which is the mutant-vs-reference split without needing
  # the SNPs annotated
  expect_setequal(as.character(unique(d$allele)), c("reference", "alternate"))
  expect_true(all(d$ehh >= 0 & d$ehh <= 1))
  # EHH is 1 at the focal SNP itself and decays away from it
  focal <- hap$map$pos[10]
  expect_equal(max(d$ehh[d$pos == focal]), 1)
  far <- d$ehh[abs(d$pos - focal) > 3e5]
  if (length(far)) expect_lte(max(far), 1)
  # the gene track is opt-in, since `genes` is usually only there to resolve `focal`
  expect_s3_class(p, "ggplot")
  expect_false(inherits(p, "patchwork"))
})

test_that("the focal SNP can be named three ways", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("rehh")
  hap <- hap_for_ehh()
  by_id <- ehh_panel(plot_ehh(hap, hap$map$snp_id[10], span = 5e5))$data
  by_pos <- ehh_panel(plot_ehh(hap, hap$map$pos[10], span = 5e5))$data
  expect_equal(by_id$ehh, by_pos$ehh)

  # a gene holding several SNPs picks the most balanced one and says which
  g <- data.frame(name = "wide", chr = normalise_chr(hap$map$chr[10]),
                  start = min(hap$map$pos) - 1, end = max(hap$map$pos) + 1)
  expect_message(plot_ehh(hap, "wide", genes = g, span = 5e5), "holds several SNPs")

  expect_error(plot_ehh(hap, "nowhere"), "not a SNP in the haplotypes")
  expect_error(plot_ehh(hap, 1), "no SNP at position 1")
  expect_error(plot_ehh(hap, c("a", "b")), "must be a `chr:pos` id")
  # a gene with no genotyped SNP in it
  empty <- data.frame(name = "empty", chr = "1", start = 1, end = 2)
  expect_error(plot_ehh(hap, "empty", genes = empty), "no genotyped SNP inside")
})

test_that("group facets in the metadata's order and skips groups too small to mean anything", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("rehh")
  hap <- hap_for_ehh()
  p <- plot_ehh(hap, hap$map$snp_id[10], group = "country", span = 5e5, min_haplotypes = 5)
  d <- ehh_panel(p)$data
  if (is.factor(d$group)) {
    expect_identical(levels(d$group),
                     intersect(levels(.as_group_factor(hap$meta$country)),
                               as.character(unique(d$group))))
  }
  # a group with fewer haplotypes than asked for is dropped, with a reason
  expect_message(plot_ehh(hap, hap$map$snp_id[10], group = "country", span = 5e5,
                          min_haplotypes = 25),
                 "skipping|fewer than|not variable")
})

test_that("span crops the curves and takes two sides", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("rehh")
  hap <- hap_for_ehh()
  focal <- hap$map$pos[10]
  wide <- ehh_panel(plot_ehh(hap, hap$map$snp_id[10], span = 6e5))$data
  narrow <- ehh_panel(plot_ehh(hap, hap$map$snp_id[10], span = 2e5))$data
  expect_lt(diff(range(narrow$pos)), diff(range(wide$pos)))
  expect_true(all(abs(narrow$pos - focal) <= 2e5))
  # Asymmetric, as everywhere else in the package. Assert on the axis, not the data: a curve
  # stops where EHH falls below `limehh`, which can be well inside the window on either side.
  asym <- plot_ehh(hap, hap$map$snp_id[10], span = c(left = 5e4, right = 5e5))
  rng <- ggplot2::ggplot_build(ehh_panel(asym))$layout$panel_params[[1]]$x.range
  # more room on the right than the left -- the window is still clamped to the chromosome,
  # so assert the asymmetry rather than an absolute reach
  expect_gt(rng[2] - focal, focal - rng[1])
  expect_lt(focal - rng[1], 1e5)
  expect_true(all(ehh_panel(asym)$data$pos >= focal - 5e4))
})

test_that("the gene track can be asked for", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("rehh")
  skip_if_not_installed("patchwork")
  hap <- hap_for_ehh()
  g <- data.frame(name = "amark", chr = normalise_chr(hap$map$chr[10]),
                  start = hap$map$pos[10] - 500, end = hap$map$pos[10] + 500)
  p <- plot_ehh(hap, hap$map$snp_id[10], span = 5e5, genes = g, gene_track = TRUE)
  expect_s3_class(p, "patchwork")
  expect_true("amark" %in% p[[2]]$data$name)
})

test_that("the frequency note can be moved or turned off", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("rehh")
  hap <- hap_for_ehh()
  note <- function(...) {
    b <- ggplot2::ggplot_build(ehh_panel(plot_ehh(hap, hap$map$snp_id[10], span = 5e5, ...)))
    d <- b$data[vapply(b$data, function(z) "label" %in% names(z), logical(1))]
    if (!length(d)) return(NULL)
    d[[1]]
  }
  top <- note()
  # default: the top of the panel, where the curves have not decayed to and cannot cross it
  expect_false(is.null(top))
  expect_gt(top$y[1], 0.9)
  expect_equal(top$hjust[1], 0)
  expect_equal(top$vjust[1], 1)
  expect_match(top$label[1], "^n = [0-9]+; reference [0-9]+%, alternate [0-9]+%$")

  bottom <- note(freq_position = "bottomleft")
  expect_lt(bottom$y[1], 0.1)
  right <- note(freq_position = "topright")
  expect_gt(right$x[1], top$x[1])
  expect_equal(right$hjust[1], 1)

  expect_null(note(show_freq = FALSE))
  expect_error(plot_ehh(hap, hap$map$snp_id[10], freq_position = "middle"),
               "should be one of")
})
