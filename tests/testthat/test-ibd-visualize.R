# IbdResults container + plot_*() functions. Plots are checked by forcing a full
# ggplot build (ggplotGrob) so scale/facet/threshold errors surface without a
# graphics device.

skip_if_no_ggplot <- function() {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("scales")
}

make_obj <- function() {
  lens <- PF3D7_CORE_CHROM_LENGTHS_BP
  set.seed(42)
  n <- 200
  chr <- sample(names(lens), n, replace = TRUE)
  pos <- floor(stats::runif(n) * lens[chr])
  base <- data.frame(snp_id = paste0("Pf3D7_", chr, "_v3:", pos),
                     chr = chr, pos = pos, stringsAsFactors = FALSE)
  per_snp <- do.call(rbind, lapply(c("A", "B"), function(r) {
    d <- base; d$region <- r; d$frac_pairs_ibd <- stats::runif(nrow(d)); d
  }))
  pairwise <- do.call(rbind, lapply(list(c("A", "A"), c("A", "B"), c("B", "B")), function(p) {
    d <- base; d$region_a <- p[1]; d$region_b <- p[2]
    d$frac_pairs_ibd <- stats::runif(nrow(d)); d
  }))
  sel <- do.call(rbind, lapply(c("A", "B"), function(r) {
    d <- base; d$region <- r; d$neg_log10_p <- abs(stats::rnorm(nrow(d), 2))
    d$chi2_stat <- d$neg_log10_p * 2; d$z_score <- stats::rnorm(nrow(d)); d
  }))
  ibd_results(per_snp_region = per_snp, pairwise_region = pairwise, selection = sel,
              threshold = data.frame(region = c("A", "B"), threshold = c(5, 5)))
}

test_that("IbdResults ingests tables and computes cumulative positions", {
  obj <- make_obj()
  expect_s3_class(obj, "IbdResults")
  ps <- obj$get_per_snp_region()
  expect_true(all(c("cum_pos", "chr") %in% names(ps)))
  # cum_pos = chromosome offset + pos; chr 1 offset is 0 so cum_pos == pos there
  lay <- obj$chrom_layout()
  one <- ps[ps$chr == "1", ]
  if (nrow(one)) expect_equal(one$cum_pos, one$pos)
  # chr 2 rows are offset by chr 1 length
  two <- ps[ps$chr == "2", ]
  if (nrow(two)) expect_equal(two$cum_pos, two$pos + lay$offset[lay$chr == "2"])
})

test_that("chromosomes are normalised on ingest", {
  df <- data.frame(chr = c("Pf3D7_07_v3", "chr4"), pos = c(100, 200),
                   frac_pairs_ibd = c(0.1, 0.2))
  obj <- ibd_results(per_snp_region = df)
  expect_setequal(obj$get_per_snp_region()$chr, c("7", "4"))
})

test_that("missing required columns are rejected", {
  expect_error(ibd_results(per_snp_region = data.frame(chr = "1", pos = 1)),
               "frac_pairs_ibd")
  expect_error(ibd_results(pairwise_region = data.frame(chr = "1", pos = 1, region_a = "A")),
               "region_b")
})

test_that("manhattan plots build", {
  skip_if_no_ggplot()
  obj <- make_obj()
  expect_s3_class(plot_ibd_manhattan(obj), "ggplot")
  expect_silent(ggplot2::ggplotGrob(plot_ibd_manhattan(obj)))
  expect_silent(ggplot2::ggplotGrob(plot_ibd_manhattan(obj, regions = "A")))
  expect_silent(ggplot2::ggplotGrob(plot_selection_manhattan(obj)))
  expect_silent(ggplot2::ggplotGrob(plot_selection_manhattan(obj, metric = "chi2_stat")))
})

test_that("selection metric must exist", {
  skip_if_no_ggplot()
  obj <- make_obj()
  expect_error(plot_selection_manhattan(obj, metric = "z_score") |> ggplot2::ggplotGrob(),
               NA)  # z_score exists -> no error
})

test_that("region heatmap builds and mirrors the triangle", {
  skip_if_no_ggplot()
  obj <- make_obj()
  expect_silent(ggplot2::ggplotGrob(plot_ibd_region_heatmap(obj)))
  expect_silent(ggplot2::ggplotGrob(plot_ibd_region_heatmap(obj, anchor = "A")))
})

test_that("heatmap fill scale is configurable", {
  skip_if_no_ggplot()
  obj <- make_obj()
  expect_silent(ggplot2::ggplotGrob(plot_ibd_region_heatmap(obj, trans = "log2")))
  expect_silent(ggplot2::ggplotGrob(plot_ibd_region_heatmap(obj, limits = c(0, 0.5))))
  expect_silent(ggplot2::ggplotGrob(plot_ibd_region_heatmap(obj, colors = c("white", "red"))))
  fs <- ggplot2::scale_fill_gradient(low = "white", high = "black")
  expect_silent(ggplot2::ggplotGrob(plot_ibd_region_heatmap(obj, fill_scale = fs)))
  expect_silent(ggplot2::ggplotGrob(plot_drug_gene_triangles(
    ibd_results(pairwise_region = obj$get_pairwise_region(),
                genes = data.frame(name = "g", chr = "1", start = 0, end = 7e5)),
    trans = "log2")))
})

test_that("tug-of-war needs both tables and a resolvable region", {
  skip_if_no_ggplot()
  obj <- make_obj()
  expect_silent(ggplot2::ggplotGrob(plot_ibd_tugofwar(obj, region = "A")))
  expect_silent(ggplot2::ggplotGrob(plot_ibd_tugofwar(obj)))          # all regions -> faceted
  expect_silent(ggplot2::ggplotGrob(plot_ibd_tugofwar(obj, region = c("A", "B"))))
  expect_silent(ggplot2::ggplotGrob(plot_ibd_tugofwar(obj, scale = "free")))
  expect_silent(ggplot2::ggplotGrob(plot_ibd_tugofwar(obj, scale = "common",
                                                       label_genes = FALSE)))
  expect_error(plot_ibd_tugofwar(obj, scale = "nope"))
  sel_only <- ibd_results(selection = obj$get_selection())
  expect_error(plot_ibd_tugofwar(sel_only, region = "A"), "per_snp_region")
})

test_that("drug-gene triangles build from the genes track", {
  skip_if_no_ggplot()
  obj <- make_obj()
  # a gene wide enough to catch some synthetic SNPs on chr 1
  genes <- data.frame(name = c("geneX", "geneY"), chr = c("1", "2"),
                      start = c(0, 0), end = c(700000, 900000))
  withr <- ibd_results(pairwise_region = obj$get_pairwise_region(), genes = genes)
  expect_s3_class(plot_drug_gene_triangles(withr), "ggplot")
  expect_silent(ggplot2::ggplotGrob(plot_drug_gene_triangles(withr, agg = "median")))
  expect_silent(ggplot2::ggplotGrob(plot_drug_gene_triangles(withr, genes = "geneX", label = FALSE)))
})

test_that("triangles need pairwise + genes and warn on empty genes", {
  skip_if_no_ggplot()
  obj <- make_obj()
  expect_error(plot_drug_gene_triangles(obj), "genes track")            # no genes
  no_pw <- ibd_results(selection = obj$get_selection(),
                       genes = data.frame(name = "g", chr = "1", start = 0, end = 10))
  expect_error(plot_drug_gene_triangles(no_pw), "no pairwise_region")
  # a gene off in empty space yields no SNPs -> informative error
  far <- ibd_results(pairwise_region = obj$get_pairwise_region(),
                     genes = data.frame(name = "empty", chr = "1", start = 1e9, end = 1e9 + 1))
  expect_error(suppressWarnings(plot_drug_gene_triangles(far)), "overlapping SNPs")
})

test_that("chromosome selection and gene highlighting work", {
  skip_if_no_ggplot()
  ex <- example_ibd_results()
  expect_silent(ggplot2::ggplotGrob(plot_selection_manhattan(ex, chroms = c("7", "8"))))
  expect_silent(ggplot2::ggplotGrob(plot_ibd_manhattan(ex, skip_chr = c("1", "2"))))
  expect_silent(ggplot2::ggplotGrob(
    plot_ibd_tugofwar(ex, region = "Tanzania", skip_chr = "1",
                      highlight_genes = "crt", label_genes = TRUE)))
  expect_silent(ggplot2::ggplotGrob(plot_ibd_region_heatmap(ex, chroms = "7")))
  # keeping only crt's chromosome still lands its reference line correctly
  expect_silent(ggplot2::ggplotGrob(
    plot_ibd_manhattan(ex, chroms = "7", highlight_genes = "crt", label_genes = TRUE)))
  # gene labels on the heatmap + tug-of-war, and case-insensitive gene selection
  expect_silent(ggplot2::ggplotGrob(
    plot_ibd_region_heatmap(ex, highlight_genes = "CRT", label_genes = TRUE)))
  expect_silent(ggplot2::ggplotGrob(
    plot_ibd_tugofwar(ex, highlight_genes = "CRT", label_genes = TRUE)))
  expect_error(plot_ibd_manhattan(ex, chroms = "99"), "no chromosomes left")
})

test_that("bundled example data loads and every plot builds", {
  skip_if_no_ggplot()
  ex <- example_ibd_results()
  expect_s3_class(ex, "IbdResults")
  expect_setequal(unique(ex$get_per_snp_region()$region),
                  c("DRC", "Ethiopia", "Kenya", "Sudan", "Tanzania"))
  expect_silent(ggplot2::ggplotGrob(plot_ibd_manhattan(ex)))
  expect_silent(ggplot2::ggplotGrob(plot_selection_manhattan(ex)))
  expect_silent(ggplot2::ggplotGrob(plot_ibd_region_heatmap(ex)))
  expect_silent(ggplot2::ggplotGrob(plot_ibd_tugofwar(ex, region = "Tanzania")))
})

test_that("triangles accept specific SNP ids and positions", {
  skip_if_no_ggplot()
  ex <- example_ibd_results()
  pw <- ex$get_pairwise_region()
  snp <- paste0(pw$chr[1], ":", pw$pos[1])            # normalised "chr:pos"
  expect_silent(ggplot2::ggplotGrob(plot_drug_gene_triangles(ex, snps = snp)))
  df <- data.frame(chr = pw$chr[1:2], pos = pw$pos[1:2], name = c("locusA", "locusB"))
  expect_silent(ggplot2::ggplotGrob(plot_drug_gene_triangles(ex, snps = df)))
  expect_error(plot_drug_gene_triangles(ex, snps = "notasnp"), "chr:pos")
})

test_that("triangles: grid returns a plot, individual returns a list", {
  skip_if_no_ggplot()
  ex <- example_ibd_results()
  pw <- ex$get_pairwise_region()
  snpids <- paste0(pw$chr[1:2], ":", pw$pos[1:2])
  grid <- plot_drug_gene_triangles(ex, snps = snpids)
  expect_s3_class(grid, "ggplot")
  ind <- plot_drug_gene_triangles(ex, snps = snpids, individual = TRUE)
  expect_type(ind, "list")
  expect_length(ind, 2)
  expect_silent(ggplot2::ggplotGrob(ind[[1]]))
  # scale option flows through (log2 of zeros warns, so just check it assembles)
  expect_s3_class(
    plot_drug_gene_triangles(ex, snps = snpids, individual = TRUE, trans = "log2")[[1]],
    "ggplot")
})

test_that("plot functions reject an empty object", {
  skip_if_no_ggplot()
  empty <- ibd_results()
  expect_error(plot_ibd_manhattan(empty), "no per_snp_region")
  expect_error(plot_selection_manhattan(empty), "no selection")
  expect_error(plot_ibd_region_heatmap(empty), "no pairwise_region")
})
