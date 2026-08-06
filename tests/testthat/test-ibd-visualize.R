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
    d <- base; d$group <- r; d$frac_pairs_ibd <- stats::runif(nrow(d)); d
  }))
  pairwise <- do.call(rbind, lapply(list(c("A", "A"), c("A", "B"), c("B", "B")), function(p) {
    d <- base; d$group_a <- p[1]; d$group_b <- p[2]
    d$frac_pairs_ibd <- stats::runif(nrow(d)); d
  }))
  sel <- do.call(rbind, lapply(c("A", "B"), function(r) {
    d <- base; d$group <- r; d$neg_log10_p <- abs(stats::rnorm(nrow(d), 2))
    d$chi2_stat <- d$neg_log10_p * 2; d$z_score <- stats::rnorm(nrow(d)); d
  }))
  ibd_results(per_snp_group = per_snp, pairwise_group = pairwise, selection = sel,
              threshold = data.frame(group = c("A", "B"), threshold = c(5, 5)))
}

test_that("IbdResults ingests tables and computes cumulative positions", {
  obj <- make_obj()
  expect_s3_class(obj, "IbdResults")
  ps <- obj$get_per_snp_group()
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
  obj <- ibd_results(per_snp_group = df)
  expect_setequal(obj$get_per_snp_group()$chr, c("7", "4"))
})

test_that("missing required columns are rejected", {
  expect_error(ibd_results(per_snp_group = data.frame(chr = "1", pos = 1)),
               "frac_pairs_ibd")
  expect_error(ibd_results(pairwise_group = data.frame(chr = "1", pos = 1, group_a = "A")),
               "group_b")
})

test_that("manhattan plots build", {
  skip_if_no_ggplot()
  obj <- make_obj()
  expect_s3_class(plot_ibd_sharing_manhattan(obj), "ggplot")
  expect_silent(ggplot2::ggplotGrob(plot_ibd_sharing_manhattan(obj)))
  expect_silent(ggplot2::ggplotGrob(plot_ibd_sharing_manhattan(obj, groups = "A")))
  expect_silent(ggplot2::ggplotGrob(plot_selection_manhattan(obj)))
  expect_silent(ggplot2::ggplotGrob(plot_selection_manhattan(obj, metric = "chi2_stat")))
})

test_that("non-finite per-group thresholds are skipped, not warned", {
  skip_if_no_ggplot()
  # A group with no valid SNPs gets a NaN Bonferroni threshold (an empty cell in
  # ibd_selection_statistic's per-group threshold file). geom_hline must not warn
  # ("Removed N rows ..."); finite thresholds still draw, NaN ones are dropped.
  base <- make_obj()$get_selection()
  base <- base[, c("chr", "pos", "group", "neg_log10_p")]
  n_hlines <- function(p) {
    b <- ggplot2::ggplot_build(p)
    sum(vapply(b$data, function(d)
      if ("yintercept" %in% names(d)) sum(is.finite(d$yintercept)) else 0L, integer(1)))
  }

  all_nan <- ibd_results(selection = base,
    threshold = data.frame(group = c("A", "B"), threshold = c(NA_real_, NA_real_)))
  p_nan <- plot_selection_manhattan(all_nan, draw_threshold = TRUE)
  expect_silent(ggplot2::ggplotGrob(p_nan))
  expect_equal(n_hlines(p_nan), 0)

  mixed <- ibd_results(selection = base,
    threshold = data.frame(group = c("A", "B"), threshold = c(5, NA_real_)))
  p_mixed <- plot_selection_manhattan(mixed, draw_threshold = TRUE)
  expect_silent(ggplot2::ggplotGrob(p_mixed))
  expect_equal(n_hlines(p_mixed), 1)
})

test_that("selection metric must exist", {
  skip_if_no_ggplot()
  obj <- make_obj()
  expect_error(plot_selection_manhattan(obj, metric = "z_score") |> ggplot2::ggplotGrob(),
               NA)  # z_score exists -> no error
})

test_that("group heatmap builds and mirrors the triangle", {
  skip_if_no_ggplot()
  obj <- make_obj()
  expect_silent(ggplot2::ggplotGrob(plot_ibd_pairwise_group_heatmap(obj)))
  expect_silent(ggplot2::ggplotGrob(plot_ibd_pairwise_group_heatmap(obj, anchor = "A")))
})

test_that("heatmap fill scale is configurable", {
  skip_if_no_ggplot()
  obj <- make_obj()
  expect_silent(ggplot2::ggplotGrob(plot_ibd_pairwise_group_heatmap(obj, trans = "log2")))
  expect_silent(ggplot2::ggplotGrob(plot_ibd_pairwise_group_heatmap(obj, limits = c(0, 0.5))))
  expect_silent(ggplot2::ggplotGrob(plot_ibd_pairwise_group_heatmap(obj, colors = c("white", "red"))))
  fs <- ggplot2::scale_fill_gradient(low = "white", high = "black")
  expect_silent(ggplot2::ggplotGrob(plot_ibd_pairwise_group_heatmap(obj, fill_scale = fs)))
  expect_silent(ggplot2::ggplotGrob(plot_pairwise_ibd_for_genes(
    ibd_results(pairwise_group = obj$get_pairwise_group(),
                genes = data.frame(name = "g", chr = "1", start = 0, end = 7e5)),
    trans = "log2")))
})

test_that("tug-of-war needs both tables and a resolvable group", {
  skip_if_no_ggplot()
  obj <- make_obj()
  expect_silent(ggplot2::ggplotGrob(plot_ibd_tugofwar(obj, group = "A")))
  expect_silent(ggplot2::ggplotGrob(plot_ibd_tugofwar(obj)))          # all groups -> faceted
  expect_silent(ggplot2::ggplotGrob(plot_ibd_tugofwar(obj, group = c("A", "B"))))
  expect_silent(ggplot2::ggplotGrob(plot_ibd_tugofwar(obj, scale = "free")))
  expect_silent(ggplot2::ggplotGrob(plot_ibd_tugofwar(obj, scale = "common",
                                                       label_genes = FALSE)))
  expect_error(plot_ibd_tugofwar(obj, scale = "nope"))
  sel_only <- ibd_results(selection = obj$get_selection())
  expect_error(plot_ibd_tugofwar(sel_only, group = "A"), "per_snp_group")
})

test_that("drug-gene triangles build from the genes track", {
  skip_if_no_ggplot()
  obj <- make_obj()
  # a gene wide enough to catch some synthetic SNPs on chr 1
  genes <- data.frame(name = c("geneX", "geneY"), chr = c("1", "2"),
                      start = c(0, 0), end = c(700000, 900000))
  withr <- ibd_results(pairwise_group = obj$get_pairwise_group(), genes = genes)
  expect_s3_class(plot_pairwise_ibd_for_genes(withr), "ggplot")
  expect_silent(ggplot2::ggplotGrob(plot_pairwise_ibd_for_genes(withr, agg = "median")))
  expect_silent(ggplot2::ggplotGrob(plot_pairwise_ibd_for_genes(withr, genes = "geneX", label = FALSE)))
})

test_that("triangles need pairwise + genes and warn on empty genes", {
  skip_if_no_ggplot()
  obj <- make_obj()
  expect_error(plot_pairwise_ibd_for_genes(obj), "genes track")            # no genes
  no_pw <- ibd_results(selection = obj$get_selection(),
                       genes = data.frame(name = "g", chr = "1", start = 0, end = 10))
  expect_error(plot_pairwise_ibd_for_genes(no_pw), "no pairwise_group")
  # a gene off in empty space yields no SNPs -> informative error
  far <- ibd_results(pairwise_group = obj$get_pairwise_group(),
                     genes = data.frame(name = "empty", chr = "1", start = 1e9, end = 1e9 + 1))
  expect_error(suppressWarnings(plot_pairwise_ibd_for_genes(far)), "overlapping SNPs")
})

test_that("chromosome selection and gene highlighting work", {
  skip_if_no_ggplot()
  ex <- example_ibd_results()
  expect_silent(ggplot2::ggplotGrob(plot_selection_manhattan(ex, chroms = c("7", "8"))))
  expect_silent(ggplot2::ggplotGrob(plot_ibd_sharing_manhattan(ex, skip_chr = c("1", "2"))))
  expect_silent(ggplot2::ggplotGrob(
    plot_ibd_tugofwar(ex, group = "Tanzania", skip_chr = "1",
                      highlight_genes = "pfcrt", label_genes = TRUE)))
  expect_silent(ggplot2::ggplotGrob(plot_ibd_pairwise_group_heatmap(ex, chroms = "7")))
  # keeping only crt's chromosome still lands its reference line correctly
  expect_silent(ggplot2::ggplotGrob(
    plot_ibd_sharing_manhattan(ex, chroms = "7", highlight_genes = "pfcrt", label_genes = TRUE)))
  # gene labels on the heatmap + tug-of-war, and case-insensitive gene selection
  expect_silent(ggplot2::ggplotGrob(
    plot_ibd_pairwise_group_heatmap(ex, highlight_genes = "pfcrt", label_genes = TRUE)))
  expect_silent(ggplot2::ggplotGrob(
    plot_ibd_tugofwar(ex, highlight_genes = "pfcrt", label_genes = TRUE)))
  expect_error(plot_ibd_sharing_manhattan(ex, chroms = "99"), "no chromosomes left")
})

test_that("bundled example data loads and every plot builds", {
  skip_if_no_ggplot()
  ex <- example_ibd_results()
  expect_s3_class(ex, "IbdResults")
  expect_setequal(unique(ex$get_per_snp_group()$group),
                  c("DRC", "Ethiopia", "Kenya", "Sudan", "Tanzania"))
  expect_silent(ggplot2::ggplotGrob(plot_ibd_sharing_manhattan(ex)))
  expect_silent(ggplot2::ggplotGrob(plot_selection_manhattan(ex)))
  expect_silent(ggplot2::ggplotGrob(plot_ibd_pairwise_group_heatmap(ex)))
  expect_silent(ggplot2::ggplotGrob(plot_ibd_tugofwar(ex, group = "Tanzania")))
})

test_that("triangles accept specific SNP ids and positions", {
  skip_if_no_ggplot()
  ex <- example_ibd_results()
  pw <- ex$get_pairwise_group()
  snp <- paste0(pw$chr[1], ":", pw$pos[1])            # normalised "chr:pos"
  expect_silent(ggplot2::ggplotGrob(plot_pairwise_ibd_for_genes(ex, snps = snp)))
  df <- data.frame(chr = pw$chr[1:2], pos = pw$pos[1:2], name = c("locusA", "locusB"))
  expect_silent(ggplot2::ggplotGrob(plot_pairwise_ibd_for_genes(ex, snps = df)))
  expect_error(plot_pairwise_ibd_for_genes(ex, snps = "notasnp"), "chr:pos")
})

test_that("triangles: grid returns a plot, individual returns a list", {
  skip_if_no_ggplot()
  ex <- example_ibd_results()
  pw <- ex$get_pairwise_group()
  snpids <- paste0(pw$chr[1:2], ":", pw$pos[1:2])
  grid <- plot_pairwise_ibd_for_genes(ex, snps = snpids)
  expect_s3_class(grid, "ggplot")
  ind <- plot_pairwise_ibd_for_genes(ex, snps = snpids, individual = TRUE)
  expect_type(ind, "list")
  expect_length(ind, 2)
  expect_silent(ggplot2::ggplotGrob(ind[[1]]))
  # scale option flows through (log2 of zeros warns, so just check it assembles)
  expect_s3_class(
    plot_pairwise_ibd_for_genes(ex, snps = snpids, individual = TRUE, trans = "log2")[[1]],
    "ggplot")
})

test_that("plot functions reject an empty object", {
  skip_if_no_ggplot()
  empty <- ibd_results()
  expect_error(plot_ibd_sharing_manhattan(empty), "no per_snp_group")
  expect_error(plot_selection_manhattan(empty), "no selection")
  expect_error(plot_ibd_pairwise_group_heatmap(empty), "no pairwise_group")
})

test_that("requesting a gene not in the track errors (not a silent no-op)", {
  skip_if_no_ggplot()
  ex <- example_ibd_results()
  expect_error(plot_ibd_sharing_manhattan(ex, highlight_genes = "notagene"),
               "not in the track")
  expect_error(plot_selection_manhattan(ex, highlight_genes = c("pfcrt", "nope")), "nope")
  expect_error(plot_ibd_pairwise_group_heatmap(ex, highlight_genes = "nope"), "not in the track")
  expect_error(plot_pairwise_ibd_for_genes(ex, genes = "nope"), "not in the track")
  # no gene track at all -> a request errors informatively
  obj <- make_obj()
  expect_error(plot_ibd_sharing_manhattan(obj, highlight_genes = "pfcrt"), "no gene track")
})

test_that("gene labels turn on automatically when highlight_genes is given", {
  skip_if_no_ggplot()
  ex <- example_ibd_results()
  has_text <- function(p) any(vapply(ggplot2::ggplot_build(p)$data,
    function(d) "label" %in% names(d), logical(1)))
  expect_true(has_text(plot_selection_manhattan(ex, highlight_genes = "pfcrt")))  # auto
  expect_false(has_text(plot_selection_manhattan(ex)))                            # none
  expect_false(has_text(plot_selection_manhattan(ex, highlight_genes = "pfcrt",  # forced off
                                                 label_genes = FALSE)))
})

test_that("IbdResults exposes plot_*() and pos_selection_genes() as methods", {
  skip_if_no_ggplot()
  ex <- example_ibd_results()
  expect_s3_class(ex$plot_ibd_sharing_manhattan(), "ggplot")
  expect_s3_class(ex$plot_selection_manhattan(), "ggplot")
  expect_s3_class(ex$plot_ibd_tugofwar(group = "Tanzania"), "ggplot")
  expect_s3_class(ex$plot_ibd_pairwise_group_heatmap(), "ggplot")
  expect_s3_class(ex$pos_selection_genes(), "data.frame")
})

test_that("pos_selection_genes maps significant SNPs onto genes within a window", {
  ex <- example_ibd_results()
  hits <- pos_selection_genes(ex)                       # 2 kb window, drug-gene track
  expect_s3_class(hits, "data.frame")
  expect_true(all(c("group", "gene_id", "name", "min_distance", "peak_pos") %in% names(hits)))
  # every metric in the selection table is reported at the peak SNP
  sel_metrics <- setdiff(names(ex$get_selection()), c("chr", "pos", "group", "cum_pos"))
  expect_true(all(paste0("peak_", sel_metrics) %in% names(hits)))
  # every reported peak is at or above its group's neg_log10_p threshold
  thr <- ex$get_thresholds()
  for (i in seq_len(nrow(hits))) {
    tv <- thr$threshold[thr$group == hits$group[i]]
    expect_gte(hits$peak_neg_log10_p[i], tv)
  }
  # widening the window can only add (or keep) hits; strict-inside is a subset
  expect_gte(nrow(pos_selection_genes(ex, within = 5000)),
             nrow(pos_selection_genes(ex, within = 0)))
})

test_that("pos_selection_genes accepts a gene-track override", {
  ex <- example_ibd_results()
  # a full-genome scan without attaching PF3D7_GENES to the object
  all_hits <- pos_selection_genes(ex, genes = PF3D7_GENES)
  expect_gte(nrow(all_hits), nrow(pos_selection_genes(ex)))    # more genes -> >= hits
  expect_true(all(c("gene_id", "peak_neg_log10_p") %in% names(all_hits)))
  expect_equal(nrow(ex$get_genes()), 8)                        # object's track untouched
  # an object with NO track still works when genes= is supplied
  bare <- ibd_results(selection = ex$get_selection(), threshold = ex$get_thresholds())
  expect_s3_class(pos_selection_genes(bare, genes = PF_EXAMPLE_DRUG_GENES), "data.frame")
  expect_error(pos_selection_genes(bare), "no gene track")
})

test_that("bundled gene datasets are well-formed and usable as a track", {
  expect_true(all(c("Pf3D7_chrom", "start", "end", "chrom", "gene_id", "name")
                  %in% names(PF3D7_GENES)))
  expect_gt(nrow(PF3D7_GENES), 5000)
  expect_setequal(PF_EXAMPLE_DRUG_GENES$name,
    c("pfcrt", "pfdhfr", "pfmdr1", "pfdhps", "pfkelch13", "pfaat1", "pfgch1", "pfpx1"))
  # the `chrom` column is accepted as a `chr` alias by the constructor
  ibd <- ibd_results(genes = PF_EXAMPLE_DRUG_GENES)
  expect_setequal(ibd$get_genes()$chr, c("7", "4", "5", "8", "13", "6", "12"))
})

test_that("triangle label colour follows fill luminance, not a value threshold", {
  skip_if_no_ggplot()
  ramp <- plasgenomicsutilsR:::.IBD_FILL_DEFAULT
  # over a 0..0.5 scale: light tiles (incl. values above a near-zero median) read dark,
  # only genuinely dark fills read white
  tc <- plasgenomicsutilsR:::.readable_text_colour(c(0, 0.05, 0.2, 0.5), ramp,
                                                   limits = c(0, 0.5), trans = "identity")
  expect_equal(tc[1:3], rep("grey15", 3))
  expect_equal(tc[4], "white")
  # a labelled triangle still builds cleanly with the identity colour scale
  ex <- example_ibd_results()
  pw <- ex$get_pairwise_group()
  snpids <- unique(paste0(pw$chr, ":", pw$pos))[1:2]
  expect_silent(ggplot2::ggplotGrob(plot_pairwise_ibd_for_genes(ex, snps = snpids)))
})

# ---- block-based gene overlap + network -----------------------------------

make_block_obj <- function() {
  genes <- data.frame(name = c("pfcrt", "pfdhps"), chr = c("7", "8"),
                      start = c(403000, 548000), end = c(406000, 551000))
  blocks <- data.frame(
    sample1 = c("s1", "s1", "s2", "s3", "s5"),
    sample2 = c("s2", "s3", "s4", "s4", "s6"),
    chr = "Pf3D7_07_v3",
    start = c(402500, 404000, 404200, 900000, 10),
    end   = c(405000, 405500, 405800, 950000, 20),
    different = c(0, 0, 0, 0, 1))            # last row non-IBD (keeps s5/s6 "analyzed")
  meta <- data.frame(sample = paste0("s", 1:6), region = c("A", "A", "B", "B", "A", "B"))
  ibd_results(genes = genes, blocks = blocks, meta = meta,
              min_block_snp = 0, min_block_kb = 0, reference = "pf3d7")
}

test_that("blocks are stored IBD-only but analyzed samples come from all rows", {
  ibd <- make_block_obj()
  expect_equal(nrow(ibd$get_blocks()), 4)                 # different==1 row dropped
  expect_setequal(ibd$get_analyzed_samples(), paste0("s", 1:6))   # s5/s6 still counted
})

test_that("gene_ibd_overlap counts block overlap over all pairs", {
  ibd <- make_block_obj()
  ov <- gene_ibd_overlap(ibd, group = "region")
  crt <- ov[ov$gene == "pfcrt", ]
  # groups A={s1,s2,s5}, B={s3,s4,s6}: within A=C(3,2)=3, within B=3, between=9
  aa <- crt[crt$group_a == "A" & crt$group_b == "A", ]
  ab <- crt[crt$group_a == "A" & crt$group_b == "B", ]
  bb <- crt[crt$group_a == "B" & crt$group_b == "B", ]
  expect_equal(aa$n_pairs_ibd, 1); expect_equal(aa$n_pairs_total, 3)   # s1-s2
  expect_equal(ab$n_pairs_ibd, 2); expect_equal(ab$n_pairs_total, 9)   # s1-s3, s2-s4
  expect_equal(bb$n_pairs_ibd, 0); expect_equal(bb$n_pairs_total, 3)
  # a gene with no overlapping block -> all zero, denominators intact
  dhps <- ov[ov$gene == "pfdhps", ]
  expect_true(all(dhps$n_pairs_ibd == 0))
  expect_true(all(dhps$n_pairs_total > 0))
})

test_that("gene triangles use block overlap when blocks are present", {
  skip_if_no_ggplot()
  ibd <- make_block_obj()
  expect_silent(ggplot2::ggplotGrob(plot_pairwise_ibd_for_genes(ibd, genes = "pfcrt", group = "region")))
  # a precomputed overlap table plots directly
  tab <- data.frame(gene = "pfcrt", group_a = c("A", "A", "B"), group_b = c("A", "B", "B"),
                    frac_pairs_ibd = c(0.3, 0.5, 0.1))
  pre <- ibd_results(gene_overlap = tab)
  expect_silent(ggplot2::ggplotGrob(plot_pairwise_ibd_for_genes(pre)))
  expect_error(plot_pairwise_ibd_for_genes(pre, genes = "nope"), "not in the gene_overlap")
})

test_that("limits = 'shared' puts every feature on one fill scale", {
  skip_if_no_ggplot()
  # two genes with very different ranges: on their own scales each would span the ramp
  tab <- data.frame(
    gene = rep(c("hi", "lo"), each = 3),
    group_a = rep(c("A", "A", "B"), 2), group_b = rep(c("A", "B", "B"), 2),
    frac_pairs_ibd = c(0.2, 0.6, 0.9, 0.01, 0.02, 0.05))
  obj <- ibd_results(gene_overlap = tab)
  fill_lims <- function(p) ggplot2::ggplot_build(p)$plot$scales$get_scales("fill")$limits

  ind <- plot_pairwise_ibd_for_genes(obj, individual = TRUE, limits = "shared")
  expect_named(ind, c("hi", "lo"))
  expect_equal(fill_lims(ind[["hi"]]), c(0.01, 0.9))          # global range, both pages
  expect_equal(fill_lims(ind[["lo"]]), c(0.01, 0.9))
  # default keeps each page on its own scale
  free <- plot_pairwise_ibd_for_genes(obj, individual = TRUE)
  expect_null(fill_lims(free[["hi"]]))
  expect_null(fill_lims(free[["lo"]]))
  # on its own scale each gene runs to the top of the ramp regardless of how small its
  # values are; shared, only the gene holding the global maximum gets there
  tile_fills <- function(p) toupper(ggplot2::ggplot_build(p)$data[[1]]$fill)
  top <- toupper(utils::tail(.IBD_FILL_DEFAULT, 1))
  expect_true(top %in% tile_fills(free[["hi"]]))
  expect_true(top %in% tile_fills(free[["lo"]]))
  expect_true(top %in% tile_fills(ind[["hi"]]))
  expect_false(top %in% tile_fills(ind[["lo"]]))
  # explicit numeric limits still win, and the grid form accepts "shared" too
  expect_equal(fill_lims(plot_pairwise_ibd_for_genes(obj, limits = c(0, 1))), c(0, 1))
  expect_equal(fill_lims(plot_pairwise_ibd_for_genes(obj, limits = "shared")), c(0.01, 0.9))
})

test_that("plot_ibd_network builds, honours include_isolated, and parses a locus", {
  skip_if_no_ggplot()
  testthat::skip_if_not_installed("igraph")
  testthat::skip_if_not_installed("ggraph")
  ibd <- make_block_obj()
  connected <- plot_ibd_network(ibd, gene = "pfcrt", color_group = "region")
  expect_s3_class(connected, "ggplot")
  expect_silent(ggplot2::ggplotGrob(connected))
  full <- plot_ibd_network(ibd, gene = "pfcrt", include_isolated = TRUE)
  expect_silent(ggplot2::ggplotGrob(full))
  # locus form ("chr:start-end") works too
  expect_s3_class(plot_ibd_network(ibd, locus = "Pf3D7_07_v3:403000-406000"), "ggplot")
  # no blocks -> informative error
  expect_error(plot_ibd_network(ibd_results(), gene = "x"), "no IBD blocks")
})

test_that("gene overlap / triangles accept a gene-table override (not just track names)", {
  skip_if_no_ggplot()
  ibd <- make_block_obj()
  # a table of genes NOT in the object's track (chrom alias accepted)
  gtab <- data.frame(name = c("geneA", "geneB"), chrom = c("7", "8"),
                     start = c(403000, 548000), end = c(406000, 551000))
  ov <- gene_ibd_overlap(ibd, genes = gtab, group = "region")
  expect_setequal(levels(ov$gene), c("geneA", "geneB"))
  expect_silent(ggplot2::ggplotGrob(plot_pairwise_ibd_for_genes(ibd, genes = gtab, group = "region")))
})

test_that("plot_ibd_network title and subtitle can be toggled off", {
  skip_if_no_ggplot()
  testthat::skip_if_not_installed("igraph")
  testthat::skip_if_not_installed("ggraph")
  ibd <- make_block_obj()
  lab <- function(p) ggplot2::ggplot_build(p)$plot$labels
  d <- lab(plot_ibd_network(ibd, gene = "pfcrt"))                          # defaults
  expect_match(d$title, "IBD network"); expect_true(!is.null(d$subtitle))
  n <- lab(plot_ibd_network(ibd, gene = "pfcrt", title = NA, subtitle = FALSE))
  expect_null(n$title); expect_null(n$subtitle)                            # both off
  cu <- lab(plot_ibd_network(ibd, gene = "pfcrt", title = "mine"))
  expect_equal(cu$title, "mine")                                          # custom
  expect_s3_class(plot_ibd_network(ibd, gene = "pfcrt", layout = "drl"), "ggplot")  # dense layout
})

test_that("group heatmap shows chromosome bands + full genome extent, log fill is quiet", {
  skip_if_no_ggplot()
  ex <- example_ibd_results()
  # log2 fill no longer warns "introduced infinite values" (0s map to na.value)
  expect_silent(ggplot2::ggplotGrob(plot_ibd_pairwise_group_heatmap(ex, trans = "log2")))
  p <- plot_ibd_pairwise_group_heatmap(ex)
  b <- ggplot2::ggplot_build(p)
  # the chromosome-band rect layers are present, in this plot's own two-grey palette
  band <- plasgenomicsutilsR:::.CHR_BAND_TILE
  has_bands <- vapply(band, function(col) any(vapply(b$data, function(d)
    all(c("xmin", "xmax") %in% names(d)) && any(d$fill == col, na.rm = TRUE), logical(1))),
    logical(1))
  expect_true(all(has_bands))
  # the x-axis spans the full chromosome lengths (empty regions shown), not just SNP positions
  xr <- b$layout$panel_params[[1]]$x.range
  expect_gte(xr[2], max(ex$chrom_layout()$xmax) * 0.99)
})

test_that("log-transformed triangle fill is quiet and keeps 0.00 labels", {
  skip_if_no_ggplot()
  ex <- example_ibd_results()
  pw <- ex$get_pairwise_group()
  ids <- unique(paste0(pw$chr, ":", pw$pos))[1:2]
  expect_silent(ggplot2::ggplotGrob(plot_pairwise_ibd_for_genes(ex, snps = ids, trans = "log2")))
})

test_that("plot_ibd_network draws edges and an isolated grid", {
  skip_if_no_ggplot()
  testthat::skip_if_not_installed("igraph")
  testthat::skip_if_not_installed("ggraph")
  set.seed(4); samps <- paste0("s", 1:40)
  pr <- t(replicate(120, sample(samps[1:25], 2)))    # s1..s25 connected, s26..s40 isolated
  blocks <- rbind(
    data.frame(sample1 = pr[, 1], sample2 = pr[, 2], chr = "Pf3D7_07_v3",
               start = 889000, end = 900000, different = 0),
    data.frame(sample1 = "s26", sample2 = "s27", chr = "Pf3D7_01_v3",
               start = 1, end = 2, different = 1))   # keeps s26..s40 analyzed, none IBD
  ibd <- ibd_results(blocks = blocks, meta = data.frame(sample = samps, region = "A"),
                     min_block_snp = 0, min_block_kb = 0, 
                     genes = data.frame(name = "g", chr = "7", start = 889883, end = 899213),
                     reference = "pf3d7")
  seg_layers <- function(p) sum(vapply(ggplot2::ggplot_build(p)$data,
    function(d) all(c("x", "xend") %in% names(d)) && any(d$x != d$xend, na.rm = TRUE),
    logical(1)))
  expect_gt(seg_layers(plot_ibd_network(ibd, gene = "g")), 0)          # edges always drawn
  # include_isolated adds the isolated samples + the "unconnected" annotation
  p_iso <- plot_ibd_network(ibd, gene = "g", include_isolated = TRUE)
  labs <- unlist(lapply(ggplot2::ggplot_build(p_iso)$data,
                        function(d) if ("label" %in% names(d)) d$label))
  expect_true(any(grepl("unconnected", labs)))
})

test_that("spread down-weights edges inside dense neighbourhoods", {
  testthat::skip_if_not_installed("igraph")
  # a 12-node clique joined by one bridge to a 3-node path: clique edges are fully
  # redundant (Jaccard ~ 1), the bridge and path edges are not
  g <- igraph::disjoint_union(igraph::make_full_graph(12),
                              igraph::make_ring(3, circular = FALSE))
  g <- igraph::add_edges(g, c(1, 13))
  w <- .clique_relax_weights(g, spread = 1)
  el <- igraph::as_edgelist(g, names = FALSE)
  in_clique <- el[, 1] <= 12 & el[, 2] <= 12
  expect_lt(max(w[in_clique]), min(w[!in_clique]))            # clique edges pull least
  expect_true(all(w > 0 & w <= 1))
  expect_null(.clique_relax_weights(g, spread = 0))           # off -> unweighted
  # weaker spread moves the clique weights back toward 1
  expect_true(all(.clique_relax_weights(g, 0.5)[in_clique] > w[in_clique]))
})

test_that("spread opens up a dense cluster relative to sparser structure", {
  skip_if_no_ggplot()
  testthat::skip_if_not_installed("igraph")
  testthat::skip_if_not_installed("ggraph")
  # a fully IBD group of 14 alongside 20 unrelated IBD pairs; the group's edges are
  # redundant, the pairs' are not, so only the group should open up
  dense <- t(utils::combn(paste0("d", 1:14), 2))
  pairs <- cbind(paste0("p", 1:20, "a"), paste0("p", 1:20, "b"))
  blocks <- data.frame(sample1 = c(dense[, 1], pairs[, 1]),   # group edges first
                       sample2 = c(dense[, 2], pairs[, 2]),
                       chr = "Pf3D7_07_v3", start = 889000, end = 900000, different = 0)
  ibd <- ibd_results(blocks = blocks, min_block_snp = 0, min_block_kb = 0, 
                     genes = data.frame(name = "g", chr = "7", start = 889883, end = 899213),
                     reference = "pf3d7")
  # nodes are laid out in first-seen order, so the first 14 points are the dense group
  group_sep <- function(p) {
    d <- ggplot2::ggplot_build(p)$data
    pts <- d[[which(vapply(d, function(x) "shape" %in% names(x), logical(1)))[1]]]
    m <- as.matrix(stats::dist(pts[1:14, c("x", "y")])); diag(m) <- Inf
    stats::median(apply(m, 1, min))
  }
  expect_gt(group_sep(plot_ibd_network(ibd, gene = "g", spread = 1)),
            group_sep(plot_ibd_network(ibd, gene = "g", spread = 0)))
})

test_that("group_col_in_meta and set_group_order drive every table's group order", {
  want <- c("Delta", "Charlie", "Alpha")          # deliberately not alphabetical
  meta <- data.frame(sample = paste0("s", 1:6),
                     region = factor(rep(want, each = 2), levels = want),
                     stringsAsFactors = FALSE)
  per_snp <- data.frame(chr = "Pf3D7_01_v3", pos = c(10, 20, 30),
                        group = c("Alpha", "Charlie", "Delta"), frac_pairs_ibd = 0.1)
  pw <- data.frame(chr = "Pf3D7_01_v3", pos = 10,
                   group_a = c("Alpha", "Alpha", "Charlie"),
                   group_b = c("Charlie", "Delta", "Delta"), frac_pairs_ibd = 0.2)

  obj <- ibd_results(per_snp_group = per_snp, pairwise_group = pw, meta = meta,
                     group_col_in_meta = "region", reference = "pf3d7")
  expect_equal(obj$get_group_col(), "region")
  expect_equal(obj$get_group_order(), want)
  expect_equal(levels(obj$get_per_snp_group()$group), want)
  expect_equal(levels(obj$get_pairwise_group()$group_a), want)
  expect_equal(levels(obj$get_pairwise_group()$group_b), want)

  # set_group_order re-stamps everything
  obj$set_group_order(c("Alpha", "Delta", "Charlie"))
  expect_equal(levels(obj$get_per_snp_group()$group), c("Alpha", "Delta", "Charlie"))

  # a group present in the results but absent from the order would become NA -> error
  expect_error(obj$set_group_order(c("Alpha", "Delta")), "missing from")
  # a level no result uses is only a warning
  expect_warning(obj$set_group_order(c(want, "Atlantis")), "no IBD result uses")
  expect_error(obj$set_group_order(c("Alpha", "Alpha", "Charlie", "Delta")), "duplicated")

  # a non-factor meta column is accepted and natural-sorted
  meta2 <- meta; meta2$region <- as.character(meta2$region)
  obj2 <- ibd_results(per_snp_group = per_snp, meta = meta2,
                      group_col_in_meta = "region", reference = "pf3d7")
  expect_equal(obj2$get_group_order(), c("Alpha", "Charlie", "Delta"))
  expect_error(ibd_results(per_snp_group = per_snp, meta = meta,
                           group_col_in_meta = "nope", reference = "pf3d7"),
               "not a column of meta")
})

test_that("a non-alphabetical group order keeps the triangle on one side of the diagonal", {
  skip_if_no_ggplot()
  want <- c("Delta", "Charlie", "Alpha")     # string order disagrees with level order
  # pairs canonicalised the way the block/python paths do it: alphabetically
  tab <- data.frame(
    gene = "g",
    group_a = c("Alpha", "Alpha", "Alpha", "Charlie", "Charlie", "Delta"),
    group_b = c("Alpha", "Charlie", "Delta", "Charlie", "Delta", "Delta"),
    frac_pairs_ibd = c(0.5, 0.1, 0.2, 0.6, 0.3, 0.7), stringsAsFactors = FALSE)
  obj <- ibd_results(gene_overlap = tab, reference = "pf3d7")
  obj$set_group_order(want)
  p <- plot_pairwise_ibd_for_genes(obj)
  d <- ggplot2::ggplot_build(p)$data[[1]]
  # y is drawn on reversed levels, so every tile must sit at or below the diagonal
  expect_true(all(d$x + d$y <= length(want) + 1))
  expect_equal(nrow(d), 6L)                  # no pair lost or duplicated by the reorient
})

test_that("per-group thresholds do not reorder the facets", {
  skip_if_no_ggplot()
  want <- c("Delta", "Charlie", "Alpha")
  sel <- data.frame(chr = "Pf3D7_01_v3", pos = c(10, 20, 30),
                    group = want, neg_log10_p = c(5, 6, 7), stringsAsFactors = FALSE)
  thr <- data.frame(group = want, threshold = c(4, 4, 4), stringsAsFactors = FALSE)
  obj <- ibd_results(selection = sel, threshold = thr,
                     genes = data.frame(name = "g", chr = "1", start = 5, end = 40),
                     reference = "pf3d7")
  obj$set_group_order(want)
  # the threshold hline layer carries the facet column too; as a character vector it used
  # to make ggplot merge the layers' facet values alphabetically
  expect_true(is.factor(obj$get_thresholds()$group))
  b <- ggplot2::ggplot_build(
    plot_selection_manhattan(obj, draw_threshold = TRUE, label_genes = TRUE))
  expect_equal(as.character(b$layout$layout$group), want)
})

test_that("plot_ibd_network maps colour and shape independently", {
  skip_if_no_ggplot()
  testthat::skip_if_not_installed("igraph")
  testthat::skip_if_not_installed("ggraph")
  set.seed(9); samps <- paste0("s", 1:12)
  pr <- t(replicate(30, sample(samps, 2)))
  blocks <- data.frame(sample1 = pr[, 1], sample2 = pr[, 2], chr = "Pf3D7_07_v3",
                       start = 889000, end = 900000, different = 0)
  meta <- data.frame(sample = samps,
                     region = factor(rep(c("Delta", "Alpha"), each = 6),
                                     levels = c("Delta", "Alpha")),
                     marker = rep(c("wt", "mut", "mix"), 4), stringsAsFactors = FALSE)
  ibd <- ibd_results(blocks = blocks, meta = meta, min_block_snp = 0, min_block_kb = 0, 
                     genes = data.frame(name = "g", chr = "7", start = 889883, end = 899213),
                     reference = "pf3d7")
  pts <- function(p) {
    d <- ggplot2::ggplot_build(p)$data
    d[[which(vapply(d, function(z) "shape" %in% names(z), logical(1)))[1]]]
  }
  # colour only
  p1 <- plot_ibd_network(ibd, gene = "g", color_group = "region")
  expect_equal(length(unique(pts(p1)$colour)), 2)
  expect_equal(length(unique(pts(p1)$shape)), 1)
  # shape only -- nodes keep a single colour
  p2 <- plot_ibd_network(ibd, gene = "g", shape_group = "marker")
  expect_equal(length(unique(pts(p2)$colour)), 1)
  expect_equal(length(unique(pts(p2)$shape)), 3)
  # both at once
  p3 <- plot_ibd_network(ibd, gene = "g", color_group = "region", shape_group = "marker")
  expect_equal(length(unique(pts(p3)$colour)), 2)
  expect_equal(length(unique(pts(p3)$shape)), 3)
  # neither
  expect_equal(length(unique(pts(plot_ibd_network(ibd, gene = "g"))$colour)), 1)
  expect_error(plot_ibd_network(ibd, gene = "g", color_group = "nope"), "no column")
})

test_that("colors/shapes accept named or positional values", {
  lv <- c("Delta", "Charlie", "Alpha")
  pal <- function(n) paste0("auto", seq_len(n))
  # named: maps by name in any order, partial mappings keep the automatic rest
  got <- .match_scale_values(c(Alpha = "red", Delta = "blue"), lv, "colors", pal)
  expect_equal(unname(got[c("Delta", "Charlie", "Alpha")]), c("blue", "auto2", "red"))
  expect_equal(names(got), lv)                       # order follows the levels
  # unnamed: positional against the level order
  expect_equal(unname(.match_scale_values(c("a", "b", "c"), lv, "colors", pal)),
               c("a", "b", "c"))
  # too few positional values is an error; a stray name is only a warning
  expect_error(.match_scale_values(c("a", "b"), lv, "colors", pal), "one per level")
  expect_warning(.match_scale_values(c(Nope = "x"), lv, "colors", pal), "not among the levels")
  # more shape levels than distinguishable shapes is an error, not a silent recycle
  expect_error(.shape_palette(length(.SHAPE_PALETTE) + 1), "distinguishable")
})

test_that("gene_ibd_pairs lists the IBD pairs over each gene with their coverage", {
  genes <- data.frame(name = c("g1", "g2"), chr = c("7", "7"),
                      start = c(1000, 5000), end = c(2000, 6000),
                      gene_id = c("PF3D7_g1", "PF3D7_g2"), stringsAsFactors = FALSE)
  blocks <- data.frame(
    sample1 = c("a", "c", "e", "z", "m"),
    sample2 = c("b", "d", "f", "a", "n"),
    chr = "Pf3D7_07_v3",
    # hmmibd end is inclusive, so subtract one from the half-open end we want
    start = c(500, 1500, 3000, 5499, 500),
    end   = c(2499, 2499, 3999, 6499, 2499) - 1,
    different = c(0, 0, 0, 0, 1), stringsAsFactors = FALSE)
  ibd <- ibd_results(blocks = blocks, genes = genes, min_block_snp = 0, min_block_kb = 0, reference = "pf3d7")

  out <- gene_ibd_pairs(ibd)
  # only pairs that are IBD (different == 0) and actually reach a gene
  expect_setequal(paste(out$sample1, out$sample2), c("a b", "c d", "a z"))

  g1 <- out[out$gene == "g1", ]
  a <- g1[g1$sample1 == "a" & g1$sample2 == "b", ]
  expect_equal(a$coverage, "complete")
  expect_equal(c(a$covered_start, a$covered_end), c(1000, 2000))   # the gene's own bounds
  expect_equal(a$percent_covered, 100)
  cd <- g1[g1$sample1 == "c", ]
  expect_equal(cd$coverage, "partial")
  expect_equal(c(cd$covered_start, cd$covered_end), c(1500, 2000))
  expect_equal(cd$percent_covered, 50)

  # pairs are order-normalised; several genes arrive in one table
  expect_true(all(out$sample1 < out$sample2))
  expect_setequal(as.character(unique(out$gene)), c("g1", "g2"))

  # `within` widens selection only -- coverage stays measured on the gene
  far <- ibd_results(genes = genes, min_block_snp = 0, min_block_kb = 0, reference = "pf3d7",
                     blocks = data.frame(sample1 = "a", sample2 = "b", chr = "Pf3D7_07_v3",
                                         start = 2100, end = 2199, different = 0))
  expect_equal(nrow(gene_ibd_pairs(far)), 0)
  w <- gene_ibd_pairs(far, within = 500)
  expect_equal(nrow(w), 1)
  expect_equal(w$covered_bp, 0)
  expect_true(is.na(w$covered_start))

  expect_error(gene_ibd_pairs(ibd_results(genes = genes, reference = "pf3d7")), "no IBD blocks")
})

test_that("within pads the SNP fallback but never an explicit snps= request", {
  skip_if_no_ggplot()
  # a gene with no SNP inside it, and SNPs 3 kb either side
  genes <- data.frame(name = "g", chr = "1", start = 10000, end = 11000,
                      stringsAsFactors = FALSE)
  pw <- data.frame(chr = "Pf3D7_01_v3", pos = c(7000, 14000),
                   group_a = "A", group_b = "B", frac_pairs_ibd = c(0.2, 0.4),
                   stringsAsFactors = FALSE)
  ibd <- ibd_results(pairwise_group = pw, genes = genes, reference = "pf3d7")

  # unpadded, the gene contains nothing -- the honest answer on a sparse panel
  expect_error(suppressWarnings(plot_pairwise_ibd_for_genes(ibd)), "no.*overlapping SNPs")
  # the warning names the padding it did use
  expect_warning(expect_error(plot_pairwise_ibd_for_genes(ibd, within = 500)),
                 "within 500 bp")
  # padding wide enough to reach the flanking SNPs makes the gene plottable
  p <- plot_pairwise_ibd_for_genes(ibd, within = 4000)
  expect_equal(levels(p$data$gene), "g")
  expect_equal(p$data$frac_pairs_ibd, mean(c(0.2, 0.4)))   # both flanking SNPs aggregated

  # an explicitly named SNP is never widened, however large `within` is
  one <- plot_pairwise_ibd_for_genes(ibd, snps = "Pf3D7_01_v3:7000", within = 1e6)
  expect_equal(nlevels(one$data$gene), 1L)
  expect_equal(one$data$frac_pairs_ibd, 0.2)               # only that SNP, not its neighbour
})

test_that("ibd_results filters short / SNP-poor IBD blocks by default", {
  # hmmibd end is inclusive, so `end` here is one less than the half-open end
  blocks <- data.frame(
    sample1 = c("a", "c", "e", "g"), sample2 = c("b", "d", "f", "h"),
    chr = "Pf3D7_07_v3",
    start = 0,
    end   = c(20000, 14999, 20000, 15000) - 1,
    Nsnp  = c(30, 30, 14, 15),
    different = 0, stringsAsFactors = FALSE)

  keep <- ibd_results(blocks = blocks, reference = "pf3d7")$get_blocks()
  expect_equal(keep$sample1, c("a", "g"))            # long+SNP-rich, and exactly at bounds

  # every compared pair is still known, so no denominator shrinks
  obj <- ibd_results(blocks = blocks, reference = "pf3d7")
  expect_setequal(obj$get_analyzed_samples(), c("a", "b", "c", "d", "e", "f", "g", "h"))

  # either criterion can be switched off
  expect_equal(ibd_results(blocks = blocks, min_block_kb = 0,
                           reference = "pf3d7")$get_blocks()$sample1, c("a", "c", "g"))
  expect_equal(ibd_results(blocks = blocks, min_block_snp = 0,
                           reference = "pf3d7")$get_blocks()$sample1, c("a", "e", "g"))
  expect_equal(nrow(ibd_results(blocks = blocks, min_block_snp = 0, min_block_kb = 0,
                                reference = "pf3d7")$get_blocks()), 4)

  # no Nsnp column: warn, and fall back to the length criterion alone
  expect_warning(
    got <- ibd_results(blocks = blocks[, setdiff(names(blocks), "Nsnp")],
                       reference = "pf3d7")$get_blocks(),
    "no 'Nsnp' column")
  expect_equal(got$sample1, c("a", "e", "g"))
})

test_that("plot_ibd_network sharing = 'complete' requires the whole interval", {
  skip_if_no_ggplot()
  testthat::skip_if_not_installed("igraph")
  testthat::skip_if_not_installed("ggraph")
  genes <- data.frame(name = "g", chr = "7", start = 1000, end = 2000,
                      stringsAsFactors = FALSE)
  # a,b span the gene; c,d cover only its back half; e,f sit outside it entirely
  blocks <- data.frame(
    sample1 = c("a", "c", "e"), sample2 = c("b", "d", "f"), chr = "Pf3D7_07_v3",
    start = c(500, 1500, 3000),
    end   = c(2500, 2500, 4000) - 1,     # hmmibd end is inclusive
    different = 0, stringsAsFactors = FALSE)
  ibd <- ibd_results(blocks = blocks, genes = genes,
                     min_block_snp = 0, min_block_kb = 0, reference = "pf3d7")
  n_edges <- function(p) {
    d <- ggplot2::ggplot_build(p)$data
    i <- which(vapply(d, function(z) all(c("x", "xend") %in% names(z)), logical(1)))[1]
    if (is.na(i)) 0L else nrow(d[[i]])
  }
  expect_equal(n_edges(plot_ibd_network(ibd, gene = "g")), 2)                    # a-b and c-d
  expect_equal(n_edges(plot_ibd_network(ibd, gene = "g", sharing = "complete")), 1)  # only a-b

  # the subtitle says which criterion produced the graph
  expect_match(plot_ibd_network(ibd, gene = "g")$labels$subtitle, "part of the interval")
  expect_match(plot_ibd_network(ibd, gene = "g", sharing = "complete")$labels$subtitle,
               "the whole interval")

  # `within` applies to both: padding makes "complete" ask for more, so a-b stops qualifying
  expect_equal(n_edges(plot_ibd_network(ibd, gene = "g", within = 1000)), 2)
  expect_error(plot_ibd_network(ibd, gene = "g", within = 1000, sharing = "complete"),
               "no samples share")
  expect_error(plot_ibd_network(ibd, gene = "g", sharing = "nope"), "'arg' should be one of")
})

test_that("both threshold plots accept the same draw_threshold values", {
  skip_if_not_installed("ggplot2")
  ibd <- example_ibd_results()
  # the FDR columns a current `ibd_selection_statistic` writes
  thr <- ibd$get_thresholds()
  sel <- ibd$get_selection()
  skip_if(is.null(thr) || is.null(sel))

  fns <- list(manhattan = plot_selection_manhattan, tugofwar = plot_ibd_tugofwar)
  for (nm in names(fns)) {
    for (w in list(TRUE, FALSE, "bonferroni", "none")) {
      p <- fns[[nm]](ibd, draw_threshold = w)
      expect_s3_class(p, "ggplot")
    }
    # a bad value is rejected the same way by both, rather than reaching `&&`
    expect_error(fns[[nm]](ibd, draw_threshold = "nope"), "should be one of|must be")
  }
})

test_that("an FDR line is drawn when the thresholds carry one", {
  skip_if_not_installed("ggplot2")
  sel <- data.frame(
    group = "a", chr = "Pf3D7_07_v3", pos = seq(1000, by = 1000, length.out = 50),
    neg_log10_p = seq(0.1, 8, length.out = 50), stringsAsFactors = FALSE)
  thr <- data.frame(group = "a", threshold = 6, neg_log10_p_fdr_threshold = 3)
  ibd <- ibd_results(selection = sel, threshold = thr, genes = PF_EXAMPLE_DRUG_GENES)

  n_lines <- function(p) sum(vapply(p$layers,
    function(l) inherits(l$geom, "GeomHline"), logical(1)))
  bonf <- plot_selection_manhattan(ibd, draw_threshold = "bonferroni")
  both <- plot_selection_manhattan(ibd, draw_threshold = "both")
  none <- plot_selection_manhattan(ibd, draw_threshold = FALSE)
  expect_gt(n_lines(both), n_lines(bonf))
  expect_gt(n_lines(bonf), n_lines(none))

  # an older run with no FDR column says so instead of drawing the wrong line
  old <- ibd_results(selection = sel, threshold = data.frame(group = "a", threshold = 6))
  expect_error(plot_selection_manhattan(old, draw_threshold = "fdr"), "no FDR threshold")
})

test_that("the permutation line is drawn, and 'all' draws every kind available", {
  skip_if_not_installed("ggplot2")
  sel <- data.frame(
    group = "a", chr = "Pf3D7_07_v3", pos = seq(1000, by = 1000, length.out = 50),
    neg_log10_p = seq(0.1, 20, length.out = 50), stringsAsFactors = FALSE)
  thr <- data.frame(group = "a", threshold = 6, neg_log10_p_fdr_threshold = 3,
                    neg_log10_p_perm_threshold = 15,
                    neg_log10_p_emp_fdr_threshold = 11)
  psg <- data.frame(group = "a", chr = sel$chr, pos = sel$pos,
                    frac_pairs_ibd = seq(0, 0.4, length.out = 50))
  ibd <- ibd_results(selection = sel, threshold = thr, per_snp_group = psg,
                     genes = PF_EXAMPLE_DRUG_GENES)

  # the tug-of-war also draws a plain hline at its mirror axis; only the threshold lines
  # carry a `threshold` column
  hlines <- function(p) Filter(
    function(l) inherits(l$geom, "GeomHline") && !is.null(l$data$threshold), p$layers)
  at <- function(p) unname(sort(vapply(hlines(p),
    function(l) l$data$threshold[1], numeric(1))))

  for (nm in c("plot_selection_manhattan", "plot_ibd_tugofwar")) {
    f <- get(nm)
    expect_equal(at(f(ibd, draw_threshold = "permutation")), 15)
    expect_equal(at(f(ibd, draw_threshold = "empirical")), 11)
    expect_equal(at(f(ibd, draw_threshold = "all")), c(3, 6, 11, 15))
    # each kind gets its own linetype, so "all" is readable
    lt <- vapply(hlines(f(ibd, draw_threshold = "all")),
                 function(l) as.character(l$aes_params$linetype), character(1))
    expect_equal(length(unique(lt)), 4L)
  }

  # the empirical-FDR line sits between the family-wise permutation line and plain BH,
  # which is what controlling FDR against a data-driven null should look like
  a <- at(plot_selection_manhattan(ibd, draw_threshold = "all"))
  expect_true(a[2] < a[3] && a[3] < a[4])   # fdr < bonferroni < empirical < permutation

  # a run without --permute says so rather than silently dropping the line
  no_perm <- ibd_results(selection = sel, per_snp_group = psg,
                         threshold = data.frame(group = "a", threshold = 6))
  expect_error(plot_selection_manhattan(no_perm, draw_threshold = "permutation"),
               "no permutation threshold")
  expect_error(plot_selection_manhattan(no_perm, draw_threshold = "empirical"),
               "no empirical threshold")
  # "all" on that run falls back to what is there
  expect_equal(at(plot_selection_manhattan(no_perm, draw_threshold = "all")), 6)
})

# --------------------------------------------------------------------------- #
#  Fill-scale breaks and labels                                                #
# --------------------------------------------------------------------------- #

test_that("fill breaks are the transform's round numbers, topped by the maximum", {
  lg <- plasgenomicsutilsR:::.fill_breaks(c(0.00054, 0.52), "log2")
  expect_equal(max(lg), 0.52)                      # the top of the bar is labelled
  # everything below it is a power of two, which is what log2 should read as
  pw <- lg[lg < 0.52]
  expect_true(all(abs(log2(pw) - round(log2(pw))) < 1e-8))

  lin <- plasgenomicsutilsR:::.fill_breaks(c(0, 0.52), "identity")
  expect_equal(max(lin), 0.52)
  expect_true(all(lin >= 0 & lin <= 0.52))
})

test_that("a round break too close to the maximum is dropped, not printed on top of it", {
  b <- plasgenomicsutilsR:::.fill_breaks(c(0, 0.52), "identity")
  expect_false(any(abs(b - 0.5) < 1e-8))           # 0.5 would collide with 0.52
  # ...but a maximum that *is* round keeps it, with no duplicate
  b2 <- plasgenomicsutilsR:::.fill_breaks(c(0, 0.5), "identity")
  expect_equal(max(b2), 0.5)
  expect_equal(anyDuplicated(b2), 0L)
})

test_that("a non-linear scale says so in its title", {
  expect_equal(plasgenomicsutilsR:::.fill_scale_name("pairs IBD", "identity"), "pairs IBD")
  expect_equal(plasgenomicsutilsR:::.fill_scale_name("pairs IBD", "log2"),
               "pairs IBD (log2)")
  expect_equal(plasgenomicsutilsR:::.fill_scale_name("pairs IBD", "sqrt"),
               "pairs IBD (sqrt)")
})

test_that("degenerate limits fall back to the default breaks", {
  expect_null(plasgenomicsutilsR:::.fill_breaks(NULL, "identity"))
  expect_null(plasgenomicsutilsR:::.fill_breaks(c(1, 1), "identity"))
  expect_null(plasgenomicsutilsR:::.fill_breaks(c(NA, 1), "identity"))
  expect_null(plasgenomicsutilsR:::.fill_breaks(c(0, 1), "log2"))   # log needs lo > 0
})

test_that("labels keep small values legible instead of rounding them to zero", {
  # the real case: a log2 scale whose bottom break is a few ten-thousandths
  lab <- plasgenomicsutilsR:::.fill_labels(c(0.00054, 0.0053, 0.053, 0.52))
  expect_equal(lab, c("0.00054", "0.0053", "0.053", "0.52"))
  expect_false(any(lab == "0"))                       # 2 decimal places would give "0"
  # no ten-decimal padding, and a true zero still reads as 0
  expect_equal(plasgenomicsutilsR:::.fill_labels(c(0, 0.25, 0.5)), c("0", "0.25", "0.5"))
  # precision is the least that keeps the breaks distinct: 3 digits separates these
  expect_equal(plasgenomicsutilsR:::.fill_labels(c(0.1234, 0.1236)), c("0.123", "0.124"))
  # ...and it grows when that is not enough
  expect_equal(plasgenomicsutilsR:::.fill_labels(c(0.12341, 0.12342)),
               c("0.12341", "0.12342"))
})

test_that("the drawn legend ends at the data maximum, shared or per-page", {
  skip_if_not_installed("ggplot2")
  # a precomputed overlap table, so the values under the legend are known exactly
  grid <- expand.grid(group_a = c("a", "b", "c"), group_b = c("a", "b", "c"),
                      gene = c("g1", "g2"), stringsAsFactors = FALSE)
  grid <- grid[grid$group_a <= grid$group_b, ]
  set.seed(1)
  grid$frac_pairs_ibd <- c(seq(0.001, 0.10, length.out = 6),      # g1 tops out at 0.10
                           seq(0.001, 0.52, length.out = 6))      # g2 tops out at 0.52
  ibd <- ibd_results(gene_overlap = grid)

  top <- function(p) {
    g <- ggplot2::ggplot_build(p)
    as.numeric(utils::tail(g$plot$scales$get_scales("fill")$get_labels(), 1))
  }
  faceted <- plot_pairwise_ibd_for_genes(ibd)
  expect_equal(top(faceted), 0.52, tolerance = 0.01)

  # per page, each scales to its own maximum...
  free <- plot_pairwise_ibd_for_genes(ibd, individual = TRUE)
  expect_equal(unname(sort(vapply(free, top, numeric(1)))), c(0.10, 0.52), tolerance = 0.01)

  # ...and limits = "shared" pins every page to the overall maximum, which is the whole
  # point: before this, a log scale's last label stopped below the strongest colour
  shared <- plot_pairwise_ibd_for_genes(ibd, individual = TRUE, limits = "shared")
  expect_equal(unname(vapply(shared, top, numeric(1))), c(0.52, 0.52), tolerance = 0.01)

  lg <- plot_pairwise_ibd_for_genes(ibd, individual = TRUE, trans = "log2",
                                    limits = "shared")
  expect_equal(unname(vapply(lg, top, numeric(1))), c(0.52, 0.52), tolerance = 0.01)
  # and the log scale announces itself
  title <- ggplot2::ggplot_build(lg[[1]])$plot$scales$get_scales("fill")$name
  expect_equal(title, "pairs IBD (log2)")
})


test_that("the heatmap's chromosome bands are both grey, so zero tiles stay visible", {
  skip_if_not_installed("ggplot2")
  ibd <- example_ibd_results()
  skip_if(is.null(ibd$get_pairwise_group()))

  rect_fills <- function(p) unique(unlist(lapply(
    Filter(function(x) inherits(x$geom, "GeomRect"), p$layers),
    function(l) l$aes_params$fill)))
  # the fill scale starts at white, so a zero tile over a white band would disappear:
  # both bands must be grey, and distinguishable from each other
  fills <- rect_fills(plot_ibd_pairwise_group_heatmap(ibd))
  expect_equal(length(fills), 2L)
  expect_false(any(tolower(fills) %in% c("white", "#ffffff")))
  expect_false(fills[1] == fills[2])
  rgb <- vapply(fills, function(f) mean(grDevices::col2rgb(f)), numeric(1))
  expect_true(all(rgb < 250))                       # both clearly off-white

  # the point/bar plots keep the plain grey-on-panel default -- nothing there is white
  expect_equal(length(rect_fills(plot_ibd_sharing_manhattan(ibd))), 1L)
})

test_that("a band layer given no second fill draws only the first band", {
  skip_if_not_installed("ggplot2")
  layout <- example_ibd_results()$chrom_layout()
  rows <- function(ls) sum(vapply(ls, function(l) nrow(l$data), numeric(1)))

  one <- plasgenomicsutilsR:::.chr_band_layer(layout)
  expect_equal(length(one), 1L)
  expect_equal(rows(one), sum(layout$band == "a"))

  # two fills cover every chromosome, and each layer carries a SCALAR fill -- a per-row
  # vector would break the moment the plot is faceted
  two <- plasgenomicsutilsR:::.chr_band_layer(layout, c("grey80", "grey95"))
  expect_equal(length(two), 2L)
  expect_equal(rows(two), nrow(layout))
  expect_true(all(vapply(two, function(l) length(l$aes_params$fill), numeric(1)) == 1))

  expect_null(plasgenomicsutilsR:::.chr_band_layer(layout, c(NA, NA)))
})
