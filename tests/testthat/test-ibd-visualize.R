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
  ibd_results(genes = genes, blocks = blocks, meta = meta, reference = "pf3d7")
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
  connected <- plot_ibd_network(ibd, gene = "pfcrt", group = "region")
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
  # a grey chromosome-band rect layer is present
  has_bands <- any(vapply(b$data, function(d)
    all(c("xmin", "xmax") %in% names(d)) && any(d$fill == "#ebebeb", na.rm = TRUE), logical(1)))
  expect_true(has_bands)
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
  ibd <- ibd_results(blocks = blocks,
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
