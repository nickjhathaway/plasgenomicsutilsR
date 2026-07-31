# Population-structure analysis: LD-pruned genotypes -> PCA / UMAP, and sNMF
# admixture, with plot_*() functions. Heavy analysis packages (SNPRelate, gdsfmt,
# uwot, LEA) are optional (Suggests), guarded at call time.

# ---- genotypes -------------------------------------------------------------

#' LD-prune a VCF and return the genotype matrix
#'
#' Converts a VCF to GDS (only when needed), LD-prunes, and returns the pruned
#' genotype matrix (samples x SNPs, coded 0/1/2, `NA` for missing) via
#' \pkg{SNPRelate}.
#'
#' @param vcf Path to a (bgzipped) VCF.
#' @param gds Optional GDS path; derived from `vcf` if `NULL`.
#' @param ld_threshold,slide_max_bp,slide_max_n,autosome_only Passed to
#'   [SNPRelate::snpgdsLDpruning()] (defaults 0.2 / 20000 / 200 / `FALSE`).
#' @param maf,missing_rate Optional MAF / per-SNP missing-rate cutoffs for pruning.
#' @param seed Random seed for the pruning.
#' @return A list with `genotype` (matrix), `sample.id`, and `snp.id`.
#' @export
run_ld_prune <- function(vcf, gds = NULL, ld_threshold = 0.2, slide_max_bp = 20000,
                         slide_max_n = 200, autosome_only = FALSE, maf = NaN,
                         missing_rate = NaN, seed = 42) {
  .need_package("SNPRelate", "run_ld_prune()")
  .need_package("gdsfmt", "run_ld_prune()")
  if (is.null(gds)) gds <- sub("\\.vcf(\\.gz)?$", ".gds", vcf, ignore.case = TRUE)
  if (identical(gds, vcf)) gds <- paste0(vcf, ".gds")
  if (!file.exists(gds) || file.mtime(gds) < file.mtime(vcf)) {
    SNPRelate::snpgdsVCF2GDS(vcf, gds, method = "biallelic.only", verbose = FALSE)
  }
  h <- SNPRelate::snpgdsOpen(gds)
  on.exit(SNPRelate::snpgdsClose(h), add = TRUE)
  set.seed(seed)
  snpset <- SNPRelate::snpgdsLDpruning(
    h, autosome.only = autosome_only, ld.threshold = ld_threshold,
    slide.max.bp = slide_max_bp, slide.max.n = slide_max_n,
    maf = maf, missing.rate = missing_rate, verbose = FALSE)
  snp_ids <- unlist(snpset, use.names = FALSE)
  geno <- SNPRelate::snpgdsGetGeno(h, snp.id = snp_ids, with.id = TRUE, verbose = FALSE)
  list(genotype = geno$genotype, sample.id = geno$sample.id, snp.id = geno$snp.id)
}

# mean-impute missing genotypes per SNP; drop all-missing columns
.impute_geno <- function(mat) {
  mat <- as.matrix(mat)
  cm <- colMeans(mat, na.rm = TRUE)
  keep <- !is.nan(cm)
  mat <- mat[, keep, drop = FALSE]
  cm <- cm[keep]
  na <- which(is.na(mat), arr.ind = TRUE)
  if (nrow(na)) mat[na] <- cm[na[, "col"]]
  mat
}

# ---- PCA + UMAP container --------------------------------------------------

#' Compute PCA and UMAP from a genotype matrix
#'
#' Mean-imputes missing genotypes, runs PCA ([stats::prcomp()]) and, optionally, a
#' UMAP embedding (\pkg{uwot}) with PCA initialisation. Returns a `pop_structure`
#' object the `plot_*()` functions read.
#'
#' @param geno A genotype matrix (samples x SNPs, 0/1/2, `NA` allowed) or the list
#'   returned by [run_ld_prune()].
#' @param samples Sample ids (defaults to the genotype list's `sample.id`, or matrix
#'   row names).
#' @param meta Optional per-sample metadata: a data frame with a `sample` column
#'   (plus e.g. `region`, `country`) used for colouring.
#' @param n_pcs Number of PCs to summarise (variance explained).
#' @param umap Compute a UMAP embedding.
#' @param umap_pca,n_neighbors,min_dist UMAP parameters (defaults 30 / 15 / 0.1);
#'   `n_neighbors` is clamped to the sample count.
#' @param seed Random seed for UMAP.
#' @return A `pop_structure` object (a list with `samples`, `pca` scores,
#'   `pca_var`, `umap`, `meta`).
#' @export
pop_structure <- function(geno, samples = NULL, meta = NULL, n_pcs = 50, umap = TRUE,
                          umap_pca = 30, n_neighbors = 15, min_dist = 0.1, seed = 42) {
  if (is.list(geno) && !is.null(geno$genotype)) {
    if (is.null(samples)) samples <- geno$sample.id
    mat <- geno$genotype
  } else {
    mat <- as.matrix(geno)
    if (is.null(samples)) samples <- rownames(mat)
  }
  if (is.null(samples)) samples <- as.character(seq_len(nrow(mat)))
  mat <- .impute_geno(mat)

  pca <- stats::prcomp(mat, center = TRUE, scale. = FALSE)
  ve <- pca$sdev^2 / sum(pca$sdev^2)
  npc <- min(n_pcs, length(ve))
  pca_var <- data.frame(
    PC = seq_len(npc),
    var_explained = round(ve[seq_len(npc)] * 100, 2),
    cumulative = round(cumsum(ve)[seq_len(npc)] * 100, 2))

  umap_df <- NULL
  if (umap) {
    .need_package("uwot", "the UMAP embedding in pop_structure()")
    nn <- max(2L, min(n_neighbors, nrow(mat) - 1L))
    # umap_pca may be a count (>= 1) or a variance fraction (0 < x < 1)
    npca <- if (umap_pca > 0 && umap_pca < 1) n_pcs_for_variance(pca, umap_pca) else as.integer(umap_pca)
    npca <- min(npca, ncol(mat), nrow(mat) - 1L)
    set.seed(seed)
    emb <- uwot::umap(mat, pca = npca, pca_center = TRUE, n_neighbors = nn,
                      min_dist = min_dist)
    umap_df <- data.frame(sample = samples, UMAP1 = emb[, 1], UMAP2 = emb[, 2],
                          stringsAsFactors = FALSE)
  }
  structure(list(samples = samples, pca = pca$x, prcomp = pca, pca_var = pca_var,
                 umap = umap_df, meta = meta), class = "pop_structure")
}

#' @export
print.pop_structure <- function(x, ...) {
  cat("<pop_structure>", length(x$samples), "samples,",
      ncol(x$pca), "PCs\n")
  cat("  PC1-2 variance:", x$pca_var$var_explained[1], "% /",
      x$pca_var$var_explained[2], "%\n")
  cat("  UMAP:", if (is.null(x$umap)) "-" else "yes",
      "  metadata:", if (is.null(x$meta)) "-" else paste(ncol(x$meta), "cols"), "\n")
  invisible(x)
}

.ps_meta_join <- function(df, meta) {
  if (is.null(meta) || !"sample" %in% names(meta)) return(df)
  merge(df, meta, by = "sample", all.x = TRUE, sort = FALSE)
}

# ---- plots -----------------------------------------------------------------

#' PCA scatter plot
#'
#' @param x A `pop_structure` object or a [PopStructure] R6 object.
#' @param pcs Which two PCs to plot (default `c(1, 2)`).
#' @param colour Metadata column to colour points by (needs `meta`).
#' @param colors Optional named `level -> colour` vector for the colour scale
#'   (e.g. from [meta_colors()]); a `PopStructure` supplies its stored map.
#' @param point_size,point_alpha Point aesthetics.
#' @param legend_point_size Size of the coloured dots in the legend (default `3`, larger
#'   than the plotted points so the key is easy to read); `NULL` leaves it as-is.
#' @return A ggplot object.
#' @export
plot_pca <- function(x, pcs = c(1, 2), colour = NULL, colors = NULL,
                     point_size = 1.6, point_alpha = 0.8, legend_point_size = 3) {
  .need_package("ggplot2", "plot_pca()")
  if (inherits(x, "PopStructure")) {
    if (is.null(colors) && !is.null(colour)) colors <- x$get_colors()[[colour]]
    x <- x$as_ps()
  }
  a <- pcs[1]; b <- pcs[2]
  df <- data.frame(sample = x$samples, .x = x$pca[, a], .y = x$pca[, b],
                   stringsAsFactors = FALSE)
  df <- .ps_meta_join(df, x$meta)
  ve <- x$pca_var$var_explained
  p <- ggplot2::ggplot(df, ggplot2::aes(.data$.x, .data$.y)) +
    ggplot2::labs(x = sprintf("PC%d (%.1f%%)", a, ve[a]),
                  y = sprintf("PC%d (%.1f%%)", b, ve[b])) +
    ggplot2::theme_minimal(base_size = 11)
  p + .scatter_points(colour, point_size, point_alpha, colors, legend_point_size)
}

#' UMAP scatter plot
#'
#' @param x A `pop_structure` object (built with `umap = TRUE`) or a [PopStructure].
#' @param colour Metadata column to colour points by (needs `meta`).
#' @param colors Optional named `level -> colour` vector for the colour scale
#'   (e.g. from [meta_colors()]); a `PopStructure` supplies its stored map.
#' @param point_size,point_alpha Point aesthetics.
#' @param legend_point_size Size of the coloured dots in the legend (default `3`, larger
#'   than the plotted points so the key is easy to read); `NULL` leaves it as-is.
#' @return A ggplot object.
#' @export
plot_umap <- function(x, colour = NULL, colors = NULL, point_size = 1.6,
                      point_alpha = 0.8, legend_point_size = 3) {
  .need_package("ggplot2", "plot_umap()")
  if (inherits(x, "PopStructure")) {
    if (is.null(colors) && !is.null(colour)) colors <- x$get_colors()[[colour]]
    x <- x$as_ps()
  }
  if (is.null(x$umap)) stop("this pop_structure has no UMAP (build with umap = TRUE)",
                            call. = FALSE)
  df <- .ps_meta_join(x$umap, x$meta)
  df$.x <- df$UMAP1; df$.y <- df$UMAP2
  ggplot2::ggplot(df, ggplot2::aes(.data$.x, .data$.y)) +
    ggplot2::labs(x = "UMAP 1", y = "UMAP 2") +
    ggplot2::coord_equal() +
    ggplot2::theme_minimal(base_size = 11) +
    .scatter_points(colour, point_size, point_alpha, colors, legend_point_size)
}

# shared point + colour layers for the scatter plots
.scatter_points <- function(colour, point_size, point_alpha, colors = NULL,
                            legend_point_size = NULL) {
  if (is.null(colour)) {
    return(list(ggplot2::geom_point(size = point_size, alpha = point_alpha)))
  }
  out <- list(
    ggplot2::geom_point(ggplot2::aes(colour = .data[[colour]]),
                        size = point_size, alpha = point_alpha),
    ggplot2::labs(colour = colour))
  if (!is.null(colors)) out <- c(out, list(ggplot2::scale_colour_manual(values = colors)))
  if (!is.null(legend_point_size))
    out <- c(out, list(ggplot2::guides(
      colour = ggplot2::guide_legend(override.aes = list(size = legend_point_size)))))
  out
}

# ---- sNMF admixture --------------------------------------------------------

#' Run sNMF (LEA) admixture on a genotype matrix
#'
#' Writes a `.geno` file (missing coded as 9) and runs [LEA::snmf()] over a range of
#' `K`, returning the sNMF project. Cross-entropy is computed so [snmf_best_k()] can
#' pick `K`.
#'
#' @param geno A genotype matrix (samples x SNPs, 0/1/2, `NA` allowed) or the list
#'   from [run_ld_prune()].
#' @param K Integer vector of ancestral-population counts to fit (default `1:10`).
#' @param rep Repetitions per `K`.
#' @param alpha Regularisation.
#' @param seed Random seed.
#' @param cpu CPU cores.
#' @param cache Reuse a previously computed project when the genotypes and all
#'   parameters match (sNMF is slow); set `FALSE` to always recompute.
#' @param cache_dir Directory for the cached `.geno`/project (default a per-session
#'   temp dir). Point it at a persistent path to reuse across sessions.
#' @param verbose Show LEA's (voluminous) console output. Default `FALSE` runs it
#'   quietly; set `log_file` to capture it instead of discarding it.
#' @param log_file Optional file to write LEA's output to when `verbose = FALSE`.
#' @return An `snmf_fit` object: a list with the LEA `project`, the fitted `K` range,
#'   and `samples`.
#' @export
run_snmf <- function(geno, K = 1:10, rep = 10, alpha = 10, seed = 42, cpu = 1,
                     cache = TRUE, cache_dir = NULL, verbose = FALSE, log_file = NULL) {
  .need_package("LEA", "run_snmf()")
  if (is.list(geno) && !is.null(geno$genotype)) {
    samples <- geno$sample.id
    mat <- geno$genotype
  } else {
    mat <- as.matrix(geno)
    samples <- rownames(mat)
  }
  K <- as.integer(K)
  if (is.null(cache_dir)) cache_dir <- file.path(tempdir(), "plas_snmf")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  key <- substr(rlang::hash(list(mat, K, rep, alpha, seed)), 1, 16)
  geno_file <- file.path(cache_dir, paste0("snmf_", key, ".geno"))
  proj_file <- sub("\\.geno$", ".snmfProject", geno_file)
  samp_file <- file.path(cache_dir, paste0("snmf_", key, ".samples.rds"))

  if (cache && file.exists(proj_file)) {
    project <- .run_quiet(function() LEA::load.snmfProject(proj_file), verbose, log_file)
    if (file.exists(samp_file)) samples <- readRDS(samp_file)
  } else {
    mat[is.na(mat)] <- 9L
    project <- .run_quiet(function() {
      LEA::write.geno(mat, output.file = geno_file)
      LEA::snmf(geno_file, K = K, repetitions = rep, alpha = alpha, entropy = TRUE,
                project = "new", seed = seed, CPU = cpu)
    }, verbose, log_file)
    if (!is.null(samples)) saveRDS(samples, samp_file)
  }
  structure(list(project = project, K = K, samples = samples, geno_file = geno_file),
            class = "snmf_fit")
}

# run a thunk quietly (LEA prints a lot): discard stdout, or send it to `log_file`.
.run_quiet <- function(fn, verbose = FALSE, log_file = NULL) {
  if (verbose) return(fn())
  con <- if (is.null(log_file)) nullfile() else log_file
  out <- NULL
  utils::capture.output(out <- fn(), file = con, type = "output")
  out
}

#' Number of PCs explaining a target cumulative variance
#'
#' Handy for setting how many principal components feed the UMAP: pass the fraction of
#' variance you want the PCA step to capture.
#'
#' @param x A [stats::prcomp()] result or a numeric vector of eigenvalues/variances.
#' @param target Cumulative variance to reach, as a proportion (`0.1`) or a percent
#'   (`10`).
#' @return The smallest number of leading PCs whose cumulative variance >= `target`.
#' @export
n_pcs_for_variance <- function(x, target = 0.8) {
  ve <- if (inherits(x, "prcomp")) x$sdev^2 / sum(x$sdev^2)
        else if (is.numeric(x)) x / sum(x)
        else stop("`x` must be a prcomp result or a numeric variance vector", call. = FALSE)
  if (target > 1) target <- target / 100
  w <- which(cumsum(ve) >= target)
  if (length(w)) w[1] else length(ve)
}

#' @export
print.snmf_fit <- function(x, ...) {
  cat("<snmf_fit>", length(x$samples), "samples,  K =",
      paste(range(x$K), collapse = "-"), "\n")
  invisible(x)
}

.snmf_project <- function(x) if (inherits(x, "snmf_fit")) x$project else x

#' Pick the best K from an sNMF fit by cross-entropy
#'
#' @param x An [run_snmf()] result (or a raw LEA project, with `K` given).
#' @param K Candidate K values (defaults to the fitted range from [run_snmf()]).
#' @param stat Combine replicates by `"mean"` (default) or `"min"` cross-entropy.
#' @return The K minimising the summarised cross-entropy.
#' @export
snmf_best_k <- function(x, K = NULL, stat = c("mean", "min")) {
  .need_package("LEA", "snmf_best_k()")
  stat <- match.arg(stat)
  project <- .snmf_project(x)
  if (is.null(K)) {
    if (inherits(x, "snmf_fit")) K <- x$K
    else stop("pass K (the fitted range) when giving a raw sNMF project", call. = FALSE)
  }
  fn <- if (stat == "min") min else mean
  vals <- vapply(K, function(k) fn(LEA::cross.entropy(project, K = k)), numeric(1))
  K[which.min(vals)]
}

#' Best-run Q (ancestry proportion) matrix for a given K
#'
#' @param x An [run_snmf()] result (or a raw LEA project).
#' @param K The number of ancestral populations.
#' @param run Which replicate; defaults to the lowest cross-entropy run.
#' @return A samples-by-K matrix of ancestry proportions, with sample ids as row names
#'   when available.
#' @export
snmf_q <- function(x, K, run = NULL) {
  .need_package("LEA", "snmf_q()")
  project <- .snmf_project(x)
  if (is.null(run)) run <- which.min(LEA::cross.entropy(project, K = K))
  q <- LEA::Q(project, K = K, run = run)
  colnames(q) <- paste0("K", seq_len(ncol(q)))
  if (inherits(x, "snmf_fit") && !is.null(x$samples)) rownames(q) <- x$samples
  q
}

#' Sample order for admixture bars
#'
#' Orders samples within each `group` by hierarchical clustering of their ancestry
#' vectors. Compute it once (e.g. at the best K) and pass it as `sample_order` to
#' [plot_admixture()] for every K, so bars stay in the same position across K.
#'
#' @param q A samples-by-K ancestry matrix.
#' @param samples Sample ids (default `rownames(q)`).
#' @param meta,group Optional metadata + the column to order within.
#' @return An ordered character vector of sample ids.
#' @export
admixture_order <- function(q, samples = NULL, meta = NULL, group = NULL) {
  q <- as.matrix(q)
  if (is.null(samples)) samples <- rownames(q)
  if (is.null(samples)) samples <- as.character(seq_len(nrow(q)))
  grp <- rep("all", length(samples))
  if (!is.null(meta) && !is.null(group) && "sample" %in% names(meta) && group %in% names(meta)) {
    grp <- as.character(meta[[group]][match(samples, meta$sample)])
  }
  ord <- character(0)
  for (g in unique(grp)) {
    idx <- which(grp == g)
    if (length(idx) > 2) {
      hc <- stats::hclust(stats::dist(q[idx, , drop = FALSE]), method = "ward.D2")
      idx <- idx[hc$order]
    }
    ord <- c(ord, samples[idx])
  }
  ord
}

#' Admixture (STRUCTURE) bar plot from a Q matrix
#'
#' Stacked ancestry-proportion bars, one per sample, optionally faceted by a grouping
#' column, with an optional group colour strip whose colours can match another plot
#' (e.g. UMAP points).
#'
#' @param q A samples-by-K ancestry matrix (e.g. from [snmf_q()]).
#' @param samples Sample ids (defaults to `rownames(q)`).
#' @param meta Optional metadata data frame with a `sample` column.
#' @param group Optional metadata column to facet by (e.g. `"region"`).
#' @param order_within Order samples within each group by clustering (ignored when
#'   `sample_order` is supplied).
#' @param sample_order Optional explicit sample order (see [admixture_order()]); keeps
#'   bars in the same place across different K.
#' @param colours Optional fill colours for the K clusters.
#' @param group_bar Draw a coloured strip above the bars keyed by `group`.
#' @param group_colours Named `level -> colour` vector for the group strip (see
#'   [meta_colors()]); pass the same mapping you use for the UMAP to match colours.
#' @param border Outline each sample's bar (default `TRUE`) so neighbours with nearly
#'   identical ancestry stay distinct; set `FALSE` for borderless bars.
#' @param border_colour,border_linewidth Colour and width of the per-sample outline.
#' @return A ggplot object.
#' @export
plot_admixture <- function(q, samples = NULL, meta = NULL, group = NULL,
                           order_within = TRUE, sample_order = NULL, colours = NULL,
                           group_bar = FALSE, group_colours = NULL,
                           border = TRUE, border_colour = "black",
                           border_linewidth = 0.15) {
  .need_package("ggplot2", "plot_admixture()")
  q <- as.matrix(q)
  if (is.null(colnames(q))) colnames(q) <- paste0("K", seq_len(ncol(q)))
  if (is.null(samples)) samples <- rownames(q)
  if (is.null(samples)) samples <- as.character(seq_len(nrow(q)))
  df <- data.frame(sample = samples, q, check.names = FALSE, stringsAsFactors = FALSE)
  df <- .ps_meta_join(df, meta)

  if (is.null(sample_order)) {
    sample_order <- if (order_within) admixture_order(q, samples, meta, group) else samples
  }
  df <- df[match(sample_order, df$sample), , drop = FALSE]
  df$sample <- factor(df$sample, levels = sample_order)

  long <- stats::reshape(
    df[, c("sample", if (!is.null(group)) group, colnames(q))],
    direction = "long", varying = colnames(q), v.names = "q",
    times = colnames(q), timevar = "cluster", idvar = "sample")
  long$sample <- factor(long$sample, levels = sample_order)

  p <- ggplot2::ggplot(long, ggplot2::aes(.data$sample, .data$q, fill = .data$cluster)) +
    ggplot2::geom_col(width = 1, position = "stack",
                      colour = if (border) border_colour else NA,
                      linewidth = border_linewidth) +
    ggplot2::scale_fill_manual(values = if (is.null(colours)) .pick_palette(ncol(q)) else colours) +
    ggplot2::labs(x = NULL, y = "ancestry", fill = "cluster") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank())

  if (group_bar && !is.null(group) && group %in% names(df)) {
    .need_package("ggnewscale", "the group colour bar")
    gdf <- data.frame(sample = factor(df$sample, levels = sample_order),
                      grp = as.character(df[[group]]), stringsAsFactors = FALSE)
    if (is.null(group_colours)) group_colours <- meta_colors(df, cols = group)[[group]]
    p <- p + ggnewscale::new_scale_fill() +
      ggplot2::geom_tile(data = gdf, ggplot2::aes(.data$sample, y = 1.06, fill = .data$grp),
                         height = 0.05, inherit.aes = FALSE) +
      ggplot2::scale_fill_manual(values = group_colours, name = group) +
      ggplot2::coord_cartesian(clip = "off")
  }
  if (!is.null(group) && group %in% names(df)) {
    p <- p + ggplot2::facet_grid(stats::as.formula(paste("~", group)),
                                 scales = "free_x", space = "free")
  }
  p
}

# ---- PopStructure R6 wrapper ----------------------------------------------

#' Population-structure workspace (PCA + UMAP + admixture)
#'
#' An R6 object that bundles a genotype matrix, its PCA (the full [stats::prcomp()]
#' result), an optional UMAP embedding, per-sample metadata, a shared metadata colour
#' map, and an sNMF admixture fit. Because it keeps the fitted objects, you can save it
#' (`saveRDS`) and re-plot without recomputing, colour PCA/UMAP/admixture consistently
#' (same `level -> colour` map), reorder admixture bars once and reuse the order across
#' K, and sub-select to a set of samples (or a metadata match) for output.
#'
#' @examples
#' ps <- example_pop_structure()
#' ps$run_umap(pca_components = 0.5)      # PCs covering 50% of variance
#' \dontrun{
#' ps$run_snmf(K = 1:6)                   # cached + quiet
#' ps$plot_admixture(K = ps$best_k(), group = "population")
#' west <- ps$subset(population = c("PopA", "PopB"))
#' }
#' @export
PopStructure <- R6::R6Class("PopStructure",
  public = list(
    #' @description Build from a genotype matrix (or a [run_ld_prune()] list).
    #' @param geno Genotype matrix (samples x SNPs, 0/1/2, `NA`) or `run_ld_prune()` list.
    #' @param samples Sample ids (default row names / the list's `sample.id`).
    #' @param meta Optional metadata (data frame with a `sample` column).
    #' @param n_pcs Number of PCs to summarise.
    #' @param colors Optional named list of colour maps (see [meta_colors()]).
    initialize = function(geno, samples = NULL, meta = NULL, n_pcs = 50, colors = NULL) {
      if (is.list(geno) && !is.null(geno$genotype)) {
        if (is.null(samples)) samples <- geno$sample.id
        geno <- geno$genotype
      }
      mat <- as.matrix(geno)
      if (is.null(samples)) samples <- rownames(mat)
      if (is.null(samples)) samples <- as.character(seq_len(nrow(mat)))
      rownames(mat) <- samples
      private$geno_mat <- mat
      private$sample_ids <- samples
      private$active_ids <- samples
      private$n_pcs <- n_pcs
      private$ps <- pop_structure(mat, samples = samples, n_pcs = n_pcs, umap = FALSE)
      private$colors <- if (is.null(colors)) list() else colors
      if (!is.null(meta)) self$add_meta(meta)
      invisible(self)
    },

    #' @description Attach/replace metadata; auto-assigns colours for new columns.
    #' @param meta Data frame with a `sample` column.
    #' @param colors Optional colour overrides (`column -> (level -> colour)`).
    add_meta = function(meta, colors = NULL) {
      meta <- as.data.frame(meta, stringsAsFactors = FALSE)
      if (!"sample" %in% names(meta)) stop("`meta` needs a `sample` column", call. = FALSE)
      private$meta_df <- meta
      private$ps$meta <- meta
      auto <- meta_colors(meta)
      for (nm in names(auto)) if (is.null(private$colors[[nm]])) private$colors[[nm]] <- auto[[nm]]
      if (!is.null(colors)) self$set_colors(colors)
      invisible(self)
    },

    #' @description Set/override colour maps for metadata columns.
    #' @param colors Named list `column -> (level -> colour)`.
    set_colors = function(colors) {
      for (nm in names(colors)) private$colors[[nm]] <- colors[[nm]]
      invisible(self)
    },

    #' @description Fix the level order of a metadata column. Because every plot reads
    #'   the shared metadata and colour map, this one order flows through the legends
    #'   (PCA / UMAP), the admixture facet order, and the colour strips. Existing colours
    #'   follow their level (only the order changes); unlisted levels are dropped.
    #' @param column Metadata column name.
    #' @param levels The desired level order.
    set_levels = function(column, levels) {
      m <- private$meta_df
      if (is.null(m) || !column %in% names(m))
        stop("no metadata column '", column, "'", call. = FALSE)
      m[[column]] <- factor(as.character(m[[column]]), levels = levels)
      private$meta_df <- m
      private$ps$meta <- m
      cur <- private$colors[[column]]
      private$colors[[column]] <-
        if (!is.null(cur) && all(levels %in% names(cur))) cur[levels]
        else meta_colors(m, cols = column)[[column]]
      invisible(self)
    },

    #' @description Compute a UMAP embedding.
    #' @param pca_components PCs feeding UMAP: a count (`>= 1`) or a variance fraction
    #'   (`0 < x < 1`, e.g. `0.1` uses enough PCs for 10% of variance, via
    #'   [n_pcs_for_variance()]).
    #' @param n_neighbors,min_dist UMAP parameters.
    #' @param seed Random seed.
    run_umap = function(pca_components = 30, n_neighbors = 15, min_dist = 0.1, seed = 42) {
      private$ps <- pop_structure(private$geno_mat, samples = private$sample_ids,
                                  meta = private$meta_df, n_pcs = private$n_pcs,
                                  umap = TRUE, umap_pca = pca_components,
                                  n_neighbors = n_neighbors, min_dist = min_dist, seed = seed)
      invisible(self)
    },

    #' @description Fit sNMF admixture (cached and quiet by default).
    #' @param K,rep,alpha,seed,cpu,cache,cache_dir,verbose,log_file Passed to [run_snmf()].
    run_snmf = function(K = 1:10, rep = 10, alpha = 10, seed = 42, cpu = 1,
                        cache = TRUE, cache_dir = NULL, verbose = FALSE, log_file = NULL) {
      private$snmf_fit <- run_snmf(
        list(genotype = private$geno_mat, sample.id = private$sample_ids),
        K = K, rep = rep, alpha = alpha, seed = seed, cpu = cpu,
        cache = cache, cache_dir = cache_dir, verbose = verbose, log_file = log_file)
      invisible(self)
    },

    #' @description Best K (cross-entropy) from the fitted sNMF.
    #' @param stat Combine replicates by `"mean"` or `"min"`.
    best_k = function(stat = c("mean", "min")) snmf_best_k(private$snmf_fit, stat = match.arg(stat)),

    #' @description Q (ancestry) matrix at K, restricted to the active samples.
    #' @param K Number of ancestral populations (default the best K).
    #' @param run Replicate (default the lowest cross-entropy run).
    q = function(K = NULL, run = NULL) {
      if (is.null(private$snmf_fit)) stop("run_snmf() first", call. = FALSE)
      if (is.null(K)) K <- self$best_k()
      qq <- snmf_q(private$snmf_fit, K = K, run = run)
      qq[rownames(qq) %in% private$active_ids, , drop = FALSE]
    },

    #' @description Restrict to a set of samples in place (used by `subset()`).
    #' @param samples Sample ids to keep.
    restrict = function(samples) {
      private$active_ids <- intersect(private$sample_ids, samples)
      invisible(self)
    },

    #' @description A new `PopStructure` limited to given samples and/or metadata
    #'   matches (does not recompute PCA/UMAP/sNMF -- the embeddings are shared and
    #'   simply filtered).
    #' @param samples Sample ids to keep.
    #' @param ... `column = value(s)` metadata filters (e.g. `region = "West Africa"`).
    subset = function(samples = NULL, ...) {
      keep <- private$active_ids
      if (!is.null(samples)) keep <- intersect(keep, samples)
      filt <- list(...)
      if (length(filt) && !is.null(private$meta_df)) {
        m <- private$meta_df
        ok <- rep(TRUE, nrow(m))
        for (nm in names(filt)) if (nm %in% names(m)) ok <- ok & (m[[nm]] %in% filt[[nm]])
        keep <- intersect(keep, m$sample[ok])
      }
      new <- self$clone(deep = FALSE)
      new$restrict(keep)
      new
    },

    #' @description The genotype matrix for the active samples (samples x SNPs).
    genotype = function() private$geno_mat[private$active_ids, , drop = FALSE],
    #' @description PCA scores for the active samples.
    pca_scores = function() private$ps$pca[private$idx(), , drop = FALSE],
    #' @description PCA variance-explained table.
    pca_variance = function() private$ps$pca_var,
    #' @description The full [stats::prcomp()] object (all samples).
    prcomp = function() private$ps$prcomp,
    #' @description UMAP data frame for the active samples (or `NULL`).
    umap_df = function() if (is.null(private$ps$umap)) NULL else
      private$ps$umap[private$ps$umap$sample %in% private$active_ids, , drop = FALSE],
    #' @description Metadata for the active samples.
    get_meta = function() if (is.null(private$meta_df)) NULL else
      private$meta_df[private$meta_df$sample %in% private$active_ids, , drop = FALSE],
    #' @description The shared colour maps.
    get_colors = function() private$colors,
    #' @description Active sample ids.
    get_samples = function() private$active_ids,

    #' @description A `pop_structure` S3 view (active samples) for the `plot_*()` fns.
    as_ps = function() {
      idx <- private$idx()
      structure(list(samples = private$active_ids,
                     pca = private$ps$pca[idx, , drop = FALSE],
                     prcomp = private$ps$prcomp,
                     pca_var = private$ps$pca_var,
                     umap = self$umap_df(),
                     meta = self$get_meta()),
                class = "pop_structure")
    },

    #' @description PCA scatter coloured by a metadata column (shared colours).
    #' @param colour Metadata column to colour by.
    #' @param pcs Which two PCs.
    #' @param ... Passed to [plot_pca()].
    plot_pca = function(colour = NULL, pcs = c(1, 2), ...)
      plot_pca(self, pcs = pcs, colour = colour, ...),

    #' @description UMAP scatter coloured by a metadata column (shared colours).
    #' @param colour Metadata column to colour by.
    #' @param ... Passed to [plot_umap()].
    plot_umap = function(colour = NULL, ...) plot_umap(self, colour = colour, ...),

    #' @description Admixture bars; the group strip reuses the shared colour map, so it
    #'   matches the UMAP/PCA colouring.
    #' @param K Number of ancestral populations (default best K).
    #' @param group Metadata column to facet by and colour the strip with.
    #' @param colour Metadata column for the strip colours (default `group`).
    #' @param sample_order Explicit sample order (see [admixture_order()]).
    #' @param group_bar Draw the group colour strip.
    #' @param ... Passed to [plot_admixture()].
    plot_admixture = function(K = NULL, group = NULL, colour = group, sample_order = NULL,
                              group_bar = !is.null(group), ...) {
      qq <- self$q(K)
      plot_admixture(qq, meta = self$get_meta(), group = group, sample_order = sample_order,
                     group_bar = group_bar,
                     group_colours = if (is.null(colour)) NULL else private$colors[[colour]], ...)
    },

    #' @description Combined UMAP + admixture figure (see [plot_structure_figure()]).
    #' @param group Metadata column to facet/colour the admixture by.
    #' @param ... Passed to [plot_structure_figure()].
    plot_figure = function(group = NULL, ...) plot_structure_figure(self, group = group, ...),

    #' @description Per-SNP population differentiation between the levels of a metadata
    #'   column (see [pop_diff()]); uses the full genotype matrix for the active samples.
    #' @param group Metadata column defining the groups.
    #' @param ... Passed to [pop_diff()] (e.g. `statistic = "fst"`).
    pop_diff = function(group = NULL, ...) pop_diff(self, group = group, ...),

    #' @description Per-SNP Jost's D between the levels of a metadata column ([jost_d()]).
    #' @param group Metadata column defining the groups.
    #' @param ... Passed to [jost_d()].
    jost_d = function(group = NULL, ...) jost_d(self, group = group, ...),

    #' @description Group-pair differentiation table across all statistics
    #'   (see [pop_diff_table()]).
    #' @param group Metadata column defining the groups.
    #' @param ... Passed to [pop_diff_table()].
    pop_diff_table = function(group = NULL, ...) pop_diff_table(self, group = group, ...),

    #' @description Group x group differentiation triangle heatmap (see
    #'   [plot_diff_heatmap()]); metadata annotation is resolved against this object's
    #'   metadata automatically.
    #' @param group Metadata column defining the groups.
    #' @param ... Passed to [plot_diff_heatmap()], plus `statistic` (`"jost_d"` default,
    #'   `"gst"`, `"gst_hedrick"`, `"fst"`) selecting the measure.
    plot_diff_heatmap = function(group = NULL, ...) {
      dots <- list(...)
      statistic <- if (is.null(dots$statistic)) "jost_d" else dots$statistic
      dots$statistic <- NULL
      do.call(plot_diff_heatmap,
              c(list(pop_diff(self, group = group, statistic = statistic),
                     meta = private$meta_df), dots))
    },

    #' @description Alias of `plot_diff_heatmap()` for Jost's D.
    #' @param group Metadata column defining the groups.
    #' @param ... Passed to [plot_diff_heatmap()].
    plot_jost_d_heatmap = function(group = NULL, ...)
      plot_diff_heatmap(jost_d(self, group = group), meta = private$meta_df, ...),

    #' @description Save the whole workspace (genotype, PCA, UMAP, metadata, colours,
    #'   and sNMF fit) to an `.rds` file so it can be reloaded without recomputing.
    #'   Note: an sNMF fit references LEA project files on disk -- run `run_snmf()` with
    #'   a persistent `cache_dir` if you want the admixture to survive a reload.
    #' @param file Destination path.
    #' @param compress Passed to [saveRDS()] (default `"xz"` for a compact file).
    save = function(file, compress = "xz") {
      saveRDS(self, file, compress = compress)
      invisible(self)
    },

    #' @description Compact summary.
    #' @param ... Ignored.
    print = function(...) {
      cat("<PopStructure>", length(private$active_ids), "of",
          length(private$sample_ids), "samples,", ncol(private$ps$pca), "PCs\n")
      cat("  UMAP:", if (is.null(private$ps$umap)) "-" else "yes",
          "  sNMF:", if (is.null(private$snmf_fit)) "-"
                     else paste0("K ", paste(range(private$snmf_fit$K), collapse = "-")),
          "  meta:", if (is.null(private$meta_df)) "-"
                     else paste(setdiff(names(private$meta_df), "sample"), collapse = ", "),
          "\n")
      invisible(self)
    }
  ),
  private = list(
    geno_mat = NULL, sample_ids = NULL, active_ids = NULL, meta_df = NULL,
    colors = NULL, ps = NULL, snmf_fit = NULL, n_pcs = NULL,
    idx = function() match(private$active_ids, private$sample_ids)
  )
)

#' Load a saved PopStructure workspace
#'
#' Reads an `.rds` written by `PopStructure$save()` (or plain [saveRDS()]).
#'
#' @param file Path to the `.rds` file.
#' @return The [PopStructure] object.
#' @export
load_pop_structure <- function(file) {
  obj <- readRDS(file)
  if (!inherits(obj, "PopStructure"))
    stop("`file` does not contain a PopStructure object", call. = FALSE)
  obj
}

# ---- example ---------------------------------------------------------------

#' Public population-structure example datasets
#'
#' Builds a [PopStructure] object from a bundled **public** genotype matrix. Two datasets
#' ship:
#' * `"ghana_cambodia"` (default) -- 60 Pf7 samples (30 Ghana, 30 Cambodia) at 49
#'   biallelic SNPs; two well-separated populations, a clean minimal demo.
#' * `"africa"` -- 258 published East/Central-African samples (`country`: DRC, Kenya,
#'   Tanzania, Uganda; finer `site` sub-regions; a macro `region`) at the 2,000 SNPs that
#'   most differentiate the sites (top Jost's D; see [top_differentiating_snps()]), so the
#'   regional structure is clear -- a richer, multi-region demo for the combined UMAP +
#'   admixture figure ([plot_structure_figure()]).
#'
#' @param dataset Which bundled dataset to load.
#' @param umap Also compute a UMAP embedding (needs \pkg{uwot}); skipped with a message
#'   if the package is missing.
#' @param seed Random seed for the UMAP.
#' @return A [PopStructure] object with a `meta` data frame.
#' @examples
#' ps <- example_pop_structure(umap = FALSE)
#' ps
#' @export
example_pop_structure <- function(dataset = c("ghana_cambodia", "africa"),
                                  umap = TRUE, seed = 42) {
  dataset <- match.arg(dataset)
  file <- switch(dataset,
                 ghana_cambodia = "pop_structure_ghana_cambodia.rds",
                 africa = "pop_structure_africa.rds")
  f <- system.file("extdata", file, package = "plasgenomicsutilsR")
  if (!nzchar(f)) stop("example genotype data not found in the installed package",
                       call. = FALSE)
  d <- readRDS(f)
  ps <- PopStructure$new(d$genotype, meta = d$meta)
  if (umap) {
    if (requireNamespace("uwot", quietly = TRUE)) {
      # params that spread each dataset nicely out of the box
      if (dataset == "africa")
        ps$run_umap(pca_components = 0.8, n_neighbors = 25, min_dist = 0.4, seed = seed)
      else ps$run_umap(seed = seed)
    } else message("install 'uwot' to add a UMAP embedding to the example")
  }
  ps
}
