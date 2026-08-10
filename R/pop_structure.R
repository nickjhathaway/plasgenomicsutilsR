# Population-structure analysis: LD-pruned genotypes -> PCA / UMAP, and sNMF
# admixture, with plot_*() functions. Heavy analysis packages (SNPRelate, gdsfmt,
# uwot, LEA) are optional (Suggests), guarded at call time.

# ---- genotypes -------------------------------------------------------------

# SNPRelate parses a VCF as text, so a binary BCF reaches it as mojibake and fails deep
# inside the parser with "invalid multibyte string". Hand it a VCF instead: reuse one
# already sitting next to the BCF if there is one, otherwise convert with bcftools into
# `vcf_dir` (default: alongside the BCF). Reusing rather than re-converting is the point
# -- these files are large and otherwise accumulate one copy per analysis.
.as_text_vcf <- function(path, vcf_dir = NULL) {
  if (!grepl("\\.bcf$", path, ignore.case = TRUE)) return(path)
  base <- sub("\\.bcf$", "", path, ignore.case = TRUE)
  out <- if (is.null(vcf_dir)) paste0(base, ".vcf.gz")
         else file.path(vcf_dir, paste0(basename(base), ".vcf.gz"))

  for (cand in unique(c(out, paste0(base, ".vcf.gz"), paste0(base, ".vcf")))) {
    if (!file.exists(cand)) next
    if (file.mtime(cand) < file.mtime(path))
      warning(sprintf("%s is older than the BCF beside it; delete it to reconvert",
                      basename(cand)), call. = FALSE)
    message("reusing the existing VCF: ", cand)
    return(cand)
  }

  bcftools <- Sys.which("bcftools")
  if (!nzchar(bcftools))
    stop(sprintf(paste0("`%s` is a BCF, which SNPRelate cannot read, and bcftools is ",
                        "not on PATH to convert it.\n  bcftools view -Oz -o %s %s"),
                 basename(path), out, path), call. = FALSE)
  if (!is.null(vcf_dir)) dir.create(vcf_dir, showWarnings = FALSE, recursive = TRUE)
  message("converting BCF to ", out, " (SNPRelate reads VCF text only)")
  log <- system2(bcftools, c("view", "-Oz", "-o", shQuote(out), shQuote(path)),
                 stdout = TRUE, stderr = TRUE)
  st <- attr(log, "status")
  if ((!is.null(st) && st != 0L) || !file.exists(out))
    stop("bcftools failed to convert ", path, ":\n  ",
         paste(utils::tail(log, 3), collapse = "\n  "), call. = FALSE)
  out
}

#' Load genotypes from a VCF, optionally LD-pruned
#'
#' Converts a VCF to GDS (only when needed) and returns the genotype matrix
#' (samples x SNPs, coded 0/1/2, `NA` for missing) via \pkg{SNPRelate}, LD-pruned by
#' default.
#'
#' Which you want depends on the question. Pruning is right for PCA, UMAP and admixture,
#' where correlated SNPs would let one locus dominate the structure. It is wrong wherever
#' the correlation between neighbouring SNPs *is* the signal -- differentiation
#' ([pop_diff()]) and haplotypes ([plot_region_haplotypes()]) -- because it keeps one SNP out
#' of each correlated run and drops the rest. Holding both is cheap: the GDS is reused, so a
#' second call with `prune = FALSE` only re-reads it.
#'
#' @param vcf Path to a (bgzipped) VCF, or a **BCF** -- SNPRelate reads VCF text only, so
#'   a BCF is converted first with `bcftools`, reusing any VCF already sitting next to it
#'   rather than making another copy.
#' @param gds Optional GDS path; derived from `vcf` if `NULL`.
#' @param prune LD-prune (default `TRUE`). `FALSE` returns **every** biallelic SNP,
#'   unpruned -- use this for the genotype matrix fed to [pop_diff()] /
#'   [pop_diff_table()], since LD-pruning removes the very SNPs that carry the
#'   differentiation signal.
#' @param ld_threshold,slide_max_bp,slide_max_n,autosome_only Passed to
#'   [SNPRelate::snpgdsLDpruning()] (defaults 0.2 / 20000 / 200 / `FALSE`); ignored
#'   when `prune = FALSE`.
#' @param maf,missing_rate Optional MAF / per-SNP missing-rate cutoffs for pruning.
#' @param seed Random seed for the pruning.
#' @param vcf_dir Where to put the VCF converted from a BCF (default: alongside the BCF).
#'   Point it somewhere scratch to keep converted copies out of the data directory.
#' @param allele Which allele the returned dosage counts. \pkg{SNPRelate} counts the
#'   **reference** allele; the default `"alt"` flips that so the matrix means what the
#'   rest of the package says it means. Only reported allele frequencies (and the
#'   arbitrary sign of a PCA axis) depend on this -- every diversity, differentiation, LD
#'   and selection statistic here is symmetric in `p` and `1 - p`.
#' @return A list with `genotype` (matrix; sample row names and `chr:pos` column names),
#'   `sample.id`, `snp.id`, and the two facts the matrix itself cannot carry: `allele` (which
#'   allele the dosages count) and `pruned`. [PopStructure] keeps both, so anything that names
#'   a call or warns about pruning can ask instead of assuming.
#' @seealso [PopStructure], [pop_structure()]
#' @export
load_genotypes <- function(vcf, gds = NULL, prune = TRUE, ld_threshold = 0.2,
                         slide_max_bp = 20000, slide_max_n = 200, autosome_only = FALSE,
                         maf = NaN, missing_rate = NaN, seed = 42, vcf_dir = NULL,
                         allele = c("alt", "ref")) {
  .need_package("SNPRelate", "load_genotypes()")
  .need_package("gdsfmt", "load_genotypes()")
  allele <- match.arg(allele)
  if (!file.exists(vcf)) stop(sprintf("no such file: %s", vcf), call. = FALSE)
  vcf <- .as_text_vcf(vcf, vcf_dir)
  if (is.null(gds)) gds <- sub("\\.vcf(\\.gz)?$", ".gds", vcf, ignore.case = TRUE)
  if (identical(gds, vcf)) gds <- paste0(vcf, ".gds")
  if (!file.exists(gds) || file.mtime(gds) < file.mtime(vcf)) {
    SNPRelate::snpgdsVCF2GDS(vcf, gds, method = "biallelic.only", verbose = FALSE)
  }
  h <- SNPRelate::snpgdsOpen(gds)
  on.exit(SNPRelate::snpgdsClose(h), add = TRUE)
  snpinfo <- SNPRelate::snpgdsSNPList(h)              # snp.id, chromosome, position
  if (prune) {
    set.seed(seed)
    snpset <- SNPRelate::snpgdsLDpruning(
      h, autosome.only = autosome_only, ld.threshold = ld_threshold,
      slide.max.bp = slide_max_bp, slide.max.n = slide_max_n,
      maf = maf, missing.rate = missing_rate, verbose = FALSE)
    snp_ids <- unlist(snpset, use.names = FALSE)
  } else {
    snp_ids <- snpinfo$snp.id
  }
  geno <- SNPRelate::snpgdsGetGeno(h, snp.id = snp_ids, with.id = TRUE, verbose = FALSE)
  idx <- match(geno$snp.id, snpinfo$snp.id)
  rownames(geno$genotype) <- geno$sample.id
  colnames(geno$genotype) <- paste0(snpinfo$chromosome[idx], ":", snpinfo$position[idx])
  # SNPRelate counts the *reference* allele, while this package documents and reports
  # alt dosage everywhere, so flip once here rather than leaving each caller to guess.
  # Every statistic downstream is symmetric in p <-> 1 - p, so this changes only reported
  # allele frequencies and the arbitrary sign of a PCA axis, never a differentiation,
  # diversity, LD or selection value.
  if (identical(allele, "alt")) geno$genotype <- 2L - geno$genotype
  # Carried through so downstream code never has to guess which allele a 2 means. Nothing in
  # a bare matrix says whether it counts reference or alternate alleles, and the two are
  # indistinguishable after the fact, so a plot that names the calls has to be told.
  list(genotype = geno$genotype, sample.id = geno$sample.id, snp.id = geno$snp.id,
       allele = allele, pruned = prune)
}

#' Deprecated name for load_genotypes()
#'
#' `run_ld_prune()` was renamed to [load_genotypes()]: the old name reads oddly for what it
#' mostly does, and reads as a contradiction with `prune = FALSE`. Kept so existing scripts
#' keep working.
#'
#' @param ... Passed to [load_genotypes()].
#' @return See [load_genotypes()].
#' @export
run_ld_prune <- function(...) {
  warning("`run_ld_prune()` has been renamed to `load_genotypes()`", call. = FALSE)
  load_genotypes(...)
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
#'   returned by [load_genotypes()].
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

# Grouping values for a set of samples as a factor, carrying the metadata column's own
# level order through. Every plot derives its group ordering from here, so setting levels
# on the metadata (or via PopStructure$set_levels()) drives them all the same way. A
# column that is not a factor is natural-sorted (see .levels_of()).
.group_factor <- function(meta, group, samples = NULL) {
  v <- meta[[group]]
  if (!is.null(samples)) v <- v[match(samples, meta$sample)]
  .as_group_factor(v)
}

# The levels actually present, in the column's level order (empty levels dropped).
.group_order <- function(v) {
  f <- .as_group_factor(v)
  levels(f)[levels(f) %in% as.character(f[!is.na(f)])]
}

# A factor keeps its own level order; anything else gets one from .levels_of().
.as_group_factor <- function(v) {
  if (is.factor(v)) v else factor(as.character(v), levels = .levels_of(v))
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
#'   from [load_genotypes()].
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

#' Cross-entropy of every sNMF replicate, summarised per K
#'
#' sNMF fits `rep` independent replicates at each K and scores each by cross-entropy
#' (lower is better). [snmf_q()] and [plot_admixture()] use the **minimum**-cross-entropy
#' replicate, so the `min` column is the one that describes the ancestry actually plotted;
#' `mean` and `max` show how much the replicates disagreed, which is worth a look before
#' trusting a K.
#'
#' Picking K is a separate judgement from picking a replicate: read the elbow of `min`
#' across K with [plot_snmf_cross_entropy()]. A curve that keeps falling or is flat means
#' the data do not support a well-defined K, whatever [snmf_best_k()] returns.
#'
#' @param x An [run_snmf()] result (or a raw LEA project, with `K` given).
#' @param K Candidate K values (defaults to the fitted range).
#' @return A tibble, one row per K: `K`, `n_runs`, `min`, `mean`, `max`, and `best_run`
#'   (the replicate index attaining `min`, i.e. the one [snmf_q()] returns).
#' @seealso [plot_snmf_cross_entropy()], [snmf_best_k()], [snmf_q()]
#' @export
snmf_cross_entropy <- function(x, K = NULL) {
  .need_package("LEA", "snmf_cross_entropy()")
  project <- .snmf_project(x)
  if (is.null(K)) {
    if (inherits(x, "snmf_fit")) K <- x$K
    else stop("pass K (the fitted range) when giving a raw sNMF project", call. = FALSE)
  }
  rows <- lapply(K, function(k) {
    ce <- as.numeric(LEA::cross.entropy(project, K = k))
    ce <- ce[is.finite(ce)]
    if (!length(ce)) {
      return(data.frame(K = k, n_runs = 0L, min = NA_real_, mean = NA_real_,
                        max = NA_real_, best_run = NA_integer_))
    }
    data.frame(K = k, n_runs = length(ce), min = min(ce), mean = mean(ce),
               max = max(ce), best_run = which.min(ce))
  })
  tibble::as_tibble(do.call(rbind, rows))
}

#' Cross-entropy elbow plot for choosing K
#'
#' Cross-entropy against K, so the elbow (or the absence of one) is visible. The line
#' follows `stat`, and `show_range = TRUE` adds the replicate min-max band -- a wide band
#' means the replicates disagreed and that K is not reproducible.
#'
#' @param x An [run_snmf()] result, or a table from [snmf_cross_entropy()].
#' @param K Candidate K values (defaults to the fitted range).
#' @param stat Which summary the line follows: `"min"` (default -- the replicate that
#'   [snmf_q()] plots) or `"mean"`.
#' @param show_range Draw the replicate min-max band (default `TRUE`).
#' @param best_k K to mark in red; `NULL` (default) marks the K minimising `stat`, `NA`
#'   marks none.
#' @param point_size,line_width Point and line sizes.
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' fit <- run_snmf(geno, K = 1:10, rep = 10)
#' snmf_cross_entropy(fit)          # the numbers
#' plot_snmf_cross_entropy(fit)     # the elbow
#' }
#' @export
plot_snmf_cross_entropy <- function(x, K = NULL, stat = c("min", "mean"),
                                    show_range = TRUE, best_k = NULL,
                                    point_size = 2.4, line_width = 0.6) {
  .need_package("ggplot2", "plot_snmf_cross_entropy()")
  stat <- match.arg(stat)
  ce <- if (is.data.frame(x) && all(c("K", "min", "mean") %in% names(x))) x
        else snmf_cross_entropy(x, K)
  ce <- ce[is.finite(ce[[stat]]), , drop = FALSE]
  if (!nrow(ce)) stop("no finite cross-entropy values to plot", call. = FALSE)
  if (is.null(best_k)) best_k <- ce$K[which.min(ce[[stat]])]
  ce$.y <- ce[[stat]]
  ce$.best <- !is.na(best_k) & ce$K == best_k

  p <- ggplot2::ggplot(ce, ggplot2::aes(.data$K, .data$.y))
  if (show_range && any(is.finite(ce$max))) {
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$min, ymax = .data$max),
                                  fill = "grey80", alpha = 0.5)
  }
  p <- p +
    ggplot2::geom_line(linewidth = line_width, colour = "grey30") +
    ggplot2::geom_point(ggplot2::aes(colour = .data$.best), size = point_size) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey25"),
                                 guide = "none") +
    ggplot2::scale_x_continuous(breaks = ce$K) +
    ggplot2::labs(x = "K", y = sprintf("cross-entropy (%s of replicates)", stat),
                  title = "sNMF cross-entropy by K") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  if (!is.na(best_k) && any(ce$.best)) {
    p <- p + ggplot2::annotate("text", x = best_k, y = ce$.y[ce$.best],
                               label = paste0("  best K = ", best_k),
                               hjust = 0, vjust = -0.9, colour = "firebrick", size = 3.4)
  }
  p
}

#' Admixture bar plots across every K, as pages
#'
#' One [plot_admixture()] per K, returned as a named list -- hand it straight to
#' [save_plot()] for a multi-page PDF, one K per page. Optionally leads with the
#' cross-entropy elbow ([plot_snmf_cross_entropy()]) so the page that tells you which K to
#' believe comes first, with the best K marked in red.
#'
#' Bars stay in the same place from page to page when a shared `sample_order` is used,
#' which is what makes the pages comparable: a sample sits at the same x on every K.
#'
#' @param x A [PopStructure] with a fitted sNMF, or an [run_snmf()] result (then give
#'   `meta` and `samples`).
#' @param K Which K values to draw (default: every fitted K with a Q matrix, K >= 2 --
#'   K = 1 is a single block and carries no information).
#' @param group Metadata column to facet by (e.g. `"region"`).
#' @param cross_entropy_first Lead with the cross-entropy elbow page (default `TRUE`).
#' @param sample_order Explicit sample order shared by every page. Overrides
#'   `sample_order_best_k`.
#' @param sample_order_best_k Derive one shared sample order from the best K's Q and use
#'   it on every page (default `TRUE`). `FALSE` lets each page cluster its own samples,
#'   so bars move between pages.
#' @param best_k The K to treat as best; `NULL` (default) uses [snmf_best_k()].
#' @param stat How [snmf_best_k()] and the elbow combine replicates: `"mean"` (default)
#'   or `"min"`.
#' @param meta,samples Metadata and sample ids, needed only when `x` is a raw
#'   [run_snmf()] result.
#' @param ... Passed to [plot_admixture()] (e.g. `group_bar`, `border`, `colours`).
#' @return A named list of ggplots: `"cross_entropy"` (when asked for) then `"K=2"`,
#'   `"K=3"`, ... The best K's page is titled as such.
#' @examples
#' \dontrun{
#' ps$run_snmf(K = 1:12)
#' pages <- ps$plot_admixture_multi_k(group = "region")
#' save_plot("admixture_all_k.pdf", pages)     # one K per page
#' }
#' @seealso [plot_snmf_cross_entropy()], [plot_admixture()], [save_plot()]
#' @export
plot_admixture_multi_k <- function(x, K = NULL, group = NULL,
                                   cross_entropy_first = TRUE,
                                   sample_order = NULL, sample_order_best_k = TRUE,
                                   best_k = NULL, stat = c("mean", "min"),
                                   meta = NULL, samples = NULL, ...) {
  .need_package("ggplot2", "plot_admixture_multi_k()")
  stat <- match.arg(stat)
  is_ps <- inherits(x, "PopStructure")
  fit <- if (is_ps) x$get_snmf_fit() else x
  if (is.null(fit)) stop("run_snmf() first", call. = FALSE)
  if (is.null(meta)) meta <- if (is_ps) x$get_meta() else NULL
  if (is.null(samples)) samples <- if (is_ps) x$get_samples() else fit$samples

  fitted_k <- if (inherits(fit, "snmf_fit")) fit$K else K
  if (is.null(K)) K <- fitted_k[fitted_k >= 2]
  K <- as.integer(K)
  if (!length(K)) stop("no K values to draw (K = 1 alone carries no structure)", call. = FALSE)
  if (is.null(best_k)) best_k <- snmf_best_k(fit, K = fitted_k, stat = stat)

  q_at <- function(k) {
    qq <- snmf_q(fit, K = k)
    if (!is.null(samples)) qq <- qq[rownames(qq) %in% samples, , drop = FALSE]
    qq
  }
  # one shared order keeps a sample at the same x on every page
  if (is.null(sample_order) && isTRUE(sample_order_best_k)) {
    sample_order <- admixture_order(q_at(best_k), meta = meta, group = group)
  }

  pages <- list()
  if (isTRUE(cross_entropy_first)) {
    pages[["cross_entropy"]] <- plot_snmf_cross_entropy(fit, K = fitted_k, stat = stat,
                                                        best_k = best_k)
  }
  for (k in K) {
    ttl <- if (!is.na(best_k) && k == best_k) sprintf("K = %d  (best by %s cross-entropy)", k, stat)
           else sprintf("K = %d", k)
    p <- plot_admixture(q_at(k), meta = meta, group = group,
                        sample_order = sample_order, ...) +
      ggplot2::labs(title = ttl)
    pages[[sprintf("K=%d", k)]] <- p
  }
  pages
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

# A legend of `n` keys, wrapped so it does not run off the canvas: down the side it grows
# into extra columns, along the bottom into extra rows. `rows` caps the keys per column
# (or per row when horizontal); NULL uses a default that keeps a right-hand legend inside a
# typical page.
.LEGEND_MAX_KEYS <- 10L

.LEGEND_MAX_ROWS_HORIZONTAL <- 2L

.legend_wrap <- function(n, rows = NULL, position = "right", order = 0) {
  if (identical(position, "none")) return("none")
  if (position %in% c("bottom", "top")) {
    # along the bottom, keep it shallow so it does not eat the panel's height
    nr <- if (is.null(rows)) .LEGEND_MAX_ROWS_HORIZONTAL else max(1L, as.integer(rows))
    return(ggplot2::guide_legend(order = order, nrow = min(n, nr)))
  }
  per <- if (is.null(rows)) .LEGEND_MAX_KEYS else max(1L, as.integer(rows))
  ggplot2::guide_legend(order = order, ncol = max(1L, ceiling(n / per)))
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
  grp <- factor(rep("all", length(samples)))
  if (!is.null(meta) && !is.null(group) && "sample" %in% names(meta) && group %in% names(meta)) {
    grp <- .group_factor(meta, group, samples)
  }
  ord <- character(0)
  # groups are swept in the column's level order, so samples sit in the requested order
  for (g in .group_order(grp)) {
    idx <- which(!is.na(grp) & grp == g)
    if (length(idx) > 2) {
      hc <- stats::hclust(stats::dist(q[idx, , drop = FALSE]), method = "ward.D2")
      idx <- idx[hc$order]
    }
    ord <- c(ord, samples[idx])
  }
  # samples whose group is not one of the levels are left out (set_levels() drops unlisted
  # levels), which silently shrinks the plot -- so say how many went
  if (anyNA(grp)) {
    warning(sum(is.na(grp)), " sample(s) have no level in the grouping column and are ",
            "not plotted; add their level with set_levels() to keep them", call. = FALSE)
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
#' @param border Outline each sample's bar (default `TRUE`, matching
#'   [plot_structure_figure()]) so neighbours with nearly identical ancestry stay distinct;
#'   `FALSE` for borderless bars. With very many samples in a narrow render the outlines
#'   can swamp the fills -- either render wider ([save_plot()] uses the attached width) or
#'   set `border = FALSE` / a thinner `border_linewidth`.
#' @param border_colour,border_linewidth Colour and width of the per-sample outline.
#' @param legend_position Where the legends go: `"right"` (default), `"bottom"`, `"top"`,
#'   `"left"`, or `"none"`. A large K makes for a tall legend stack, so `"bottom"` often
#'   fits better on a wide, short admixture panel.
#' @param legend_rows Keys per legend column (or per row when the legend is horizontal)
#'   before wrapping into another column/row. `NULL` (default) wraps a side legend every
#'   `r plasgenomicsutilsR:::.LEGEND_MAX_KEYS` keys and splits a horizontal one over two
#'   rows, which keeps `K` = 15 plus a group strip on the page. The suggested output height
#'   accounts for whatever this works out to.
#' @return A ggplot object.
#' @export
plot_admixture <- function(q, samples = NULL, meta = NULL, group = NULL,
                           order_within = TRUE, sample_order = NULL, colours = NULL,
                           group_bar = FALSE, group_colours = NULL,
                           border = TRUE, border_colour = "black",
                           border_linewidth = 0.15,
                           legend_position = c("right", "bottom", "top", "left", "none"),
                           legend_rows = NULL) {
  .need_package("ggplot2", "plot_admixture()")
  legend_position <- match.arg(legend_position)
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
  # reshape() leaves `cluster` a character, so K10 would sort before K2. Keep the Q
  # matrix's own column order, which is K1..K<K>.
  long$cluster <- factor(long$cluster, levels = colnames(q))

  cl_cols <- if (is.null(colours)) .pick_palette(ncol(q)) else colours
  p <- ggplot2::ggplot(long, ggplot2::aes(.data$sample, .data$q, fill = .data$cluster)) +
    ggplot2::geom_col(width = 1, position = "stack",
                      colour = if (border) border_colour else NA,
                      linewidth = border_linewidth) +
    ggplot2::scale_fill_manual(values = cl_cols, drop = FALSE,
                               guide = .legend_wrap(ncol(q), legend_rows, legend_position,
                                                    order = 1)) +
    ggplot2::labs(x = NULL, y = "ancestry", fill = "cluster") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank())

  if (group_bar && !is.null(group) && group %in% names(df)) {
    .need_package("ggnewscale", "the group colour bar")
    # the strip's facet column has to stay the same factor as the bars': a character copy
    # here makes ggplot combine the two layers' facet values alphabetically
    gdf <- data.frame(sample = factor(df$sample, levels = sample_order),
                      grp = .as_group_factor(df[[group]]), stringsAsFactors = FALSE)
    gdf[[group]] <- gdf$grp   # carry the facet variable so facet_grid subsets the strip
    if (is.null(group_colours)) group_colours <- meta_colors(df, cols = group)[[group]]
    p <- p + ggnewscale::new_scale_fill() +
      ggplot2::geom_tile(data = gdf, ggplot2::aes(.data$sample, y = 1.06, fill = .data$grp),
                         height = 0.05, inherit.aes = FALSE) +
      ggplot2::scale_fill_manual(values = group_colours, name = group,
                                 guide = .legend_wrap(length(group_colours), legend_rows,
                                                      legend_position, order = 2)) +
      ggplot2::coord_cartesian(clip = "off")
  }
  p <- p + ggplot2::theme(legend.position = legend_position,
                          legend.box = if (legend_position %in% c("bottom", "top"))
                            "horizontal" else "vertical")
  n_groups <- if (!is.null(group) && group %in% names(df)) length(unique(df[[group]])) else 1L
  if (!is.null(group) && group %in% names(df)) {
    p <- p + ggplot2::facet_grid(stats::as.formula(paste("~", group)),
                                 scales = "free_x", space = "free")
  }
  # Width scales with the number of bars so save_plot() makes it wide enough to read.
  # Height has to clear a side legend as well, or a large K clips its own keys -- measure
  # the built legend rather than estimating from the key count, since key size, text size
  # and how many columns the guide wrapped into all contribute.
  # Measured on a throwaway device: measuring on the current one leaves a stray Rplots.pdf
  # in a script and an empty figure in a notebook chunk (see .with_null_device()).
  panel_h <- 4
  height <- .with_null_device({
    gt <- try(ggplot2::ggplotGrob(p), silent = TRUE)
    if (inherits(gt, "try-error")) panel_h
    else if (legend_position %in% c("right", "left")) {
      # a side legend spans the canvas, so the canvas has to be at least that tall
      max(panel_h, .guide_box_height(gt, c("guide-box-right", "guide-box-left")) + 0.25)
    } else {
      # a legend along the bottom/top is a layout row: it adds to the height instead of
      # competing with it, so give the panel its space on top of the legend's
      panel_h + .guide_box_height(gt, c("guide-box-bottom", "guide-box-top"))
    }
  })
  attr(p, "plasgenomics_dims") <- c(
    width = round(max(6, min(24, nrow(q) * 0.06 + n_groups * 0.3)), 1),
    height = round(height, 1))
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
    #' @description Build from a genotype matrix (or a [load_genotypes()] list).
    #' @param geno Genotype matrix (samples x SNPs, 0/1/2, `NA`) or `load_genotypes()` list.
    #' @param samples Sample ids (default row names / the list's `sample.id`).
    #' @param meta Optional metadata (data frame with a `sample` column).
    #' @param n_pcs Number of PCs to summarise.
    #' @param colors Optional named list of colour maps (see [meta_colors()]).
    #' @param pruned Whether the SNPs were LD-pruned. Taken from a [load_genotypes()] list
    #'   when it says so. Pruning is right for PCA / admixture and wrong for looking at
    #'   haplotypes, since it drops SNPs precisely for being correlated with a neighbour.
    #' @param allele Which allele the dosages count, `"alt"` or `"ref"`. Taken from a
    #'   [load_genotypes()] list when it says so; a bare matrix cannot say, and the two codings
    #'   are indistinguishable afterwards, so anything that names the calls has to be told.
    initialize = function(geno, samples = NULL, meta = NULL, n_pcs = 50, colors = NULL,
                          allele = NULL, pruned = NULL, full = NULL) {
      if (is.list(geno) && !is.null(geno$genotype)) {
        if (is.null(samples)) samples <- geno$sample.id
        if (is.null(allele)) allele <- geno$allele
        if (is.null(pruned)) pruned <- geno$pruned
        geno <- geno$genotype
      }
      private$was_pruned <- if (is.null(pruned)) NULL else isTRUE(pruned)
      private$allele_counted <- if (is.null(allele)) NULL else
        match.arg(allele, c("alt", "ref"))
      mat <- as.matrix(geno)
      if (is.null(samples)) samples <- rownames(mat)
      if (is.null(samples)) samples <- as.character(seq_len(nrow(mat)))
      rownames(mat) <- samples
      private$geno_mat <- mat
      private$sample_ids <- samples
      private$active_ids <- samples
      # The primary panel is named for what it is, so `$genotype("full")` finds it when the
      # object was built unpruned and nothing else has to know which call made it.
      private$primary <- if (isTRUE(pruned)) "pruned" else
        if (isFALSE(pruned)) "full" else "genotypes"
      private$panel_list <- stats::setNames(
        list(list(genotype = mat, allele = private$allele_counted,
                  pruned = private$was_pruned)), private$primary)
      private$n_pcs <- n_pcs
      private$ps <- pop_structure(mat, samples = samples, n_pcs = n_pcs, umap = FALSE)
      private$colors <- if (is.null(colors)) list() else colors
      if (!is.null(full)) self$add_panel("full", full, pruned = FALSE)
      if (!is.null(meta)) self$add_meta(meta)
      invisible(self)
    },

    #' @description Register another genotype panel under a name, so one object can hold the
    #'   pruned SNPs for PCA / admixture and the full set for the analyses where the
    #'   correlation between neighbouring SNPs is the signal.
    #' @param name Panel name (`"full"` and `"pruned"` are the ones other functions ask for).
    #' @param geno Genotype matrix, or a [load_genotypes()] list.
    #' @param allele,pruned What this panel is; taken from a `load_genotypes()` list when it
    #'   says so.
    add_panel = function(name, geno, allele = NULL, pruned = NULL) {
      private$ensure_panels()
      if (!is.character(name) || length(name) != 1 || !nzchar(name))
        stop("`name` must be a single non-empty string", call. = FALSE)
      if (is.list(geno) && !is.null(geno$genotype)) {
        if (is.null(allele)) allele <- geno$allele
        if (is.null(pruned)) pruned <- geno$pruned
        if (is.null(rownames(geno$genotype)) && !is.null(geno$sample.id))
          rownames(geno$genotype) <- geno$sample.id
        geno <- geno$genotype
      }
      mat <- as.matrix(geno)
      if (is.null(rownames(mat)) || is.null(colnames(mat)))
        stop("a panel needs sample row names and `chr:pos` column names", call. = FALSE)
      gone <- setdiff(private$sample_ids, rownames(mat))
      if (length(gone))
        stop("panel \"", name, "\" is missing ", length(gone), " of this object's ",
             length(private$sample_ids), " samples", call. = FALSE)
      private$panel_list[[name]] <- list(
        genotype = mat[private$sample_ids, , drop = FALSE],
        allele = if (is.null(allele)) NULL else match.arg(allele, c("alt", "ref")),
        pruned = if (is.null(pruned)) NULL else isTRUE(pruned))
      invisible(self)
    },

    #' @description The panels this object holds, primary first.
    panels = function() {
      private$ensure_panels()
      unique(c(private$primary, names(private$panel_list)))
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

    #' @description Per-K cross-entropy summary of the sNMF replicates
    #'   (see [snmf_cross_entropy()]).
    #' @param ... Passed to [snmf_cross_entropy()].
    cross_entropy = function(...) {
      private$require_snmf()
      snmf_cross_entropy(private$snmf_fit, ...)
    },

    #' @description Cross-entropy elbow plot for choosing K
    #'   (see [plot_snmf_cross_entropy()]).
    #' @param ... Passed to [plot_snmf_cross_entropy()].
    plot_cross_entropy = function(...) {
      private$require_snmf()
      plot_snmf_cross_entropy(private$snmf_fit, ...)
    },

    #' @description The fitted sNMF result (`NULL` before `run_snmf()`).
    get_snmf_fit = function() private$snmf_fit,

    #' @description One admixture plot per K as pages, for a multi-page PDF
    #'   (see [plot_admixture_multi_k()]).
    #' @param ... Passed to [plot_admixture_multi_k()].
    plot_admixture_multi_k = function(...) {
      private$require_snmf()
      plot_admixture_multi_k(self, ...)
    },

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
    #' @param panel Panel to return by name; the primary one when `NULL`.
    #' @param prefer Panel to use *if the object has it*, falling back to the primary one --
    #'   how an analysis asks for the panel it wants without requiring it.
    genotype = function(panel = NULL, prefer = NULL) {
      nm <- private$pick_panel(panel, prefer)
      private$panel_list[[nm]]$genotype[private$active_ids, , drop = FALSE]
    },
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
    #' @description Which allele the dosages count (`"alt"` / `"ref"`), or `NULL` when the
    #'   object does not record it (built from a bare matrix, or saved by an older version).
    #' @param panel Which panel to report on; the primary one when `NULL`.
    allele = function(panel = NULL) private$panel_list[[private$pick_panel(panel)]]$allele,
    #' @description Whether the SNPs were LD-pruned (`TRUE` / `FALSE`), or `NULL` when the
    #'   object does not record it.
    #' @param panel Which panel to report on; the primary one when `NULL`.
    pruned = function(panel = NULL) private$panel_list[[private$pick_panel(panel)]]$pruned,
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
    #'   column (see [pop_diff()]); uses the object's genotype matrix (pass
    #'   `genotype = load_genotypes(vcf, prune = FALSE)` to run on the full unpruned set).
    #' @param group Metadata column defining the groups.
    #' @param ... Passed to [pop_diff()] (e.g. `statistic = "fst"`, `genotype = `).
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
    #'   `"gst_hedrick"`, `"fst"`) selecting the measure, and `genotype` (a full/unpruned
    #'   override, see [pop_diff()]).
    plot_diff_heatmap = function(group = NULL, ...) {
      dots <- list(...)
      statistic <- if (is.null(dots$statistic)) "jost_d" else dots$statistic
      genotype  <- dots$genotype
      dots$statistic <- NULL; dots$genotype <- NULL
      do.call(plot_diff_heatmap,
              c(list(pop_diff(self, group = group, statistic = statistic, genotype = genotype),
                     meta = private$meta_df), dots))
    },

    #' @description Alias of `plot_diff_heatmap()` for Jost's D.
    #' @param group Metadata column defining the groups.
    #' @param ... Passed to [plot_diff_heatmap()].
    plot_jost_d_heatmap = function(group = NULL, ...)
      plot_diff_heatmap(jost_d(self, group = group), meta = private$meta_df, ...),

    #' @description Per-SNP differentiation in long form (see [pop_diff_snps()]).
    #' @param group Metadata column defining the groups.
    #' @param ... Passed to [pop_diff()] (e.g. `statistic = `, `genotype = `).
    pop_diff_snps = function(group = NULL, ...) pop_diff_snps(pop_diff(self, group = group, ...)),

    #' @description Genome-wide differentiation Manhattan (see [plot_diff_manhattan()]).
    #' @param group Metadata column defining the groups.
    #' @param ... Passed to [plot_diff_manhattan()]; `statistic` / `genotype` go to [pop_diff()].
    plot_diff_manhattan = function(group = NULL, ...) {
      dots <- list(...)
      statistic <- if (is.null(dots$statistic)) "jost_d" else dots$statistic
      genotype  <- dots$genotype
      dots$statistic <- NULL; dots$genotype <- NULL
      do.call(plot_diff_manhattan,
              c(list(pop_diff(self, group = group, statistic = statistic, genotype = genotype)),
                dots))
    },

    #' @description Within-group diversity (see [pop_diversity()]).
    #' @param group Metadata column defining the groups.
    #' @param ... Passed to [pop_diversity()] (e.g. `by = `, `accessible = `).
    diversity = function(group = NULL, ...) pop_diversity(self, group = group, ...),

    #' @description Multilocus index of association (see [ld_index()]).
    #' @param group Metadata column defining the groups.
    #' @param ... Passed to [ld_index()].
    ld_index = function(group = NULL, ...) ld_index(self, group = group, ...),

    #' @description Beta scores for balancing selection (see [beta_score()]).
    #' @param group Metadata column defining the groups.
    #' @param ... Passed to [beta_score()].
    beta_score = function(group = NULL, ...) beta_score(self, group = group, ...),

    #' @description Phased haplotypes for a haplotype-homozygosity scan
    #'   (see [parasite_haplotypes()]).
    #' @param ... Passed to [parasite_haplotypes()] (e.g. `fws = `, `maf = `).
    haplotypes = function(...) parasite_haplotypes(self, ...),

    #' @description Genotype heatmap over one region, samples clustered
    #'   (see [plot_region_haplotypes()]).
    #' @param region The interval to draw.
    #' @param ... Passed to [plot_region_haplotypes()] (e.g. `split = `, `spacing = `).
    plot_region_haplotypes = function(region, ...)
      plot_region_haplotypes(self, region, ...),

    #' @description Integrated haplotype score (see [run_ihs()]); builds the haplotypes
    #'   first unless one is supplied.
    #' @param group Metadata column defining the groups.
    #' @param hap Optional [parasite_haplotypes()] object to reuse.
    #' @param ... Passed to [run_ihs()].
    ihs = function(group = NULL, hap = NULL, ...) {
      run_ihs(if (is.null(hap)) parasite_haplotypes(self) else hap, group = group, ...)
    },

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
    allele_counted = NULL, was_pruned = NULL, panel_list = NULL, primary = NULL,
    told = character(0),

    # Resolve a panel name. `panel` is a requirement (absent -> error); `prefer` is a wish
    # (absent -> the primary panel, with one note when that means handing pruned SNPs to an
    # analysis that asked for the full set).
    # An object saved before panels existed has `geno_mat` but no panel list, so build one on
    # first use rather than letting every accessor fail on a loaded object.
    ensure_panels = function() {
      if (is.null(private$panel_list) || !length(private$panel_list)) {
        private$primary <- if (isTRUE(private$was_pruned)) "pruned" else
          if (isFALSE(private$was_pruned)) "full" else "genotypes"
        private$panel_list <- stats::setNames(
          list(list(genotype = private$geno_mat, allele = private$allele_counted,
                    pruned = private$was_pruned)), private$primary)
      }
      invisible(TRUE)
    },

    pick_panel = function(panel = NULL, prefer = NULL) {
      private$ensure_panels()
      if (!is.null(panel)) {
        if (!panel %in% names(private$panel_list))
          stop("no panel called \"", panel, "\"; this object has: ",
               paste(names(private$panel_list), collapse = ", "), call. = FALSE)
        return(panel)
      }
      if (!is.null(prefer) && prefer %in% names(private$panel_list)) return(prefer)
      if (identical(prefer, "full") &&
          isTRUE(private$panel_list[[private$primary]]$pruned) &&
          !"full" %in% private$told) {
        private$told <- c(private$told, "full")
        message("this reads best on the full SNP set, but the object only holds a pruned ",
                "panel; add one with ",
                "`$add_panel(\"full\", load_genotypes(vcf, prune = FALSE))`")
      }
      private$primary
    },
    colors = NULL, ps = NULL, snmf_fit = NULL, n_pcs = NULL,
    idx = function() match(private$active_ids, private$sample_ids),
    require_snmf = function() {
      if (is.null(private$snmf_fit)) stop("run_snmf() first", call. = FALSE)
      invisible(NULL)
    }
  )
)

# An R6 object carries its own copy of the methods it was built with, so one saved by an
# earlier version of the package comes back without anything added since. Move the saved state
# into an instance built by the class as it stands now; a field the saved object never had
# keeps the current default.
.refresh_pop_structure <- function(obj) {
  old <- obj$.__enclos_env__$private
  fresh <- PopStructure$new(matrix(c(0L, 1L, 1L, 0L), nrow = 2,
                                   dimnames = list(c("a", "b"), NULL)))
  new <- fresh$.__enclos_env__$private
  for (nm in names(PopStructure$private_fields)) {
    v <- old[[nm]]
    if (!is.null(v)) new[[nm]] <- v
  }
  # Only non-NULL fields are copied, so an object saved before panels existed would keep the
  # placeholder's panel list -- which describes the two-sample dummy, not the genotypes just
  # copied over. Drop it and let it be rebuilt from those on first use.
  if (is.null(old$panel_list)) { new$panel_list <- NULL; new$primary <- NULL }
  fresh
}

#' Load a saved PopStructure workspace
#'
#' Reads an `.rds` written by `PopStructure$save()` (or plain [saveRDS()]). The workspace is
#' re-bound to the installed version of the class, so a file written by an older version of
#' the package gains the methods added since.
#'
#' @param file Path to the `.rds` file.
#' @return The [PopStructure] object.
#' @export
load_pop_structure <- function(file) {
  obj <- readRDS(file)
  if (!inherits(obj, "PopStructure"))
    stop("`file` does not contain a PopStructure object", call. = FALSE)
  .refresh_pop_structure(obj)
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
  # The shipped fixtures are alt dosage (load_genotypes()'s default). `ghana_cambodia` also
  # carries a `full` panel: the same sparse genome-wide SNPs plus every biallelic SNP around
  # pfcrt / pfdhps / pfkelch13, so the locus, haplotype and EHH examples have real density
  # while PCA / UMAP / admixture keep reading the thinned set they were tuned on.
  ps <- PopStructure$new(d$genotype, meta = d$meta, allele = d$allele %||% "alt",
                         pruned = d$pruned, full = d$full)
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
