# plasgenomicsutilsR

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/nickjhathaway/plasgenomicsutilsR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/nickjhathaway/plasgenomicsutilsR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

> **Version 0.1.1** — early development; APIs, defaults, and outputs may change
> between versions.

R utilities for **visualizing and analyzing Plasmodium genomics data** — the
R-side companion to the Python package
[`plasgenomicsutils`](https://github.com/nickjhathaway/plasgenomicsutils), which
does the heavy post-processing compute (VCF filtering/harmonization, IBD
analysis) and writes plain tables that this package reads and plots.

Reference-genome facts are namespaced by species (see `get_reference()`), so the
tools generalize beyond *Plasmodium falciparum*.

## Install

```r
# install.packages("remotes")
remotes::install_github("nickjhathaway/plasgenomicsutilsR")
```

That base install pulls only the required (CRAN) dependencies. Every plotting and
analysis dependency is a **`Suggests`** — so it is optional and *not* installed
automatically — and each function checks for what it needs and prints an install hint if
it is missing. You only install what you actually use:

- **CRAN** — the plots use `ggplot2`, `scales`, `patchwork`, `ggnewscale`, `ggtext`; the
  UMAP uses `uwot`.
- **Bioconductor** — the population-structure tools use **`SNPRelate`**, **`gdsfmt`**
  (LD-pruning a VCF → genotypes) and **`LEA`** (sNMF admixture).

Because these come from two repositories, install them Bioconductor-aware. The simplest
one-liner uses `pak`, which reads this package's `biocViews` and resolves CRAN **and**
Bioconductor deps together:

```r
# install.packages("pak")
pak::pak("nickjhathaway/plasgenomicsutilsR", dependencies = TRUE)
```

Or add the Bioconductor repositories first (your usual `setRepositories` trick), then let
`remotes` install the suggested packages:

```r
setRepositories(ind = 1:3)   # CRAN + Bioconductor software + annotation
remotes::install_github("nickjhathaway/plasgenomicsutilsR", dependencies = TRUE)
```

Or install the Bioconductor packages explicitly and the CRAN ones normally:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("SNPRelate", "gdsfmt", "LEA"))
install.packages(c("ggplot2", "scales", "patchwork", "ggnewscale", "ggtext", "uwot"))
```

(CI installs the full set the same way — `setup-r-dependencies` is `biocViews`-aware — so
`R CMD check` exercises the SNPRelate/LEA code paths, not just the CRAN ones.)

## What's here

- **`IbdResults` + `plot_*()`** — an R6 container over the tables the Python
  `plasgenomicsutils ibd` tools emit, and genome-wide plots reading from it:
  - `plot_ibd_manhattan()` — per-SNP fraction of pairs IBD along the genome
  - `plot_selection_manhattan()` — the IBD selection statistic, with the
    Bonferroni threshold line
  - `plot_ibd_region_heatmap()` — region × region IBD as tiles along the genome
  - `plot_ibd_tugofwar()` — selection (top) vs IBD (bottom) mirror for one region
  - `plot_drug_gene_triangles()` — region × region IBD-sharing triangles, per gene
    (SNPs strictly inside the gene, no flanking) or per specific locus via
    `snps = "Pf3D7_07_v3:403222"`. A gene's cell aggregates **all** its in-gene SNPs
    (`agg = "mean"`/`"median"`/`"max"` — nothing is picked or dropped). Returns a
    faceted grid, or `individual = TRUE` for a list of one-per-feature plots (e.g. a
    multi-page PDF) with the legend tucked into the empty upper triangle

  ggplot2 / scales are optional (Suggests); each plot returns a ggplot object.
- **Population structure** — `PopStructure` is an R6 workspace bundling a genotype
  matrix, its PCA (full `prcomp`), an optional UMAP, per-sample metadata, a **shared
  colour map**, and an sNMF admixture fit, so PCA / UMAP / admixture colour and order
  consistently and the whole thing can be `saveRDS`-ed and re-plotted without
  recomputing. It offers `run_umap()` (PC count *or* a fraction of variance to capture,
  via `n_pcs_for_variance()`), `run_snmf()` (**cached and quiet** — sNMF is slow and
  noisy), `subset()` (down to given samples or a metadata match), `best_k()` / `q()`,
  and `plot_pca()` / `plot_umap()` / `plot_admixture()` (bars ordered once via
  `admixture_order()`, per-sample borders so near-identical neighbours stay distinct, and
  a group strip coloured to match the UMAP). `set_levels()` fixes a metadata column's
  order once and it flows through every legend, facet, and strip; `save()` /
  `load_pop_structure()` persist the whole workspace. `run_ld_prune()` builds the genotype
  matrix from a VCF (SNPRelate LD-pruning). SNPRelate / uwot / LEA / ggnewscale / patchwork
  are optional (Suggests).
- **`plot_structure_figure()`** — the UMAP and sNMF admixture as **one** figure: a shared
  theme (matching fonts), one region colour map (UMAP points match the admixture colour
  strips), collected legends, region-faceted bars (colour strip, no text) laid out over
  configurable `rows`, and `orientation = "vertical"`/`"horizontal"`, auto-sized.
  `example_pop_structure()` ships two **public** demos: `"ghana_cambodia"` (minimal) and
  `"africa"` (258 East-African samples across DRC / Kenya / Tanzania / Uganda sites).
- **Population differentiation** — `pop_diff()` computes a per-SNP differentiation
  statistic for every pair of metadata groups from the **full, unpruned** genotypes:
  **Jost's D** (`jost_d()`), **Nei's Gst**, **Hedrick's standardized G′st**, and
  **Hudson's Fst**. `pop_diff_matrix()` collapses it to a group × group summary and
  `plot_diff_heatmap()` draws a triangle heatmap (legend in the empty corner, optional
  clustering **dendrogram** and one or more metadata **annotation** strips with custom
  colours), and `pop_diff_table()` returns every statistic × summary per group-pair in one
  data frame. Because most of the *P. falciparum* genome is barely differentiated, a
  genome-wide **mean looks near-zero** — use `stat = "top_mean"` (mean of the top few % of
  SNPs per pair), `"max"`, or a `trans = "sqrt"` fill to see the signal.
  `top_differentiating_snps()` picks the most differentiating markers (round-robin across
  pairs); the `"africa"` example fixture is itself built from them, which sharpens its
  UMAP/admixture structure. Estimators follow Jost (2008), Nei & Chesser (1983), Hedrick
  (2005), and Hudson et al. (1992) / Bhatia et al. (2013).
- Reference registry: `get_reference()`, `available_references()`,
  `normalise_chr()`, `PF3D7_CORE_CHROM_LENGTHS_BP`.

Fws is computed by the companion Python package (`plasgenomicsutils calculate_fws`),
a reimplementation of `moimix::getFws`.

A small **public** example dataset (five African countries) ships with the package:

```r
library(plasgenomicsutilsR)

ibd <- example_ibd_results()          # bundled public example
plot_ibd_manhattan(ibd)
plot_selection_manhattan(ibd, metric = "neg_log10_p")
plot_ibd_tugofwar(ibd, region = "Tanzania")
plot_ibd_region_heatmap(ibd, trans = "log2")   # log2 fill reads best for IBD
```

The genome-wide plots share a few options:

- **Colour scale** (`plot_ibd_region_heatmap`, `plot_drug_gene_triangles`):
  `trans` (`"log2"`, `"sqrt"`, …), a custom `colors` ramp, `limits = c(lo, hi)`
  (extremes squished so a few high values don't crush the scale), or a full
  `fill_scale` override. The default is a light single-hue ramp.
- **Skip chromosomes**: `chroms = c("7", "13")` keeps only those, or
  `skip_chr = "1"` drops them — the rest are re-laid-out contiguously.
- **Highlight genes**: gene positions *and* display names come from the `genes`
  track you pass to `ibd_results()` — a data frame of `name`, `chr`, `start`, `end`
  (the bundled [`EXAMPLE_DRUG_GENES`] is one). The Manhattan, tug-of-war, and region
  heatmap draw a reference line at each gene; `highlight_genes = c("crt", "dhps")`
  selects which to show (case-insensitive) and `label_genes = TRUE` labels them (top
  panel only, just outside the plot). The label text is the track's `name`, so set it
  to whatever you want displayed (`"CRT"` vs `"crt"`) by passing your own track:

```r
my_genes <- data.frame(name = c("CRT", "DHPS"), chr = c("7", "8"),
                       start = c(403222, 548200), end = c(406317, 550616))
ibd <- ibd_results(per_snp_region = "...", selection = "...", genes = my_genes)
plot_selection_manhattan(ibd, chroms = c("7", "8"),
                         highlight_genes = c("crt", "dhps"), label_genes = TRUE)
```

On your own data, point `ibd_results()` at the tables the Python `plasgenomicsutils
ibd` tools write:

```r
ibd <- ibd_results(
  per_snp_region  = "ibd_analysis.per_snp_per_region.tsv.gz",
  pairwise_region = "ibd_analysis.per_snp_pairwise_region.tsv.gz",
  selection       = "ibd_selection_analysis.per_region.selection_stats.tsv.gz",
  threshold       = "ibd_selection_analysis.per_region.threshold.txt",
  reference       = "pf3d7"
)
save_plot("ibd_manhattan.pdf", plot_ibd_manhattan(ibd), width = 9, height = 4)
```

`save_plot()` wraps `ggplot2::ggsave()` and, for `.pdf` output, defaults to the cairo
PDF device (better font embedding), falling back to the standard `pdf` device where
cairo is unavailable or unreliable (e.g. Windows). Force a device with `device =`, or
use `pdf_device()` directly with `ggsave()`.

Drug-gene triangles read a `genes` track (`name`, `chr`, `start`, `end`) passed
to `ibd_results()`; a SNP belongs to a gene when its position falls in the gene
interval:

```r
ibd <- ibd_results(pairwise_region = "...pairwise_region.tsv.gz",
                   genes = "drug_resistance_genes.tsv")
plot_drug_gene_triangles(ibd)
```

## Development

```r
devtools::load_all()
devtools::test()
devtools::document()
```

## License

GPL-3.
