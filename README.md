# plasgenomicsutilsR

R utilities for **visualizing and analyzing Plasmodium genomics data** — the
R-side companion to the Python package
[`plasgenomicsutils`](https://github.com/nickjhathaway/plasgenomicsutils), which
does the heavy post-processing compute (VCF filtering/harmonization, IBD
analysis) and writes plain tables that this package reads and plots.

Reference-genome facts are namespaced by species (see `get_reference()`), so the
tools generalize beyond *Plasmodium falciparum*.

## Install

```r
# install.packages("devtools")
devtools::install_github("nickjhathaway/plasgenomicsutilsR")
```

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
