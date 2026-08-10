# Manhattan plot of a haplotype-homozygosity scan

Draws `neg_log10_p` (or the raw statistic) along the genome, one panel
per group or population pair, over the chromosome bands and optional
gene markers used by the other genome-wide plots in the package.

## Usage

``` r
plot_ihs(
  scan,
  metric = c("neg_log10_p", "ihs", "value"),
  threshold = NULL,
  genes = NULL,
  highlight_genes = NULL,
  label_genes = NULL,
  chroms = NULL,
  skip_chr = NULL,
  zoom = NULL,
  zoom_pad = 0.05,
  genes_for_track = NULL,
  gene_label_angle = 0,
  reference = DEFAULT_REFERENCE,
  point_size = 0.6,
  point_alpha = 0.7,
  colours = NULL
)
```

## Arguments

- scan:

  The tibble from
  [`run_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_ihs.md),
  [`run_rsb()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_rsb.md)
  or
  [`run_xpehh()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_xpehh.md).

- metric:

  `"neg_log10_p"` (default) or the statistic itself (`"ihs"` /
  `"value"`).

- threshold:

  Draw a dashed significance line at this height. `NULL` uses the
  `-log10(p)` of a 1% two-sided tail when plotting `neg_log10_p`, the
  convention these scans are usually read at; `NA` draws none.

- genes:

  Gene table to mark (e.g.
  [PF_EXAMPLE_DRUG_GENES](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF_EXAMPLE_DRUG_GENES.md));
  `NULL` for none.

- highlight_genes, label_genes:

  Which genes to mark and whether to name them.

- chroms, skip_chr:

  Chromosomes to keep or drop.

- zoom:

  Optional single interval to crop to, keeping the same data and the
  same coordinates as the genome-wide plot: a chromosome (`"7"`), a
  range (`"7:728,081-988,719"`), a gene name from `genes`, or a one-row
  data frame with chr/start/end. Every gene in the window is drawn and
  named unless `label_genes = FALSE`.

- zoom_pad:

  Context to add around `zoom`, clamped to the chromosome. One value
  pads both sides (default 5%); two pad the left and the right, either
  in that order or named – `c(left = 5000, right = 40000)`, and naming
  only one side pads only that side. Each side is read on its own: below
  1 it is a fraction of the interval's span, at or above 1 it is base
  pairs, so `c(0.1, 20000)` is legal.

- genes_for_track:

  Optional gene table for the track drawn beneath a zoomed plot (e.g.
  [PF3D7_GENES](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF3D7_GENES.md)),
  so every gene in the window is shown and named while the plot's own
  short track still supplies the marked positions inside the panel.
  Without it the track and the marks come from the same genes, which
  means marking a whole annotation just to see the neighbours.

- gene_label_angle:

  Rotation for the gene names in that track, in degrees. `0` (default)
  centres each name under its gene; `45` or `90` runs it down to the
  left, which is what keeps long systematic ids from colliding over a
  dense annotation.

- reference:

  Reference id for the chromosome layout.

- point_size, point_alpha, colours:

  Point and band aesthetics.

## Value

A ggplot object.

## Examples

``` r
ps <- example_pop_structure(umap = FALSE)
hap <- parasite_haplotypes(ps, maf = 0.05)
plot_ihs(run_ihs(hap, group = "country"), genes = PF_EXAMPLE_DRUG_GENES)
```
