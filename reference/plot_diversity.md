# Windowed diversity along the genome

One panel per group, with the statistic plotted at the window
mid-points.

## Usage

``` r
plot_diversity(
  div,
  metric = "pi",
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
  colours = NULL,
  colors = NULL
)
```

## Arguments

- div:

  A
  [`pop_diversity()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diversity.md)
  result computed with `by = "window"`.

- metric:

  Which column to draw: `"pi"` (default), `"he"`, `"theta_w"`,
  `"tajima_d"`, `"hap_div"`, `"seg_sites"`, ...

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

- point_size, point_alpha, colours, colors:

  Point and band aesthetics.

## Value

A ggplot object.

## Examples

``` r
ps <- example_pop_structure(umap = FALSE)
d <- pop_diversity(ps, group = "country", by = "window", window = 500000)
plot_diversity(d, metric = "he")
```
