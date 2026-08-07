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
  reference = DEFAULT_REFERENCE,
  point_size = 0.6,
  point_alpha = 0.7,
  colours = NULL
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

- reference:

  Reference id for the chromosome layout.

- point_size, point_alpha, colours:

  Point and band aesthetics.

## Value

A ggplot object.

## Examples

``` r
ps <- example_pop_structure(umap = FALSE)
d <- pop_diversity(ps, group = "country", by = "window", window = 500000)
plot_diversity(d, metric = "he")
```
