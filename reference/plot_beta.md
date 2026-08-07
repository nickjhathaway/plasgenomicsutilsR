# Manhattan plot of beta scores

Manhattan plot of beta scores

## Usage

``` r
plot_beta(
  b,
  threshold = NULL,
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

- b:

  The tibble from
  [`beta_score()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/beta_score.md).

- threshold:

  Draw a dashed line at this beta; `NULL` uses the 99th percentile of
  the scores actually plotted, the empirical tail these scans are read
  at.

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
b <- beta_score(ps, group = "country", window = 300000, min_window_snps = 1)
plot_beta(b)
```
