# Coverage per chromosome, per sample

A sample-by-chromosome tile of mean depth (or any other column in the
table). A chromosome that is systematically low across samples is a
reference or amplification problem, not a sample problem – the two read
very differently here.

## Usage

``` r
plot_coverage_by_chrom(
  cov,
  metric = "mean",
  relative = TRUE,
  sample_order = NULL
)
```

## Arguments

- cov:

  A coverage table from
  [`read_coverage()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/read_coverage.md).

- metric:

  Column to fill by (default `"mean"`).

- relative:

  Divide each sample's value by its own genome-wide value, so the tile
  shows relative rather than absolute depth and deep samples do not
  dominate.

- sample_order:

  Optional sample order; defaults to increasing genome-wide mean.

## Value

A ggplot object.
