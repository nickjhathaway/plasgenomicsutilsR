# Per-sample coverage overview

Mean depth against the fraction of the genome reaching a usable depth,
one point per sample, with the QC floors drawn. The failures separate
along whichever axis they failed on, which is the quickest way to tell a
shallow run from an uneven one.

## Usage

``` r
plot_coverage_summary(
  cov,
  threshold = COVERAGE_QC_THRESHOLD,
  min_mean = COVERAGE_MIN_MEAN,
  min_breadth = COVERAGE_MIN_BREADTH,
  label_failures = TRUE
)
```

## Arguments

- cov:

  A coverage table from
  [`read_coverage()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/read_coverage.md).

- threshold:

  Depth whose breadth column is used (default 10; the table must have
  been produced with that threshold).

- min_mean, min_breadth:

  Floors for mean depth and percent of bases at `threshold`.

- label_failures:

  Name the samples that fail.

## Value

A ggplot object.

## Examples

``` r
cov <- data.frame(sample = paste0("s", 1:4), chrom = "genome",
                  mean = c(60, 45, 30, 3), pct_ge_10x = c(96, 91, 41, 7))
plot_coverage_summary(cov)
```
