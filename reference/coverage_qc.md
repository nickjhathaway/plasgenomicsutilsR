# Per-sample coverage QC verdict

Reduces a coverage table to one row per sample – the genome-wide row –
and applies the two floors a sample has to clear: enough average depth,
and enough of the genome actually reaching a usable depth. The second
matters more: selective whole-genome amplification can give a
respectable mean while leaving much of the genome at zero, and only the
breadth column shows that.

## Usage

``` r
coverage_qc(
  cov,
  threshold = COVERAGE_QC_THRESHOLD,
  min_mean = COVERAGE_MIN_MEAN,
  min_breadth = COVERAGE_MIN_BREADTH
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

## Value

A tibble with one row per sample: `sample`, `mean`, `median`, `sd`,
`pct_ge_<threshold>x`, `pct_zero`, `pass`, and `fail_reason`.

## Examples

``` r
cov <- data.frame(sample = c("a", "b"), chrom = "genome", mean = c(45, 4),
                  median = c(44, 0), sd = c(9, 7), pct_zero = c(1, 62),
                  pct_ge_10x = c(96, 30))
coverage_qc(cov)
#> # A tibble: 2 × 8
#>   sample  mean median    sd pct_zero pct_ge_10x pass  fail_reason          
#>   <chr>  <dbl>  <dbl> <dbl>    <dbl>      <dbl> <lgl> <chr>                
#> 1 b          4      0     7       62         30 FALSE low depth and breadth
#> 2 a         45     44     9        1         96 TRUE  NA                   
```
