# Read a coverage table

Identifier columns are read as text. Sample names are very often bare
digits – a sequencing id, or a BAM named after one – and guessing a type
would turn those into doubles, printing `4089106922` as `4.09e+09` and
silently failing to join against the metadata.

## Usage

``` r
read_coverage(path)
```

## Arguments

- path:

  TSV(.gz) written by `plasgenomicsutils coverage_depth_stats`
  (per-sample/per-chromosome), or by `coverage_dropout_regions`.

## Value

A tibble.

## Examples

``` r
f <- tempfile(fileext = ".tsv")
write.table(data.frame(sample = c("s1", "s2"), chrom = "genome",
                       mean = c(48, 4), pct_ge_10x = c(95, 12)),
            f, sep = "\t", quote = FALSE, row.names = FALSE)
read_coverage(f)
#> # A tibble: 2 × 4
#>   sample chrom   mean pct_ge_10x
#>   <chr>  <chr>  <dbl>      <dbl>
#> 1 s1     genome    48         95
#> 2 s2     genome     4         12
```
