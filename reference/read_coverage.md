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
