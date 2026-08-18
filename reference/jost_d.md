# Per-SNP Jost's D between metadata groups

Convenience wrapper for `pop_diff(..., statistic = "jost_d")`. See
[`pop_diff()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff.md).

## Usage

``` r
jost_d(x, group = NULL, meta = NULL, clamp = TRUE)
```

## Arguments

- x:

  A
  [PopStructure](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PopStructure.md)
  (uses its full genotype matrix for the active samples) or a genotype
  matrix (samples x SNPs, 0/1/2, `NA`).

- group:

  Metadata column (for a `PopStructure`) or a per-sample grouping vector
  (for a matrix).

- meta:

  When `x` is a matrix, an optional data frame with a `sample` column
  plus `group`; otherwise `group` is a vector aligned to the matrix
  rows.

- clamp:

  Clamp small negative estimates to 0 (default `TRUE`).

## Value

A `pop_diff` object (of Jost's D).

## Examples

``` r
ps <- example_pop_structure("africa", umap = FALSE)
jost_d(ps, group = "region")$groups
#> [1] "Central Africa" "East Africa"   
```
