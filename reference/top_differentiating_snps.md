# The SNPs that most differentiate groups

Picks markers from a
[`pop_diff()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff.md)
/
[`jost_d()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/jost_d.md)
result. `"roundrobin"` (default) walks the pairwise comparisons taking
each one's next-highest SNP in turn until `n` unique SNPs are collected
(so every pair contributes its top differentiators); `"max"` ranks SNPs
by their largest value across all pairs.

## Usage

``` r
top_differentiating_snps(jd, n, method = c("roundrobin", "max"))
```

## Arguments

- jd:

  A
  [`pop_diff()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff.md)
  /
  [`jost_d()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/jost_d.md)
  result.

- n:

  How many SNPs to return.

- method:

  `"roundrobin"` or `"max"`.

## Value

A character vector of SNP ids (a subset of the result's SNPs).

## Examples

``` r
ps <- example_pop_structure("africa", umap = FALSE)
top_differentiating_snps(pop_diff(ps, group = "site"), 5)
#> [1] "Pf3D7_08_v3:549992"  "Pf3D7_08_v3:544528"  "Pf3D7_08_v3:587930" 
#> [4] "Pf3D7_04_v3:1083335" "Pf3D7_07_v3:895873" 
```
