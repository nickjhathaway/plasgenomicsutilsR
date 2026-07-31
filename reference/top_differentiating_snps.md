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
