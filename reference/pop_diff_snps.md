# Per-SNP differentiation in long form

Unpacks a
[`pop_diff()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff.md)
result (a SNP x pair matrix) into a tidy long table, parsing the
`chr:pos` SNP ids into genomic coordinates – so you can look at
differentiation SNP by SNP or feed it to
[`plot_diff_manhattan()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_diff_manhattan.md).

## Usage

``` r
pop_diff_snps(pd)
```

## Arguments

- pd:

  A
  [`pop_diff()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff.md)
  /
  [`jost_d()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/jost_d.md)
  result (its SNPs must be `chr:pos` ids, as from
  [`run_ld_prune()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_ld_prune.md)).

## Value

A tibble with `snp`, `chr`, `pos`, `a`, `b`, `pair`, `statistic`, and
`value` (one row per SNP x pair).
