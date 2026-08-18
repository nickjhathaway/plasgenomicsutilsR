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
  [`load_genotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/load_genotypes.md)).

## Value

A tibble with `snp`, `chr`, `pos`, `a`, `b`, `pair`, `statistic`, and
`value` (one row per SNP x pair).

## Examples

``` r
ps <- example_pop_structure("africa", umap = FALSE)
sn <- pop_diff_snps(pop_diff(ps, group = "region"))
head(sn[order(-sn$value), c("snp", "chr", "pos", "value")])
#> # A tibble: 6 × 4
#>   snp                chr            pos value
#>   <chr>              <chr>        <dbl> <dbl>
#> 1 Pf3D7_08_v3:549992 Pf3D7_08_v3 549992 0.582
#> 2 Pf3D7_08_v3:544528 Pf3D7_08_v3 544528 0.565
#> 3 Pf3D7_08_v3:546840 Pf3D7_08_v3 546840 0.539
#> 4 Pf3D7_08_v3:536890 Pf3D7_08_v3 536890 0.437
#> 5 Pf3D7_04_v3:242142 Pf3D7_04_v3 242142 0.404
#> 6 Pf3D7_01_v3:469121 Pf3D7_01_v3 469121 0.335
```
