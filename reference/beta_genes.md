# Summarise beta scores per gene

Mean and maximum beta over the SNPs inside each gene, giving the
per-gene view used to rank candidates for balancing selection.

## Usage

``` r
beta_genes(b, genes = NULL, min_snps = 3)
```

## Arguments

- b:

  The tibble from
  [`beta_score()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/beta_score.md).

- genes:

  Gene table (`name`, `chr` or `chrom`, `start`, `end`); defaults to
  [PF3D7_GENES](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF3D7_GENES.md).
  Coordinates are 0-based half-open.

- min_snps:

  Genes with fewer scored SNPs than this are dropped.

## Value

A tibble with `group`, `gene`, `chr`, `start`, `end`, `n_snps`,
`beta_mean`, `beta_max`, sorted by `beta_mean` within each group.

## Examples

``` r
ps <- example_pop_structure(umap = FALSE)
b <- beta_score(ps, group = "country", window = 200000, min_window_snps = 1)
beta_genes(b, genes = PF_EXAMPLE_DRUG_GENES, min_snps = 1)
#> # A tibble: 5 × 8
#>   group    gene      chr     start     end n_snps beta_mean beta_max
#>   <fct>    <chr>     <chr>   <dbl>   <dbl>  <int>     <dbl>    <dbl>
#> 1 Cambodia pfdhps    8      548199  550616      4     8.67      8.67
#> 2 Cambodia pfcrt     7      403221  406317      9    -0.732     3.73
#> 3 Ghana    pfcrt     7      403221  406317      5    -0.520    16.8 
#> 4 Ghana    pfdhps    8      548199  550616      4    -3.08      4.49
#> 5 Ghana    pfkelch13 13    1724816 1726997      1    -3.94     -3.94
```
