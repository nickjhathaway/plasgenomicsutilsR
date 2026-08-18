# Group x group differentiation summary matrix

Collapses the per-SNP pairwise values into a symmetric group-by-group
matrix. Because most of the genome is barely differentiated in *P.
falciparum*, a genome-wide `"mean"` looks near-zero; `"top_mean"` (mean
of the top `top` fraction of SNPs per pair) or `"max"` surface the
differentiation that lives in a minority of loci.

## Usage

``` r
pop_diff_matrix(pd, stat = c("mean", "median", "top_mean", "max"), top = 0.05)

jost_d_matrix(pd, stat = "mean", top = 0.05)
```

## Arguments

- pd:

  A
  [`pop_diff()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff.md)
  /
  [`jost_d()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/jost_d.md)
  result.

- stat:

  `"mean"` (default), `"median"`, `"top_mean"`, or `"max"`.

- top:

  Fraction of highest-value SNPs to average for `stat = "top_mean"`
  (default `0.05`, i.e. the top 5%).

## Value

A symmetric numeric matrix (0 on the diagonal).

## Examples

``` r
ps <- example_pop_structure("africa", umap = FALSE)
pd <- pop_diff(ps, group = "site")
# the genome-wide mean is near zero between neighbours; the tail separates them
round(pop_diff_matrix(pd, "top_mean", top = 0.05)[1:3, 1:3], 4)
#>               DRC Kenya_East Kenya_West
#> DRC        0.0000     0.2486     0.2781
#> Kenya_East 0.2486     0.0000     0.2483
#> Kenya_West 0.2781     0.2483     0.0000
```
