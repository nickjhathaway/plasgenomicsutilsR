# Number of PCs explaining a target cumulative variance

Handy for setting how many principal components feed the UMAP: pass the
fraction of variance you want the PCA step to capture.

## Usage

``` r
n_pcs_for_variance(x, target = 0.8)
```

## Arguments

- x:

  A [`stats::prcomp()`](https://rdrr.io/r/stats/prcomp.html) result or a
  numeric vector of eigenvalues/variances.

- target:

  Cumulative variance to reach, as a proportion (`0.1`) or a percent
  (`10`).

## Value

The smallest number of leading PCs whose cumulative variance \>=
`target`.

## Examples

``` r
ps <- example_pop_structure(umap = FALSE)
# how many PCs to feed UMAP for a given share of the variance
n_pcs_for_variance(ps$prcomp(), 0.5)
#> [1] 5
n_pcs_for_variance(ps$prcomp(), 50)     # percent is accepted too
#> [1] 5
```
