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
