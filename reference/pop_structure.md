# Compute PCA and UMAP from a genotype matrix

Mean-imputes missing genotypes, runs PCA
([`stats::prcomp()`](https://rdrr.io/r/stats/prcomp.html)) and,
optionally, a UMAP embedding (uwot) with PCA initialisation. Returns a
`pop_structure` object the `plot_*()` functions read.

## Usage

``` r
pop_structure(
  geno,
  samples = NULL,
  meta = NULL,
  n_pcs = 50,
  umap = TRUE,
  umap_pca = 30,
  n_neighbors = 15,
  min_dist = 0.1,
  seed = 42
)
```

## Arguments

- geno:

  A genotype matrix (samples x SNPs, 0/1/2, `NA` allowed) or the list
  returned by
  [`load_genotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/load_genotypes.md).

- samples:

  Sample ids (defaults to the genotype list's `sample.id`, or matrix row
  names).

- meta:

  Optional per-sample metadata: a data frame with a `sample` column
  (plus e.g. `region`, `country`) used for colouring.

- n_pcs:

  Number of PCs to summarise (variance explained).

- umap:

  Compute a UMAP embedding.

- umap_pca, n_neighbors, min_dist:

  UMAP parameters (defaults 30 / 15 / 0.1); `n_neighbors` is clamped to
  the sample count.

- seed:

  Random seed for UMAP.

## Value

A `pop_structure` object (a list with `samples`, `pca` scores,
`pca_var`, `umap`, `meta`).

## Examples

``` r
# PCA (and optionally UMAP) from a genotype matrix
G <- matrix(stats::rbinom(30 * 80, 2, 0.4), 30, 80)
rownames(G) <- paste0("s", 1:30)
ps <- pop_structure(G, umap = FALSE)
dim(ps$pca)
#> [1] 30 30
head(ps$pca_var)
#>   PC var_explained cumulative
#> 1  1          7.86       7.86
#> 2  2          7.45      15.31
#> 3  3          6.91      22.22
#> 4  4          6.55      28.76
#> 5  5          5.90      34.66
#> 6  6          5.44      40.10
```
