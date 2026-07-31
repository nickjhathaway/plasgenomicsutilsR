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
  [`run_ld_prune()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_ld_prune.md).

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
