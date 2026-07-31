# Sample order for admixture bars

Orders samples within each `group` by hierarchical clustering of their
ancestry vectors. Compute it once (e.g. at the best K) and pass it as
`sample_order` to
[`plot_admixture()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_admixture.md)
for every K, so bars stay in the same position across K.

## Usage

``` r
admixture_order(q, samples = NULL, meta = NULL, group = NULL)
```

## Arguments

- q:

  A samples-by-K ancestry matrix.

- samples:

  Sample ids (default `rownames(q)`).

- meta, group:

  Optional metadata + the column to order within.

## Value

An ordered character vector of sample ids.
