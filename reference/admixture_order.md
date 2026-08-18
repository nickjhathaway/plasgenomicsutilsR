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

## Examples

``` r
q <- matrix(c(0.9, 0.1, 0.2, 0.8, 0.85, 0.15, 0.1, 0.9), ncol = 2, byrow = TRUE)
rownames(q) <- paste0("s", 1:4)
meta <- data.frame(sample = rownames(q), region = c("A", "A", "B", "B"))
admixture_order(q, meta = meta, group = "region")
#> [1] "s1" "s2" "s3" "s4"
```
