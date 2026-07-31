# Pick the best K from an sNMF fit by cross-entropy

Pick the best K from an sNMF fit by cross-entropy

## Usage

``` r
snmf_best_k(x, K = NULL, stat = c("mean", "min"))
```

## Arguments

- x:

  An
  [`run_snmf()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_snmf.md)
  result (or a raw LEA project, with `K` given).

- K:

  Candidate K values (defaults to the fitted range from
  [`run_snmf()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_snmf.md)).

- stat:

  Combine replicates by `"mean"` (default) or `"min"` cross-entropy.

## Value

The K minimising the summarised cross-entropy.
