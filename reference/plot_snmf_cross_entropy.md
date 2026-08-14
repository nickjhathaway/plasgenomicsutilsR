# Cross-entropy elbow plot for choosing K

Cross-entropy against K, so the elbow (or the absence of one) is
visible. The line follows `stat`, and `show_range = TRUE` adds the
replicate min-max band – a wide band means the replicates disagreed and
that K is not reproducible.

## Usage

``` r
plot_snmf_cross_entropy(
  x,
  K = NULL,
  stat = c("min", "mean"),
  show_range = TRUE,
  best_k = NULL,
  point_size = 2.4,
  line_width = 0.6
)
```

## Arguments

- x:

  An
  [`run_snmf()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_snmf.md)
  result, or a table from
  [`snmf_cross_entropy()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/snmf_cross_entropy.md).

- K:

  Candidate K values (defaults to the fitted range).

- stat:

  Which summary the line follows, and which the red marker minimises:
  `"min"` (default, as in
  [`snmf_best_k()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/snmf_best_k.md))
  or `"mean"`. See
  [`snmf_best_k()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/snmf_best_k.md)
  for why `"min"` is the default and when `"mean"` is worth asking for.

- show_range:

  Draw the replicate min-max band (default `TRUE`).

- best_k:

  K to mark in red; `NULL` (default) marks the K minimising `stat`, `NA`
  marks none.

- point_size, line_width:

  Point and line sizes.

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- run_snmf(geno, K = 1:10, rep = 10)
snmf_cross_entropy(fit)          # the numbers
plot_snmf_cross_entropy(fit)     # the elbow
} # }
```
