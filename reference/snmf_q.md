# Best-run Q (ancestry proportion) matrix for a given K

Best-run Q (ancestry proportion) matrix for a given K

## Usage

``` r
snmf_q(x, K, run = NULL)
```

## Arguments

- x:

  An
  [`run_snmf()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_snmf.md)
  result (or a raw LEA project).

- K:

  The number of ancestral populations.

- run:

  Which replicate; defaults to the lowest cross-entropy run.

## Value

A samples-by-K matrix of ancestry proportions, with sample ids as row
names when available.

## Examples

``` r
if (FALSE) { # \dontrun{
q <- snmf_q(fit, K = 5)           # the best-fitting replicate at K = 5
dim(q)
} # }
```
