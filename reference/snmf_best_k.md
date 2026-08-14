# Pick the best K from an sNMF fit by cross-entropy

sNMF fits `rep` replicates at each K, so a K has a spread of
cross-entropies rather than one. They are summarised by their
**minimum** across replicates, which is what LEA's own
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) of an sNMF
project shows, and the criterion its vignette reads a best K off. It is
also the replicate
[`snmf_q()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/snmf_q.md)
returns: sNMF's objective is non-convex, replicates land in different
local optima, and the best-fitting one is the model whose ancestry you
go on to plot – so comparing minima compares the models actually used at
each K.

## Usage

``` r
snmf_best_k(x, K = NULL, stat = c("min", "mean"))
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

  Combine replicates by `"min"` (default) or `"mean"` cross-entropy.

## Value

The K minimising the summarised cross-entropy.

## Details

`stat = "mean"` averages the replicates instead. Reach for it when the
replicate count differs between K values (`n_runs` in
[`snmf_cross_entropy()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/snmf_cross_entropy.md)),
since a minimum over more replicates is expected to be smaller whether
or not that K fits better; averaging is not sensitive to how many were
run.

## See also

[`snmf_cross_entropy()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/snmf_cross_entropy.md)
for the per-K spread,
[`plot_snmf_cross_entropy()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_snmf_cross_entropy.md)
for the elbow.
