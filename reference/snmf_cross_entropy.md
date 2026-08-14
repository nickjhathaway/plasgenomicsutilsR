# Cross-entropy of every sNMF replicate, summarised per K

sNMF fits `rep` independent replicates at each K and scores each by
cross-entropy (lower is better).
[`snmf_q()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/snmf_q.md)
and
[`plot_admixture()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_admixture.md)
use the **minimum**-cross-entropy replicate, so the `min` column is the
one that describes the ancestry actually plotted; `mean` and `max` show
how much the replicates disagreed, which is worth a look before trusting
a K.

## Usage

``` r
snmf_cross_entropy(x, K = NULL)
```

## Arguments

- x:

  An
  [`run_snmf()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_snmf.md)
  result (or a raw LEA project, with `K` given).

- K:

  Candidate K values (defaults to the fitted range).

## Value

A tibble, one row per K: `K`, `n_runs`, `min`, `mean`, `max`, and
`best_run` (the replicate index attaining `min`, i.e. the one
[`snmf_q()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/snmf_q.md)
returns).

## Details

`min` is also what K is chosen on, by
[`snmf_best_k()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/snmf_best_k.md)
and by the
[`plot_snmf_cross_entropy()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_snmf_cross_entropy.md)
elbow, so the K you settle on and the ancestry you draw at it are scored
by the same number. A curve that keeps falling or is flat means the data
do not support a well-defined K, whatever
[`snmf_best_k()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/snmf_best_k.md)
returns; `n_runs` is worth a glance first, since a K whose replicates
mostly failed has its minimum taken over fewer of them.

## See also

[`plot_snmf_cross_entropy()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_snmf_cross_entropy.md),
[`snmf_best_k()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/snmf_best_k.md),
[`snmf_q()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/snmf_q.md)
