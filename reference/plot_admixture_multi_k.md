# Admixture bar plots across every K, as pages

One
[`plot_admixture()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_admixture.md)
per K, returned as a named list – hand it straight to
[`save_plot()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/save_plot.md)
for a multi-page PDF, one K per page. Optionally leads with the
cross-entropy elbow
([`plot_snmf_cross_entropy()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_snmf_cross_entropy.md))
so the page that tells you which K to believe comes first, with the best
K marked in red.

## Usage

``` r
plot_admixture_multi_k(
  x,
  K = NULL,
  group = NULL,
  cross_entropy_first = TRUE,
  sample_order = NULL,
  sample_order_best_k = TRUE,
  best_k = NULL,
  stat = c("min", "mean"),
  meta = NULL,
  samples = NULL,
  ...
)
```

## Arguments

- x:

  A
  [PopStructure](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PopStructure.md)
  with a fitted sNMF, or an
  [`run_snmf()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_snmf.md)
  result (then give `meta` and `samples`).

- K:

  Which K values to draw (default: every fitted K with a Q matrix, K \>=
  2 – K = 1 is a single block and carries no information).

- group:

  Metadata column to facet by (e.g. `"region"`).

- cross_entropy_first:

  Lead with the cross-entropy elbow page (default `TRUE`).

- sample_order:

  Explicit sample order shared by every page. Overrides
  `sample_order_best_k`.

- sample_order_best_k:

  Derive one shared sample order from the best K's Q and use it on every
  page (default `TRUE`). `FALSE` lets each page cluster its own samples,
  so bars move between pages.

- best_k:

  The K to treat as best; `NULL` (default) uses
  [`snmf_best_k()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/snmf_best_k.md).

- stat:

  How
  [`snmf_best_k()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/snmf_best_k.md)
  and the elbow combine replicates: `"min"` (default) or `"mean"`.

- meta, samples:

  Metadata and sample ids, needed only when `x` is a raw
  [`run_snmf()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_snmf.md)
  result.

- ...:

  Passed to
  [`plot_admixture()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_admixture.md)
  (e.g. `group_bar`, `border`, `colours`).

## Value

A named list of ggplots: `"cross_entropy"` (when asked for) then
`"K=2"`, `"K=3"`, ... The best K's page is titled as such.

## Details

Bars stay in the same place from page to page when a shared
`sample_order` is used, which is what makes the pages comparable: a
sample sits at the same x on every K.

## See also

[`plot_snmf_cross_entropy()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_snmf_cross_entropy.md),
[`plot_admixture()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_admixture.md),
[`save_plot()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/save_plot.md)

## Examples

``` r
if (FALSE) { # \dontrun{
ps$run_snmf(K = 1:12)
pages <- ps$plot_admixture_multi_k(group = "region")
save_plot("admixture_all_k.pdf", pages)     # one K per page
} # }
```
