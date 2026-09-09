# The IBD pairs a network draws

The edge list behind
[`plot_ibd_pair_network()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_pair_network.md):
one row per sample pair sharing more than `min_ibd` of the genome.
[`ibd_pair_clusters()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_pair_clusters.md)
answers which samples those edges connect; this is the sharing itself.

## Usage

``` r
ibd_pair_links(
  pairs,
  weight = NULL,
  min_ibd = 0.01,
  samples = NULL,
  add_meta_cols = NULL,
  meta = NULL
)
```

## Arguments

- pairs:

  The per-pair table from
  `plasgenomicsutils ibd_fraction_and_snp_density`
  (`*.pair_ibd_fraction.tsv.gz`): a path, a data frame, or an
  [IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
  carrying one.

- weight:

  Column holding the fraction. Defaults to `ibd_fraction_accessible`,
  the callable-genome denominator; `ibd_fraction_full_genome` divides by
  the whole genome instead, so it reads lower for the same pair.

- min_ibd:

  Count a pair as linked only above this fraction (default `0.01`). Pass
  what you plotted with, or the clusters will not match the components.

- samples:

  Optional subset of samples to keep.

- add_meta_cols:

  Metadata columns to attach to both endpoints. Each `col` becomes
  `sample1_col` and `sample2_col`, in the order asked for. A sample
  missing from `meta` gets `NA`; factor columns keep their level order,
  so a region stays in map order rather than turning alphabetical.

- meta:

  Sample metadata for `add_meta_cols`; taken from the
  [IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
  when `pairs` is one, so it only needs giving alongside a plain table.

## Value

A tibble of `sample1`, `sample2` and `ibd_fraction` (the `weight`
column's value), highest sharing first, plus a pair of columns for each
of `add_meta_cols`.

## Details

`add_meta_cols` carries metadata onto both ends of each edge, which is
what turns the list into something to summarise over – whether a pair is
within a site or between two, say:

    e <- ibd_all$ibd_pair_links(min_ibd = 0.03, add_meta_cols = "region")
    table(within_region = e$sample1_region == e$sample2_region)

## See also

[`ibd_pair_clusters()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_pair_clusters.md),
[`plot_ibd_pair_network()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_pair_network.md),
[`pair_fraction_summary()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pair_fraction_summary.md)
for the sharing already summarised by group.

## Examples

``` r
if (FALSE) { # \dontrun{
ibd_pair_links(ibd_all, min_ibd = 0.03)

# the same edges, labelled with where each end came from
ibd_all$ibd_pair_links(min_ibd = 0.03, add_meta_cols = c("region", "country"))
} # }
```
