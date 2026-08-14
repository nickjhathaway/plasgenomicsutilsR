# Summarise genome-wide IBD sharing between metadata groups

Reduces the per-pair IBD fraction table to one row per pair of `group`
levels, over every sample pair that spans them: two groups' row covers
all `n_a x n_b` cross-group pairs, and a group's row against itself
covers its `choose(n, 2)` within-group pairs. Within-group rows are the
ones that say whether a group is internally related at all, so they are
included by default and the comparison to the between-group rows is the
point.

## Usage

``` r
pair_fraction_summary(
  x,
  group = NULL,
  meta = NULL,
  value = .PAIR_FRACTION_COL,
  probs = c(0.25, 0.75),
  within = TRUE,
  min_pairs = 1L
)
```

## Arguments

- x:

  An
  [IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
  holding a pair table (`ibd_results(pair_fraction = )` or
  `$set_pair_fraction()`), or the pair table itself as a data frame or
  path.

- group:

  Metadata column defining the groups. Defaults to the object's declared
  group column, then to the first non-`sample` column of `meta`.

- meta:

  Sample metadata (`sample` plus `group`). Taken from `x` when it is an
  [IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md).

- value:

  Column holding the sharing measure (default
  `"ibd_fraction_accessible"`, what
  `plasgenomicsutils ibd_fraction_and_snp_density` writes). Pass
  `"ibd_fraction_full_genome"` to summarise against the full genome
  instead.

- probs:

  Quantiles to report, as proportions (default the quartiles,
  `c(0.25, 0.75)`). Each becomes a `q<pct>` column.

- within:

  Include each group's row against itself (default `TRUE`).

- min_pairs:

  Groups contributing fewer than this many pairs are dropped, with a
  note (default `1`, so only empty combinations go).

## Value

A tibble, one row per group pair, ordered by the object's group order
(or the column's factor levels):

- `group_a`, `group_b`:

  the two groups, `group_a == group_b` on a within-group row.

- `n_samples_a`, `n_samples_b`:

  samples of each group present in the pair table.

- `n_pairs`:

  sample pairs the summary is over. On a complete table this is
  `n_samples_a * n_samples_b` between groups and
  `choose(n_samples_a, 2)` within one; short of that, the table is
  missing pairs.

- `mean`, `median`, `sd`, `min`, `max`:

  over those pairs.

- `q25`, `q75`:

  the requested quantiles (type 7,
  [`stats::quantile()`](https://rdrr.io/r/stats/quantile.html)'s
  default).

## Details

The summary is over *pairs*, not samples, so a highly related cluster
inside one group pulls its mean up through every pair it takes part in.
The median and the quartiles are there because relatedness is skewed:
most pairs share almost nothing and a few share a great deal, so a mean
on its own reads as a much more related cohort than any pair actually
is.

## See also

[`plot_ibd_pair_network()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_pair_network.md),
which draws the same table pair by pair.

## Examples

``` r
pairs <- data.frame(
  sample1 = c("s1", "s1", "s1", "s2", "s2", "s3"),
  sample2 = c("s2", "s3", "s4", "s3", "s4", "s4"),
  ibd_fraction_accessible = c(0.02, 0.31, 0.04, 0.28, 0.03, 0.05))
meta <- data.frame(sample = paste0("s", 1:4),
                   region = c("north", "north", "south", "south"))
pair_fraction_summary(pairs, group = "region", meta = meta)
#> # A tibble: 3 × 12
#>   group_a group_b n_samples_a n_samples_b n_pairs  mean median     sd   min
#>   <fct>   <fct>         <int>       <int>   <int> <dbl>  <dbl>  <dbl> <dbl>
#> 1 north   north             2           2       1 0.02    0.02 NA      0.02
#> 2 north   south             2           2       4 0.165   0.16  0.151  0.03
#> 3 south   south             2           2       1 0.05    0.05 NA      0.05
#> # ℹ 3 more variables: max <dbl>, q25 <dbl>, q75 <dbl>
```
