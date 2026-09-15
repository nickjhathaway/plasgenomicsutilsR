# Block extension between every pair of groups, as a triangle per locus

The length counterpart of
[`plot_pairwise_ibd_for_genes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_pairwise_ibd_for_genes.md).
That plot fills each group-pair cell with the *fraction of pairs
sharing* a gene; this one fills it with how much **longer** those shared
segments are than the sharing pairs' own genome-wide background, the
`paired_ratio` of
[`ibd_block_extension_test()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_block_extension_test.md).

## Usage

``` r
plot_pairwise_block_extension(
  x,
  loci = NULL,
  min_pairs = 5L,
  group = NULL,
  within = 0,
  sharing = c("overlap", "complete"),
  min_ref_blocks = 1L,
  individual = FALSE,
  label = TRUE,
  digits = 2,
  ncol = NULL,
  limits = "shared",
  colors = NULL,
  fill_scale = NULL,
  meta = NULL,
  colours = NULL
)
```

## Arguments

- x:

  An
  [IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
  built with `blocks =` and `meta =`, or a data frame already returned
  by `ibd_block_extension_test(..., pairs = "all")`.

- loci:

  Loci to draw, as
  [`ibd_block_extension_test()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_block_extension_test.md)
  takes them. Ignored when `x` is already a result table.

- min_pairs:

  Group pairs with fewer sharing pairs are left grey (default `5`).

- group, within, sharing, min_ref_blocks, meta:

  Passed to
  [`ibd_block_extension_test()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_block_extension_test.md).

- individual:

  Return a named list of one plot per locus instead of a faceted grid.

- label:

  Draw the ratio and pair count on each tile (default `TRUE`).

- digits:

  Decimal places for the ratio (default `2`).

- ncol:

  Columns in the facet grid.

- limits:

  Fill limits, or `"shared"` (default) to pin every locus to one
  symmetric scale so the pages are comparable.

- colors, colours:

  Fill ramp, low to high through the midpoint.

- fill_scale:

  A complete `ggplot2` fill scale, replacing the built-in one.

## Value

A `ggplot`, or a named list of them when `individual = TRUE`.

## Details

The diagonal is the within-group result. The off-diagonal is what the
diagonal cannot say: whether the long segments continue across a group
boundary, which is the difference between a haplotype that is spreading
between regions and one that is expanding inside each.

## Reading the tiles

The fill diverges around **1** on a log2 scale, so 2 and 0.5 sit the
same distance from the middle. Note that 1 is *not* the null – every
locus is inflated by length-biased sampling, so a whole triangle sitting
above 1 is the expected picture, not a finding. Use
[`ibd_block_extension_scan()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_block_extension_scan.md)
to say what a given ratio is worth in a given group.

Cells are grey where the group pair had fewer than `min_pairs` sharing
pairs, which is common and is not the same as a ratio near 1. Every
drawn cell carries its pair count under the ratio, because a cell
resting on five pairs and one resting on five hundred look identical
otherwise and routinely differ by more than the colour does.

## See also

[`ibd_block_extension_test()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_block_extension_test.md),
[`ibd_block_extension_scan()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_block_extension_scan.md),
[`plot_pairwise_ibd_for_genes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_pairwise_ibd_for_genes.md)
for the sharing-fraction counterpart.

## Examples

``` r
if (FALSE) { # \dontrun{
ibd <- ibd_results(blocks = "hmm.txt", meta = meta, group_col_in_meta = "region")
plot_pairwise_block_extension(ibd, c("pfcrt", "pfgch1"))
} # }
```
