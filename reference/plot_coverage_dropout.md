# Coverage dropouts along the genome

The fraction of samples below depth in each window, across the genome,
with the merged dropout regions marked. Regions where the line sits near
1 are amplified in almost nobody – they will read as invariant rather
than as missing unless they are excluded.

## Usage

``` r
plot_coverage_dropout(
  windows,
  min_depth = 5,
  min_frac_samples = 0.9,
  genes = NULL,
  highlight_genes = NULL,
  label_genes = NULL,
  chroms = NULL,
  skip_chr = NULL,
  reference = DEFAULT_REFERENCE
)
```

## Arguments

- windows:

  Per-window table from `coverage_depth_stats --windows-output`, or the
  already-merged regions from `coverage_dropout_regions`.

- min_depth:

  A sample counts as uncovered in a window below this mean depth
  (ignored when `windows` is already merged).

- min_frac_samples:

  Draw the flag line at this fraction.

- genes, highlight_genes, label_genes:

  Optional gene track, as in
  [`plot_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ihs.md).

- chroms, skip_chr, reference:

  As in
  [`plot_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ihs.md).

## Value

A ggplot object.
