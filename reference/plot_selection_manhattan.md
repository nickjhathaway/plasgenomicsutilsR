# IBD selection-statistic Manhattan plot

The selection statistic along the genome, with the Bonferroni threshold
drawn as a dashed line when plotting `neg_log10_p` and thresholds are
available.

## Usage

``` r
plot_selection_manhattan(
  x,
  metric = c("neg_log10_p", "chi2_stat", "z_score"),
  regions = NULL,
  chroms = NULL,
  skip_chr = NULL,
  highlight_genes = NULL,
  label_genes = FALSE,
  draw_threshold = TRUE,
  point_size = 0.5,
  point_alpha = 0.6,
  colours = NULL
)
```

## Arguments

- x:

  An
  [IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
  object.

- metric:

  Which column to plot: `"neg_log10_p"` (default), `"chi2_stat"`, or
  `"z_score"`.

- regions:

  Optional character vector to keep (needs a `region` column).

- chroms:

  Optional chromosomes to keep (any spelling); others are dropped and
  the remaining ones re-laid-out contiguously.

- skip_chr:

  Optional chromosomes to drop (complement of `chroms`).

- highlight_genes:

  Optional gene names (from the object's `genes` track) to draw as
  reference lines; default all genes in the track.

- label_genes:

  Label the highlighted genes (top panel only, outside the plot).

- draw_threshold:

  Draw the significance threshold line (only for `neg_log10_p`).

- point_size, point_alpha:

  Point aesthetics.

- colours:

  Optional length-2 colour vector for the alternating chromosome bands.

## Value

A ggplot object.
