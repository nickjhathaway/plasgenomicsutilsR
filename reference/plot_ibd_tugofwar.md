# IBD / selection "tug-of-war" mirror plot

The selection statistic hangs from the top (bars descending) while the
per-SNP IBD fraction rises from the bottom (bars ascending), sharing one
centred axis so peaks of each line up. A single left axis carries both
halves, its tick labels tinted to their track (selection on top, IBD on
the bottom).

## Usage

``` r
plot_ibd_tugofwar(
  x,
  region = NULL,
  metric = "neg_log10_p",
  scale = c("common", "free"),
  chroms = NULL,
  skip_chr = NULL,
  highlight_genes = NULL,
  label_genes = FALSE,
  draw_threshold = TRUE,
  selection_colour = "#fd8d3c",
  ibd_colour = "#2166ac"
)
```

## Arguments

- x:

  An
  [IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
  object (needs both selection and per_snp_region tables).

- region:

  Region(s) to plot. `NULL` (default) plots every region, faceted one
  per row; a single region gives one panel; a vector facets those
  regions.

- metric:

  Selection metric for the top track (default `"neg_log10_p"`).

- scale:

  `"common"` (default) scales every panel to a shared maximum so they
  are directly comparable; `"free"` scales each region to its own
  maximum for more per-region detail (the axis then reads as
  within-region percentages).

- chroms:

  Optional chromosomes to keep (any spelling); others are dropped and
  the rest re-laid-out contiguously.

- skip_chr:

  Optional chromosomes to drop (complement of `chroms`).

- highlight_genes:

  Optional gene names from the `genes` track to mark with lines.

- label_genes:

  Label the highlighted genes (top panel only, outside the plot).

- draw_threshold:

  Draw the per-region significance threshold (only for `neg_log10_p`
  with `scale = "common"`, when thresholds are available).

- selection_colour, ibd_colour:

  Bar and axis colours for the two tracks.

## Value

A ggplot object.
