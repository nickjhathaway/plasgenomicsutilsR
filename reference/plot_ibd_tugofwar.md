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
  group = NULL,
  metric = "neg_log10_p",
  scale = c("common", "free"),
  chroms = NULL,
  skip_chr = NULL,
  highlight_genes = NULL,
  label_genes = NULL,
  draw_threshold = TRUE,
  selection_colour = "#fd8d3c",
  ibd_colour = "#2166ac"
)
```

## Arguments

- x:

  An
  [IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
  object (needs both selection and per_snp_group tables).

- group:

  Group(s) to plot. `NULL` (default) plots every group, faceted one per
  row; a single group gives one panel; a vector facets those groups.

- metric:

  Selection metric for the top track (default `"neg_log10_p"`).

- scale:

  `"common"` (default) scales every panel to a shared maximum so they
  are directly comparable; `"free"` scales each group to its own maximum
  for more per-group detail (the axis then reads as within-group
  percentages).

- chroms:

  Optional chromosomes to keep (any spelling); others are dropped and
  the rest re-laid-out contiguously.

- skip_chr:

  Optional chromosomes to drop (complement of `chroms`).

- highlight_genes:

  Optional gene names from the `genes` track to mark with lines;
  requesting a name not in the track is an error.

- label_genes:

  Label the genes. `NULL` (default) labels them only when
  `highlight_genes` is given; `TRUE`/`FALSE` forces it.

- draw_threshold:

  Which significance line(s) to draw (only for `neg_log10_p`, and only
  when not `normalized`): `TRUE` / `"bonferroni"`, `"fdr"`,
  `"permutation"`, `"empirical"`, `"both"`, `"all"` or `FALSE` – the
  same set
  [`plot_selection_manhattan()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_selection_manhattan.md)
  takes, in the same colours, resolved by the same helper. Each line is
  mapped through the same transform as the mirrored selection half, so
  it lands where the data does.

- selection_colour, ibd_colour:

  Bar and axis colours for the two tracks.

## Value

A ggplot object.
