# IBD Manhattan plot

Per-SNP fraction of pairs IBD along the genome. If the table carries a
`region` column the plot is faceted one row per region.

## Usage

``` r
plot_ibd_manhattan(
  x,
  regions = NULL,
  chroms = NULL,
  skip_chr = NULL,
  highlight_genes = NULL,
  label_genes = FALSE,
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

  Label the highlighted genes with their names.

- point_size, point_alpha:

  Point aesthetics.

- colours:

  Optional length-2 colour vector for the alternating chromosome bands.

## Value

A ggplot object.
