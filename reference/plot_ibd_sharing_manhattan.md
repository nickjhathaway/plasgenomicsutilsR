# IBD Manhattan plot

Per-SNP fraction of pairs IBD along the genome. If the table carries a
`group` column the plot is faceted one row per group.

## Usage

``` r
plot_ibd_sharing_manhattan(
  x,
  groups = NULL,
  chroms = NULL,
  skip_chr = NULL,
  highlight_genes = NULL,
  label_genes = NULL,
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

- groups:

  Optional character vector to keep (needs a `group` column).

- chroms:

  Optional chromosomes to keep (any spelling); others are dropped and
  the remaining ones re-laid-out contiguously.

- skip_chr:

  Optional chromosomes to drop (complement of `chroms`).

- highlight_genes:

  Optional gene names (from the object's `genes` track) to draw as
  reference lines; default all genes in the track. Requesting a name not
  in the track is an error (rather than silently drawing nothing).

- label_genes:

  Label the genes with their names. `NULL` (default) labels them only
  when `highlight_genes` is given; `TRUE`/`FALSE` forces it.

- point_size, point_alpha:

  Point aesthetics.

- colours:

  Optional length-2 colour vector for the alternating chromosome bands.

## Value

A ggplot object.
