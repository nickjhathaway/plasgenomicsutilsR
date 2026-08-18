# Group-by-group IBD heatmap along the genome

Per-SNP IBD sharing between group pairs as tiles along the genome, one
facet per anchor group. The stored upper-triangle (`group_a <= group_b`)
is mirrored so every anchor shows all partner groups. Alternating grey
chromosome bands and thin boundary lines mark where each chromosome
starts and ends, and the x-axis spans the **full chromosome lengths** so
un-genotyped (empty) regions are visible rather than collapsed out.

## Usage

``` r
plot_ibd_pairwise_group_heatmap(
  x,
  anchor = NULL,
  chroms = NULL,
  skip_chr = NULL,
  trans = "identity",
  colors = NULL,
  limits = NULL,
  fill_scale = NULL,
  highlight_genes = NULL,
  label_genes = NULL,
  colours = NULL
)
```

## Arguments

- x:

  An
  [IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
  object.

- anchor:

  Optional single group to show (one panel) instead of all.

- chroms:

  Optional chromosomes to keep (any spelling); others are dropped and
  the remaining ones re-laid-out contiguously.

- skip_chr:

  Optional chromosomes to drop (complement of `chroms`).

- trans:

  Fill-scale transform, e.g. `"identity"` (default), `"log2"`, `"sqrt"`.

- colors, colours:

  Optional colour ramp for the fill (defaults to a single-hue
  light-to-dark sequential scale that stays readable when most values
  are near 0).

- limits:

  Optional `c(lo, hi)` fill limits; values outside are squished into
  range, so a few extremes near 1 don't crush the rest of the scale.

- fill_scale:

  Optional ggplot2 fill scale that fully overrides the above.

- highlight_genes:

  Optional gene names from the `genes` track to mark with lines;
  requesting a name not in the track is an error.

- label_genes:

  Label the genes above the top panel. `NULL` (default) labels them only
  when `highlight_genes` is given; `TRUE`/`FALSE` forces it.

## Value

A ggplot object.

## Examples

``` r
plot_ibd_pairwise_group_heatmap(example_ibd_results())
```
