# Region-by-region IBD heatmap along the genome

Per-SNP IBD sharing between region pairs as tiles along the genome, one
facet per anchor region. The stored upper-triangle
(`region_a <= region_b`) is mirrored so every anchor shows all partner
regions.

## Usage

``` r
plot_ibd_region_heatmap(
  x,
  anchor = NULL,
  chroms = NULL,
  skip_chr = NULL,
  trans = "identity",
  colors = NULL,
  limits = NULL,
  fill_scale = NULL,
  highlight_genes = NULL,
  label_genes = FALSE
)
```

## Arguments

- x:

  An
  [IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
  object.

- anchor:

  Optional single region to show (one panel) instead of all.

- chroms:

  Optional chromosomes to keep (any spelling); others are dropped and
  the remaining ones re-laid-out contiguously.

- skip_chr:

  Optional chromosomes to drop (complement of `chroms`).

- trans:

  Fill-scale transform, e.g. `"identity"` (default), `"log2"`, `"sqrt"`.

- colors:

  Optional colour ramp for the fill (defaults to a single-hue
  light-to-dark sequential scale that stays readable when most values
  are near 0).

- limits:

  Optional `c(lo, hi)` fill limits; values outside are squished into
  range, so a few extremes near 1 don't crush the rest of the scale.

- fill_scale:

  Optional ggplot2 fill scale that fully overrides the above.

- highlight_genes:

  Optional gene names from the `genes` track to mark with lines.

- label_genes:

  Label the highlighted genes (top panel only, outside the plot).

## Value

A ggplot object.
