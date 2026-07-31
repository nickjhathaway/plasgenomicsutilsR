# IBD "triangle" panels for genes or specific SNPs

Draws a region-by-region IBD-sharing triangle per feature: a **gene**
(per-SNP pairwise IBD for SNPs falling strictly inside the gene
interval, aggregated) or a **specific SNP** (the sharing at that single
position). One facet per feature – use it to ask whether a gene or locus
is itself shared between regions.

## Usage

``` r
plot_drug_gene_triangles(
  x,
  genes = NULL,
  snps = NULL,
  agg = c("mean", "median", "max"),
  individual = FALSE,
  label = TRUE,
  digits = 2,
  ncol = NULL,
  trans = "identity",
  colors = NULL,
  limits = NULL,
  fill_scale = NULL
)
```

## Arguments

- x:

  An
  [IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
  object (needs `pairwise_region`; a `genes` track for `genes`).

- genes:

  Gene names to include from the track (default: all in the track;
  case-insensitive). Ignored if only `snps` is given.

- snps:

  Specific SNPs to draw, as `"chr:pos"` ids (e.g.
  `"Pf3D7_07_v3:403222"`) or a data frame with `chr`, `pos` (and
  optional `name`). Each is one facet.

- agg:

  How a gene's value in each region-pair cell is computed from **all**
  its in-gene SNPs (no SNP is picked or dropped): `"mean"` (default),
  `"median"`, or `"max"` (the strongest sharing at that pair). A
  single-SNP feature is unaffected.

- individual:

  If `TRUE`, return a **named list** of one-triangle plots (one per
  feature, e.g. to write a multi-page PDF) with the legend tucked into
  the empty upper triangle; if `FALSE` (default), one faceted grid plot.

- label:

  Draw the value in each tile.

- digits:

  Decimal places for the tile labels.

- ncol:

  Facet columns for the grid (default: ggplot2 chooses).

- trans:

  Fill-scale transform, e.g. `"identity"` (default), `"log2"`, `"sqrt"`.

- colors:

  Optional colour ramp for the fill (default: the pairwise-sharing
  ramp).

- limits:

  Optional `c(lo, hi)` fill limits (values outside are squished).

- fill_scale:

  Optional ggplot2 fill scale that fully overrides the above.

## Value

A ggplot object (grid), or a named list of ggplot objects when
`individual`.

## Details

Genes come from the `genes` track on the
[IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
object and membership is strict: a SNP belongs to a gene only when its
position is inside the interval (no flanking). Target individual loci
with `snps` when a specific variant matters more than a whole gene.
