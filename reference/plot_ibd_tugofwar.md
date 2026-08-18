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
  top = NULL,
  top_label = NULL,
  metric = "neg_log10_p",
  scan_abs = NULL,
  scale = c("common", "free"),
  chroms = NULL,
  skip_chr = NULL,
  zoom = NULL,
  zoom_pad = 0.05,
  genes_for_track = NULL,
  gene_label_angle = 0,
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
  object (needs a per_snp_group table, and a selection table unless
  `top` supplies the upper track).

- group:

  Group(s) to plot. `NULL` (default) plots every group, faceted one per
  row; a single group gives one panel; a vector facets those groups.

- top:

  Optional per-SNP scan to hang from the top instead of the object's own
  IBD selection statistic: a
  [`run_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_ihs.md)
  result, a
  [`beta_score()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/beta_score.md)
  table, or any table with chr, pos, the `metric` column and (to face
  the IBD groups) a matching `group` column. Two-population scans
  ([`run_rsb()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_rsb.md),
  [`run_xpehh()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_xpehh.md))
  are indexed by population pair rather than group and cannot be
  mirrored against per-group IBD.

- top_label:

  Name for the upper track on the shared axis; defaults to `"selection"`
  for the object's statistic and `"scan"` for a `top` table.

- metric:

  Selection metric for the top track (default `"neg_log10_p"`).

- scan_abs:

  Plot the magnitude of a signed `metric` (iHS, Rsb, a z-score). `NULL`
  (default) does so whenever the metric has negative values, labelling
  the axis `|metric|`. The mirror cannot show the sign – its top half
  only spans zero to the maximum – so `FALSE` on signed values is an
  error rather than a silently clipped plot; use
  [`plot_ibd_locus()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_locus.md)
  when the sign matters.

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

- zoom:

  Optional single interval to crop to, keeping the same data and the
  same coordinates as the genome-wide plot: a chromosome (`"7"`), a
  range (`"7:728,081-988,719"`), a gene name from the object's track, or
  a one-row data frame with chr/start/end. Every gene in the window is
  drawn and named unless `label_genes = FALSE`.

- zoom_pad:

  Context to add around `zoom`, clamped to the chromosome. One value
  pads both sides (default 5%); two pad the left and the right, either
  in that order or named – `c(left = 5000, right = 40000)`, and naming
  only one side pads only that side. Each side is read on its own: below
  1 it is a fraction of the interval's span, at or above 1 it is base
  pairs, so `c(0.1, 20000)` is legal.

- genes_for_track:

  Optional gene table for the track drawn beneath a zoomed plot (e.g.
  [PF3D7_GENES](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF3D7_GENES.md)),
  so every gene in the window is shown and named while the plot's own
  short track still supplies the marked positions inside the panel.
  Without it the track and the marks come from the same genes, which
  means marking a whole annotation just to see the neighbours.

- gene_label_angle:

  Rotation for the gene names in that track, in degrees. `0` (default)
  centres each name under its gene; `45` or `90` runs it down to the
  left, which is what keeps long systematic ids from colliding over a
  dense annotation.

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

## Examples

``` r
# IBD sharing above, the selection statistic mirrored below
plot_ibd_tugofwar(example_ibd_results())
```
