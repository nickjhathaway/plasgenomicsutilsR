# IBD selection-statistic Manhattan plot

The selection statistic along the genome, with the Bonferroni threshold
drawn as a dashed line when plotting `neg_log10_p` and thresholds are
available.

## Usage

``` r
plot_selection_manhattan(
  x,
  metric = c("neg_log10_p", "chi2_stat", "z_score"),
  groups = NULL,
  chroms = NULL,
  skip_chr = NULL,
  zoom = NULL,
  zoom_pad = 0.05,
  genes_for_track = NULL,
  gene_label_angle = 0,
  highlight_genes = NULL,
  label_genes = NULL,
  draw_threshold = TRUE,
  point_size = 0.5,
  point_alpha = 0.6,
  colours = NULL,
  colors = NULL
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

- groups:

  Optional character vector to keep (needs a `group` column).

- chroms:

  Optional chromosomes to keep (any spelling); others are dropped and
  the remaining ones re-laid-out contiguously.

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

  Optional gene names (from the object's `genes` track) to draw as
  reference lines; default all genes in the track. Requesting a name not
  in the track is an error.

- label_genes:

  Label the genes. `NULL` (default) labels them only when
  `highlight_genes` is given; `TRUE`/`FALSE` forces it.

- draw_threshold:

  Which significance line(s) to draw (only for `neg_log10_p`). The two
  read off a chi-square(1) are warm, the two built by permutation cool:
  `TRUE` / `"bonferroni"` (red dashed), `"fdr"` (orange dot-dash),
  `"permutation"` (family-wise, green solid), `"empirical"` (FDR over
  the permutation's own p-values, blue long-dash). Also `"both"` for the
  first two, `"all"` for every kind present, or `FALSE` for none. Prefer
  a permutation line where it disagrees with a parametric one: it is
  built by re-shuffling the IBD segments themselves, so it needs no
  chi-square(1) reference and it accounts for one segment spanning many
  SNPs, and it can land well above Bonferroni's. `lambda_gc` says how
  far that reference is from fitting; the further from 1, the less the
  parametric lines mean. Naming a kind the run did not write is an
  error; `"all"` draws whichever kinds the threshold table carries.

- point_size, point_alpha:

  Point aesthetics.

- colours, colors:

  Optional length-2 colour vector for the alternating chromosome bands.

## Value

A ggplot object.

## Examples

``` r
plot_selection_manhattan(example_ibd_results(), genes = PF_EXAMPLE_DRUG_GENES)
```
