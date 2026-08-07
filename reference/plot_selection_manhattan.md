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
  highlight_genes = NULL,
  label_genes = NULL,
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
  SNPs. On real data it lands several times higher than Bonferroni's.
  `lambda_gc` says how far that reference is from fitting; the further
  from 1, the less the parametric lines mean. Naming a kind the run did
  not write is an error; `"all"` draws whichever kinds the threshold
  table carries.

- point_size, point_alpha:

  Point aesthetics.

- colours:

  Optional length-2 colour vector for the alternating chromosome bands.

## Value

A ggplot object.
