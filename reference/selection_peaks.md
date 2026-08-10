# Merge a selection scan's significant SNPs into peaks

Turns a per-SNP scan into one row per implicated locus. A sweep leaves a
run of neighbouring significant SNPs, so the SNP count answers "how
dense is my panel here" while the peak count answers "how many things
did I find" – and only the second belongs in a result.

## Usage

``` r
selection_peaks(
  x,
  criterion = c("bonferroni", "fdr", "permutation", "empirical", "top", "value"),
  metric = NULL,
  cutoff = NULL,
  top = 0.01,
  gap = PEAK_GAP_BP,
  min_snps = 1L,
  pad = 0,
  genes = NULL,
  thresholds = NULL
)
```

## Arguments

- x:

  An
  [IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md),
  or a scan data frame with `chr`, `pos` and a value column.

- criterion:

  How a SNP qualifies:

  - `"bonferroni"` (default) / `"fdr"` — the scan's own `significant` /
    `significant_fdr` flag, or the matching stored threshold line.

  - `"permutation"` — the `significant_perm` flag from
    `ibd_selection_statistic --permute`: family-wise control against a
    null drawn from the data rather than a chi-square(1). Prefer it when
    `lambda_gc` is far from 1.

  - `"empirical"` — the `significant_fdr_perm` flag from the same run:
    Benjamini-Hochberg over the permutation's own per-SNP p-values. The
    FDR counterpart of `"permutation"`, and the only FDR here resting on
    calibrated p-values.

  - `"top"` — the top `top` fraction within each group, the convention
    for scans whose null is not trustworthy.

  - `"value"` — at or above `cutoff`.

- metric:

  Column to rank by; guessed from the scan when `NULL` (`neg_log10_p`,
  then `beta`, `value`, `ihs`).

- cutoff:

  Value cutoff for `criterion = "value"`.

- top:

  Tail fraction for `criterion = "top"` (default 0.01).

- gap:

  Merge hits separated by at most this many bp (default 20000).

- min_snps:

  Drop peaks supported by fewer than this many significant SNPs. Raising
  it to 2 or 3 is the cheapest way to suppress isolated single-SNP hits.

- pad:

  Widen each peak by this many bp on both sides before reporting and
  before matching genes.

- genes:

  Optional gene table (`name`, `chr` or `chrom`, `start`, `end`) to
  annotate each peak with the genes it overlaps;
  [PF3D7_GENES](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF3D7_GENES.md)
  is a sensible argument.

- thresholds:

  Threshold table, when `x` is a bare data frame carrying no flag
  column.

## Value

A tibble with one row per peak: the grouping column, `chr`, `start`,
`end`, `width`, `n_snps` (significant SNPs in the peak), `peak_pos` (the
single highest-scoring SNP), `peak_value`, `mean_value`, and – when a
gene table was given – five annotation columns. Sorted by `peak_value`.

Everything is anchored on `peak_pos`, not the interval's midpoint: the
midpoint is an artefact of where merging started and stopped, and can
land in a gap between genes.

- `peak_interval_genes` — a **list-column** of every gene the interval
  spans, in genomic order. A peak of a few hundred kb, as real IBD
  sharing regions are, can cross dozens, so this is the set to take into
  a follow-up rather than something to read off the screen. Unnest it
  for one row per peak x gene: `tidyr::unnest(pk, peak_interval_genes)`;
  [`lengths()`](https://rdrr.io/r/base/lengths.html) of it is `n_genes`.

- `gene_at_peak` — the gene(s) whose span **contains** the peak SNP,
  comma-separated. Usually one, since gene spans rarely overlap, and
  **empty when the peak SNP is intergenic** — about a third of peaks on
  a real cohort.

- `nearest_gene`, `distance_to_gene` — of the genes the peak **covers**,
  the closest one to the peak SNP and the gap to it in bp, `0` when the
  SNP is inside it. This is what to read when `gene_at_peak` is empty;
  those intergenic peaks are typically within a kb or two of a gene in
  the same peak. Candidates are restricted to the interval, so a peak
  that covers no gene reports none rather than pointing at something far
  outside it — which matters when `genes` is a short list, where most
  peaks cover nothing from it.

- `n_genes` — how many genes the interval spans, i.e.
  `lengths(peak_interval_genes)`.

They agree by construction: `n_genes == 0` implies an empty
`peak_interval_genes` and empty name columns,
`n_genes == lengths(peak_interval_genes)`, and a non-empty
`gene_at_peak` implies `distance_to_gene == 0`.

## Details

Works on any of the per-SNP scans in the package: an
[IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
(its selection table), or the tibble from
[`run_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_ihs.md),
[`run_rsb()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_rsb.md),
[`run_xpehh()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_xpehh.md)
or
[`beta_score()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/beta_score.md).

## See also

[`ihs_genes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ihs_genes.md)
and
[`beta_genes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/beta_genes.md)
for the per-gene view; this is the per-locus one, which does not need a
gene to exist where the signal is.

## Examples

``` r
ps <- example_pop_structure(umap = FALSE)
b <- beta_score(ps, group = "country", window = 300000, min_window_snps = 1)
selection_peaks(b, criterion = "top", top = 0.2, genes = PF_EXAMPLE_DRUG_GENES)
#> # A tibble: 5 × 14
#>   group    chr          start    end width n_snps peak_pos peak_value mean_value
#>   <chr>    <chr>        <dbl>  <dbl> <dbl>  <int>    <dbl>      <dbl>      <dbl>
#> 1 Cambodia Pf3D7_01_v3 5.31e5 5.31e5     1      1   531105      0.536      0.536
#> 2 Cambodia Pf3D7_01_v3 5.77e5 5.77e5     1      1   577066      0.484      0.484
#> 3 Ghana    Pf3D7_07_v3 5.99e5 5.99e5     1      1   598578      0.275      0.275
#> 4 Ghana    Pf3D7_09_v3 1.48e6 1.48e6     1      1  1475529      0.169      0.169
#> 5 Ghana    Pf3D7_07_v3 5.37e5 5.37e5     1      1   536889      0.142      0.142
#> # ℹ 5 more variables: peak_interval_genes <list>, gene_at_peak <chr>,
#> #   nearest_gene <chr>, distance_to_gene <dbl>, n_genes <int>
```
