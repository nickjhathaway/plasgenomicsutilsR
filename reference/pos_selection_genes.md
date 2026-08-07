# Genes under positive selection

Intersects the significant SNPs of the IBD selection statistic (those at
or above the per-group Bonferroni threshold on `neg_log10_p`) with a
gene track, returning the genes hit by a selection signal. Because VCF
filtering can leave the peak SNP just outside a gene, a SNP counts for a
gene when it falls within `within` bp of the gene interval (default 2
kb), not only strictly inside it.

## Usage

``` r
pos_selection_genes(
  x,
  within = 2000,
  genes = NULL,
  groups = NULL,
  threshold = NULL
)
```

## Arguments

- x:

  An
  [IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
  object (needs a `selection` table; per-group thresholds come from the
  object, or pass `threshold`).

- within:

  Maximum distance in bp between a significant SNP and the gene interval
  for the SNP to count (default `2000`). `0` requires the SNP to be
  strictly inside.

- genes:

  Optional gene track to scan (path or data frame with `chr`/`chrom`,
  `start`, `end`, and optionally `name`, `gene_id`), overriding the
  object's track. Pass
  [PF3D7_GENES](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF3D7_GENES.md)
  here to scan every gene without attaching it to the object (which
  would make the plot functions draw a line per gene). Default: the
  object's `genes`.

- groups:

  Optional subset of groups to scan (needs a `group` column); default
  all.

- threshold:

  Optional threshold override on the `neg_log10_p` scale: a single
  number for every group, or a named vector (`group -> threshold`).
  Default uses the object's thresholds.

## Value

A tibble, one row per (group, gene) hit, sorted by group then descending
peak `neg_log10_p`: `group` (dropped for a group-less selection table),
`gene_id`, `name`, `chr`, `gene_start`, `gene_end`, `n_snps`
(significant SNPs in the window), `min_distance` (0 when a SNP is inside
the gene), `peak_pos`, and a `peak_<metric>` column for every metric in
the selection table (e.g. `peak_neg_log10_p`, `peak_maf`).

## Details

Thresholding is always on `neg_log10_p` (the only metric the Bonferroni
threshold is defined for); the returned row reports **every** metric
column of the selection table at that gene's peak SNP (`peak_*`), so you
still see `maf` / `z_score` / `chi2_stat`, etc.

## Examples

``` r
ibd <- example_ibd_results()
pos_selection_genes(ibd)                          # object's drug-gene track, 2 kb window
#> # A tibble: 0 × 14
#> # ℹ 14 variables: group <chr>, gene_id <chr>, name <chr>, chr <chr>,
#> #   gene_start <int>, gene_end <int>, n_snps <int>, min_distance <dbl>,
#> #   peak_pos <int>, peak_maf <dbl>, peak_z_score <dbl>, peak_chi2_stat <dbl>,
#> #   peak_neg_log10_p <dbl>, peak_significant <dbl>
pos_selection_genes(ibd, within = 0)              # strictly-inside only
#> # A tibble: 0 × 14
#> # ℹ 14 variables: group <chr>, gene_id <chr>, name <chr>, chr <chr>,
#> #   gene_start <int>, gene_end <int>, n_snps <int>, min_distance <dbl>,
#> #   peak_pos <int>, peak_maf <dbl>, peak_z_score <dbl>, peak_chi2_stat <dbl>,
#> #   peak_neg_log10_p <dbl>, peak_significant <dbl>
pos_selection_genes(ibd, genes = PF3D7_GENES)     # scan every gene, track left untouched
#> # A tibble: 56 × 14
#>    group gene_id    name  chr   gene_start gene_end n_snps min_distance peak_pos
#>    <chr> <chr>      <chr> <chr>      <int>    <int>  <int>        <dbl>    <dbl>
#>  1 DRC   PF3D7_070… pfhs… 7         381591   384614      1         1677   386290
#>  2 DRC   PF3D7_070… PF3D… 7         385582   388321      1            0   386290
#>  3 DRC   PF3D7_081… pfab… 8         521186   524009      1          898   524906
#>  4 DRC   PF3D7_081… pfpp… 8         525057   527604      1          151   524906
#>  5 DRC   PF3D7_070… pfme… 7         413559   421749      1            0   417012
#>  6 DRC   PF3D7_081… PF3D… 8         593648   595412      1          573   595984
#>  7 DRC   PF3D7_081… PF3D… 8         596189   600328      1          205   595984
#>  8 DRC   PF3D7_081… pfca… 8         567823   573148      1            0   568929
#>  9 DRC   PF3D7_062… PF3D… 6        1183895  1185617      1          493  1186109
#> 10 DRC   PF3D7_062… pfga… 6        1185995  1188644      1            0  1186109
#> # ℹ 46 more rows
#> # ℹ 5 more variables: peak_maf <dbl>, peak_z_score <dbl>, peak_chi2_stat <dbl>,
#> #   peak_neg_log10_p <dbl>, peak_significant <lgl>
```
