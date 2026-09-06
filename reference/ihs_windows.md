# Windowed summary of an iHS scan

The fraction of SNPs in each window whose `abs(ihs)` exceeds `threshold`
– the summary iHS is normally read through, rather than SNP by SNP.

## Usage

``` r
ihs_windows(
  scan,
  window = 50000,
  step = NULL,
  threshold = 2,
  min_snps = 10,
  metric = "ihs"
)
```

## Arguments

- scan:

  A
  [`run_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_ihs.md)
  result, or any table with `chr`, `pos`, the `metric` column and
  optionally `group`.

- window:

  Window width in base pairs.

- step:

  Distance between window starts; defaults to `window` (windows that
  tile without overlapping). A smaller `step` slides the window and
  smooths the track, at the cost of neighbouring windows sharing SNPs.

- threshold:

  The `abs(metric)` a SNP must exceed to be counted as extreme. The
  conventional 2 for iHS.

- min_snps:

  Drop windows holding fewer scored SNPs than this. A window with three
  SNPs can read 100% and mean nothing; this is what keeps the sparse
  edges of the data from becoming the tallest peaks.

- metric:

  Column to summarise (default `"ihs"`); its magnitude is used, so an
  unpolarized scan needs no other handling.

## Value

A tibble with `group` (when the scan has one), `chr`, `start`, `end`,
`pos` (the window midpoint), `n_snps`, `n_extreme`, `frac_extreme` and
`max_abs`.

## Details

Per-SNP iHS is a high-variance statistic: each score rests on one focal
SNP's EHH decay, so neighbouring SNPs in the same haplotype disagree
freely and a genome-wide plot of them looks like grass. A sweep does not
raise one SNP, it raises a *run* of them, and counting the run is both
what the original description of the statistic proposed and what makes
the result comparable to any other windowed track – an IBD fraction,
say.

It also sidesteps the per-SNP p-value, which rehh computes from a normal
approximation and does not correct for multiple testing: a nominal
`p < 0.01` line is crossed by 1% of a neutral genome, so on a
whole-genome scan it marks a tail, not a finding.

## References

Voight, B. F., Kudaravalli, S., Wen, X. & Pritchard, J. K. (2006) A map
of recent positive selection in the human genome. *PLoS Biology* 4, e72.
[doi:10.1371/journal.pbio.0040072](https://doi.org/10.1371/journal.pbio.0040072)

## See also

[`run_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_ihs.md),
[`plot_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ihs.md),
[`plot_ibd_tugofwar()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_tugofwar.md)

## Examples

``` r
ps <- example_pop_structure(umap = FALSE)
hap <- parasite_haplotypes(ps, maf = 0.05)
ihs_windows(run_ihs(hap, group = "country"), window = 1e5, min_snps = 3)
#> # A tibble: 8 × 9
#>   group    chr         start    end    pos n_snps n_extreme frac_extreme max_abs
#>   <fct>    <chr>       <dbl>  <dbl>  <dbl>  <int>     <dbl>        <dbl>   <dbl>
#> 1 Cambodia Pf3D7_07_v3 3  e5 4   e5 3.5 e5      4         2       0.5      3.09 
#> 2 Ghana    Pf3D7_07_v3 3  e5 4   e5 3.5 e5     48         8       0.167    3.38 
#> 3 Cambodia Pf3D7_07_v3 4  e5 5   e5 4.35e5     30         4       0.133    2.59 
#> 4 Ghana    Pf3D7_07_v3 4  e5 5   e5 4.35e5    154        10       0.0649   3.22 
#> 5 Cambodia Pf3D7_08_v3 5  e5 6   e5 5.50e5     19         0       0        1.39 
#> 6 Ghana    Pf3D7_08_v3 5  e5 6   e5 5.50e5     47         1       0.0213   2.35 
#> 7 Ghana    Pf3D7_13_v3 1.6e6 1.7 e6 1.65e6      7         0       0        0.604
#> 8 Ghana    Pf3D7_13_v3 1.7e6 1.80e6 1.75e6     28         4       0.143    3.33 
```
