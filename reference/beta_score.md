# Beta: balancing selection from clustered allele frequencies

For each core SNP, compares an estimate of theta weighted towards
neighbours whose folded allele frequency matches the core's against
Watterson's theta from the same window. Large positive values mark
neighbourhoods where frequencies are piled up around an
intermediate-frequency variant – the footprint of long-term balancing
selection, and in *P. falciparum* typically an antigen rather than a
drug target.

## Usage

``` r
beta_score(
  x,
  group = NULL,
  window = BETA_WINDOW,
  p = BETA_P,
  min_freq = 0,
  min_window_snps = 2,
  het = c("missing", "dosage"),
  min_samples = 4,
  meta = NULL,
  genotype = NULL
)
```

## Arguments

- x:

  A
  [PopStructure](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PopStructure.md)
  or a genotype matrix (samples x SNPs, alt dosage, `chr:pos` column
  names).

- group:

  Metadata column name, or a vector aligned to the rows; `NULL` pools
  every sample.

- window:

  Total window width in bp; neighbours are taken within `window / 2` on
  each side (default 1000).

- p:

  Sharpness of the frequency-similarity weighting (default 2).

- min_freq:

  Skip core SNPs whose folded frequency is at or below this.

- min_window_snps:

  Skip core SNPs with fewer neighbours than this in the window.

- het:

  How a heterozygous call is read; see
  [`pop_diversity()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diversity.md).

- min_samples:

  Skip groups smaller than this.

- meta, genotype:

  As in
  [`pop_diversity()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diversity.md).

## Value

A tibble with `group`, `chr`, `pos`, `snp_id`, `freq`, `n_called`,
`n_window_snps` and `beta`, one row per core SNP per group.

**Reading `beta`.** This is Siewert & Voight's folded Beta1: how much
more the SNPs around this one share its allele frequency than drift
alone would produce. Long-term balancing selection holds a variant at
intermediate frequency for many generations, and its neighbours
hitch-hike to *similar* frequencies – so a cluster of matched
frequencies is the signature, which is what beta measures.

The scale is not a z-score and has no p-value here: it depends on the
window, the SNP density and the sample size, so a value meaningful in
one scan is not comparable to another. **Rank rather than threshold.**
Larger is more evidence of balancing selection; around zero is what
neutrality gives; negative means the neighbourhood is *less*
frequency-matched than expected, which is a directional-sweep pattern,
not balancing selection – for that read
[`run_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_ihs.md)
instead.

In practice take the top tail –
`selection_peaks(b, criterion = "top", top = 0.01, metric = "beta")` –
and annotate it with
[`annotate_snps()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/annotate_snps.md).
In *P. falciparum* the expected hits are the antigens under long-term
frequency-dependent selection from host immunity – on a real 249-sample
cohort the top 1% is led by *pfama1* and *pfdblmsp*, which is the
positive control to look for; a scan that surfaces none of them is a
reason to check the input before believing anything else in it. Note the
classic examples *var* and *rifin* will **not** appear: they are
subtelomeric, so a core-genome filter removes them before the scan ever
sees them.

## References

Siewert, K. M. & Voight, B. F. (2017) Detecting long-term balancing
selection using allele frequency correlation. *Molecular Biology and
Evolution* 34, 2996-3005.
[doi:10.1093/molbev/msx209](https://doi.org/10.1093/molbev/msx209)

## See also

[`beta_genes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/beta_genes.md)
to summarise per gene,
[`run_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_ihs.md)
for directional selection.

## Examples

``` r
ps <- example_pop_structure(umap = FALSE)
beta_score(ps, group = "country", window = 200000, min_window_snps = 1)
#> # A tibble: 661 × 8
#>    group    chr            pos snp_id        freq n_called n_window_snps    beta
#>    <fct>    <chr>        <dbl> <chr>        <dbl>    <dbl>         <int>   <dbl>
#>  1 Cambodia Pf3D7_01_v3 531104 Pf3D7_01_v… 0.6          25             1   0.536
#>  2 Cambodia Pf3D7_01_v3 577065 Pf3D7_01_v… 0.36         25             1   0.484
#>  3 Cambodia Pf3D7_07_v3 375090 Pf3D7_07_v… 0.0714       28            95   2.05 
#>  4 Cambodia Pf3D7_07_v3 375130 Pf3D7_07_v… 0.615        13            95 -12.7  
#>  5 Cambodia Pf3D7_07_v3 375470 Pf3D7_07_v… 0.933        30            95   2.18 
#>  6 Cambodia Pf3D7_07_v3 375817 Pf3D7_07_v… 0.933        30            95   2.18 
#>  7 Cambodia Pf3D7_07_v3 376527 Pf3D7_07_v… 0.857        14            95  -3.38 
#>  8 Cambodia Pf3D7_07_v3 377476 Pf3D7_07_v… 0.967        30            95   3.53 
#>  9 Cambodia Pf3D7_07_v3 378903 Pf3D7_07_v… 0.0333       30            95   3.53 
#> 10 Cambodia Pf3D7_07_v3 379066 Pf3D7_07_v… 0.05         20            95   4.47 
#> # ℹ 651 more rows
```
