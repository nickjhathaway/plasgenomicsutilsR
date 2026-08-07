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
#> # A tibble: 16 × 8
#>    group    chr             pos snp_id       freq n_called n_window_snps    beta
#>    <fct>    <chr>         <dbl> <chr>       <dbl>    <dbl>         <int>   <dbl>
#>  1 Cambodia Pf3D7_01_v3  531105 Pf3D7_01_… 0.6          25             1  0.536 
#>  2 Cambodia Pf3D7_01_v3  577066 Pf3D7_01_… 0.36         25             1  0.484 
#>  3 Cambodia Pf3D7_07_v3  536889 Pf3D7_07_… 0.0476       21             1 -0.277 
#>  4 Cambodia Pf3D7_07_v3  598578 Pf3D7_07_… 0.524        21             1 -0.266 
#>  5 Cambodia Pf3D7_08_v3 1271193 Pf3D7_08_… 0.767        30             1  0.228 
#>  6 Cambodia Pf3D7_08_v3 1356248 Pf3D7_08_… 0.833        30             1  0.0928
#>  7 Ghana    Pf3D7_04_v3 1006055 Pf3D7_04_… 0.467        30             1 -0.158 
#>  8 Ghana    Pf3D7_04_v3 1071662 Pf3D7_04_… 0.133        30             1 -0.249 
#>  9 Ghana    Pf3D7_07_v3   80466 Pf3D7_07_… 0.0333       30             1 -0.0749
#> 10 Ghana    Pf3D7_07_v3  161320 Pf3D7_07_… 0.2          30             1 -0.123 
#> 11 Ghana    Pf3D7_07_v3  536889 Pf3D7_07_… 0.2          20             1  0.142 
#> 12 Ghana    Pf3D7_07_v3  598578 Pf3D7_07_… 0.737        19             1  0.275 
#> 13 Ghana    Pf3D7_08_v3 1356248 Pf3D7_08_… 0.533        30             1 -0.242 
#> 14 Ghana    Pf3D7_08_v3 1388102 Pf3D7_08_… 0.0435       23             1 -0.268 
#> 15 Ghana    Pf3D7_09_v3 1416040 Pf3D7_09_… 0.0333       30             1  0.137 
#> 16 Ghana    Pf3D7_09_v3 1475529 Pf3D7_09_… 0.0556       18             1  0.169 
```
