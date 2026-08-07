# Integrated haplotype score (iHS)

Scans each group for recent positive directional selection,
standardising the integrated EHH ratio within allele-frequency bins so
scores are comparable along the genome. Large `abs(ihs)` marks a SNP
whose haplotype background is unusually long for its frequency.

## Usage

``` r
run_ihs(
  hap,
  group = NULL,
  meta = NULL,
  polarized = FALSE,
  freqbin = 0.05,
  min_maf = 0.05,
  min_samples = 4,
  threads = 1
)
```

## Arguments

- hap:

  A
  [`parasite_haplotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/parasite_haplotypes.md)
  object.

- group:

  Metadata column naming the grouping, a vector aligned to the haplotype
  rows, or `NULL` to scan every sample as one population.

- meta:

  Metadata (defaults to the one carried by `hap`).

- polarized:

  Treat allele 1 as derived (needs a real ancestral state).

- freqbin:

  Width of the allele-frequency bins used for standardisation.

- min_maf:

  Minor-allele frequency floor applied at standardisation.

- min_samples:

  Skip groups smaller than this.

- threads:

  Threads for rehh.

## Value

A tibble with `group`, `chr`, `pos`, `snp_id`, `freq_minor`, `ihs` and
`neg_log10_p`.

## Details

Without an outgroup there is no ancestral state to polarise by, so
`polarized = FALSE` by default and the comparison is major versus minor
allele rather than ancestral versus derived. That is the standard
treatment for *P. falciparum* and it means the *sign* of `ihs` should
not be read as "selection on the derived allele" – use `abs(ihs)` and
`neg_log10_p`.

## References

Voight, B. F., Kudaravalli, S., Wen, X. & Pritchard, J. K. (2006) A map
of recent positive selection in the human genome. *PLoS Biology* 4, e72.
[doi:10.1371/journal.pbio.0040072](https://doi.org/10.1371/journal.pbio.0040072)

Gautier, M., Klassmann, A. & Vitalis, R. (2017) rehh 2.0: a
reimplementation of the R package rehh to detect positive selection from
haplotype structure. *Molecular Ecology Resources* 17, 78-90.
[doi:10.1111/1755-0998.12634](https://doi.org/10.1111/1755-0998.12634)

## See also

[`ihs_genes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ihs_genes.md),
[`plot_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ihs.md),
[`run_rsb()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_rsb.md),
[`beta_score()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/beta_score.md)

## Examples

``` r
ps <- example_pop_structure(umap = FALSE)
hap <- parasite_haplotypes(ps, maf = 0.05)
run_ihs(hap, group = "country")
#> Warning: If alleles are unpolarized, 'freqbin' should be set to 1 (one bin).
#> Warning: The number of markers with allele frequencies in bin [0.05,0.1) is less than 10: you should probably increase bin width.
#> Warning: The number of markers with allele frequencies in bin [0.1,0.15) is less than 10: you should probably increase bin width.
#> Warning: The number of markers with allele frequencies in bin [0.15,0.2) is less than 10: you should probably increase bin width.
#> Warning: The number of markers with allele frequencies in bin [0.2,0.25) is less than 10: you should probably increase bin width.
#> Warning: The number of markers with allele frequencies in bin [0.25,0.3) is less than 10: you should probably increase bin width.
#> Warning: The number of markers with allele frequencies in bin [0.3,0.35) is less than 10: you should probably increase bin width.
#> Warning: The number of markers with allele frequencies in bin [0.35,0.4) is less than 10: you should probably increase bin width.
#> Warning: The number of markers with allele frequencies in bin [0.4,0.45) is less than 10: you should probably increase bin width.
#> Warning: The number of markers with allele frequencies in bin [0.45,0.5) is less than 10: you should probably increase bin width.
#> Warning: The number of markers with allele frequencies in bin [0.5,0.55) is less than 10: you should probably increase bin width.
#> Warning: The number of markers with allele frequencies in bin [0.55,0.6) is less than 10: you should probably increase bin width.
#> Warning: The number of markers with allele frequencies in bin [0.6,0.65) is less than 10: you should probably increase bin width.
#> Warning: The number of markers with allele frequencies in bin [0.65,0.7) is less than 10: you should probably increase bin width.
#> Warning: The number of markers with allele frequencies in bin [0.7,0.75) is less than 10: you should probably increase bin width.
#> Warning: The number of markers with allele frequencies in bin [0.75,0.8) is less than 10: you should probably increase bin width.
#> Warning: The number of markers with allele frequencies in bin [0.8,0.85) is less than 10: you should probably increase bin width.
#> Warning: The number of markers with allele frequencies in bin [0.85,0.9) is less than 10: you should probably increase bin width.
#> Warning: The number of markers with allele frequencies in bin [0.9,0.95) is less than 10: you should probably increase bin width.
#> Warning: If alleles are unpolarized, 'freqbin' should be set to 1 (one bin).
#> Warning: The number of markers with allele frequencies in bin [0.05,0.1) is less than 10: you should probably increase bin width.
#> Warning: The number of markers with allele frequencies in bin [0.1,0.15) is less than 10: you should probably increase bin width.
#> Warning: The number of markers with allele frequencies in bin [0.15,0.2) is less than 10: you should probably increase bin width.
#> Warning: The number of markers with allele frequencies in bin [0.2,0.25) is less than 10: you should probably increase bin width.
#> Warning: The number of markers with allele frequencies in bin [0.25,0.3) is less than 10: you should probably increase bin width.
#> Warning: The number of markers with allele frequencies in bin [0.3,0.35) is less than 10: you should probably increase bin width.
#> Warning: The number of markers with allele frequencies in bin [0.35,0.4) is less than 10: you should probably increase bin width.
#> Warning: The number of markers with allele frequencies in bin [0.4,0.45) is less than 10: you should probably increase bin width.
#> Warning: The number of markers with allele frequencies in bin [0.45,0.5) is less than 10: you should probably increase bin width.
#> Warning: The number of markers with allele frequencies in bin [0.5,0.55) is less than 10: you should probably increase bin width.
#> Warning: The number of markers with allele frequencies in bin [0.55,0.6) is less than 10: you should probably increase bin width.
#> Warning: The number of markers with allele frequencies in bin [0.6,0.65) is less than 10: you should probably increase bin width.
#> Warning: The number of markers with allele frequencies in bin [0.65,0.7) is less than 10: you should probably increase bin width.
#> Warning: The number of markers with allele frequencies in bin [0.7,0.75) is less than 10: you should probably increase bin width.
#> Warning: The number of markers with allele frequencies in bin [0.75,0.8) is less than 10: you should probably increase bin width.
#> Warning: The number of markers with allele frequencies in bin [0.8,0.85) is less than 10: you should probably increase bin width.
#> Warning: The number of markers with allele frequencies in bin [0.85,0.9) is less than 10: you should probably increase bin width.
#> Warning: The number of markers with allele frequencies in bin [0.9,0.95) is less than 10: you should probably increase bin width.
#> # A tibble: 24 × 7
#>    group    chr             pos snp_id             freq_minor    ihs neg_log10_p
#>    <fct>    <chr>         <dbl> <chr>                   <dbl>  <dbl>       <dbl>
#>  1 Cambodia Pf3D7_04_v3   92597 Pf3D7_04_v3:92597      0.1    -1.20       0.639 
#>  2 Cambodia Pf3D7_04_v3  544673 Pf3D7_04_v3:544673     0.467  NA         NA     
#>  3 Cambodia Pf3D7_04_v3  898668 Pf3D7_04_v3:898668     0.0667 NA         NA     
#>  4 Cambodia Pf3D7_04_v3 1006055 Pf3D7_04_v3:10060…     0.133  -0.133      0.0486
#>  5 Cambodia Pf3D7_08_v3 1271193 Pf3D7_08_v3:12711…     0.233   0.707      0.319 
#>  6 Cambodia Pf3D7_08_v3 1356248 Pf3D7_08_v3:13562…     0.167  -0.707      0.319 
#>  7 Cambodia Pf3D7_13_v3  173108 Pf3D7_13_v3:173108     0.167   0.707      0.319 
#>  8 Cambodia Pf3D7_13_v3 1523439 Pf3D7_13_v3:15234…     0.1     0.491      0.205 
#>  9 Cambodia Pf3D7_13_v3 1872427 Pf3D7_13_v3:18724…     0.1     1.41       0.797 
#> 10 Cambodia Pf3D7_14_v3  524416 Pf3D7_14_v3:524416     0.233  -0.707      0.319 
#> # ℹ 14 more rows
```
