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
  freqbin = NULL,
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

  Width of the allele-frequency bins iHS is standardised within, or
  `NULL` (default) to pick one from `polarized`: **1** (a single bin)
  when unpolarized, `0.05` when polarized. The binning exists to control
  for *derived* allele frequency, and an unpolarized scan has no
  ancestral state – only `FREQ_MAJ`/`FREQ_MIN` – so major/minor is not
  derived/ancestral and binning by it controls nothing. rehh warns about
  this and about the resulting sparse bins above 0.5; a single bin
  silences both because it is the right answer, not because the warning
  was noise. Not cosmetic: on a real 249-sample cohort the two settings
  share only 25 of their top 50 \|iHS\| hits, and 4-12% of SNPs change
  sign (most in the 0.05-0.1 and \>0.3 MAF bands).

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
#> # A tibble: 24 × 7
#>    group    chr             pos snp_id            freq_minor     ihs neg_log10_p
#>    <fct>    <chr>         <dbl> <chr>                  <dbl>   <dbl>       <dbl>
#>  1 Cambodia Pf3D7_04_v3   92597 Pf3D7_04_v3:92597     0.1    -1.50        0.878 
#>  2 Cambodia Pf3D7_04_v3  544673 Pf3D7_04_v3:5446…     0.467  -0.0946      0.0340
#>  3 Cambodia Pf3D7_04_v3  898668 Pf3D7_04_v3:8986…     0.0667 -1.38        0.775 
#>  4 Cambodia Pf3D7_04_v3 1006055 Pf3D7_04_v3:1006…     0.133  -0.105       0.0381
#>  5 Cambodia Pf3D7_08_v3 1271193 Pf3D7_08_v3:1271…     0.233   0.279       0.108 
#>  6 Cambodia Pf3D7_08_v3 1356248 Pf3D7_08_v3:1356…     0.167   0.355       0.141 
#>  7 Cambodia Pf3D7_13_v3  173108 Pf3D7_13_v3:1731…     0.167   0.916       0.444 
#>  8 Cambodia Pf3D7_13_v3 1523439 Pf3D7_13_v3:1523…     0.1     0.711       0.321 
#>  9 Cambodia Pf3D7_13_v3 1872427 Pf3D7_13_v3:1872…     0.1     1.91        1.25  
#> 10 Cambodia Pf3D7_14_v3  524416 Pf3D7_14_v3:5244…     0.233  -0.420       0.171 
#> # ℹ 14 more rows
```
