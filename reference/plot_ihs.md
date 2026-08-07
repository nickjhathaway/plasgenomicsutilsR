# Manhattan plot of a haplotype-homozygosity scan

Draws `neg_log10_p` (or the raw statistic) along the genome, one panel
per group or population pair, over the chromosome bands and optional
gene markers used by the other genome-wide plots in the package.

## Usage

``` r
plot_ihs(
  scan,
  metric = c("neg_log10_p", "ihs", "value"),
  threshold = NULL,
  genes = NULL,
  highlight_genes = NULL,
  label_genes = NULL,
  chroms = NULL,
  skip_chr = NULL,
  reference = DEFAULT_REFERENCE,
  point_size = 0.6,
  point_alpha = 0.7,
  colours = NULL
)
```

## Arguments

- scan:

  The tibble from
  [`run_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_ihs.md),
  [`run_rsb()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_rsb.md)
  or
  [`run_xpehh()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_xpehh.md).

- metric:

  `"neg_log10_p"` (default) or the statistic itself (`"ihs"` /
  `"value"`).

- threshold:

  Draw a dashed significance line at this height. `NULL` uses the
  `-log10(p)` of a 1% two-sided tail when plotting `neg_log10_p`, the
  convention these scans are usually read at; `NA` draws none.

- genes:

  Gene table to mark (e.g.
  [PF_EXAMPLE_DRUG_GENES](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF_EXAMPLE_DRUG_GENES.md));
  `NULL` for none.

- highlight_genes, label_genes:

  Which genes to mark and whether to name them.

- chroms, skip_chr:

  Chromosomes to keep or drop.

- reference:

  Reference id for the chromosome layout.

- point_size, point_alpha, colours:

  Point and band aesthetics.

## Value

A ggplot object.

## Examples

``` r
ps <- example_pop_structure(umap = FALSE)
hap <- parasite_haplotypes(ps, maf = 0.05)
plot_ihs(run_ihs(hap, group = "country"), genes = PF_EXAMPLE_DRUG_GENES)
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
```
