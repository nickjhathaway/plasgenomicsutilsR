# Build phased haplotypes for a haplotype-homozygosity scan

Turns a genotype matrix into the complete, unambiguous 0/1 haplotypes
that
[`run_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_ihs.md)
and its cross-population relatives need.

## Usage

``` r
parasite_haplotypes(
  x,
  samples = NULL,
  fws = NULL,
  min_fws = IHS_MIN_FWS,
  het = c("sample", "missing"),
  maf = IHS_MIN_MAF,
  max_snp_missing = 0.1,
  max_sample_missing = 0.2,
  impute = TRUE,
  seed = 42,
  meta = NULL,
  genotype = NULL
)
```

## Arguments

- x:

  A
  [PopStructure](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PopStructure.md)
  or a genotype matrix (samples x SNPs, alt dosage 0/1/2, `NA` missing,
  `chr:pos` column names).

- samples:

  Restrict to these samples before anything else.

- fws:

  Per-sample Fws: a named numeric vector, or a data frame with `sample`
  and `fws` columns (e.g. read from `plasgenomicsutils calculate_fws`).
  `NULL` skips the monoclonal gate.

- min_fws:

  Fws floor for the gate (default 0.95).

- het:

  `"sample"` draws the allele at the population frequency; `"missing"`
  leaves the call to imputation.

- maf:

  Minor-allele frequency floor (default 0.03).

- max_snp_missing, max_sample_missing:

  Missingness ceilings for SNPs and samples.

- impute:

  Fill remaining gaps by drawing at the SNP's allele frequency. `FALSE`
  instead drops every SNP that still has a gap.

- seed:

  Random seed.

- meta, genotype:

  As in
  [`pop_diversity()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diversity.md).

## Value

A `parasite_haplotypes` object: `hap` (samples x SNPs, 0/1), `map`
(`chr`, `pos`, `snp_id`), `meta`, and a `filtering` record of what was
dropped.

## Details

The steps, in order, each reported in the result:

1.  **Monoclonal gate.** With `fws` supplied, samples below `min_fws`
    are dropped: a polyclonal infection is a mixture of haplotypes, not
    one haplotype.

2.  **Mixed calls.** Whatever heterozygous calls remain are either
    resolved by drawing the allele at the population frequency
    (`het = "sample"`, the usual choice for *Plasmodium*) or set to
    missing.

3.  **Filtering.** SNPs below `maf` or above `max_snp_missing`, then
    samples above `max_sample_missing`.

4.  **Imputation.** Remaining gaps are filled by drawing at the SNP's
    allele frequency, so the matrix is complete.

Both the allele draw and the imputation are random; `seed` makes a run
reproducible. Because the result depends on that draw, a signal worth
reporting should survive repeating the whole thing under a different
seed.

## See also

[`run_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_ihs.md),
[`run_rsb()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_rsb.md),
[`run_xpehh()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_xpehh.md)

## Examples

``` r
ps <- example_pop_structure(umap = FALSE)
parasite_haplotypes(ps, maf = 0.05)
#> <parasite_haplotypes> 60 haplotypes x 27 SNPs
#>   from            : 60 samples x 49 SNPs
#>   mixed calls     : 79 resolved by allele draw 
#>   SNPs dropped    : 12 missing, 10 MAF
#>   samples dropped : 0 missing
#>   imputed calls   : 16  seed: 42 
```
