# LD-prune a VCF and return the genotype matrix

Converts a VCF to GDS (only when needed), LD-prunes, and returns the
pruned genotype matrix (samples x SNPs, coded 0/1/2, `NA` for missing)
via SNPRelate.

## Usage

``` r
run_ld_prune(
  vcf,
  gds = NULL,
  ld_threshold = 0.2,
  slide_max_bp = 20000,
  slide_max_n = 200,
  autosome_only = FALSE,
  maf = NaN,
  missing_rate = NaN,
  seed = 42
)
```

## Arguments

- vcf:

  Path to a (bgzipped) VCF.

- gds:

  Optional GDS path; derived from `vcf` if `NULL`.

- ld_threshold, slide_max_bp, slide_max_n, autosome_only:

  Passed to
  [`SNPRelate::snpgdsLDpruning()`](https://rdrr.io/pkg/SNPRelate/man/snpgdsLDpruning.html)
  (defaults 0.2 / 20000 / 200 / `FALSE`).

- maf, missing_rate:

  Optional MAF / per-SNP missing-rate cutoffs for pruning.

- seed:

  Random seed for the pruning.

## Value

A list with `genotype` (matrix), `sample.id`, and `snp.id`.
