# LD-prune a VCF and return the genotype matrix

Converts a VCF to GDS (only when needed), LD-prunes, and returns the
pruned genotype matrix (samples x SNPs, coded 0/1/2, `NA` for missing)
via SNPRelate.

## Usage

``` r
run_ld_prune(
  vcf,
  gds = NULL,
  prune = TRUE,
  ld_threshold = 0.2,
  slide_max_bp = 20000,
  slide_max_n = 200,
  autosome_only = FALSE,
  maf = NaN,
  missing_rate = NaN,
  seed = 42,
  vcf_dir = NULL,
  allele = c("alt", "ref")
)
```

## Arguments

- vcf:

  Path to a (bgzipped) VCF, or a **BCF** – SNPRelate reads VCF text
  only, so a BCF is converted first with `bcftools`, reusing any VCF
  already sitting next to it rather than making another copy.

- gds:

  Optional GDS path; derived from `vcf` if `NULL`.

- prune:

  LD-prune (default `TRUE`). `FALSE` returns **every** biallelic SNP,
  unpruned – use this for the genotype matrix fed to
  [`pop_diff()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff.md)
  /
  [`pop_diff_table()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff_table.md),
  since LD-pruning removes the very SNPs that carry the differentiation
  signal.

- ld_threshold, slide_max_bp, slide_max_n, autosome_only:

  Passed to
  [`SNPRelate::snpgdsLDpruning()`](https://rdrr.io/pkg/SNPRelate/man/snpgdsLDpruning.html)
  (defaults 0.2 / 20000 / 200 / `FALSE`); ignored when `prune = FALSE`.

- maf, missing_rate:

  Optional MAF / per-SNP missing-rate cutoffs for pruning.

- seed:

  Random seed for the pruning.

- vcf_dir:

  Where to put the VCF converted from a BCF (default: alongside the
  BCF). Point it somewhere scratch to keep converted copies out of the
  data directory.

- allele:

  Which allele the returned dosage counts. SNPRelate counts the
  **reference** allele; the default `"alt"` flips that so the matrix
  means what the rest of the package says it means. Only reported allele
  frequencies (and the arbitrary sign of a PCA axis) depend on this –
  every diversity, differentiation, LD and selection statistic here is
  symmetric in `p` and `1 - p`.

## Value

A list with `genotype` (matrix; sample row names and `chr:pos` column
names), `sample.id`, and `snp.id`.
