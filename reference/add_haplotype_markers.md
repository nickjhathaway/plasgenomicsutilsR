# Add markers to a haplotype set from a callset of their own

Splices one or more markers into a
[`parasite_haplotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/parasite_haplotypes.md)
object, reading them as **allele indices** so a site with more than two
alleles keeps them apart.

## Usage

``` r
add_haplotype_markers(
  hap,
  vcf,
  het = c("missing", "draw"),
  impute = TRUE,
  seed = 42
)
```

## Arguments

- hap:

  A
  [`parasite_haplotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/parasite_haplotypes.md)
  object.

- vcf:

  A VCF/BCF holding the marker(s) to add. Needs `bcftools` on `PATH`.

- het:

  What to do with a call carrying two different alleles, which cannot be
  one haplotype: `"missing"` (default) sets it aside for the object's
  imputation, `"draw"` picks between the alleles that sample carries,
  weighted by how common each is among the unambiguous samples.

- impute:

  Fill the resulting gaps by drawing at the marker's allele frequencies
  (default `TRUE`). rehh cannot read a missing call, so a marker left
  with gaps would drop those haplotypes from every curve.

- seed:

  Seed for the draws.

## Value

The `parasite_haplotypes` object with the markers added.

## Details

This exists because the dosage route cannot carry the information. A
genotype matrix counts copies of one allele, so at a triallelic site
every non-reference call collapses to the same number no matter which
alternate it carries –
[`load_genotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/load_genotypes.md)
says as much when asked for `variants = "all"`. A marker that was
dropped from the main callset *for being multiallelic*, then called
separately to get it back, therefore has to enter at the haplotype
layer, where an allele is an index rather than a count.

Samples are matched by name, not position. The marker is inserted in
coordinate order, so the haplotypes still read along the genome, and a
position already present is an error rather than a silent replacement.

## See also

[`parasite_haplotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/parasite_haplotypes.md),
[`plot_ehh()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ehh.md)

## Examples

``` r
if (FALSE) { # \dontrun{
hap <- parasite_haplotypes(ps, fws = fws)
hap <- add_haplotype_markers(hap, "PF3D7_1343700.1-AA469.bcf", het = "draw")
plot_ehh(hap, "Pf3D7_13_v3:1725591")      # three curves, one per allele
} # }
```
