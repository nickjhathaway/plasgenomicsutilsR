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

## Allele count

not a reduction, and not a problem here: `scan_hh()` reports only two
allele frequencies whatever a marker carries, which looks like a
reduction and is not one for this statistic. Rsb's input is `iNES`, a
**site-level** homozygosity computed over every haplotype: at a
four-allele marker `scan_hh()`'s `iNES` is identical to
[`rehh::calc_ehhs()`](https://rdrr.io/pkg/rehh/man/calc_ehhs.html)'s
over the full data. No haplotype is dropped and no allele is ignored;
only the reported `FREQ_` columns reduce, and this function does not
return them.

`iNES` normalises by the focal site's own homozygosity, so splitting the
same haplotypes into more allele classes moves it by 1% or less. That is
why Rsb carries no multiallelic warning while
[`run_xpehh()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_xpehh.md)
does – see there.

iES falls as the focal splits, and XP-EHH divides two of them: No
haplotype is dropped and no allele ignored – `scan_hh()`'s `iES` at a
four-allele marker is identical to
[`rehh::calc_ehhs()`](https://rdrr.io/pkg/rehh/man/calc_ehhs.html)'s
over the full data, and only the reported `FREQ_` columns reduce. But
`iES` is a site-level homozygosity **pooled across allele classes**, so
it falls as the focal is split into more of them. Measured on identical
haplotypes:

|                      |              |              |
|----------------------|--------------|--------------|
| alleles at the focal | `iES`        | `iNES`       |
| 2                    | 991.8        | 2406.7       |
| 3                    | 631.5 (-36%) | 2388.6 (-1%) |
| 4                    | 558.5 (-44%) | 2395.8 (0%)  |

XP-EHH divides one population's `iES` by the other's, so a marker
carrying **different** numbers of alleles in the two populations puts
them on different footings, and the score partly measures the diversity
difference rather than the haplotype-length difference. It warns when
that happens, and reports `n_alleles_pop1` / `n_alleles_pop2` so those
markers can be dropped.
[`run_rsb()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_rsb.md)
reads `iNES` instead and is not affected.

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
