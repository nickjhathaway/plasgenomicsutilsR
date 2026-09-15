# What each sample carries at one variant

The set of alleles each sample's genotype names at one position, as a
named character vector – exactly the `sample -> state` shape
[`ibd_block_extension_by_allele()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_block_extension_by_allele.md)
takes for its `allele =` argument.

## Usage

``` r
allele_states(vcf, position, names = c("base", "index"), one_based = FALSE)
```

## Arguments

- vcf:

  Path to a VCF/BCF. Needs `bcftools` on `PATH`.

- position:

  `"chr:pos"`. **0-based**, like every position in this package
  ([`?"plasgenomicsutilsR-coordinates"`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plasgenomicsutilsR-coordinates.md));
  pass `one_based = TRUE` to give a VCF `POS` instead. The chromosome
  spelling is normalised, so `13` and `Pf3D7_13_v3` both work.

- names:

  How to label a state. `"base"` (default) uses the alleles themselves,
  which is what a carrier contrast wants – `carrier = "C"` names
  something the reader can check against the callset. `"index"` uses the
  positional wording a legend wants (`"reference"`, `"alternate 1"`),
  which is what
  [`plot_region_haplotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_region_haplotypes.md)
  shows.

- one_based:

  Read `position` as a 1-based VCF `POS`.

## Value

A named character vector, one entry per sample in the callset, `NA`
where the genotype is missing.

## Details

This closes the gap that made a carrier contrast awkward to reach.
`allele =` has always accepted a named vector, but nothing in the
package produced one: the only code that read a per-sample allele set
was private and wired into
[`plot_region_haplotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_region_haplotypes.md),
so "the D384A carriers" was reachable from a hand-built metadata column
and from nowhere else.

## Why a set and not a dosage

A dosage says how many alternate copies a sample has, not **which**
alternate. At a codon carrying three independently arisen changes that
is the whole question, so the state is the allele set: `"C"`, `"G"`, or
`"C + G"` for a mixed infection carrying both. A sample with no call is
`NA` and is left out of every stratum downstream.

## See also

[`ibd_block_extension_by_allele()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_block_extension_by_allele.md),
[`plot_region_haplotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_region_haplotypes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
st <- allele_states("px1.bcf", "Pf3D7_13_v3:1725591")
table(st)
ibd_block_extension_by_allele(ibd, loci, allele = st, carrier = "C", reference = "A")
} # }
```
