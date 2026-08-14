# The SNPs an EHH plot had to choose between

[`plot_ehh()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ehh.md)
measures decay from a single focal SNP, and a gene usually holds many.
It takes the one with the most balanced alleles and reports which – this
returns the whole shortlist it chose from, chosen SNP first, so the
decision is visible and can be overridden by passing a `chr:pos` back to
`focal`.

## Usage

``` r
ehh_candidates(
  x,
  focal,
  group = NULL,
  genes = NULL,
  min_haplotypes = 10,
  reference = DEFAULT_REFERENCE
)
```

## Arguments

- x:

  A
  [`parasite_haplotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/parasite_haplotypes.md)
  object, or a
  [PopStructure](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PopStructure.md).

- focal:

  As in
  [`plot_ehh()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ehh.md):
  a `chr:pos` id, a bare position, or a gene name in `genes`.

- group:

  Optional metadata column, for the per-group columns.

- genes:

  Gene table used to resolve `focal` (e.g.
  [PF3D7_GENES](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF3D7_GENES.md)).

- min_haplotypes:

  Groups smaller than this are left out, as in
  [`plot_ehh()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ehh.md).

- reference:

  Reference id, used when `focal` names a whole chromosome.

## Value

A tibble of `snp_id`, `chr`, `pos`, `maf`, `n_hap` and `chosen`, plus
the per-group columns when `group` is given. The chosen SNP is the first
row; the rest follow by descending `maf`.

## Details

`maf` is the metric that decides it: minor-allele frequency across the
haplotypes, which
[`plot_ehh()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ehh.md)
maximises because EHH measured from a near-singleton is a flat line at 1
that says nothing about a sweep.

With `group`, each group gets a `maf_<group>` column and
`n_groups_variable` counts the groups the SNP actually varies in. That
is the other way a panel comes back empty: a SNP can be well balanced
overall and monomorphic inside one group, which drops that group's curve
with "the focal SNP is not variable in it". Sorting on
`n_groups_variable` finds a SNP that works everywhere, when one exists.

## See also

[`plot_ehh()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ehh.md)

## Examples

``` r
ps <- example_pop_structure(umap = FALSE)
hap <- parasite_haplotypes(ps, maf = 0.05)
ehh_candidates(hap, "pfcrt", genes = PF_EXAMPLE_DRUG_GENES)
#> # A tibble: 11 × 6
#>    snp_id             chr            pos    maf n_hap chosen
#>    <chr>              <chr>        <dbl>  <dbl> <int> <lgl> 
#>  1 Pf3D7_07_v3:403624 Pf3D7_07_v3 403624 0.5       60 TRUE  
#>  2 Pf3D7_07_v3:404835 Pf3D7_07_v3 404835 0.5       60 FALSE 
#>  3 Pf3D7_07_v3:405837 Pf3D7_07_v3 405837 0.5       60 FALSE 
#>  4 Pf3D7_07_v3:405361 Pf3D7_07_v3 405361 0.483     60 FALSE 
#>  5 Pf3D7_07_v3:405599 Pf3D7_07_v3 405599 0.483     60 FALSE 
#>  6 Pf3D7_07_v3:406230 Pf3D7_07_v3 406230 0.433     60 FALSE 
#>  7 Pf3D7_07_v3:403336 Pf3D7_07_v3 403336 0.333     60 FALSE 
#>  8 Pf3D7_07_v3:403674 Pf3D7_07_v3 403674 0.117     60 FALSE 
#>  9 Pf3D7_07_v3:403661 Pf3D7_07_v3 403661 0.1       60 FALSE 
#> 10 Pf3D7_07_v3:403686 Pf3D7_07_v3 403686 0.0667    60 FALSE 
#> 11 Pf3D7_07_v3:405561 Pf3D7_07_v3 405561 0.05      60 FALSE 
```
