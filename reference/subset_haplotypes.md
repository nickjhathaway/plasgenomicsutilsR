# Keep only some of the haplotypes

A
[`parasite_haplotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/parasite_haplotypes.md)
object restricted to certain samples or metadata groups, for when a scan
or an EHH curve is only worth reading within one part of the cohort – a
mutation that segregates in a single region, say, where pooling
everything buries it.

## Usage

``` r
subset_haplotypes(x, samples = NULL, meta = NULL, ...)
```

## Arguments

- x:

  A
  [`parasite_haplotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/parasite_haplotypes.md)
  object.

- samples:

  Optional sample ids to keep.

- meta:

  Metadata to match `...` against; defaults to the object's own.

- ...:

  Metadata filters as `column = values`, e.g. `region = "Southwest"`. A
  value that no sample has is an error rather than a silently empty
  result.

## Value

`x` with fewer haplotypes;
[`print()`](https://rdrr.io/r/base/print.html) reports the restriction.

## Details

Metadata columns are matched the way
[PopStructure](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PopStructure.md)'s
`$subset()` does, `column = values`, and a column takes several values:

    subset_haplotypes(hap, region = c("North", "Southwest"))

The **SNP panel is left alone**: the same columns, the same MAF and
missingness filtering that built them. That is deliberate, so two
subsets stay comparable to each other and to the whole – and it costs
nothing for haplotype work, since a SNP that is monomorphic within the
subset is dropped when the scan or curve is computed anyway. Rebuild
with
[`parasite_haplotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/parasite_haplotypes.md)
on a subset
[PopStructure](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PopStructure.md)
instead when you want the filtering itself redone within the group.

## See also

[`parasite_haplotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/parasite_haplotypes.md),
[`plot_ehh()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ehh.md),
[`run_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_ihs.md)

## Examples

``` r
ps <- example_pop_structure(umap = FALSE)
hap <- parasite_haplotypes(ps, maf = 0.05)
subset_haplotypes(hap, country = "Ghana")
#> <parasite_haplotypes> 30 haplotypes x 351 SNPs
#>   from            : 60 samples x 719 SNPs
#>   mixed calls     : 814 resolved by allele draw 
#>   SNPs dropped    : 221 missing, 147 MAF
#>   samples dropped : 0 missing
#>   imputed calls   : 178  seed: 42 
#>   subset          : 30 of 60 haplotypes (country: Ghana) 
```
