# The samples a haplotype set kept

[`parasite_haplotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/parasite_haplotypes.md)
drops samples twice over: the polyclonal ones the Fws gate removes, and
any left too incomplete by `max_sample_missing`. This is what came
through both – the cohort the EHH and iHS scans actually ran on.

## Usage

``` r
haplotype_samples(hap)
```

## Arguments

- hap:

  A
  [`parasite_haplotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/parasite_haplotypes.md)
  object.

## Value

A character vector of sample ids, in the order the haplotype matrix
holds them.

## Details

It is the set to reuse when another analysis should describe the same
samples. `hap$filtering` says how many went to each cause, so it is
visible whether the Fws gate was the whole story.

## See also

[`subset_genotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/subset_genotypes.md)
to narrow a panel to them,
[`parasite_haplotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/parasite_haplotypes.md).

## Examples

``` r
ps <- example_pop_structure(umap = FALSE)
hap <- parasite_haplotypes(ps, maf = 0.05)
length(haplotype_samples(hap))
#> [1] 60
hap$filtering[c("n_dropped_polyclonal", "n_dropped_sample_missing")]
#> $n_dropped_polyclonal
#> [1] 0
#> 
#> $n_dropped_sample_missing
#> [1] 0
#> 
```
