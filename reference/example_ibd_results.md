# Load the bundled example IBD results

Returns an
[IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
built from the small public example dataset shipped in `inst/extdata`:
per-SNP-per-region IBD, pairwise-region IBD, and the per-region
selection statistic (with thresholds) for five African countries, plus
the
[EXAMPLE_DRUG_GENES](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/EXAMPLE_DRUG_GENES.md)
track. Use it to try the `plot_*()` functions without your own data.

## Usage

``` r
example_ibd_results()
```

## Value

An
[IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
object.

## Details

The data are derived from publicly available *P. falciparum*
whole-genome samples (five countries), run through the
`plasgenomicsutils` IBD tools and downsampled to a small SNP panel for a
compact fixture.

## Examples

``` r
ibd <- example_ibd_results()
ibd
#> <IbdResults>  reference: pf3d7 
#>   per_snp_region : 6785 rows 
#>   pairwise_region: 20355 rows 
#>   selection      : 6785 rows 
#>   thresholds     : 5 
#>   genes          : 5 
```
