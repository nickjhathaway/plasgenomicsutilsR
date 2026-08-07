# Load the bundled example IBD results

Returns an
[IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
built from the small public example dataset shipped in `inst/extdata`:
per-SNP-per-group IBD, pairwise-group IBD, and the per-group selection
statistic (with thresholds) for five African countries, plus the
[PF_EXAMPLE_DRUG_GENES](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF_EXAMPLE_DRUG_GENES.md)
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
#>   per_snp_group : 6785 rows 
#>   pairwise_group: 20355 rows 
#>   selection      : 6785 rows 
#>   thresholds     : 5 
#>   genes          : 8 
```
