# Create an [IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md) object

Convenience wrapper for `IbdResults$new()`.

## Usage

``` r
ibd_results(...)
```

## Arguments

- ...:

  Passed to
  [IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)'s
  constructor; see there for the arguments.

## Value

An
[IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
object.

## Examples

``` r
# every argument is optional -- hold only what you plan to plot
ibd <- example_ibd_results()
ibd
#> <IbdResults>  reference: pf3d7 
#>   per_snp_group : 6785 rows 
#>   pairwise_group: 20355 rows 
#>   selection      : 6785 rows 
#>   thresholds     : 5 
#>   genes          : 8 
ibd$get_group_order()
#> NULL
```
