# Look up a reference genome's facts by id

Look up a reference genome's facts by id

## Usage

``` r
get_reference(ref_id = DEFAULT_REFERENCE)
```

## Arguments

- ref_id:

  Reference id (case-insensitive), e.g. "pf3d7".

## Value

A list with `ref_id`, `species`, `assembly`, `core_chrom_lengths_bp`,
and `bp_per_cm`.

## Examples

``` r
ref <- get_reference("pf3d7")
ref$bp_per_cm
#> [1] 15000
head(ref$core_chrom_lengths_bp)
#>       1       2       3       4       5       6 
#>  640851  947102 1067971 1200490 1343557 1418242 
```
