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
