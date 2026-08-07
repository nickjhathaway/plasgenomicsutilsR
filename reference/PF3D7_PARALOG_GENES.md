# Pf3D7 paralogous / hypervariable gene families

Members of the *P. falciparum* multi-gene / hypervariable families (var,
rifin, stevor, surfin, ...) whose short reads mismap; commonly masked in
population-genetic analyses. Intersect against this with
[`bed_intersect()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/bed_intersect.md)
to flag genes that fall in these families.

## Usage

``` r
PF3D7_PARALOG_GENES
```

## Format

A data frame with columns `Pf3D7_chrom`, `start`, `end`, `chrom`,
`gene_id`, `description`; coordinates are 0-based half-open (see
[plasgenomicsutilsR-coordinates](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plasgenomicsutilsR-coordinates.md)).

## Source

The companion Python package (`plasgenomicsutils`
`builtin:pf3d7_paralog_genes`).

## See also

[PF3D7_CORE_REGIONS](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF3D7_CORE_REGIONS.md),
[`bed_intersect()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/bed_intersect.md)
