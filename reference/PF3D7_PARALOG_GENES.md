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
`gene_id`, `description` and `strand` (`"+"` / `"-"`, as in
[PF3D7_GENES](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF3D7_GENES.md));
coordinates are 0-based half-open (see
[plasgenomicsutilsR-coordinates](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plasgenomicsutilsR-coordinates.md)).

## Source

The companion Python package (`plasgenomicsutils`
`builtin:pf3d7_paralog_genes`), with `strand` joined on `gene_id` from
the same VEuPathDB / PlasmoDB GFF
[PF3D7_GENES](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF3D7_GENES.md)
is built from. About a quarter of these are pseudogenes, so they are
absent from
[PF3D7_GENES](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF3D7_GENES.md),
which holds protein-coding genes only.

## See also

[PF3D7_CORE_REGIONS](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF3D7_CORE_REGIONS.md),
[`bed_intersect()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/bed_intersect.md)
