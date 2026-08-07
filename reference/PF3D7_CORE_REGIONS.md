# Pf3D7 core genome regions

The **core** (non-subtelomeric, non-hypervariable) intervals of the *P.
falciparum* 3D7 core chromosomes. A locus is "core" when it overlaps one
of these intervals and "subtelomeric / internally hypervariable"
otherwise – intersect a gene or SNP table against this with
[`bed_intersect()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/bed_intersect.md).

## Usage

``` r
PF3D7_CORE_REGIONS
```

## Format

A data frame with columns `Pf3D7_chrom`, `start`, `end`, `chrom`;
coordinates are 0-based half-open (see
[plasgenomicsutilsR-coordinates](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plasgenomicsutilsR-coordinates.md)).

## Source

MalariaGEN `regions-20130225` core/non-core boundaries (with a few core
extensions bringing `hrp2` / `hrp3` into core), as shipped by the
companion Python package (`plasgenomicsutils`
`builtin:pf3d7_core_regions`).

## See also

[PF3D7_PARALOG_GENES](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF3D7_PARALOG_GENES.md),
[`bed_intersect()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/bed_intersect.md)
