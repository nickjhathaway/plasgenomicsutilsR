# Intersect two sets of genomic intervals

A lightweight, dependency-free "bedtools intersect"-style overlap
between two interval tables (e.g. genes vs. core regions, SNPs vs.
paralog masks). Two intervals overlap when they are on the same
chromosome and share at least one base
(`start1 < end2 & end1 > start2`); chromosome names are matched via
[`normalise_chr()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/normalise_chr.md)
so `"Pf3D7_07_v3"` and `"7"` agree.

## Usage

``` r
bed_intersect(
  locs1,
  locs2,
  chrom1 = "chr",
  start1 = "start",
  end1 = "end",
  chrom2 = "chr",
  start2 = "start",
  end2 = "end"
)
```

## Arguments

- locs1, locs2:

  Interval tables (data frames / tibbles).

- chrom1, start1, end1:

  Column names for the chromosome / start / end in `locs1` (defaults
  `"chr"`, `"start"`, `"end"`; `"chrom"` is accepted as a chromosome
  alias).

- chrom2, start2, end2:

  As above for `locs2`.

## Value

A list of three tibbles:

- `overlap`:

  one row per overlapping `locs1`x`locs2` pair – every `locs1` column,
  then every `locs2` column suffixed `.2`, then `overlap_start`,
  `overlap_end`, `overlap_bp`.

- `only1`:

  `locs1` rows overlapping nothing in `locs2`.

- `only2`:

  `locs2` rows overlapping nothing in `locs1`.

## Coordinates

Intervals are **0-based half-open** `[start, end)`, as in BED and
throughout this package – so abutting intervals (`end1 == start2`) do
not overlap. See
[plasgenomicsutilsR-coordinates](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plasgenomicsutilsR-coordinates.md).

## Examples

``` r
genes <- data.frame(name = c("a", "b"), chr = c("7", "7"),
                    start = c(100, 5000), end = c(200, 5200))
core  <- data.frame(chr = "7", start = 50, end = 1000)
hit <- bed_intersect(genes, core)
hit$overlap$name   # "a"  (inside core)
#> [1] "a"
hit$only1$name     # "b"  (subtelomeric)
#> [1] "b"
```
