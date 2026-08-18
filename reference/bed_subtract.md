# Subtract one set of genomic intervals from another

The complement of
[`bed_intersect()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/bed_intersect.md),
and a dependency-free `bedtools subtract`: what is left of `locs1` once
everything in `locs2` is removed, cut at **base** resolution. An
interval partly covered comes back as the pieces that are not, one row
each, carrying its original columns; one covered end to end disappears.

## Usage

``` r
bed_subtract(
  locs1,
  locs2,
  chrom1 = "chr",
  start1 = "start",
  end1 = "end",
  chrom2 = "chr",
  start2 = "start",
  end2 = "end",
  min_width = 1
)
```

## Arguments

- locs1:

  Interval table to subtract from (a data frame / tibble).

- locs2:

  What to remove. An interval table, a character vector of `chr:pos` SNP
  ids, a genotype matrix or
  [`load_genotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/load_genotypes.md)
  list (its column names are the ids), or a
  [PopStructure](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PopStructure.md)
  (its genotype panel).

- chrom1, start1, end1:

  Column names in `locs1` (defaults `"chr"`, `"start"`, `"end"`;
  `"chrom"` is accepted as a chromosome alias).

- chrom2, start2, end2:

  As above for `locs2` when it is a table.

- min_width:

  Drop leftover pieces narrower than this (default `1`, i.e. keep every
  base). Raise it to ignore slivers.

## Value

A tibble of the uncovered pieces: every `locs1` column, with
`start`/`end` replaced by the piece's own bounds, plus `piece` (which
piece of that row this is) and `width`. Rows of `locs1` covered
completely are absent. Input row order is kept.

## Details

The use it was written for is filling gaps in a callset.
[`aa_intervals()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/aa_intervals.md)
gives the genomic span of each codon you care about; subtracting the
SNPs a panel already has leaves the bases that were never called, which
is what to hand `bcftools mpileup -R`. Because the cut is per base, a
codon with one of its three bases in the panel still returns the other
two – which is the point, since a residue cannot be read from one base.

## Coordinates

Intervals are **0-based half-open** `[start, end)`, as in BED and
throughout this package, so an interval abutting another
(`end1 == start2`) loses nothing. See
[plasgenomicsutilsR-coordinates](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plasgenomicsutilsR-coordinates.md).

## See also

[`bed_intersect()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/bed_intersect.md)
for the overlap,
[`write_bed()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/write_bed.md)
to write the result out.

## Examples

``` r
cds <- read_gff_cds(system.file("extdata", "pf3d7_drug_gene_cds.gff",
                                package = "plasgenomicsutilsR"))
want <- aa_intervals(data.frame(transcript_id = c("pfcrt", "pfcrt", "pfkelch13"),
                                aa_position = c(72, 76, 580)), cds)

# a panel that happens to carry only the middle base of the K76T codon
have <- paste0("Pf3D7_07_v3:", want$start[want$aa_position == 76] + 1)
gaps <- bed_subtract(want, have)
gaps[, c("name", "aa_position", "start", "end", "width")]
#> # A tibble: 4 × 5
#>   name                  aa_position   start     end width
#>   <chr>                       <int>   <dbl>   <dbl> <dbl>
#> 1 PF3D7_0709000.1-AA72           72  403611  403614     3
#> 2 PF3D7_0709000.1-AA76           76  403623  403624     1
#> 3 PF3D7_0709000.1-AA76           76  403625  403626     1
#> 4 PF3D7_1343700.1-AA580         580 1725257 1725260     3
```
