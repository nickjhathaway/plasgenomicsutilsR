# Read a table of short tandem repeats

Turns the BED a repeat finder writes into a table that says, for each
repeat, what its unit is, how long that unit really is (its *period*),
and how many bases the run covers – the three things slippage depends
on.
[`flag_tandem_repeats()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/flag_tandem_repeats.md)
applies the thresholds,
[`merge_tandem_repeats()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/merge_tandem_repeats.md)
joins runs that flow into one another, and
[`tandem_repeats_to_avoid()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/tandem_repeats_to_avoid.md)
does both.

## Usage

``` r
tandem_repeats(x)
```

## Arguments

- x:

  Path to a BED file (plain or gzipped), or a data frame holding the
  first four BED fields – by name (`chrom`/`chr`, `start`, `end`,
  `name`) or, failing that, by position.

## Value

A tibble with one row per repeat: `chrom`, `chr`, `start`, `end`,
`width`, `name`, `repeat_unit`, `unit_size`, `period`, `copies`. Input
order is kept.

## Details

Only the first four BED fields are read. The name field must be
`UNIT_xCOPIES` (`AACCCTA_x50.5714`), optionally preceded by a
`chrom-start-end__` location, which is what
`elucidator TandemRepeatFinderOutputToBed` and
`findSimpleTandemRepeatLocations` write
(`Pf3D7_01_v3-2-356__AACCCTA_x50.5714`); the location only restates the
coordinates and is ignored. A name of any other form gives `NA` for the
unit, period and copies, and such a row can still be flagged on its
total width alone.

**Period, not unit size.** A finder spells the same run as `AT`, `ATAT`
or `ATATAT` depending on where its alignment started. `period` is the
smallest number of bases the unit is a repetition of (2 for all three
spellings; 1 for `AAAA`; 4 for `TATT`), so a dinucleotide threshold
reaches every dinucleotide run however it was reported. `unit_size`
keeps the reported spelling's length. The same span is sometimes
reported twice with different units (the two finders disagree on the
alignment); both rows are kept, and merging collapses them.

## Coordinates

0-based half-open `[start, end)`, as in BED and throughout the package;
see
[plasgenomicsutilsR-coordinates](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plasgenomicsutilsR-coordinates.md).
`chr` is the chromosome normalised for matching
([`normalise_chr()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/normalise_chr.md)),
`chrom` as the file spells it.

## See also

[`pf3d7_tandem_repeats()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pf3d7_tandem_repeats.md)
for the bundled Pf3D7 table,
[`tandem_repeats_to_avoid()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/tandem_repeats_to_avoid.md)
for the whole recipe.

## Examples

``` r
bed <- data.frame(chrom = "Pf3D7_07_v3", start = c(100, 130, 500),
                  end = c(120, 160, 560),
                  name = c("Pf3D7_07_v3-100-120__A_x20",
                           "Pf3D7_07_v3-130-160__ATAT_x7.5",
                           "Pf3D7_07_v3-500-560__TATT_x15"))
tandem_repeats(bed)[, c("start", "end", "width", "repeat_unit", "period")]
#> # A tibble: 3 × 5
#>   start   end width repeat_unit period
#>   <dbl> <dbl> <dbl> <chr>        <int>
#> 1   100   120    20 A                1
#> 2   130   160    30 ATAT             2
#> 3   500   560    60 TATT             4
```
