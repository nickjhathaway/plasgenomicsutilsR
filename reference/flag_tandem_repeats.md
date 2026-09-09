# Flag the tandem repeats a polymerase is likely to slip on

Adds an `avoid` column to a
[`tandem_repeats()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/tandem_repeats.md)
table. A repeat is flagged when it is at least `min_width` bases long
whatever its unit, or when its **period** has its own, shorter,
threshold and the run reaches that. Short units slip at shorter lengths,
so the defaults descend with the period:

## Usage

``` r
flag_tandem_repeats(
  x,
  min_width_by_period = c(`1` = 11, `2` = 12, `3` = 21),
  min_width = 50
)
```

## Arguments

- x:

  A
  [`tandem_repeats()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/tandem_repeats.md)
  table, or anything
  [`tandem_repeats()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/tandem_repeats.md)
  accepts; or the block table from
  [`merge_tandem_repeats()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/merge_tandem_repeats.md).

- min_width_by_period:

  Minimum total width (bp) to flag, per period. A named numeric vector
  whose names are periods, or an unnamed one read as periods 1, 2, 3,
  ... in order. `NULL` for no per-period rule.

- min_width:

  Minimum total width (bp) to flag a repeat of any period. Also the
  length from which a *merged* run is avoided in
  [`tandem_repeats_to_avoid()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/tandem_repeats_to_avoid.md).

## Value

`x` as a tibble with a logical `avoid` column.

## Details

|                   |              |                       |
|-------------------|--------------|-----------------------|
| period            | flagged from | i.e.                  |
| 1 (homopolymer)   | 11 bp        | `AAAAAAAAAAA`         |
| 2 (dinucleotide)  | 12 bp        | six copies of `AT`    |
| 3 (trinucleotide) | 21 bp        | seven copies of `AAT` |
| any               | 50 bp        |                       |

These are the thresholds the HEOME Pf3D7 design used. Tighten or loosen
them per period, add periods
(`c("1" = 11, "2" = 12, "3" = 21, "4" = 30)`), or drop the period rules
altogether (`min_width_by_period = NULL`) to flag on total length only.

**Records or blocks.** Given a
[`tandem_repeats()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/tandem_repeats.md)
table, each repeat is judged on its own width. Given the blocks from
[`merge_tandem_repeats()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/merge_tandem_repeats.md),
each block is judged on its *combined* width, held to the lowest
threshold among the periods it contains: a block holding a homopolymer
is flagged from 11 bp however much of it is something else, so an 8 bp
`A` run flowing into a 9 bp `ATA` run is avoided as one stretch although
neither would be on its own.
[`tandem_repeats_to_avoid()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/tandem_repeats_to_avoid.md)
offers both as `rule`.

## See also

[`tandem_repeats_to_avoid()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/tandem_repeats_to_avoid.md),
[`merge_tandem_repeats()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/merge_tandem_repeats.md)

## Examples

``` r
bed <- data.frame(chrom = "Pf3D7_07_v3", start = c(100, 200, 300, 400),
                  end = c(110, 212, 320, 460),
                  name = c("Pf3D7_07_v3-100-110__A_x10",      # 10 bp homopolymer: kept
                           "Pf3D7_07_v3-200-212__AT_x6",      # 12 bp dinucleotide: flagged
                           "Pf3D7_07_v3-300-320__AAT_x6.67",  # 20 bp trinucleotide: kept
                           "Pf3D7_07_v3-400-460__TATTG_x12")) # 60 bp: flagged
flag_tandem_repeats(bed)[, c("repeat_unit", "period", "width", "avoid")]
#> # A tibble: 4 × 4
#>   repeat_unit period width avoid
#>   <chr>        <int> <dbl> <lgl>
#> 1 A                1    10 FALSE
#> 2 AT               2    12 TRUE 
#> 3 AAT              3    20 FALSE
#> 4 TATTG            5    60 TRUE 

# merged first: an 8 bp A run and a 9 bp ATA run 3 bp apart make one 20 bp block, which
# the homopolymer rule (11 bp) catches although neither repeat is flagged alone
close <- data.frame(chrom = "Pf3D7_07_v3", start = c(100, 111), end = c(108, 120),
                    name = c("A_x8", "ATA_x3"))
flag_tandem_repeats(close)$avoid
#> [1] FALSE FALSE
flag_tandem_repeats(merge_tandem_repeats(close, gap = 10))[, c("width", "periods", "avoid")]
#> # A tibble: 1 × 3
#>   width periods avoid
#>   <dbl> <chr>   <lgl>
#> 1    20 1,3     TRUE 
```
