# The tandem repeats a target design should avoid

The whole recipe in one call: merge repeats that run into one another
([`merge_tandem_repeats()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/merge_tandem_repeats.md)),
apply the period and length thresholds
([`flag_tandem_repeats()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/flag_tandem_repeats.md)),
keep what is flagged. Subtract the result from your targets with
[`bed_subtract()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/bed_subtract.md),
padding by a few bases so a primer cannot end right at a repeat's edge.

## Usage

``` r
tandem_repeats_to_avoid(
  x,
  min_width_by_period = c(`1` = 11, `2` = 12, `3` = 21),
  min_width = 50,
  gap = 0,
  dedupe = FALSE,
  rule = c("record", "block")
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
  `tandem_repeats_to_avoid()`.

- gap:

  Records separated by at most this many bases are merged too (default
  `0`: only overlapping or abutting ones).

- dedupe:

  Drop a repeat that lies entirely inside another on the same chromosome
  whose unit is the same up to rotation or reverse complement – `TAT`
  within a longer `ATA` run, `AT` within `ATAT` – before merging
  (default `FALSE`). The two finders often report one run twice that
  way. Such a record adds nothing: the container has the same period and
  at least the width, so the blocks and every flag come out the same;
  only `n_repeats` and the `repeat_units` list, which no longer names
  both spellings, change.

- rule:

  `"record"` or `"block"`; see above.

## Value

The
[`merge_tandem_repeats()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/merge_tandem_repeats.md)
table restricted to the blocks to avoid, without the `avoid` column.

## Details

`rule` says what the thresholds are applied to:

- `"record"` (the default): each repeat is judged on its own width, then
  merged; a merged run is kept when any repeat in it was flagged, or
  when the run as a whole reaches `min_width`, since a chain of
  individually harmless short repeats adds up to a stretch that slips
  like a long one. This reproduces the mask the HEOME Pf3D7 design used
  (per-record thresholds; every repeat merged and runs of 50 bp or more
  kept; the two sets merged), except that a short unflagged repeat
  abutting a flagged one is absorbed into its block, extending it by at
  most those few bases.

- `"block"`: repeats are merged first (with `gap`) and each block is
  judged on its *combined* width, held to the lowest threshold among the
  periods it contains. An 8 bp `A` run within 10 bp of a 9 bp `ATA` run
  is neither flagged alone, but together with the bases between them
  they are a 20 bp stretch holding a homopolymer, and 20 is past the
  homopolymer's 11. The combined width is the block's span, gap
  included.

## See also

[`pf3d7_tandem_repeats()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pf3d7_tandem_repeats.md)
for the bundled input,
[`bed_subtract()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/bed_subtract.md)
for what to do with the output,
[`write_bed()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/write_bed.md)
to write it.

## Examples

``` r
avoid <- tandem_repeats_to_avoid(pf3d7_tandem_repeats())
nrow(avoid); sum(avoid$width)                # blocks, and bases masked
#> [1] 75931
#> [1] 4156065

# stricter: merge repeats within 10 bp and judge each block on its combined width,
# with homopolymers and dinucleotides held to the trinucleotide length
strict <- tandem_repeats_to_avoid(pf3d7_tandem_repeats(), rule = "block", gap = 10,
                                  min_width_by_period = c("1" = 21, "2" = 21, "3" = 21))
nrow(strict); sum(strict$width)
#> [1] 56086
#> [1] 4169744

# the targets: a few genes minus the repeats, with 10 bp of clearance around each
genes <- PF3D7_GENES[PF3D7_GENES$name %in% c("pfcrt", "pfdhfr", "pfkelch13"), ]
targets <- bed_subtract(genes, avoid, pad = 10)
targets[, c("name", "start", "end", "piece", "width")]
#> # A tibble: 18 × 5
#>    name        start     end piece width
#>    <chr>       <dbl>   <dbl> <int> <dbl>
#>  1 pfdhfr     748087  749914     1  1827
#>  2 pfcrt      403221  403363     1   142
#>  3 pfcrt      403480  403823     2   343
#>  4 pfcrt      403899  404123     3   224
#>  5 pfcrt      404281  404418     4   137
#>  6 pfcrt      404478  404524     5    46
#>  7 pfcrt      404561  404692     6   131
#>  8 pfcrt      404725  404877     7   152
#>  9 pfcrt      404929  405103     8   174
#> 10 pfcrt      405146  405263     9   117
#> 11 pfcrt      405342  405454    10   112
#> 12 pfcrt      405510  405635    11   125
#> 13 pfcrt      405736  405883    12   147
#> 14 pfcrt      405921  405961    13    40
#> 15 pfcrt      406022  406129    14   107
#> 16 pfcrt      406175  406317    15   142
#> 17 pfkelch13 1724816 1726561     1  1745
#> 18 pfkelch13 1726617 1726997     2   380
if (FALSE) { # \dontrun{
write_bed(targets, "targets_no_repeats.bed")
} # }
```
