# Merge tandem repeats that run into one another

Repeat finders report a stretch like `AAAAAAAAAAATATATATAT` as two
records, an `A` homopolymer and an `AT` run, that overlap or abut; one
trinucleotide often flows into another the same way. For avoiding
slippage the whole stretch is one region, so this joins records that
overlap, abut, or lie within `gap` bases, and keeps a record of what
went into each block.

## Usage

``` r
merge_tandem_repeats(x, gap = 0, dedupe = FALSE)
```

## Arguments

- x:

  A
  [`tandem_repeats()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/tandem_repeats.md)
  table, or anything
  [`tandem_repeats()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/tandem_repeats.md)
  accepts. An `avoid` column
  ([`flag_tandem_repeats()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/flag_tandem_repeats.md))
  is carried through as *any member flagged*.

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

## Value

A tibble of blocks in chromosome / start order: `chrom`, `chr`, `start`,
`end`, `width`, `name` (`chrom-start-end`), `n_repeats` (records
merged), `repeat_units` (the distinct units, comma-separated, left to
right along the block), `periods` (the distinct periods, ascending),
`min_period`, and `avoid` when `x` had it.

## See also

[`tandem_repeats_to_avoid()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/tandem_repeats_to_avoid.md),
[`bed_merge()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/bed_merge.md)
for plain intervals

## Examples

``` r
bed <- data.frame(chrom = "Pf3D7_07_v3", start = c(100, 115, 300),
                  end = c(116, 130, 330),
                  name = c("Pf3D7_07_v3-100-116__A_x16",
                           "Pf3D7_07_v3-115-130__AT_x7.5",
                           "Pf3D7_07_v3-300-330__AAT_x10"))
merge_tandem_repeats(bed)[, c("start", "end", "n_repeats", "repeat_units", "periods")]
#> # A tibble: 2 × 5
#>   start   end n_repeats repeat_units periods
#>   <dbl> <dbl>     <int> <chr>        <chr>  
#> 1   100   130         2 A,AT         1,2    
#> 2   300   330         1 AAT          3      

# the same run reported as `A` and `AA`: dedupe keeps one spelling
twice <- rbind(bed, data.frame(chrom = "Pf3D7_07_v3", start = 100, end = 116,
                               name = "Pf3D7_07_v3-100-116__AA_x8"))
merge_tandem_repeats(twice)$repeat_units[1]
#> [1] "A,AA,AT"
merge_tandem_repeats(twice, dedupe = TRUE)$repeat_units[1]
#> [1] "A,AT"
```
