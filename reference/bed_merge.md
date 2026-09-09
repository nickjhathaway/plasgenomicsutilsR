# Merge overlapping or nearby intervals

A dependency-free `bedtools merge`: intervals on one chromosome that
overlap, abut, or lie within `gap` bases of each other become one. The
complement of what
[`bed_subtract()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/bed_subtract.md)
does to a mask internally, exposed for when the merged blocks are what
you want to look at or write out.

## Usage

``` r
bed_merge(x, gap = 0, chrom = "chr", start = "start", end = "end")
```

## Arguments

- x:

  An interval table.

- gap:

  Intervals separated by at most this many bases are merged too (default
  `0`: only overlapping or abutting ones).

- chrom, start, end:

  Column names in `x` (`"chrom"` is accepted as an alias of the default
  `"chr"`).

## Value

A tibble of merged blocks in chromosome / start order: `chrom` (as `x`
spells it, taken from the block's first interval), `chr` (normalised),
`start`, `end`, `width`, and `n`, how many input intervals the block
absorbed.

## Coordinates

0-based half-open `[start, end)`; abutting intervals (`end1 == start2`)
are merged, as `bedtools merge` does by default. See
[plasgenomicsutilsR-coordinates](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plasgenomicsutilsR-coordinates.md).

## See also

[`merge_tandem_repeats()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/merge_tandem_repeats.md),
which also summarises what was merged;
[`bed_subtract()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/bed_subtract.md),
[`bed_intersect()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/bed_intersect.md)

## Examples

``` r
iv <- data.frame(chr = "7", start = c(100, 150, 300, 320), end = c(200, 250, 310, 330))
bed_merge(iv)              # 100-250 and two singletons
#> # A tibble: 3 × 6
#>   chrom chr   start   end width     n
#>   <chr> <chr> <dbl> <dbl> <dbl> <int>
#> 1 7     7       100   250   150     2
#> 2 7     7       300   310    10     1
#> 3 7     7       320   330    10     1
bed_merge(iv, gap = 10)    # 300-330 as well
#> # A tibble: 2 × 6
#>   chrom chr   start   end width     n
#>   <chr> <chr> <dbl> <dbl> <dbl> <int>
#> 1 7     7       100   250   150     2
#> 2 7     7       300   330    30     2
```
