# Group-pair differentiation summary across statistics

One tidy table of every pairwise comparison with, for each statistic and
summary, its value – so Jost's D, Nei's Gst, Hedrick's G'st and Hudson's
Fst can be read together.

## Usage

``` r
pop_diff_table(
  x,
  group = NULL,
  statistics = c("jost_d", "gst", "gst_hedrick", "fst"),
  stats = c("mean", "top_mean", "max"),
  top = 0.05,
  meta = NULL,
  clamp = TRUE
)
```

## Arguments

- x:

  A
  [PopStructure](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PopStructure.md)
  (uses its full genotype matrix for the active samples) or a genotype
  matrix (samples x SNPs, 0/1/2, `NA`).

- group:

  Metadata column (for a `PopStructure`) or a per-sample grouping vector
  (for a matrix).

- statistics:

  Measures to include (any of `"jost_d"`, `"gst"`, `"gst_hedrick"`,
  `"fst"`; default all four).

- stats:

  Per-pair summaries to include (`"mean"`, `"median"`, `"top_mean"`,
  `"max"`; default mean + top_mean + max).

- top:

  Fraction of top SNPs for `"top_mean"` (default `0.05`).

- meta:

  When `x` is a matrix, an optional data frame with a `sample` column
  plus `group`; otherwise `group` is a vector aligned to the matrix
  rows.

- clamp:

  Clamp small negative estimates to 0 (default `TRUE`).

## Value

A data frame: `a`, `b`, `n_snps`, then one `statistic_stat` column each.

## Examples

``` r
if (FALSE) { # \dontrun{
ps <- example_pop_structure("africa")
pop_diff_table(ps, group = "site")
} # }
```
