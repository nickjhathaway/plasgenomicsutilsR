# Which intervals each SNP falls in

Takes any per-SNP table –
[`run_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_ihs.md),
[`run_rsb()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_rsb.md),
[`run_xpehh()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_xpehh.md),
[`beta_score()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/beta_score.md),
[`pop_diff_snps()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff_snps.md),
an `IbdResults` selection table, anything with a `snp_id` or a
`chr`/`pos` pair – and reports the intervals covering each SNP. The
usual use is "which gene is this hit in", which otherwise means writing
the overlap join by hand every time.

## Usage

``` r
annotate_snps(
  scan,
  intervals,
  within = 0,
  keep = c("all", "hits"),
  collapse = FALSE,
  prefix = "",
  one_based_snps = FALSE
)
```

## Arguments

- scan:

  A data frame with `snp_id` (`"chr:pos0"`) or with `chr` and `pos`
  columns.

- intervals:

  A data frame with `name`, `chr` (or `chrom`), `start` and `end` –
  [PF3D7_GENES](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF3D7_GENES.md),
  [PF_EXAMPLE_DRUG_GENES](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF_EXAMPLE_DRUG_GENES.md),
  [PF3D7_CORE_REGIONS](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF3D7_CORE_REGIONS.md)
  or your own BED.

- within:

  Widen every interval by this many bp on both sides before testing, for
  catching a hit just outside a short gene (default `0`).

- keep:

  `"all"` (default) keeps SNPs matching nothing, with `NA` annotations;
  `"hits"` keeps only SNPs inside an interval.

- collapse:

  Return one row per SNP, with multiple hits pasted into `name` and
  counted in `n_intervals`, instead of one row per SNP x interval
  (default `FALSE`).

- prefix:

  Prepend this to the added column names, to keep two annotations side
  by side (e.g. `prefix = "core_"`).

- one_based_snps:

  The positions in `scan` are 1-based (VCF `POS`), so shift them before
  testing. Genotype-matrix column names from
  [`load_genotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/load_genotypes.md)
  are 1-based, and joining those against these 0-based intervals without
  saying so moves every SNP one base and can turn a near-miss into a
  hit. Default `FALSE`, the package convention.

## Value

`scan` with `name`, `interval_start`, `interval_end`,
`distance_to_midpoint` and any `gene_id` carried through – plus
`n_intervals` when `collapse = TRUE`. The original columns and their
order are preserved.

## Details

Positions are 0-based and intervals half-open `[start, end)`, the
package convention
([`?"plasgenomicsutilsR-coordinates"`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plasgenomicsutilsR-coordinates.md)),
so a SNP at the interval's `end` is *outside* it. Chromosome names are
normalised on both sides, so `Pf3D7_07_v3`, `chr7` and `7` match.

A SNP in no interval is kept with `NA` (or dropped by `keep = "hits"`),
and a SNP in several – overlapping gene spans, or nested intervals –
yields one row per interval rather than being silently collapsed. Use
`collapse = TRUE` for one row per SNP with the names pasted together
instead.

## See also

[`bed_intersect()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/bed_intersect.md)
for interval-to-interval overlap;
[`selection_peaks()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/selection_peaks.md)
to merge neighbouring significant SNPs into loci before annotating.

## Examples

``` r
ihs <- data.frame(snp_id = c("Pf3D7_07_v3:403500", "Pf3D7_07_v3:1"),
                  ihs = c(4.2, 0.1))
annotate_snps(ihs, PF_EXAMPLE_DRUG_GENES)
#> # A tibble: 2 × 7
#>   snp_id      ihs name  interval_start interval_end distance_to_midpoint gene_id
#>   <chr>     <dbl> <chr>          <int>        <int>                <dbl> <chr>  
#> 1 Pf3D7_07…   4.2 pfcrt         403221       406317                -1269 PF3D7_…
#> 2 Pf3D7_07…   0.1 NA                NA           NA                   NA NA     
annotate_snps(ihs, PF_EXAMPLE_DRUG_GENES, keep = "hits")
#> # A tibble: 1 × 7
#>   snp_id      ihs name  interval_start interval_end distance_to_midpoint gene_id
#>   <chr>     <dbl> <chr>          <int>        <int>                <dbl> <chr>  
#> 1 Pf3D7_07…   4.2 pfcrt         403221       406317                -1269 PF3D7_…
```
