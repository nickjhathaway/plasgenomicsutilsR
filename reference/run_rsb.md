# Cross-population extended haplotype homozygosity (Rsb)

Compares the site-specific integrated EHH of two groups. Unlike
[`run_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_ihs.md),
Rsb finds sweeps that have gone to completion in one population, because
the comparison is against the other population rather than against the
other allele.

## Usage

``` r
run_rsb(
  hap,
  group,
  meta = NULL,
  pairs = NULL,
  polarized = FALSE,
  min_samples = 4,
  threads = 1
)
```

## Arguments

- hap:

  A
  [`parasite_haplotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/parasite_haplotypes.md)
  object.

- group:

  Metadata column naming the grouping (required – there must be at least
  two groups to compare).

- meta:

  Metadata (defaults to the one carried by `hap`).

- pairs:

  Optional list of `c(group_a, group_b)` pairs; defaults to all pairs.

- polarized:

  Treat allele 1 as derived (needs a real ancestral state).

- min_samples:

  Skip groups smaller than this.

- threads:

  Threads for rehh.

## Value

A tibble with `pair`, `pop1`, `pop2`, `chr`, `pos`, `snp_id`, `value`
(Rsb) and `neg_log10_p`.

## References

Tang, K., Thornton, K. R. & Stoneking, M. (2007) A new approach for
using genome scans to detect recent positive selection in the human
genome. *PLoS Biology* 5, e171.
[doi:10.1371/journal.pbio.0050171](https://doi.org/10.1371/journal.pbio.0050171)
