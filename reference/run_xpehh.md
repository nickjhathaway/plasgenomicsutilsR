# Cross-population extended haplotype homozygosity (XP-EHH)

The allele-aware sibling of
[`run_rsb()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_rsb.md):
it contrasts the integrated EHH of the same allele between two
populations.

## Usage

``` r
run_xpehh(
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

A tibble shaped like
[`run_rsb()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_rsb.md)'s,
with `value` holding XP-EHH. Read it exactly as Rsb – a standardised log
ratio, positive when the extended haplotype is in `pop1` – the
difference being that XP-EHH integrates to a fixed distance while Rsb
uses the site-specific EHH, so XP-EHH is the more sensitive of the two
to a sweep that has gone nearly to fixation, where within-population
statistics like iHS lose power.

## References

Sabeti, P. C. et al. (2007) Genome-wide detection and characterization
of positive selection in human populations. *Nature* 449, 913-918.
[doi:10.1038/nature06250](https://doi.org/10.1038/nature06250)
