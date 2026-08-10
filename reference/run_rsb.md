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

**Reading `value`.** Rsb is a log ratio of site-specific extended
haplotype homozygosity between the two populations, standardised to
roughly a standard normal under neutrality. So it is a z-score, and its
**sign says which population**: positive means haplotype homozygosity
extends further in `pop1` than `pop2`, i.e. the sweep is in `pop1`;
negative points at `pop2`. Swapping the pair flips the sign.

Magnitude reads like any z: \|Rsb\| above ~2 is unremarkable in a genome
scan, above ~4 is worth a look, and real sweeps in *P. falciparum* run
higher still. Rather than picking a cutoff by eye, use `neg_log10_p` –
the two-sided normal p-value rehh derives from `value` – and correct it:
at ~20k SNPs a Bonferroni line sits near 5.6, which is the convention in
the literature.
[`selection_peaks()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/selection_peaks.md)
will merge the survivors into loci, since one sweep spans many SNPs, and
[`annotate_snps()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/annotate_snps.md)
says which genes they are in.

Two cautions. The standardisation assumes most of the genome is neutral,
so a genome-wide p-value is relative to *this* comparison and not
comparable across pairs with different sample sizes. And Rsb contrasts
two populations, so a high score means they differ – it cannot by itself
tell a sweep in one from a loss of variation in the other; that is what
the EHH decay curves are for.

## References

Tang, K., Thornton, K. R. & Stoneking, M. (2007) A new approach for
using genome scans to detect recent positive selection in the human
genome. *PLoS Biology* 5, e171.
[doi:10.1371/journal.pbio.0050171](https://doi.org/10.1371/journal.pbio.0050171)
