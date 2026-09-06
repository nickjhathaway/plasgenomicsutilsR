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
  maxgap = NA,
  scalegap = NA,
  discard_at_border = NULL,
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

- maxgap:

  Largest gap between consecutive SNPs, in base pairs, that the EHH
  integration may cross; `NA` (the default, and rehh's) lets it cross
  any gap. This matters more than its default suggests. A region with no
  SNPs – a centromere, a masked hypervariable block – has nothing to
  break the haplotype, so EHH runs flat across it and the integral
  accumulates `EHH x gap length`. The SNPs flanking such a hole then
  score on the width of the hole rather than on their haplotypes, in
  either direction: the ratio is diluted towards zero when both alleles
  carry EHH into the gap, and inflated when only one does. Pick a value
  from the data's own spacing (several dozen times the median gap leaves
  ordinary density untouched while stopping at a real hole) rather than
  from a round number.

- scalegap:

  Gaps wider than this are counted as being exactly this wide, rather
  than stopping the integration outright; `NA` (default) does not
  rescale. A softer form of `maxgap` – it caps a hole's contribution
  instead of refusing to cross it.

- discard_at_border:

  Return `NA` instead of a truncated integral when the integration runs
  into the end of a chromosome or a gap wider than `maxgap`. `NULL`
  (default) ties it to `maxgap`: off when no `maxgap` is set (so the
  markers nearest the telomeres are still scored), on when one is. On
  sparse markers this can empty the scan – if EHH never decays before
  the data runs out, every marker is at a border – so a scan that comes
  back mostly `NA` says so.

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

## Examples

``` r
if (FALSE) { # \dontrun{
hap <- parasite_haplotypes(ps, fws = fws, min_fws = 0.95)
xp <- run_xpehh(hap, group = "region", pop1 = "north", pop2 = "south")
} # }
```
