# Are a variant's carriers more related genome-wide to begin with?

The genome-wide background control for
[`ibd_block_extension_by_allele()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_block_extension_by_allele.md).
Schaffner et al.'s caution about reading IBD inside a sweep is that
isolates carrying the swept haplotype are typically more related to each
other across the *rest* of the genome as well, so a locus statistic
normalised only at the population level counts that background as locus
signal. This measures exactly that background: within each group it
splits the pairs by their state at a core variant and reports each
stratum's mean per-pair IBD off the **focal chromosome**, which is
dropped in full so the sweep cannot contribute to its own baseline.

## Usage

``` r
carrier_genome_wide_relatedness(
  x,
  allele,
  focal_chr,
  carrier = NULL,
  reference = NULL,
  group = NULL,
  min_pairs = 5L,
  meta = NULL
)
```

## Arguments

- x:

  An
  [IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
  built with `blocks =` and `meta =`.

- allele:

  The variant to split on: the name of a metadata column, or a named
  vector of `sample -> state` (as
  [`allele_states()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/allele_states.md)
  returns). Samples with `NA` are left out.

- focal_chr:

  The chromosome to exclude from the background, as the locus sits on
  it. Any spelling
  [`normalise_chr()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/normalise_chr.md)
  accepts (`"Pf3D7_07_v3"` or `"7"`).

- carrier, reference:

  Which states count as carrying and as reference; each may be a vector.
  With exactly two observed states either may be left out and the other
  is implied (`reference` the first level), with a message. With three
  or more states both must be named, since there is nothing to imply the
  second from.

- group:

  Metadata column defining the groups. Defaults to the object's declared
  group column, then to the first non-`sample` column of `meta`.

- min_pairs:

  A stratum needs this many pairs before it is reported, and the
  reference stratum needs it before a fold or a test is anchored on it
  (default `5`).

- meta:

  Sample metadata; taken from `x` when not given.

## Value

A tibble, one row per group x class, ordered by group then class:

- `group`, `class`:

  the group and one of `carrier/carrier`, `reference/reference`,
  `discordant`.

- `n_pairs`:

  every within-group pair of the class, including the pairs that share
  nothing – taking the mean over only the sharing pairs would compare
  two different denominators and hide the effect being tested.

- `mean_ibd_mb`, `median_ibd_mb`:

  megabases of IBD per pair, summed over that pair's segments off the
  focal chromosome.

- `mean_ibd_frac`:

  the mean as a share of the callable genome outside the focal
  chromosome (its core length), the scale the rest of an IBD document
  uses.

- `fold_vs_ref`:

  `mean_ibd_mb` over the group's `reference/reference` mean.

- `p_vs_ref`:

  two-sided rank-sum of the stratum against that same group's
  `reference/reference` pairs; `NA` on the reference row itself and
  where either side is under `min_pairs`.

## Details

Where `fold_vs_ref` is above 1 and `p_vs_ref` is small, the carriers are
a more related set genome-wide and a population-level normalisation at
the locus would be reading part of that as locus signal;
[`ibd_block_extension_test()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_block_extension_test.md)'s
`paired_ratio` and everything
[`ibd_block_extension_by_allele()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_block_extension_by_allele.md)
reports already divide each pair by its own baseline and are immune to
it, so the point of the table is to say how much that per-pair control
is removing, and where. Where the fold sits at 1, carriage says nothing
about background relatedness and the two normalisations agree.

## What counts as carrier and reference

Unlike
[`ibd_block_extension_by_allele()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_block_extension_by_allele.md),
which contrasts exactly two named states, `carrier` and `reference` are
**sets** here, matched with `%in%`: at *pfkelch13* every propeller
allele can count as carrier against a single `"REF"` reference, and a
mixed call naming both a carrier and a reference allele lands in
neither. A pair is `carrier/carrier` when both ends are in `carrier`,
`reference/reference` when both are in `reference`, and `discordant`
when one is a carrier and the other a reference. A pair with an end in
neither set (a third allele, or `NA`) is left out of every stratum.

## See also

[`ibd_block_extension_by_allele()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_block_extension_by_allele.md)
for the locus contrast this backstops,
[`pair_fraction_summary()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pair_fraction_summary.md)
for the whole-genome per-pair sharing by group,
[`ibd_block_extension_test()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_block_extension_test.md)
for the per-pair normalisation that absorbs this background.

## Examples

``` r
if (FALSE) { # \dontrun{
ibd <- ibd_results(blocks = "hmm.txt", meta = meta, group_col_in_meta = "region")
carrier_genome_wide_relatedness(ibd, allele = "PIN_Haplotype", focal_chr = "Pf3D7_07_v3",
                                carrier = "Present", reference = "Absent")
} # }
```
