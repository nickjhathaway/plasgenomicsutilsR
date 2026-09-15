# Split IBD block extension at a locus by carriage of a variant

The allele-level form of
[`ibd_block_extension_test()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_block_extension_test.md).
Within each group, the pairs are split by their state at a core variant
into carrier/carrier, reference/reference and discordant, the same
per-pair extension statistic is computed in each stratum, and the two
matching strata are compared. Extension measured over all pairs says a
locus is unusual in a group; this says whether it is the *haplotype*
that is unusual, which is a different and stronger claim.

## Usage

``` r
ibd_block_extension_by_allele(
  x,
  loci,
  allele,
  carrier = NULL,
  reference = NULL,
  group = NULL,
  within = 0,
  sharing = c("overlap", "complete"),
  min_ref_blocks = 1L,
  min_pairs = 5L,
  adjust = c("per_locus", "per_group", "all", "none"),
  meta = NULL
)
```

## Arguments

- x:

  An
  [IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
  built with `blocks =` and `meta =`.

- loci:

  Loci to test, as
  [`ibd_block_extension_test()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_block_extension_test.md)
  takes them.

- allele:

  The variant to split on: the name of a metadata column, or a named
  vector of `sample -> state`. Samples with `NA` are left out of every
  stratum.

- carrier, reference:

  Which states count as carrying and as reference. With exactly two
  states either may be left out and the other is implied; leaving both
  out reads the column's first two levels, `reference` first, with a
  message saying which way round it was read. With **three or more
  states both must be named** – there is nothing to imply the second one
  from, and guessing would contrast one alternate against another. Name
  them explicitly when the labels are not self-evident.

- group, within, sharing, min_ref_blocks, meta:

  Passed to
  [`ibd_block_extension_test()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_block_extension_test.md).

- min_pairs:

  Pairs a stratum needs at a locus before its length statistic is
  computed (default `5`). The fraction statistic is reported whatever
  the count, since it has a denominator either way.

- adjust:

  Benjamini-Hochberg scope for `p_length` and `p_fraction`, as in
  [`ibd_block_extension_test()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_block_extension_test.md).

## Value

A tibble, one row per locus x group, ordered by locus then by
`p_fraction`:

- `locus`, `name`, `gene_id`, `chr`, `start`, `end`, `span_bp`, `group`:

  as in
  [`ibd_block_extension_test()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_block_extension_test.md).

- `n_carrier`, `n_reference`, `n_discordant`:

  pairs IBD across the locus in each stratum. Report all three; see the
  note on the discordant one.

- `n_excluded_other_allele`:

  pairs IBD across the locus that were left out because an end carries a
  state that is neither `carrier` nor `reference`. Zero at a two-state
  locus. These are not discordant and not unknown, so they would
  otherwise be invisible – and at a multiallelic locus they can be the
  pairs that would show the interval is wider than the haplotype.

- `ratio_carrier`, `ratio_reference`:

  `paired_ratio` within each stratum.

- `ratio_contrast`:

  `ratio_carrier / ratio_reference` on the per-pair scale, so `1` means
  carriage makes no difference to segment length.

- `p_length`:

  two-sided rank-sum test of the two strata's per-pair log ratios.

- `n_carrier_possible`, `n_reference_possible`:

  pairs that *could* share, from the analyzed samples of each stratum.

- `frac_carrier`, `frac_reference`:

  sharing pairs over possible pairs.

- `odds_ratio`, `p_fraction`:

  Fisher's exact test of those two fractions. The better-powered
  comparison.

- `q_length`, `q_fraction`:

  adjusted within the family set by `adjust`.

The per-stratum full statistics are on the `strata` attribute, and the
per-pair ratios, labelled by stratum, on `pair_ratios`.

## Details

Two statistics come back per locus and group, and they answer different
questions:

- **length**, `ratio_carrier` against `ratio_reference`, tested by
  `p_length`. Among pairs that share the locus at all, are the carriers'
  segments the longer ones?

- **fraction**, `frac_carrier` against `frac_reference`, tested by
  `p_fraction`. Of all the pairs that *could* share, how many do? This
  is usually the better powered of the two, because it uses every pair
  rather than only the sharing ones, so read it first.

## Which states are contrasted

The non-carrier class is the **named reference state only**: a sample
carrying some third state is excluded from both strata, not pooled into
the reference. That is deliberate and is the whole point at a
multiallelic site. At *pfpx1* codon 384, D384A, D384G and D384Y arose
independently, so pooling the alternates into one "not reference" class
would merge origins the analysis exists to separate. To contrast two
alternates directly, name them as `carrier` and `reference`; that is a
different question and it should look different in the call. The pairs
left out either way are counted in `n_excluded_other_allele`.

## The discordant stratum is a check, not a result

A pair sharing an interval by descent shares whatever allele sits in it,
so pairs that are IBD across the locus and discordant at the variant
should be rare. A large `n_discordant` points at genotyping error, at a
recombination inside the interval, or at the interval being wider than
the haplotype – look there before reading anything else.

## See also

[`ibd_block_extension_test()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_block_extension_test.md),
[`ibd_block_extension_scan()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_block_extension_scan.md)
for the genome-wide background,
[`run_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_ihs.md)
for the haplotype-length statistic this mirrors.

## Examples

``` r
if (FALSE) { # \dontrun{
ibd <- ibd_results(blocks = "hmm.txt", meta = meta, group_col_in_meta = "region")
ibd_block_extension_by_allele(ibd, loci = "pfpx1", allele = "PIN_Haplotype",
                              carrier = "Present", reference = "Absent")
} # }
```
