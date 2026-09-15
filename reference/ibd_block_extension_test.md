# Test whether IBD blocks are longer at a locus than the sharing pairs' own background

For each locus and each group of samples, compares the length of the IBD
segment a pair shares **across the locus** with the length of that same
pair's segments elsewhere in the genome. Long segments at a locus are
the footprint of a recent sweep, but raw length is not comparable
between groups: a group with less outcrossing carries longer segments
genome-wide, so the naive comparison can run backwards. Making each pair
its own control removes that, and makes the pair – not the segment – the
unit of observation, which is what the signed-rank test needs.

## Usage

``` r
ibd_block_extension_test(
  x,
  loci = NULL,
  group = NULL,
  within = 0,
  sharing = c("overlap", "complete"),
  pairs = c("within", "all", "between"),
  min_ref_blocks = 1L,
  min_pairs = 5L,
  adjust = c("per_locus", "per_group", "all", "none"),
  keep_pair_ratios = TRUE,
  meta = NULL
)
```

## Arguments

- x:

  An
  [IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
  built with `blocks =` and `meta =`.

- loci:

  Loci to test: gene names from the object's track, `NULL` for every
  gene in it, or an interval data frame (`name`, `chr`/`chrom`, `start`,
  `end`, optionally `gene_id`). Rows of an interval frame that share a
  name, a chromosome and a `gene_id` are merged into one interval
  spanning them, so an
  [`aa_intervals()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/aa_intervals.md)
  table for a marker set can be passed as-is once it carries a `name`,
  while a gene track passed this way stays one locus per gene even where
  a name repeats across chromosomes.

- group:

  Metadata column defining the groups. Defaults to the object's declared
  group column, then to the first non-`sample` column of `meta`.

- within:

  Pad each locus by this many bp on both sides when deciding whether a
  segment counts (default `0`). It applies under either `sharing`, so
  with padding `"complete"` asks the segment to cover the padded
  interval.

- sharing:

  What a pair's segment must do at the locus to be counted, the same
  choice
  [`plot_ibd_network()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_network.md)
  offers:

  `"overlap"`

  :   (default) the segment touches the locus anywhere, which is
      [`gene_ibd_pairs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/gene_ibd_pairs.md)'s
      rule: `start < end_locus & end > start_locus`. The pair shares
      *some* of the locus.

  `"complete"`

  :   the segment spans the whole locus. Far stricter, and it changes
      what the statistic measures – see the note below on why locus
      width then matters.

- pairs:

  Which sample pairs to count. `"within"` (default) uses only pairs
  whose two samples share a group and reports a `group` column.
  `"between"` uses only cross-group pairs and `"all"` uses both; those
  two report `group_a` and `group_b` **in place of** `group`, one row
  per unordered pair of groups, which is the shape
  [`plot_pairwise_block_extension()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_pairwise_block_extension.md)
  draws. Nothing about the statistic changes – a cross-group pair is
  still measured against its own segments off the locus chromosome – so
  the diagonal of an `"all"` result is exactly the `"within"` result.

- min_ref_blocks:

  Segments a pair needs off the locus chromosome before its baseline is
  usable (default `1`). One is deliberate rather than lax: a pair with a
  single segment genome-wide is a distant pair, and that segment is
  length-biased long because short ones fell under the floor, so a thin
  baseline reads long and pushes the ratio *down*. Raising this is
  conservative in the wrong direction – it drops whole groups while
  rarely moving a conclusion.

- min_pairs:

  Pairs a group needs at a locus before it is tested (default `5`).
  Thinner strata are dropped rather than reported untested.

- adjust:

  Scope of the Benjamini-Hochberg correction, which is a real choice and
  not a detail. `"per_locus"` (default) makes one family of the groups
  tested at each locus, which is right for a handful of hand-picked
  loci. `"per_group"` makes one family of the loci tested in each group,
  which is what a genome-wide scan wants. `"all"` is a single family
  over every row, `"none"` leaves `q_paired` as `NA`.

- keep_pair_ratios:

  Attach the per-pair log ratios as the `pair_ratios` attribute (default
  `TRUE`). Set `FALSE` for a genome-wide scan: on one real cohort that
  attribute was 288 MB of a 330 MB result, and nothing downstream of a
  scan reads it.

- meta:

  Sample metadata; taken from `x` when not given.

## Value

A tibble, one row per locus x group, ordered by locus then descending
`paired_ratio`, carrying the per-pair ratios it was built from on the
`pair_ratios` attribute (`locus`, `group`, `sample1`, `sample2`,
`locus_len`, `ref_n`, `ref_median`, `log2_ratio`):

- `locus`, `name`, `gene_id`:

  locus labels, as
  [`gene_ibd_overlap()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/gene_ibd_overlap.md)
  uses them.

- `chr`, `start`, `end`, `span_bp`:

  the interval tested, 0-based half-open. Under `sharing = "complete"`,
  `span_bp` is a floor on `locus_median`, so compare loci of like width.

- `group`:

  the group both samples of every counted pair belong to.

- `n_pairs`:

  pairs IBD across the locus with a usable baseline.

- `gw_median`:

  that group's genome-wide median segment length, over all its
  within-group pairs.

- `locus_median`:

  median segment length across the locus.

- `naive_ratio`:

  `locus_median / gw_median`. Descriptive only.

- `ref_median`:

  median over the counted pairs of each pair's own baseline. Well above
  `gw_median` means the sharing pairs are unusually related.

- `paired_ratio`:

  the statistic. `1` means a pair's segment at the locus is no longer
  than its segments elsewhere.

- `z_paired`:

  `median(log2_ratio)` over the standard error of their mean. A rough
  effect size for ranking loci, not the test.

- `p_paired`, `p_two_sided`:

  one-sided (`greater`) and two-sided signed-rank p-values. Read the
  two-sided column before calling anything shorter: a one-sided p near
  `1` means "not longer", never "significantly shorter". Ties use the
  normal approximation.

- `q_paired`:

  `p_paired` adjusted within the family set by `adjust`.

## Details

Two statistics come back, and only one of them is testable:

- `naive_ratio` is the locus median segment length over the group's
  genome-wide median. Numerator and denominator come from different sets
  of pairs, because pairs IBD at a locus are selected for being related
  enough to share it. There is no error model behind it, so it gets no
  p-value. It is reported to be compared against the next one.

- `paired_ratio` is
  `2^median(log2(locus segment / that same pair's median segment on every other chromosome))`.
  Background relatedness cancels pair by pair. **This is the one to
  quote.** `p_paired` tests it.

Where a locus's `ref_median` sits well above `gw_median`, the pairs
sharing it are unusually related and `naive_ratio` was counting that
background rather than anything local.

Each pair's baseline excludes the **whole locus chromosome**, so a
sweep's own footprint cannot leak into the reference and flatten the
ratio. Segments come from `x$get_blocks()`, which has already dropped
non-IBD rows and applied the short-segment floor; taking a filtered
numerator against an unfiltered denominator inflates these ratios
several-fold, so do not re-read the raw `hmm.txt` for one side of it.

Only pairs whose two samples fall in the same group are used, since a
cross-group pair has no single group to report under.

This measures segment *length*. The count analogue, XiR,s (Henden et al.
2018,
[doi:10.1371/journal.pgen.1007279](https://doi.org/10.1371/journal.pgen.1007279)
), is computed by the companion Python pipeline and is generally the
better-powered statistic; treat this as its complement rather than its
replacement. For background on excess IBD as a selection signal see
Albrechtsen et al. (2010,
[doi:10.1534/genetics.110.113977](https://doi.org/10.1534/genetics.110.113977)
); for the principled coalescent-time version of what this ratio
approximates, ASMC (Palamara et al. 2018,
[doi:10.1038/s41588-018-0177-x](https://doi.org/10.1038/s41588-018-0177-x)
).

## Reading it

Do not run the test over segments instead of pairs. One pair contributes
several correlated segments and the p-value comes out far too small.
Splitting the pairs by carriage of a core variant, which is the IBD
analogue of iHS (see
[`run_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_ihs.md)),
turns a statement about a group into one about an allele.

A `paired_ratio` of 1 is **not** the null. A locus is a fixed target,
and a segment covers it with probability rising in the segment's own
length, so the numerator is drawn from a length-biased distribution
while each pair's baseline is not. Every locus is inflated by this,
whether or not anything happened there. Read a locus against the
distribution from scanning every gene in the same group, not against 1.

## Why `sharing` changes what locus width means

Under `"overlap"` a segment counts if it touches the locus, so the
selection weight rises with segment length *plus* locus width. Since a
gene is small against a typical IBD segment, every gene behaves like a
point and takes the same length bias: on one real cohort, locus width
explained an R-squared below 0.002 of `log2(paired_ratio)`, and binning
the scan by width only added noise.

Under `"complete"` the segment must cover the locus, so it cannot be
shorter than the locus is wide. Locus width becomes a hard floor under
`locus_median`, and the wider the locus the more of the short-segment
end of the distribution is cut away. The ratio then rises with width for
reasons that have nothing to do with selection, and loci of different
widths are no longer on one scale. So with `"complete"`, compare like
widths – bin the scan by `span_bp` – or stay with `"overlap"` when
ranking genes against each other.

## See also

[`gene_ibd_pairs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/gene_ibd_pairs.md)
for the pairs themselves,
[`gene_ibd_overlap()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/gene_ibd_overlap.md)
for the fraction of pairs sharing at a locus (the count side of the same
question),
[plasgenomicsutilsR-coordinates](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plasgenomicsutilsR-coordinates.md)
for the interval convention.

## Examples

``` r
# five related pairs, each sharing 80 kb across the locus and 20 kb elsewhere
pairs <- data.frame(sample1 = paste0("s", seq(1, 9, by = 2)),
                    sample2 = paste0("s", seq(2, 10, by = 2)))
blocks <- rbind(
  data.frame(pairs, chr = "7", start = 3.6e5, end = 4.4e5),   # over the locus
  data.frame(pairs, chr = "1", start = 1e5,   end = 1.2e5),   # baseline
  data.frame(pairs, chr = "2", start = 1e5,   end = 1.2e5))
blocks$different <- 0
blocks$Nsnp <- 40
meta <- data.frame(sample = paste0("s", 1:10), region = "north")
ibd <- ibd_results(blocks = blocks, meta = meta, group_col_in_meta = "region")
locus <- data.frame(name = "sweep", chr = "7", start = 4e5, end = 4.1e5)
ibd_block_extension_test(ibd, locus)[, c("locus", "group", "n_pairs", "paired_ratio")]
#> # A tibble: 1 × 4
#>   locus group n_pairs paired_ratio
#>   <fct> <chr>   <int>        <dbl>
#> 1 sweep north       5         4.00
```
