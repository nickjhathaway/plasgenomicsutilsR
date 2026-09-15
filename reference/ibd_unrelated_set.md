# The largest mutually unrelated set of samples

The complement of
[`ibd_pair_clusters()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_pair_clusters.md):
instead of who is related to whom, the largest group found in which
**no** pair shares more than `max_ibd` of the genome. Use it to pick
founder haplotypes for simulation, to de-duplicate a cohort before
statistics that assume independent samples (allele frequencies,
\\F\_{ST}\\, PCA, admixture), or to report how much independent signal a
cohort really carries.

## Usage

``` r
ibd_unrelated_set(
  pairs,
  weight = NULL,
  max_ibd = 0.01,
  samples = NULL,
  restarts = 200L,
  seed = 1L,
  add_meta_cols = NULL,
  meta = NULL
)
```

## Arguments

- pairs:

  The per-pair table from
  `plasgenomicsutils ibd_fraction_and_snp_density`
  (`*.pair_ibd_fraction.tsv.gz`): a path, a data frame, or an
  [IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
  carrying one.

- weight:

  Column holding the fraction. Defaults to `ibd_fraction_accessible`,
  the callable-genome denominator; `ibd_fraction_full_genome` divides by
  the whole genome instead, so it reads lower for the same pair.

- max_ibd:

  Allow no pair in the set to share more than this fraction (default
  `0.01`).

- samples:

  Optional subset of samples to keep.

- restarts:

  Randomised restarts to try (default `200L`); the first is
  deterministic.

- seed:

  Seed for tie-breaking, so a run reproduces. The caller's random state
  is restored afterwards.

- add_meta_cols:

  Metadata columns to carry onto each row, by sample. A sample missing
  from `meta` gets `NA`; factor columns keep their level order, so a
  region stays in map order rather than turning alphabetical.

- meta:

  Sample metadata for `add_meta_cols`; taken from the
  [IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
  when `pairs` is one, so it only needs giving alongside a plain table.

## Value

A tibble with one row per analysed sample, selected first then by name:

- `sample`:

  the sample.

- `selected`:

  whether it is in the unrelated set.

- `n_links`:

  how many samples it shares more than `max_ibd` with.

- `max_ibd`:

  the most it shares with any of them, `NA` when none.

- `blocked_by`:

  for an unselected sample, a selected sample it conflicts with, so the
  result can be audited; `NA` when selected.

plus one column for each of `add_meta_cols`. Carries a `set_size`
attribute and a `max_ibd_within_set` attribute – the most any two
selected samples actually share, which is at most `max_ibd` and says how
unrelated the set really is.

## Details

Keeping one sample per single-linkage cluster is **not** the same thing,
and is usually much worse. The two coincide only when clusters are
cliques; at a low cutoff the relatedness graph is typically one
sprawling component, so one-per-cluster returns a handful of samples
where this returns many.

Finding the true maximum is NP-hard, so this is a heuristic: greedy
minimum-degree with randomised restarts, then a pass that adds back any
sample that turns out to fit. The result is always **valid** (no pair
inside it exceeds `max_ibd`) and **maximal** (no sample could be added),
but is not guaranteed to be the largest such set that exists. More
`restarts` searches harder. Samples sharing with nobody are always
selected.

`weight` should hold a fraction computed from **length-filtered** IBD
segments. A raw posterior site fraction carries a background floor –
routinely a median near 0.015 across unrelated pairs – so a cutoff near
1% on one of those counts almost every pair as related and collapses the
set to nearly nothing.

Set size falls away steeply as `max_ibd` tightens, so scan a few cutoffs
rather than trusting one. The cutoff is an inclusive ceiling: a pair
sharing exactly `max_ibd` may stay, matching the strictly-greater
comparison
[`plot_ibd_pair_network()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_pair_network.md)
draws with.

## See also

[`ibd_pair_clusters()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_pair_clusters.md)
for the related-to-whom view,
[`ibd_pair_links()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_pair_links.md)
for the edges themselves,
[`plot_ibd_pair_network()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_pair_network.md)
for the picture.

## Examples

``` r
if (FALSE) { # \dontrun{
u <- ibd_unrelated_set(ibd_all, max_ibd = 0.01)
sum(u$selected)
attr(u, "max_ibd_within_set")   # <= 0.01 by construction

# how the set shrinks as the cutoff tightens
sapply(c(0.005, 0.01, 0.02, 0.05),
       function(t) sum(ibd_unrelated_set(ibd_all, max_ibd = t)$selected))

# founder haplotypes for a simulation
writeLines(u$sample[u$selected], "founders.txt")

# does the unrelated set still span every site?
u <- ibd_all$ibd_unrelated_set(max_ibd = 0.01, add_meta_cols = "region")
table(subset(u, selected)$region)
} # }
```
