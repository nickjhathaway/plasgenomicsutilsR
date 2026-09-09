# Genome-wide IBD clusters, and which samples are connected

Who is related to whom across the genome as a table: one row per sample,
saying whether it shares more than `min_ibd` with anyone, and which
single-linkage cluster it falls in. These are the components
[`plot_ibd_pair_network()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_pair_network.md)
draws – the connected blobs and the grid of unconnected samples beneath
them – so the same `weight` and `min_ibd` give the same answer the
picture shows.

## Usage

``` r
ibd_pair_clusters(
  pairs,
  weight = NULL,
  min_ibd = 0.01,
  samples = NULL,
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

- min_ibd:

  Count a pair as linked only above this fraction (default `0.01`). Pass
  what you plotted with, or the clusters will not match the components.

- samples:

  Optional subset of samples to keep.

- add_meta_cols:

  Metadata columns to carry onto each row, by sample. A sample missing
  from `meta` gets `NA`; factor columns keep their level order, so a
  region stays in map order rather than turning alphabetical.

- meta:

  Sample metadata for `add_meta_cols`; taken from the
  [IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
  when `pairs` is one, so it only needs giving alongside a plain table.

## Value

A tibble with one row per sample, connected samples first (largest
cluster first, then by name) and the unconnected after them:

- `sample`:

  the sample.

- `connected`:

  whether it shares more than `min_ibd` with any other sample.

- `cluster_id`, `cluster_size`:

  its single-linkage cluster and how many samples are in it; `NA` when
  it is connected to nobody.

- `n_links`:

  how many samples it is linked to.

- `max_ibd`:

  the largest fraction it shares with any of them, `NA` when none.

plus one column for each of `add_meta_cols`.

## Details

Single linkage means a sample joins a cluster if it shares with **any**
member, so a chain of pairs is one cluster even where its ends never
share directly. Ids run largest cluster first, so `cluster_id == 1` is
the biggest group; a sample sharing with nobody gets `NA`.
[`add_ibd_clusters()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/add_ibd_clusters.md)
and
[`gene_ibd_pairs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/gene_ibd_pairs.md)
number their per-gene clusters the same way.

## See also

[`plot_ibd_pair_network()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_pair_network.md)
for the picture,
[`ibd_pair_links()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_pair_links.md)
for the edges themselves,
[`add_ibd_clusters()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/add_ibd_clusters.md)
and
[`gene_ibd_pairs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/gene_ibd_pairs.md)
for the per-gene equivalents.

## Examples

``` r
if (FALSE) { # \dontrun{
cl <- ibd_pair_clusters(ibd_all, min_ibd = 0.03)

# the connected samples, and the unconnected ones
cl$sample[cl$connected]
cl$sample[!cl$connected]

# the biggest cluster
subset(cl, cluster_id == 1)

# with where each sample came from, to see whether a cluster spans sites
cl <- ibd_all$ibd_pair_clusters(min_ibd = 0.03, add_meta_cols = c("region", "country"))
with(subset(cl, connected), table(cluster_id, region))

# colour the network by cluster, the way add_ibd_clusters() does per gene
ibd_all$set_meta(merge(ibd_all$get_meta(), cl[, c("sample", "cluster_id")], by = "sample"))
ibd_all$plot_ibd_pair_network(min_ibd = 0.03, color_group = "cluster_id")
} # }
```
