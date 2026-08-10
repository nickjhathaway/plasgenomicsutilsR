# Add IBD cluster ids to the stored metadata

For each gene or locus, works out the single-linkage clusters of samples
sharing IBD over it and writes them into the object's metadata as
`<label>_cluster_id` (and `<label>_cluster_size` when `size = TRUE`).
Any plot that reads a metadata column can then use them – notably
[`plot_ibd_network()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_network.md),
where the clusters *are* the connected components being drawn:

## Usage

``` r
add_ibd_clusters(
  x,
  genes = NULL,
  locus = NULL,
  within = 0,
  sharing = c("overlap", "complete"),
  size = FALSE,
  prefix = ""
)
```

## Arguments

- x:

  An
  [IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
  with `blocks` and `meta` loaded.

- genes:

  Gene names in the object's track, or `NULL` for every gene in it.
  Ignored when `locus` is given.

- locus:

  A locus instead of a gene: `"chr:pos"`, `"chr:start-end"`, or a data
  frame with `chr`/`start`/`end`. Its `name` (or the coordinate string)
  labels the column.

- within:

  Pad the interval by this many bp on both sides before deciding whether
  a block overlaps it (default `0`).

- sharing:

  `"overlap"` (default) counts a pair when any IBD block touches the
  interval; `"complete"` only when a block spans the whole of it.

- size:

  Also add `<label>_cluster_size`, the number of samples in the cluster
  (default `FALSE`).

- prefix:

  Optional string put before each column name, for keeping several
  settings side by side (e.g. `prefix = "complete_"`). Re-running
  replaces any column of the same name, so changing `within` or
  `sharing` and calling again re-clusters rather than accumulating stale
  columns. Use `prefix` when you want two settings side by side.

## Value

Invisibly `x`, with the metadata extended. `x$get_meta()` shows the new
columns.

## Details

    ibd$add_ibd_clusters(genes = "pfcrt")
    plot_ibd_network(ibd, gene = "pfcrt", color_group = "pfcrt_cluster_id")

Single linkage means a sample joins a cluster if it shares with **any**
member, so a chain of pairs is one cluster even where its ends never
share directly. Ids run largest cluster first and are numbered per
interval, so cluster 1 at one gene is unrelated to cluster 1 at another.
A sample that shares with nobody over the interval gets `NA`.

`within` and `sharing` mean exactly what they do in
[`plot_ibd_network()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_network.md)
and
[`gene_ibd_pairs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/gene_ibd_pairs.md),
and they decide which edges exist – so pass the same values you plot
with, or the colours will not match the components.

## Examples

``` r
if (FALSE) { # \dontrun{
ibd <- ibd_results(blocks = "blocks.hmm.txt", meta = meta,
                   genes = PF_EXAMPLE_DRUG_GENES)
add_ibd_clusters(ibd)                       # one column per gene in the track
plot_ibd_network(ibd, gene = "pfcrt", color_group = "pfcrt_cluster_id")
} # }
```
