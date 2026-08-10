# Sample pairs sharing IBD over each gene

An adjacency list of the pairs that are IBD across one or more genes:
one row per sample pair x IBD block x gene, saying how much of the gene
the block covers. Pairs with no IBD over a gene are simply absent – use
[`gene_ibd_overlap()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/gene_ibd_overlap.md)
when you need the denominator (all pairs compared) rather than just the
sharing ones.

## Usage

``` r
gene_ibd_pairs(x, genes = NULL, within = 0)
```

## Arguments

- x:

  An
  [IbdResults](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
  built with `blocks =`.

- genes:

  Gene names from the object's track (default all), or a gene-interval
  data frame (`name`, `chr`/`chrom`, `start`, `end`). Several genes come
  back in one table.

- within:

  Pad each gene interval by this many bp on both sides when deciding
  whether a block overlaps it (default `0`). Coverage is always measured
  against the gene's own span, so a block that reaches only into the
  padding covers `0`.

## Value

A tibble with one row per pair x block x gene:

- `sample1`, `sample2`:

  the IBD pair, ordered so `sample1 < sample2`.

- `chr`:

  chromosome of the block and gene.

- `block_start`, `block_end`:

  the IBD segment, 0-based half-open.

- `gene`, `name`, `gene_id`:

  gene labels (`gene` is unique, see
  [`gene_ibd_overlap()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/gene_ibd_overlap.md)).

- `gene_start`, `gene_end`:

  the gene interval, 0-based half-open.

- `coverage`:

  `"complete"` when the block spans the whole gene, else `"partial"`.

- `covered_start`, `covered_end`:

  the covered portion of the gene – the gene's own bounds when
  `coverage` is `"complete"`.

- `gene_cluster_id`, `gene_cluster_size`:

  single-linkage cluster of samples sharing at this gene, and how many
  samples are in it. A sample joins a cluster if it shares with **any**
  member, so a chain of pairs is one cluster even where its ends never
  share directly – which is what
  [`plot_ibd_network()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_network.md)
  draws as a connected component. Ids run largest first, so `1` is the
  biggest group at that gene; they are per gene, so cluster 1 at `pfcrt`
  and cluster 1 at `pfdhps` are unrelated.

- `covered_bp`, `percent_covered`:

  width of that portion, and it as a percentage of the gene's length.

## Details

A pair appears more than once for a gene only if it has several separate
IBD segments spanning it.

## See also

[`gene_ibd_overlap()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/gene_ibd_overlap.md)
for the per-group-pair fractions,
[plasgenomicsutilsR-coordinates](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plasgenomicsutilsR-coordinates.md)
for the interval convention.

## Examples

``` r
if (FALSE) { # \dontrun{
ibd <- ibd_results(blocks = "hmm.txt", genes = PF_EXAMPLE_DRUG_GENES)
pairs <- gene_ibd_pairs(ibd, genes = c("pfcrt", "pfdhps"))
subset(pairs, coverage == "complete")
} # }
```
