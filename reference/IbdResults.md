# IBD post-analysis results

A container for the per-SNP, pairwise-group, and selection-statistic
tables produced by the Python `plasgenomicsutils ibd` tools, plus the
shared cumulative-genome coordinate the genome-wide plots use. Pass file
paths or data frames; each argument is optional so you can hold only
what you plan to plot. The `plot_*()` functions read from an object of
this class.

Expected columns (superset; extras are kept):

- `per_snp_group`: `chr`, `pos`, `frac_pairs_ibd`, optionally `group`
  (output of `analyze_ibd_matrix` per-SNP / per-SNP-per-group).

- `pairwise_group`: `chr`, `pos`, `group_a`, `group_b`, `frac_pairs_ibd`
  (per-SNP pairwise-group output).

- `selection`: `chr`, `pos`, a metric column (`neg_log10_p`,
  `chi2_stat`, or `z_score`), optionally `group` and `significant`
  (output of `ibd_selection_statistic`).

## Methods

### Public methods

- [`IbdResults$new()`](#method-IbdResults-initialize)

- [`IbdResults$set_group_order()`](#method-IbdResults-set_group_order)

- [`IbdResults$restrict_groups()`](#method-IbdResults-restrict_groups)

- [`IbdResults$subset_groups()`](#method-IbdResults-subset_groups)

- [`IbdResults$get_group_order()`](#method-IbdResults-get_group_order)

- [`IbdResults$get_group_col()`](#method-IbdResults-get_group_col)

- [`IbdResults$get_per_snp_group()`](#method-IbdResults-get_per_snp_group)

- [`IbdResults$get_pairwise_group()`](#method-IbdResults-get_pairwise_group)

- [`IbdResults$get_selection()`](#method-IbdResults-get_selection)

- [`IbdResults$get_thresholds()`](#method-IbdResults-get_thresholds)

- [`IbdResults$get_genes()`](#method-IbdResults-get_genes)

- [`IbdResults$get_blocks()`](#method-IbdResults-get_blocks)

- [`IbdResults$get_analyzed_samples()`](#method-IbdResults-get_analyzed_samples)

- [`IbdResults$get_meta()`](#method-IbdResults-get_meta)

- [`IbdResults$set_meta()`](#method-IbdResults-set_meta)

- [`IbdResults$get_gene_overlap()`](#method-IbdResults-get_gene_overlap)

- [`IbdResults$chrom_layout()`](#method-IbdResults-chrom_layout)

- [`IbdResults$reference_id()`](#method-IbdResults-reference_id)

- [`IbdResults$print()`](#method-IbdResults-print)

- [`IbdResults$plot_ibd_sharing_manhattan()`](#method-IbdResults-plot_ibd_sharing_manhattan)

- [`IbdResults$plot_selection_manhattan()`](#method-IbdResults-plot_selection_manhattan)

- [`IbdResults$plot_ibd_tugofwar()`](#method-IbdResults-plot_ibd_tugofwar)

- [`IbdResults$plot_ibd_pairwise_group_heatmap()`](#method-IbdResults-plot_ibd_pairwise_group_heatmap)

- [`IbdResults$plot_pairwise_ibd_for_genes()`](#method-IbdResults-plot_pairwise_ibd_for_genes)

- [`IbdResults$gene_ibd_overlap()`](#method-IbdResults-gene_ibd_overlap)

- [`IbdResults$gene_ibd_pairs()`](#method-IbdResults-gene_ibd_pairs)

- [`IbdResults$add_ibd_clusters()`](#method-IbdResults-add_ibd_clusters)

- [`IbdResults$plot_ibd_network()`](#method-IbdResults-plot_ibd_network)

- [`IbdResults$pos_selection_genes()`](#method-IbdResults-pos_selection_genes)

- [`IbdResults$clone()`](#method-IbdResults-clone)

------------------------------------------------------------------------

### `IbdResults$new()`

Create an IbdResults object.

#### Usage

    IbdResults$new(
      per_snp_group = NULL,
      pairwise_group = NULL,
      selection = NULL,
      threshold = NULL,
      genes = NULL,
      blocks = NULL,
      meta = NULL,
      gene_overlap = NULL,
      group_col_in_meta = NULL,
      min_block_snp = IBD_MIN_BLOCK_SNP,
      min_block_kb = IBD_MIN_BLOCK_KB,
      reference = "pf3d7"
    )

#### Arguments

- `per_snp_group`:

  Per-SNP (optionally per-group) IBD table: path or data frame.

- `pairwise_group`:

  Per-SNP pairwise-group IBD table: path or data frame.

- `selection`:

  IBD selection-statistic table: path or data frame.

- `threshold`:

  Significance threshold(s): a scalar, a named vector or
  `group`/`threshold` data frame (per group), or a path to a one-number
  file.

- `genes`:

  Optional gene-annotation track (`name`, `chr`, `start`, `end`) drawn
  as reference lines on the Manhattan plots.

- `blocks`:

  Optional hmmibd-rs IBD segments (path or data frame: `sample1`,
  `sample2`, `chr`, `start`, `end`, optional `different`). Enables
  block-based gene triangles
  ([`gene_ibd_overlap()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/gene_ibd_overlap.md)
  /
  [`plot_pairwise_ibd_for_genes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_pairwise_ibd_for_genes.md)):
  only IBD segments (`different == 0`) are kept, but the analyzed-sample
  set (the denominator) is taken from every row first, so pairs that are
  never IBD still count.

- `meta`:

  Optional sample metadata (path or data frame with a `sample` column
  plus grouping columns), used to group `blocks` pairs.

- `gene_overlap`:

  Optional precomputed per-gene per-group-pair block-overlap table (from
  `plasgenomicsutils ibd_gene_overlap`: `gene`, `group_a`, `group_b`,
  `frac_pairs_ibd`, ...), used directly by the gene triangles.

- `group_col_in_meta`:

  Name of the `meta` column that defines the grouping. It becomes the
  default `group` for the block-based tools, and if the column is a
  factor its levels set the group order for every loaded table
  (equivalent to calling
  `$set_group_order(levels(meta[[group_col_in_meta]]))`).

- `min_block_snp, min_block_kb`:

  Discard IBD segments carrying fewer than `min_block_snp` SNPs or
  shorter than `min_block_kb` kb (defaults `15` and `15`). Short,
  SNP-poor segments are commonly spurious, and this filter is
  conventionally applied before any IBD summary – it is built in so the
  blocks need not be pre-filtered. `0` disables either criterion. Only
  the IBD evidence is filtered: the analyzed-sample set (the denominator
  behind every fraction) still comes from every row of the blocks file,
  so a pair whose only segment is short still counts as compared. The
  SNP criterion needs the `Nsnp` column hmmibd-rs writes.

- `reference`:

  Reference id for chromosome lengths (default `"pf3d7"`).

#### Returns

An `IbdResults` object (invisibly self).

------------------------------------------------------------------------

### `IbdResults$set_group_order()`

Set the order of the groups for every loaded table, so facets, legends
and axes follow it. Errors if a group present in the results is missing
from `levels` (it would silently become `NA`); warns about levels no
result uses.

#### Usage

    IbdResults$set_group_order(levels)

#### Arguments

- `levels`:

  Group names in the desired order.

#### Returns

Invisibly self.

------------------------------------------------------------------------

### `IbdResults$restrict_groups()`

Drop or keep groups, in place. `subset_groups()` is the copying form.

#### Usage

    IbdResults$restrict_groups(keep = NULL, drop = NULL)

#### Arguments

- `keep`:

  Group labels to keep (`NULL` keeps all).

- `drop`:

  Group labels to remove. Applied after `keep`.

#### Returns

Invisibly self.

------------------------------------------------------------------------

### `IbdResults$subset_groups()`

A new `IbdResults` holding only some of the groups.

Everything this object carries is either summarised per group or per
group pair, so groups are the natural unit to cut on. All of the
per-SNP, group-pair, selection and threshold tables are filtered; a
group-pair row survives only when **both** of its groups do, a cell
against a dropped group having no meaning. `meta`, `blocks` and the
analysed-sample set narrow to the samples belonging to the surviving
groups, so block-derived output —
[`gene_ibd_overlap()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/gene_ibd_overlap.md),
[`gene_ibd_pairs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/gene_ibd_pairs.md),
[`plot_ibd_network()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_network.md)
— follows as well. The group order is trimmed to what survives.

Narrowing groups cannot recompute a summary, so every number kept is
still the one computed over that group's full sample set. Groups are
dropped, never re-derived.

#### Usage

    IbdResults$subset_groups(keep = NULL, drop = NULL)

#### Arguments

- `keep`:

  Group labels to keep (`NULL` keeps all).

- `drop`:

  Group labels to remove. Applied after `keep`, so passing only `drop`
  is the usual "everything except these" form.

#### Returns

A new `IbdResults`; this object is unchanged.

------------------------------------------------------------------------

### `IbdResults$get_group_order()`

The current group order (`NULL` when none has been set).

#### Usage

    IbdResults$get_group_order()

------------------------------------------------------------------------

### `IbdResults$get_group_col()`

The `meta` column naming the grouping, if one was declared.

#### Usage

    IbdResults$get_group_col()

------------------------------------------------------------------------

### `IbdResults$get_per_snp_group()`

Per-SNP (optionally per-group) IBD table with `cum_pos`.

#### Usage

    IbdResults$get_per_snp_group()

------------------------------------------------------------------------

### `IbdResults$get_pairwise_group()`

Per-SNP pairwise-group IBD table with `cum_pos`.

#### Usage

    IbdResults$get_pairwise_group()

------------------------------------------------------------------------

### `IbdResults$get_selection()`

Selection-statistic table with `cum_pos`.

#### Usage

    IbdResults$get_selection()

------------------------------------------------------------------------

### `IbdResults$get_thresholds()`

Threshold tibble (`group`, `threshold`) or `NULL`.

#### Usage

    IbdResults$get_thresholds()

------------------------------------------------------------------------

### `IbdResults$get_genes()`

Gene-annotation track with `cum_mid`, or `NULL`.

#### Usage

    IbdResults$get_genes()

------------------------------------------------------------------------

### `IbdResults$get_blocks()`

IBD segment table (`sample1`, `sample2`, `chr`, `start`, `end`), or
`NULL`.

#### Usage

    IbdResults$get_blocks()

------------------------------------------------------------------------

### `IbdResults$get_analyzed_samples()`

Analyzed-sample ids (from every block row, pre IBD filter), or `NULL`.

#### Usage

    IbdResults$get_analyzed_samples()

------------------------------------------------------------------------

### `IbdResults$get_meta()`

Sample metadata data frame, or `NULL`.

#### Usage

    IbdResults$get_meta()

------------------------------------------------------------------------

### `IbdResults$set_meta()`

Replace the metadata, keeping the declared group order applied. Used by
[`add_ibd_clusters()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/add_ibd_clusters.md)
to add derived columns; the `sample` column must survive.

#### Usage

    IbdResults$set_meta(meta)

#### Arguments

- `meta`:

  The new metadata data frame.

#### Returns

Invisibly self.

------------------------------------------------------------------------

### `IbdResults$get_gene_overlap()`

Precomputed per-gene block-overlap table, or `NULL`.

#### Usage

    IbdResults$get_gene_overlap()

------------------------------------------------------------------------

### `IbdResults$chrom_layout()`

Chromosome layout tibble (offsets, bands, axis mid-points).

#### Usage

    IbdResults$chrom_layout()

------------------------------------------------------------------------

### `IbdResults$reference_id()`

The reference id used for chromosome lengths.

#### Usage

    IbdResults$reference_id()

------------------------------------------------------------------------

### `IbdResults$print()`

Compact summary of what the object holds.

#### Usage

    IbdResults$print(...)

#### Arguments

- `...`:

  Ignored; present for the
  [`print()`](https://rdrr.io/r/base/print.html) generic.

------------------------------------------------------------------------

### `IbdResults$plot_ibd_sharing_manhattan()`

Genome-wide per-SNP IBD-sharing Manhattan. See
[`plot_ibd_sharing_manhattan()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_sharing_manhattan.md).

#### Usage

    IbdResults$plot_ibd_sharing_manhattan(...)

#### Arguments

- `...`:

  Passed to
  [`plot_ibd_sharing_manhattan()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_sharing_manhattan.md).

------------------------------------------------------------------------

### `IbdResults$plot_selection_manhattan()`

IBD selection-statistic Manhattan. See
[`plot_selection_manhattan()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_selection_manhattan.md).

#### Usage

    IbdResults$plot_selection_manhattan(...)

#### Arguments

- `...`:

  Passed to
  [`plot_selection_manhattan()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_selection_manhattan.md).

------------------------------------------------------------------------

### `IbdResults$plot_ibd_tugofwar()`

Selection/IBD "tug-of-war" mirror. See
[`plot_ibd_tugofwar()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_tugofwar.md).

#### Usage

    IbdResults$plot_ibd_tugofwar(...)

#### Arguments

- `...`:

  Passed to
  [`plot_ibd_tugofwar()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_tugofwar.md).

------------------------------------------------------------------------

### `IbdResults$plot_ibd_pairwise_group_heatmap()`

Group x group IBD heatmap along the genome. See
[`plot_ibd_pairwise_group_heatmap()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_pairwise_group_heatmap.md).

#### Usage

    IbdResults$plot_ibd_pairwise_group_heatmap(...)

#### Arguments

- `...`:

  Passed to
  [`plot_ibd_pairwise_group_heatmap()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_pairwise_group_heatmap.md).

------------------------------------------------------------------------

### `IbdResults$plot_pairwise_ibd_for_genes()`

Per-gene (or per-SNP) group x group IBD triangles. See
[`plot_pairwise_ibd_for_genes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_pairwise_ibd_for_genes.md).

#### Usage

    IbdResults$plot_pairwise_ibd_for_genes(...)

#### Arguments

- `...`:

  Passed to
  [`plot_pairwise_ibd_for_genes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_pairwise_ibd_for_genes.md).

------------------------------------------------------------------------

### `IbdResults$gene_ibd_overlap()`

Per-gene IBD-block overlap between groups (see
[`gene_ibd_overlap()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/gene_ibd_overlap.md)).

#### Usage

    IbdResults$gene_ibd_overlap(...)

#### Arguments

- `...`:

  Passed to
  [`gene_ibd_overlap()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/gene_ibd_overlap.md).

------------------------------------------------------------------------

### `IbdResults$gene_ibd_pairs()`

Sample pairs sharing IBD over each gene (see
[`gene_ibd_pairs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/gene_ibd_pairs.md)).

#### Usage

    IbdResults$gene_ibd_pairs(...)

#### Arguments

- `...`:

  Passed to
  [`gene_ibd_pairs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/gene_ibd_pairs.md).

------------------------------------------------------------------------

### `IbdResults$add_ibd_clusters()`

Add single-linkage IBD cluster ids to the metadata.

#### Usage

    IbdResults$add_ibd_clusters(...)

#### Arguments

- `...`:

  Passed to
  [`add_ibd_clusters()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/add_ibd_clusters.md).

------------------------------------------------------------------------

### `IbdResults$plot_ibd_network()`

Sample-level IBD network at a gene / locus (see
[`plot_ibd_network()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_network.md)).

#### Usage

    IbdResults$plot_ibd_network(...)

#### Arguments

- `...`:

  Passed to
  [`plot_ibd_network()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_network.md).

------------------------------------------------------------------------

### `IbdResults$pos_selection_genes()`

Genes overlapping (or within `within` bp of) above-threshold selection
SNPs. See
[`pos_selection_genes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pos_selection_genes.md).

#### Usage

    IbdResults$pos_selection_genes(...)

#### Arguments

- `...`:

  Passed to
  [`pos_selection_genes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pos_selection_genes.md).

------------------------------------------------------------------------

### `IbdResults$clone()`

The objects of this class are cloneable with this method.

#### Usage

    IbdResults$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
ibd <- example_ibd_results()
groups <- levels(factor(ibd$get_selection()$group))

# everything except one group, or only the ones named
ibd$subset_groups(drop = groups[1])
#> <IbdResults>  reference: pf3d7 
#>   per_snp_group : 5428 rows 
#>   pairwise_group: 13570 rows 
#>   selection      : 5428 rows 
#>   thresholds     : 4 
#>   genes          : 8 
pair <- ibd$subset_groups(keep = groups[1:2])
pair$plot_selection_manhattan()


# `restrict_groups()` is the same thing in place
ibd$clone(deep = TRUE)$restrict_groups(drop = groups[1])
```
