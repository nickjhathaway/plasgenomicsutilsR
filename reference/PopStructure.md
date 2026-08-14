# Population-structure workspace (PCA + UMAP + admixture)

An R6 object that bundles a genotype matrix, its PCA (the full
[`stats::prcomp()`](https://rdrr.io/r/stats/prcomp.html) result), an
optional UMAP embedding, per-sample metadata, a shared metadata colour
map, and an sNMF admixture fit. Because it keeps the fitted objects, you
can save it (`saveRDS`) and re-plot without recomputing, colour
PCA/UMAP/admixture consistently (same `level -> colour` map), reorder
admixture bars once and reuse the order across K, and sub-select to a
set of samples (or a metadata match) for output.

## Methods

### Public methods

- [`PopStructure$new()`](#method-PopStructure-initialize)

- [`PopStructure$add_panel()`](#method-PopStructure-add_panel)

- [`PopStructure$panels()`](#method-PopStructure-panels)

- [`PopStructure$add_meta()`](#method-PopStructure-add_meta)

- [`PopStructure$set_colors()`](#method-PopStructure-set_colors)

- [`PopStructure$set_levels()`](#method-PopStructure-set_levels)

- [`PopStructure$run_umap()`](#method-PopStructure-run_umap)

- [`PopStructure$run_snmf()`](#method-PopStructure-run_snmf)

- [`PopStructure$best_k()`](#method-PopStructure-best_k)

- [`PopStructure$cross_entropy()`](#method-PopStructure-cross_entropy)

- [`PopStructure$plot_cross_entropy()`](#method-PopStructure-plot_cross_entropy)

- [`PopStructure$get_snmf_fit()`](#method-PopStructure-get_snmf_fit)

- [`PopStructure$plot_admixture_multi_k()`](#method-PopStructure-plot_admixture_multi_k)

- [`PopStructure$q()`](#method-PopStructure-q)

- [`PopStructure$restrict()`](#method-PopStructure-restrict)

- [`PopStructure$subset()`](#method-PopStructure-subset)

- [`PopStructure$genotype()`](#method-PopStructure-genotype)

- [`PopStructure$pca_scores()`](#method-PopStructure-pca_scores)

- [`PopStructure$pca_variance()`](#method-PopStructure-pca_variance)

- [`PopStructure$prcomp()`](#method-PopStructure-prcomp)

- [`PopStructure$umap_df()`](#method-PopStructure-umap_df)

- [`PopStructure$get_meta()`](#method-PopStructure-get_meta)

- [`PopStructure$get_colors()`](#method-PopStructure-get_colors)

- [`PopStructure$get_colours()`](#method-PopStructure-get_colours)

- [`PopStructure$allele()`](#method-PopStructure-allele)

- [`PopStructure$positions()`](#method-PopStructure-positions)

- [`PopStructure$pruned()`](#method-PopStructure-pruned)

- [`PopStructure$get_samples()`](#method-PopStructure-get_samples)

- [`PopStructure$as_ps()`](#method-PopStructure-as_ps)

- [`PopStructure$plot_pca()`](#method-PopStructure-plot_pca)

- [`PopStructure$plot_umap()`](#method-PopStructure-plot_umap)

- [`PopStructure$plot_admixture()`](#method-PopStructure-plot_admixture)

- [`PopStructure$plot_figure()`](#method-PopStructure-plot_figure)

- [`PopStructure$pop_diff()`](#method-PopStructure-pop_diff)

- [`PopStructure$jost_d()`](#method-PopStructure-jost_d)

- [`PopStructure$pop_diff_table()`](#method-PopStructure-pop_diff_table)

- [`PopStructure$plot_diff_heatmap()`](#method-PopStructure-plot_diff_heatmap)

- [`PopStructure$plot_jost_d_heatmap()`](#method-PopStructure-plot_jost_d_heatmap)

- [`PopStructure$pop_diff_snps()`](#method-PopStructure-pop_diff_snps)

- [`PopStructure$plot_diff_manhattan()`](#method-PopStructure-plot_diff_manhattan)

- [`PopStructure$diversity()`](#method-PopStructure-diversity)

- [`PopStructure$ld_index()`](#method-PopStructure-ld_index)

- [`PopStructure$beta_score()`](#method-PopStructure-beta_score)

- [`PopStructure$haplotypes()`](#method-PopStructure-haplotypes)

- [`PopStructure$plot_region_haplotypes()`](#method-PopStructure-plot_region_haplotypes)

- [`PopStructure$ihs()`](#method-PopStructure-ihs)

- [`PopStructure$save()`](#method-PopStructure-save)

- [`PopStructure$print()`](#method-PopStructure-print)

- [`PopStructure$clone()`](#method-PopStructure-clone)

------------------------------------------------------------------------

### `PopStructure$new()`

Build from a genotype matrix (or a
[`load_genotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/load_genotypes.md)
list).

#### Usage

    PopStructure$new(
      geno,
      samples = NULL,
      meta = NULL,
      n_pcs = 50,
      colors = NULL,
      allele = NULL,
      pruned = NULL,
      full = NULL,
      one_based = FALSE,
      colours = NULL
    )

#### Arguments

- `geno`:

  Genotype matrix (samples x SNPs, 0/1/2, `NA`) or
  [`load_genotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/load_genotypes.md)
  list.

- `samples`:

  Sample ids (default row names / the list's `sample.id`).

- `meta`:

  Optional metadata (data frame with a `sample` column).

- `n_pcs`:

  Number of PCs to summarise.

- `colors`:

  Optional named list of colour maps (see
  [`meta_colors()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/meta_colors.md)).

- `allele`:

  Which allele the dosages count, `"alt"` or `"ref"`. Taken from a
  [`load_genotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/load_genotypes.md)
  list when it says so; a bare matrix cannot say, and the two codings
  are indistinguishable afterwards, so anything that names the calls has
  to be told.

- `pruned`:

  Whether the SNPs were LD-pruned. Taken from a
  [`load_genotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/load_genotypes.md)
  list when it says so. Pruning is right for PCA / admixture and wrong
  for looking at haplotypes, since it drops SNPs precisely for being
  correlated with a neighbour.

- `one_based`:

  The genotype column names carry 1-based VCF positions (a matrix built
  by an older
  [`load_genotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/load_genotypes.md),
  or read straight off a VCF); shift them to the 0-based convention on
  the way in.

------------------------------------------------------------------------

### `PopStructure$add_panel()`

Register another genotype panel under a name, so one object can hold the
pruned SNPs for PCA / admixture and the full set for the analyses where
the correlation between neighbouring SNPs is the signal.

#### Usage

    PopStructure$add_panel(name, geno, allele = NULL, pruned = NULL)

#### Arguments

- `name`:

  Panel name (`"full"` and `"pruned"` are the ones other functions ask
  for).

- `geno`:

  Genotype matrix, or a
  [`load_genotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/load_genotypes.md)
  list.

- `allele, pruned`:

  What this panel is; taken from a
  [`load_genotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/load_genotypes.md)
  list when it says so.

------------------------------------------------------------------------

### `PopStructure$panels()`

The panels this object holds, primary first.

#### Usage

    PopStructure$panels()

------------------------------------------------------------------------

### `PopStructure$add_meta()`

Attach/replace metadata; auto-assigns colours for new columns.

#### Usage

    PopStructure$add_meta(meta, colors = NULL, colours = NULL)

#### Arguments

- `meta`:

  Data frame with a `sample` column.

- `colors`:

  Optional colour overrides (`column -> (level -> colour)`).

------------------------------------------------------------------------

### `PopStructure$set_colors()`

Set/override colour maps for metadata columns.

#### Usage

    PopStructure$set_colors(colors, colours = NULL)

#### Arguments

- `colors`:

  Named list `column -> (level -> colour)`.

------------------------------------------------------------------------

### `PopStructure$set_levels()`

Fix the level order of a metadata column. Because every plot reads the
shared metadata and colour map, this one order flows through the legends
(PCA / UMAP), the admixture facet order, and the colour strips. Existing
colours follow their level (only the order changes); unlisted levels are
dropped.

#### Usage

    PopStructure$set_levels(column, levels)

#### Arguments

- `column`:

  Metadata column name.

- `levels`:

  The desired level order.

------------------------------------------------------------------------

### `PopStructure$run_umap()`

Compute a UMAP embedding.

#### Usage

    PopStructure$run_umap(
      pca_components = 30,
      n_neighbors = 15,
      min_dist = 0.1,
      seed = 42
    )

#### Arguments

- `pca_components`:

  PCs feeding UMAP: a count (`>= 1`) or a variance fraction
  (`0 < x < 1`, e.g. `0.1` uses enough PCs for 10% of variance, via
  [`n_pcs_for_variance()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/n_pcs_for_variance.md)).

- `n_neighbors, min_dist`:

  UMAP parameters.

- `seed`:

  Random seed.

------------------------------------------------------------------------

### `PopStructure$run_snmf()`

Fit sNMF admixture (cached and quiet by default).

#### Usage

    PopStructure$run_snmf(
      K = 1:10,
      rep = 10,
      alpha = 10,
      seed = 42,
      cpu = 1,
      cache = TRUE,
      cache_dir = NULL,
      verbose = FALSE,
      log_file = NULL
    )

#### Arguments

- `K, rep, alpha, seed, cpu, cache, cache_dir, verbose, log_file`:

  Passed to
  [`run_snmf()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_snmf.md).

------------------------------------------------------------------------

### `PopStructure$best_k()`

Best K (cross-entropy) from the fitted sNMF.

#### Usage

    PopStructure$best_k(stat = c("min", "mean"))

#### Arguments

- `stat`:

  Combine replicates by `"min"` or `"mean"`; see
  [`snmf_best_k()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/snmf_best_k.md).

------------------------------------------------------------------------

### `PopStructure$cross_entropy()`

Per-K cross-entropy summary of the sNMF replicates (see
[`snmf_cross_entropy()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/snmf_cross_entropy.md)).

#### Usage

    PopStructure$cross_entropy(...)

#### Arguments

- `...`:

  Passed to
  [`snmf_cross_entropy()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/snmf_cross_entropy.md).

------------------------------------------------------------------------

### `PopStructure$plot_cross_entropy()`

Cross-entropy elbow plot for choosing K (see
[`plot_snmf_cross_entropy()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_snmf_cross_entropy.md)).

#### Usage

    PopStructure$plot_cross_entropy(...)

#### Arguments

- `...`:

  Passed to
  [`plot_snmf_cross_entropy()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_snmf_cross_entropy.md).

------------------------------------------------------------------------

### `PopStructure$get_snmf_fit()`

The fitted sNMF result (`NULL` before
[`run_snmf()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_snmf.md)).

#### Usage

    PopStructure$get_snmf_fit()

------------------------------------------------------------------------

### `PopStructure$plot_admixture_multi_k()`

One admixture plot per K as pages, for a multi-page PDF (see
[`plot_admixture_multi_k()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_admixture_multi_k.md)).

#### Usage

    PopStructure$plot_admixture_multi_k(...)

#### Arguments

- `...`:

  Passed to
  [`plot_admixture_multi_k()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_admixture_multi_k.md).

------------------------------------------------------------------------

### `PopStructure$q()`

Q (ancestry) matrix at K, restricted to the active samples.

#### Usage

    PopStructure$q(K = NULL, run = NULL)

#### Arguments

- `K`:

  Number of ancestral populations (default the best K).

- `run`:

  Replicate (default the lowest cross-entropy run).

------------------------------------------------------------------------

### `PopStructure$restrict()`

Restrict to a set of samples in place (used by
[`subset()`](https://rdrr.io/r/base/subset.html)).

#### Usage

    PopStructure$restrict(samples)

#### Arguments

- `samples`:

  Sample ids to keep.

------------------------------------------------------------------------

### `PopStructure$subset()`

A new `PopStructure` limited to given samples and/or metadata matches
(does not recompute PCA/UMAP/sNMF – the embeddings are shared and simply
filtered).

#### Usage

    PopStructure$subset(samples = NULL, drop_samples = NULL, ...)

#### Arguments

- `samples`:

  Sample ids to keep.

- `drop_samples`:

  Sample ids to remove, applied after `samples` and the metadata
  filters. For taking out the one or two samples an otherwise good
  subset should not include – a contaminant, a duplicate – without
  having to spell out the keep list. Ids that are not in the object are
  reported rather than ignored, since a typo would otherwise look
  exactly like a sample that was already gone.

- `...`:

  `column = value(s)` metadata filters (e.g. `region = "West Africa"`).

------------------------------------------------------------------------

### `PopStructure$genotype()`

The genotype matrix for the active samples (samples x SNPs).

#### Usage

    PopStructure$genotype(panel = NULL, prefer = NULL)

#### Arguments

- `panel`:

  Panel to return by name; the primary one when `NULL`.

- `prefer`:

  Panel to use *if the object has it*, falling back to the primary one –
  how an analysis asks for the panel it wants without requiring it.

------------------------------------------------------------------------

### `PopStructure$pca_scores()`

PCA scores for the active samples.

#### Usage

    PopStructure$pca_scores()

------------------------------------------------------------------------

### `PopStructure$pca_variance()`

PCA variance-explained table.

#### Usage

    PopStructure$pca_variance()

------------------------------------------------------------------------

### `PopStructure$prcomp()`

The full [`stats::prcomp()`](https://rdrr.io/r/stats/prcomp.html) object
(all samples).

#### Usage

    PopStructure$prcomp()

------------------------------------------------------------------------

### `PopStructure$umap_df()`

UMAP data frame for the active samples (or `NULL`).

#### Usage

    PopStructure$umap_df()

------------------------------------------------------------------------

### `PopStructure$get_meta()`

Metadata for the active samples.

#### Usage

    PopStructure$get_meta()

------------------------------------------------------------------------

### `PopStructure$get_colors()`

The shared colour maps.

#### Usage

    PopStructure$get_colors()

------------------------------------------------------------------------

### `PopStructure$get_colours()`

Alias for `get_colors()`.

#### Usage

    PopStructure$get_colours()

------------------------------------------------------------------------

### `PopStructure$allele()`

Which allele the dosages count (`"alt"` / `"ref"`), or `NULL` when the
object does not record it (built from a bare matrix, or saved by an
older version).

#### Usage

    PopStructure$allele(panel = NULL)

#### Arguments

- `panel`:

  Which panel to report on; the primary one when `NULL`.

------------------------------------------------------------------------

### `PopStructure$positions()`

Which convention the SNP positions follow, or `NULL` when the object
does not record it (built before this was tracked – rebuild it, or pass
`one_based`).

#### Usage

    PopStructure$positions()

------------------------------------------------------------------------

### `PopStructure$pruned()`

Whether the SNPs were LD-pruned (`TRUE` / `FALSE`), or `NULL` when the
object does not record it.

#### Usage

    PopStructure$pruned(panel = NULL)

#### Arguments

- `panel`:

  Which panel to report on; the primary one when `NULL`.

------------------------------------------------------------------------

### `PopStructure$get_samples()`

Active sample ids.

#### Usage

    PopStructure$get_samples()

------------------------------------------------------------------------

### `PopStructure$as_ps()`

A `pop_structure` S3 view (active samples) for the `plot_*()` fns.

#### Usage

    PopStructure$as_ps()

------------------------------------------------------------------------

### `PopStructure$plot_pca()`

PCA scatter coloured by a metadata column (shared colours).

#### Usage

    PopStructure$plot_pca(colour = NULL, pcs = c(1, 2), ..., color = NULL)

#### Arguments

- `colour, color`:

  Metadata column to colour by (either spelling).

- `pcs`:

  Which two PCs.

- `...`:

  Passed to
  [`plot_pca()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_pca.md).

------------------------------------------------------------------------

### `PopStructure$plot_umap()`

UMAP scatter coloured by a metadata column (shared colours).

#### Usage

    PopStructure$plot_umap(colour = NULL, ..., color = NULL)

#### Arguments

- `colour, color`:

  Metadata column to colour by (either spelling).

- `...`:

  Passed to
  [`plot_umap()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_umap.md).

------------------------------------------------------------------------

### `PopStructure$plot_admixture()`

Admixture bars; the group strip reuses the shared colour map, so it
matches the UMAP/PCA colouring.

#### Usage

    PopStructure$plot_admixture(
      K = NULL,
      group = NULL,
      colour = group,
      sample_order = NULL,
      group_bar = !is.null(group),
      ...,
      color = NULL
    )

#### Arguments

- `K`:

  Number of ancestral populations (default best K).

- `group`:

  Metadata column to facet by and colour the strip with.

- `colour, color`:

  Metadata column for the strip colours (default `group`, either
  spelling). To recolour the ancestry clusters themselves, pass
  `colours` / `colors` (a palette), which goes to
  [`plot_admixture()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_admixture.md).

- `sample_order`:

  Explicit sample order (see
  [`admixture_order()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/admixture_order.md)).

- `group_bar`:

  Draw the group colour strip.

- `...`:

  Passed to
  [`plot_admixture()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_admixture.md).

------------------------------------------------------------------------

### `PopStructure$plot_figure()`

Combined UMAP + admixture figure (see
[`plot_structure_figure()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_structure_figure.md)).

#### Usage

    PopStructure$plot_figure(group = NULL, ...)

#### Arguments

- `group`:

  Metadata column to facet/colour the admixture by.

- `...`:

  Passed to
  [`plot_structure_figure()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_structure_figure.md).

------------------------------------------------------------------------

### `PopStructure$pop_diff()`

Per-SNP population differentiation between the levels of a metadata
column (see
[`pop_diff()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff.md));
uses the object's genotype matrix (pass
`genotype = load_genotypes(vcf, prune = FALSE)` to run on the full
unpruned set).

#### Usage

    PopStructure$pop_diff(group = NULL, ...)

#### Arguments

- `group`:

  Metadata column defining the groups.

- `...`:

  Passed to
  [`pop_diff()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff.md)
  (e.g. `statistic = "fst"`, `genotype = `).

------------------------------------------------------------------------

### `PopStructure$jost_d()`

Per-SNP Jost's D between the levels of a metadata column
([`jost_d()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/jost_d.md)).

#### Usage

    PopStructure$jost_d(group = NULL, ...)

#### Arguments

- `group`:

  Metadata column defining the groups.

- `...`:

  Passed to
  [`jost_d()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/jost_d.md).

------------------------------------------------------------------------

### `PopStructure$pop_diff_table()`

Group-pair differentiation table across all statistics (see
[`pop_diff_table()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff_table.md)).

#### Usage

    PopStructure$pop_diff_table(group = NULL, ...)

#### Arguments

- `group`:

  Metadata column defining the groups.

- `...`:

  Passed to
  [`pop_diff_table()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff_table.md).

------------------------------------------------------------------------

### `PopStructure$plot_diff_heatmap()`

Group x group differentiation triangle heatmap (see
[`plot_diff_heatmap()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_diff_heatmap.md));
metadata annotation is resolved against this object's metadata
automatically.

#### Usage

    PopStructure$plot_diff_heatmap(group = NULL, ...)

#### Arguments

- `group`:

  Metadata column defining the groups.

- `...`:

  Passed to
  [`plot_diff_heatmap()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_diff_heatmap.md),
  plus `statistic` (`"jost_d"` default, `"gst_hedrick"`, `"fst"`)
  selecting the measure, and `genotype` (a full/unpruned override, see
  [`pop_diff()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff.md)).

------------------------------------------------------------------------

### `PopStructure$plot_jost_d_heatmap()`

Alias of
[`plot_diff_heatmap()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_diff_heatmap.md)
for Jost's D.

#### Usage

    PopStructure$plot_jost_d_heatmap(group = NULL, ...)

#### Arguments

- `group`:

  Metadata column defining the groups.

- `...`:

  Passed to
  [`plot_diff_heatmap()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_diff_heatmap.md).

------------------------------------------------------------------------

### `PopStructure$pop_diff_snps()`

Per-SNP differentiation in long form (see
[`pop_diff_snps()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff_snps.md)).

#### Usage

    PopStructure$pop_diff_snps(group = NULL, ...)

#### Arguments

- `group`:

  Metadata column defining the groups.

- `...`:

  Passed to
  [`pop_diff()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff.md)
  (e.g. `statistic = `, `genotype = `).

------------------------------------------------------------------------

### `PopStructure$plot_diff_manhattan()`

Genome-wide differentiation Manhattan (see
[`plot_diff_manhattan()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_diff_manhattan.md)).

#### Usage

    PopStructure$plot_diff_manhattan(group = NULL, ...)

#### Arguments

- `group`:

  Metadata column defining the groups.

- `...`:

  Passed to
  [`plot_diff_manhattan()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_diff_manhattan.md);
  `statistic` / `genotype` go to
  [`pop_diff()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff.md).

------------------------------------------------------------------------

### `PopStructure$diversity()`

Within-group diversity (see
[`pop_diversity()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diversity.md)).

#### Usage

    PopStructure$diversity(group = NULL, ...)

#### Arguments

- `group`:

  Metadata column defining the groups.

- `...`:

  Passed to
  [`pop_diversity()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diversity.md)
  (e.g. `by = `, `accessible = `).

------------------------------------------------------------------------

### `PopStructure$ld_index()`

Multilocus index of association (see
[`ld_index()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ld_index.md)).

#### Usage

    PopStructure$ld_index(group = NULL, ...)

#### Arguments

- `group`:

  Metadata column defining the groups.

- `...`:

  Passed to
  [`ld_index()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ld_index.md).

------------------------------------------------------------------------

### `PopStructure$beta_score()`

Beta scores for balancing selection (see
[`beta_score()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/beta_score.md)).

#### Usage

    PopStructure$beta_score(group = NULL, ...)

#### Arguments

- `group`:

  Metadata column defining the groups.

- `...`:

  Passed to
  [`beta_score()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/beta_score.md).

------------------------------------------------------------------------

### `PopStructure$haplotypes()`

Phased haplotypes for a haplotype-homozygosity scan (see
[`parasite_haplotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/parasite_haplotypes.md)).

#### Usage

    PopStructure$haplotypes(...)

#### Arguments

- `...`:

  Passed to
  [`parasite_haplotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/parasite_haplotypes.md)
  (e.g. `fws = `, `maf = `).

------------------------------------------------------------------------

### `PopStructure$plot_region_haplotypes()`

Genotype heatmap over one region, samples clustered (see
[`plot_region_haplotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_region_haplotypes.md)).

#### Usage

    PopStructure$plot_region_haplotypes(region, ...)

#### Arguments

- `region`:

  The interval to draw.

- `...`:

  Passed to
  [`plot_region_haplotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_region_haplotypes.md)
  (e.g. `split = `, `spacing = `).

------------------------------------------------------------------------

### `PopStructure$ihs()`

Integrated haplotype score (see
[`run_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_ihs.md));
builds the haplotypes first unless one is supplied.

#### Usage

    PopStructure$ihs(group = NULL, hap = NULL, ...)

#### Arguments

- `group`:

  Metadata column defining the groups.

- `hap`:

  Optional
  [`parasite_haplotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/parasite_haplotypes.md)
  object to reuse.

- `...`:

  Passed to
  [`run_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_ihs.md).

------------------------------------------------------------------------

### `PopStructure$save()`

Save the whole workspace (genotype, PCA, UMAP, metadata, colours, and
sNMF fit) to an `.rds` file so it can be reloaded without recomputing.
Note: an sNMF fit references LEA project files on disk – run
[`run_snmf()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_snmf.md)
with a persistent `cache_dir` if you want the admixture to survive a
reload.

#### Usage

    PopStructure$save(file, compress = "xz")

#### Arguments

- `file`:

  Destination path.

- `compress`:

  Passed to [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) (default
  `"xz"` for a compact file).

------------------------------------------------------------------------

### `PopStructure$print()`

Compact summary.

#### Usage

    PopStructure$print(...)

#### Arguments

- `...`:

  Ignored.

------------------------------------------------------------------------

### `PopStructure$clone()`

The objects of this class are cloneable with this method.

#### Usage

    PopStructure$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
ps <- example_pop_structure()
ps$run_umap(pca_components = 0.5)      # PCs covering 50% of variance
if (FALSE) { # \dontrun{
ps$run_snmf(K = 1:6)                   # cached + quiet
ps$plot_admixture(K = ps$best_k(), group = "population")
west <- ps$subset(population = c("PopA", "PopB"))
} # }
```
