# plasgenomicsutilsR

> **Version 0.3.1** — early development; APIs, defaults, and outputs may
> change between versions.

R utilities for **visualizing and analyzing Plasmodium genomics data** —
the R-side companion to the Python package
[`plasgenomicsutils`](https://github.com/nickjhathaway/plasgenomicsutils),
which does the heavy post-processing compute (VCF
filtering/harmonization, IBD analysis) and writes plain tables that this
package reads and plots.

Reference-genome facts are namespaced by species (see
[`get_reference()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/get_reference.md)),
so the tools generalize beyond *Plasmodium falciparum*.

## Install

``` r

# install.packages("remotes")
remotes::install_github("nickjhathaway/plasgenomicsutilsR")
```

That base install pulls only the required (CRAN) dependencies. Every
plotting and analysis dependency is a **`Suggests`** — so it is optional
and *not* installed automatically — and each function checks for what it
needs and prints an install hint if it is missing. You only install what
you actually use:

- **CRAN** — the plots use `ggplot2`, `scales`, `patchwork`,
  `ggnewscale`, `ggtext`; the UMAP uses `uwot`.
- **Bioconductor** — the population-structure tools use **`SNPRelate`**,
  **`gdsfmt`** (LD-pruning a VCF → genotypes) and **`LEA`** (sNMF
  admixture).

Because these come from two repositories, install them
Bioconductor-aware. The simplest one-liner uses `pak`, which reads this
package’s `biocViews` and resolves CRAN **and** Bioconductor deps
together:

``` r

# install.packages("pak")
pak::pak("nickjhathaway/plasgenomicsutilsR", dependencies = TRUE)
```

Or add the Bioconductor repositories first (your usual `setRepositories`
trick), then let `remotes` install the suggested packages:

``` r

setRepositories(ind = 1:3)   # CRAN + Bioconductor software + annotation
remotes::install_github("nickjhathaway/plasgenomicsutilsR", dependencies = TRUE)
```

Or install the Bioconductor packages explicitly and the CRAN ones
normally:

``` r

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("SNPRelate", "gdsfmt", "LEA"))
install.packages(c("ggplot2", "scales", "patchwork", "ggnewscale", "ggtext", "uwot"))
```

(CI installs the full set the same way — `setup-r-dependencies` is
`biocViews`-aware — so `R CMD check` exercises the SNPRelate/LEA code
paths, not just the CRAN ones.)

## What’s here

- **`IbdResults` + `plot_*()`** — an R6 container over the tables the
  Python `plasgenomicsutils ibd` tools emit, and genome-wide plots
  reading from it:

  - [`plot_ibd_sharing_manhattan()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_sharing_manhattan.md)
    — per-SNP fraction of pairs IBD along the genome
  - [`plot_selection_manhattan()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_selection_manhattan.md)
    — the IBD selection statistic, with a `draw_threshold =` line:
    `"bonferroni"`, `"fdr"`, `"permutation"`, `"empirical"` (whichever
    the Python side wrote), or `"all"`
  - [`plot_ibd_pairwise_group_heatmap()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_pairwise_group_heatmap.md)
    — group × group IBD as tiles along the genome
  - [`plot_ibd_tugofwar()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_tugofwar.md)
    — selection (top) vs IBD (bottom) mirror for one group. `top =`
    hangs a different per-SNP scan from the upper half (an
    [`run_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_ihs.md)
    result, say), so the mirror can ask whether the IBD signal coincides
    with a *haplotype* signal rather than with a statistic derived from
    the same segments
  - [`plot_ibd_locus()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_locus.md)
    — one window on two axes: IBD sharing as a step curve on the left, a
    selection scan as points on the right, genes underneath
  - [`plot_pairwise_ibd_for_genes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_pairwise_ibd_for_genes.md)
    — group × group IBD-sharing triangles, per gene. With IBD blocks
    loaded (`ibd_results(blocks=, meta=)`) or a precomputed overlap
    table (`gene_overlap=`, from `plasgenomicsutils ibd_gene_overlap`),
    a gene’s cell is the fraction of pairs whose IBD **block overlaps
    the gene**
    ([`gene_ibd_overlap()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/gene_ibd_overlap.md)),
    so a pair counts when it shares a segment spanning the gene even
    with no genotyped SNP inside it. Without blocks it falls back to
    aggregating the pairwise IBD of SNPs inside the gene; `within`
    widens the window on either path, which matters because gene spans
    are CDS and can be short relative to a sparse panel.
    `snps = "Pf3D7_07_v3:403222"` draws a specific locus. Returns a
    faceted grid, or `individual = TRUE` for one-per-feature plots —
    pair that with `limits = "shared"` so the pages are
    colour-comparable.
  - [`plot_ibd_network()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_network.md)
    — a sample-level IBD network at one gene/locus: nodes are samples,
    an edge joins two whose pair shares an IBD block over the interval.
    Nodes take their colour from `color_group` and their **shape** from
    `shape_group` — two independent metadata columns, so two variables
    read off one plot; `colors` / `shapes` accept a named
    `level -> value` vector (mapped by name, may be partial) or an
    unnamed one (positional). `sharing` decides what an edge requires:
    `"overlap"` (default) if the pair’s segment touches the gene/locus
    anywhere, or `"complete"` if it must span the whole thing — the same
    gene can give very different graphs (on one real cohort 90% of
    pfcrt-sharing pairs share it completely but only 45% for pfdhps).
    Isolated (unconnected) nodes are optionally kept (to show the
    total N) or dropped. `spread` (default 1.5) weights edge attraction
    by `(1 - J)^spread`, `J` being the Jaccard overlap of two samples’
    IBD neighbourhoods, so a densely inter-connected group opens into a
    readable disc instead of collapsing to a blob while genuinely
    separate clusters keep their separation. Needs `igraph` + `ggraph`.
  - [`gene_ibd_pairs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/gene_ibd_pairs.md)
    — the adjacency list behind the triangles: one row per sample pair ×
    IBD block × gene, with whether the block covers the gene
    `"complete"`ly or `"partial"`ly, the covered span, and
    `percent_covered`. Pairs with no IBD over a gene are simply absent.
    Mirrors `plasgenomicsutils ibd_gene_pairs`, which writes the same
    table as a TSV.
  - [`pos_selection_genes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pos_selection_genes.md)
    — the genes hit by an above-threshold selection signal: intersects
    the significant SNPs with the `genes` track, counting a SNP for a
    gene when it lands within `within` bp of it (default 2 kb, since
    filtering can push the peak just outside a gene). Pass the full
    `PF3D7_GENES` track to scan every gene.

  Blocks are filtered on ingest: segments with fewer than **15 SNPs** or
  shorter than **15 kb** are discarded (`min_block_snp` /
  `min_block_kb`, `0` to disable), since small IBD blocks are commonly
  spurious and this filter is conventionally applied before any summary.
  Only the IBD evidence is filtered — the analyzed-sample set behind
  every denominator still comes from every row, so a pair whose only
  segment was short still counts as compared.
  [`print()`](https://rdrr.io/r/base/print.html) reports how many
  segments went.

  The gene track (`highlight_genes` / the `genes =` argument) defaults
  to every gene in the object’s track; naming a gene that is not in the
  track is an error (not a silent no-op), and passing `highlight_genes`
  labels those genes by default. ggplot2 / scales are optional
  (Suggests); each plot returns a ggplot object. Every `plot_*()` and
  [`pos_selection_genes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pos_selection_genes.md)
  is also a method, so `ibd$plot_selection_manhattan()` and
  `plot_selection_manhattan(ibd)` are interchangeable.

- **`zoom =` on every genome-wide plot** (both Manhattans, the
  tug-of-war,
  [`plot_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ihs.md),
  [`plot_beta()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_beta.md),
  [`plot_diversity()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_diversity.md))
  crops to one interval — a chromosome, a `"chr:start-end"` range, a
  gene name, or a data frame — keeping the same data and the same
  coordinates, so a locus sits at the same x as in the whole-genome
  figure. `zoom_pad` adds context, one value for both sides or two for
  the left and the right. Inside a window every gene in view is drawn at
  its real extent and named in a track underneath; `genes_for_track =`
  fills that track from a separate table (e.g. `PF3D7_GENES`) while the
  plot’s own short track still supplies the marked positions, and
  `gene_label_angle` turns long systematic ids out of each other’s way.

- **[`add_ibd_clusters()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/add_ibd_clusters.md)**
  writes the single-linkage IBD clusters at a gene or locus into the
  object’s metadata as `<gene>_cluster_id`, so any plot that colours or
  shapes by a metadata column can name the blobs
  [`plot_ibd_network()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_network.md)
  draws.

- **[`annotate_snps()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/annotate_snps.md)**
  labels any table with a `snp_id` (or chr/pos) with the intervals it
  falls in — genes, core regions, anything BED-shaped — with `within =`
  for a flanking window.

- **[`aa_intervals()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/aa_intervals.md)**
  turns “codon 76 of *pfcrt*” into a genomic interval, walking the
  protein back through the transcript’s coding exons from a GFF
  ([`read_gff_cds()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/read_gff_cds.md)
  parses it once). Resistance markers are named by amino acid while
  every plot here works in genomic coordinates, and the conversion
  depends on the exon structure, the strand and the CDS phase — so it is
  not something to do by eye. The result is shaped like the package’s
  other interval tables (`chr`, `start`, `end`,
  `name = <transcript>-AA<pos>`), so it drops straight into `genes =`,
  `mark_snps =`,
  [`annotate_snps()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/annotate_snps.md)
  or
  [`bed_intersect()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/bed_intersect.md).
  A codon straddling an intron is flagged (`spans_intron`) since its
  interval then spans the intron too.

- **[`snp_aa_positions()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/snp_aa_positions.md)**
  is the other direction: given SNP positions, which codon of which
  transcript each one sits in — `aa_position` 1-based as the literature
  numbers residues, `codon_base` 1/2/3 in transcript orientation (so on
  a minus-strand gene base 1 is the highest coordinate), `NA` for
  anything non-coding. Take the genotyped SNPs in a gene and get the
  residues they cover, without leaving R.

- **Population structure** — `PopStructure` is an R6 workspace bundling
  a genotype matrix, its PCA (full `prcomp`), an optional UMAP,
  per-sample metadata, a **shared colour map**, and an sNMF admixture
  fit, so PCA / UMAP / admixture colour and order consistently and the
  whole thing can be `saveRDS`-ed and re-plotted without recomputing. It
  offers `run_umap()` (PC count *or* a fraction of variance to capture,
  via
  [`n_pcs_for_variance()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/n_pcs_for_variance.md)),
  [`run_snmf()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_snmf.md)
  (**cached and quiet** — sNMF is slow and noisy),
  [`subset()`](https://rdrr.io/r/base/subset.html) (down to given
  samples or a metadata match), `best_k()` /
  [`q()`](https://rdrr.io/r/base/quit.html), and
  [`plot_pca()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_pca.md)
  /
  [`plot_umap()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_umap.md)
  /
  [`plot_admixture()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_admixture.md)
  (bars ordered once via
  [`admixture_order()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/admixture_order.md),
  per-sample borders so near-identical neighbours stay distinct, and a
  group strip coloured to match the UMAP). `cross_entropy()` summarises
  the sNMF replicates per K (`min` / `mean` / `max`, plus the `best_run`
  index that [`q()`](https://rdrr.io/r/base/quit.html) returns) and
  `plot_cross_entropy()` draws the elbow with the replicate spread as a
  band — a flat stretch or a wide band means the data do not pin K down,
  whatever `best_k()` says.
  [`plot_admixture_multi_k()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_admixture_multi_k.md)
  returns the elbow plus one page per K, sharing one sample order taken
  from the best K, ready for
  [`save_plot()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/save_plot.md)
  as a multi-page PDF. A long cluster legend wraps into columns (or a
  shallow row block under `legend_position = "bottom"`) and the
  suggested canvas grows to hold it, so a K of 15 keys does not run off
  the page. `set_levels()` fixes a metadata column’s order once and it
  flows through every legend, facet, and strip;
  [`save()`](https://rdrr.io/r/base/save.html) /
  [`load_pop_structure()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/load_pop_structure.md)
  persist the whole workspace.
  [`load_genotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/load_genotypes.md)
  builds the genotype matrix from a VCF (SNPRelate, LD-pruned by
  default), and records the two facts a matrix cannot carry: which
  allele its dosages count and whether it was pruned. SNPRelate / uwot /
  LEA / ggnewscale / patchwork are optional (Suggests).

- **Both SNP panels in one object** — LD pruning keeps one SNP out of
  each correlated run, which is right for PCA / UMAP / admixture and the
  opposite of what you want wherever that correlation *is* the signal:
  differentiation, diversity, LD, haplotype scans, haplotype plots.
  Register both with
  `ps$add_panel("full", load_genotypes(vcf, prune = FALSE))` and each
  analysis takes the one it needs — `ps$plot_pca()` the pruned SNPs,
  `pop_diff(ps)` the full set — instead of building a second object each
  time. `$genotype(panel =)` requires a panel by name,
  `$genotype(prefer =)` takes one if it exists, and any other name works
  too. Without a full panel the analyses that want one say so once
  rather than quietly using pruned SNPs.

- **[`plot_structure_figure()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_structure_figure.md)**
  — the UMAP and sNMF admixture as **one** figure: a shared theme
  (matching fonts), one region colour map (UMAP points match the
  admixture colour strips), collected legends, region-faceted bars
  (colour strip, no text) laid out over configurable `rows`, and
  `orientation = "vertical"`/`"horizontal"`, auto-sized.
  [`example_pop_structure()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/example_pop_structure.md)
  ships two **public** demos: `"ghana_cambodia"` (minimal) and
  `"africa"` (258 East-African samples across DRC / Kenya / Tanzania /
  Uganda sites).

- **Population differentiation** —
  [`pop_diff()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff.md)
  computes a per-SNP differentiation statistic for every pair of
  metadata groups: **Jost’s D**
  ([`jost_d()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/jost_d.md)),
  **Hedrick’s standardized G′st**, and **Hudson’s Fst**. It wants the
  **full, unpruned** genotypes (LD-pruning removes the differentiating
  SNPs), which a `PopStructure` supplies when it holds a full panel
  alongside the pruned one —
  `ps$add_panel("full", load_genotypes(vcf, prune = FALSE))` — and every
  analysis whose signal is SNP-to-SNP correlation then reads it
  automatically. Only one Gst is offered: plain Nei’s Gst is strongly
  deflated when within-group diversity is high (typical in *P.
  falciparum*), so its standardized form G′st is kept.
  [`pop_diff_matrix()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff_matrix.md)
  collapses to a group × group summary;
  [`plot_diff_heatmap()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_diff_heatmap.md)
  draws a triangle heatmap (legend in the empty corner, optional
  clustering **dendrogram** with leaf tips coloured to match the
  UMAP/admixture, and metadata **annotation** strips);
  [`pop_diff_table()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff_table.md)
  returns every statistic × summary per group-pair;
  [`pop_diff_snps()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff_snps.md)
  unpacks the per-SNP values (with `chr:pos` coordinates) and
  [`plot_diff_manhattan()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_diff_manhattan.md)
  draws them along the genome. Because most of the genome is barely
  differentiated, a genome-wide **mean looks near-zero** — use
  `stat = "top_mean"`, `"max"`, or a `trans = "sqrt"` fill to see the
  signal.
  [`top_differentiating_snps()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/top_differentiating_snps.md)
  picks the most differentiating markers (round-robin across pairs); the
  `"africa"` example fixture is itself built from them. Estimators
  follow Jost (2008), Nei & Chesser (1983), Hedrick (2005), and Hudson
  et al. (1992) / Bhatia et al. (2013).

- **Within-population diversity** —
  [`pop_diversity()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diversity.md)
  reports, for each metadata group and genome-wide / per gene / per
  window: nucleotide diversity, expected heterozygosity, Watterson’s
  theta, Tajima’s D, segregating sites, and the haplotype /
  multilocus-genotype summaries (`hap_div`, Shannon `H`, Simpson’s
  `lambda`, evenness). **`pi` is per accessible base pair** — the
  denominator is `accessible` callable sites
  (e.g. `PF3D7_CORE_REGIONS`), not the SNP count, so windows of
  different SNP density stay comparable; the per-SNP average is reported
  separately as `he` and the two are never interchangeable. The parasite
  is haploid, so a heterozygous call is read as a mixed infection and
  dropped at that site by default (`het = "dosage"` splits it instead).
  Tajima’s D keeps every site and uses the mean number of calls as *n*
  rather than discarding samples with gaps.
  [`plot_diversity()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_diversity.md)
  draws a windowed track. Tajima’s D, haplotype diversity and pi are
  verified against to floating-point.

- **Linkage disequilibrium** —
  [`ld_index()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ld_index.md)
  gives the multilocus index of association `Ia` and its standardized
  `rbarD` (verified against ), the genome-wide read on clonality. The r²
  **decay** curve is quadratic in SNP count and so runs in the Python
  package (`plasgenomicsutils ld_decay`, ~7 s on a 249-sample callset
  against ~2 min here);
  [`read_ld_decay()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/read_ld_decay.md)
  reads it back — half-decay distance and scan settings attached — and
  [`plot_ld_decay()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ld_decay.md)
  draws it.

- **Selection scans** —
  [`parasite_haplotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/parasite_haplotypes.md)
  builds the complete, phased 0/1 haplotypes needs from a callset that
  has neither: gate to monoclonal infections on Fws, resolve remaining
  mixed calls by drawing at the population frequency, filter, impute —
  reporting every sample and SNP it removed, because how the haplotypes
  were made determines what the scan can claim. Then
  [`run_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_ihs.md)
  (within a population),
  [`run_rsb()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_rsb.md)
  /
  [`run_xpehh()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_xpehh.md)
  (between two),
  [`ihs_genes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ihs_genes.md)
  for the per-gene peak, and
  [`plot_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ihs.md).
  Without an outgroup the scan is unpolarised, so read `abs(ihs)`, not
  its sign.
  [`beta_score()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/beta_score.md)
  covers the other half: long-term **balancing** selection from allele
  frequencies clustered around an intermediate-frequency core (Siewert &
  Voight’s folded Beta1), with
  [`beta_genes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/beta_genes.md)
  and
  [`plot_beta()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_beta.md)
  — the antigen counterpart to iHS’s sweeps.
  `subset_haplotypes(hap, region = c("North", "Southwest"))` keeps only
  some of the haplotypes — by sample or by metadata group, several
  values allowed — for a mutation that segregates in one part of the
  cohort and is buried by pooling. The SNP panel is left as built, so
  subsets stay comparable with each other and with the whole.
  [`plot_ehh()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ehh.md)
  is the picture behind one point on a scan: the haplotype homozygosity
  decaying either side of a focal SNP, one curve per allele **at that
  SNP** — the mutant-versus- reference comparison without needing the
  SNPs annotated.

- **[`selection_peaks()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/selection_peaks.md)**
  turns any per-SNP scan into a peak list: contiguous runs of
  above-threshold SNPs, with every gene the peak interval covers
  (`peak_interval_genes`, a list-column in genomic order), the gene at
  the peak SNP itself, and the nearest gene and its distance for a peak
  in an intergenic stretch. `criterion` picks the bar — `"bonferroni"`,
  `"fdr"`, `"permutation"`, `"empirical"`, `"top"` or `"value"` — so the
  list is built at a threshold the run actually supports.

- **[`plot_region_haplotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_region_haplotypes.md)**
  shows the genotypes themselves over a window: one row per sample, one
  column per SNP, samples clustered with a dendrogram beside them,
  coloured metadata strips down the right and a gene track underneath.
  `split =` blocks the rows by a metadata column and clusters *within*
  each block (`ComplexHeatmap`’s `row_split` semantics), so a haplotype
  shared across a group reads as a solid band rather than being
  scattered by one global ordering. `spacing = "genomic"` moves every
  mark to its real coordinate, keeping them all the same width, so the
  gaps between SNPs are what you see.

- **Coverage QC** —
  [`read_coverage()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/read_coverage.md),
  [`coverage_qc()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/coverage_qc.md),
  [`plot_coverage_summary()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_coverage_summary.md),
  [`plot_coverage_by_chrom()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_coverage_by_chrom.md)
  and
  [`plot_coverage_dropout()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_coverage_dropout.md)
  read and plot the depth tables from
  `plasgenomicsutils coverage_depth_stats` / `coverage_dropout_regions`.
  Breadth matters more than mean depth: sWGA can give a respectable
  average while leaving much of the genome at zero, and only the breadth
  column shows it.

- **Genomic intervals** —
  [`bed_intersect()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/bed_intersect.md)
  overlaps two BED-style interval tables (configurable
  `chr`/`start`/`end` columns, chromosome spellings reconciled) and
  returns `overlap` / `only1` / `only2`. Bundled region tracks
  `PF3D7_CORE_REGIONS` (core vs. subtelomeric/hypervariable) and
  `PF3D7_PARALOG_GENES` let you classify genes, e.g.
  `bed_intersect(PF3D7_GENES, PF3D7_CORE_REGIONS)$only1` are the
  subtelomeric genes.

- **Coordinates are 0-based throughout**
  ([`?"plasgenomicsutilsR-coordinates"`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plasgenomicsutilsR-coordinates.md))
  — intervals half-open `[start, end)` as in BED, and variant positions
  0-based too, so there is one rule and no part of the package to
  remember an exception for. Sources that number differently are
  converted once at the boundary: the PlasmoDB GFF where the gene
  datasets are built, `hmmibd-rs` block ends when an `IbdResults` reads
  them, and VCF `POS` in the Python package before any table reaches R.
  `PF3D7_GENES` gives each gene’s **CDS** span (the translated extent,
  introns included, UTRs excluded).

- Reference registry:
  [`get_reference()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/get_reference.md),
  [`available_references()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/available_references.md),
  [`normalise_chr()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/normalise_chr.md),
  `PF3D7_CORE_CHROM_LENGTHS_BP`.

Fws is computed by the companion Python package
(`plasgenomicsutils calculate_fws`), a reimplementation of
`moimix::getFws`.

A small **public** example dataset (five African countries) ships with
the package:

``` r

library(plasgenomicsutilsR)

ibd <- example_ibd_results()          # bundled public example
plot_ibd_sharing_manhattan(ibd)
plot_selection_manhattan(ibd, metric = "neg_log10_p")
plot_ibd_tugofwar(ibd, group = "Tanzania")
plot_ibd_pairwise_group_heatmap(ibd, trans = "log2")   # log2 fill reads best for IBD
```

The genome-wide plots share a few options:

- **Colour scale** (`plot_ibd_pairwise_group_heatmap`,
  `plot_pairwise_ibd_for_genes`): `trans` (`"log2"`, `"sqrt"`, …), a
  custom `colors` ramp, `limits = c(lo, hi)` (extremes squished so a few
  high values don’t crush the scale), or a full `fill_scale` override.
  The default is a light single-hue ramp.
- **Skip chromosomes**: `chroms = c("7", "13")` keeps only those, or
  `skip_chr = "1"` drops them — the rest are re-laid-out contiguously.
- **Highlight genes**: gene positions *and* display names come from the
  `genes` track you pass to
  [`ibd_results()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_results.md)
  — a data frame of `name`, `chr`, `start`, `end` (the bundled
  \[`PF_EXAMPLE_DRUG_GENES`\] is one). The Manhattan, tug-of-war, and
  region heatmap draw a reference line at each gene;
  `highlight_genes = c("pfcrt", "pfdhps")` selects which to show
  (case-insensitive) and `label_genes = TRUE` labels them (top panel
  only, just outside the plot). The label text is the track’s `name`, so
  set it to whatever you want displayed (`"CRT"` vs `"pfcrt"`) by
  passing your own track:

``` r

my_genes <- data.frame(name = c("CRT", "DHPS"), chr = c("7", "8"),
                       start = c(403222, 548200), end = c(406317, 550616))
ibd <- ibd_results(per_snp_group = "...", selection = "...", genes = my_genes)
plot_selection_manhattan(ibd, chroms = c("7", "8"),
                         highlight_genes = c("CRT", "DHPS"), label_genes = TRUE)
```

On your own data, point
[`ibd_results()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_results.md)
at the tables the Python `plasgenomicsutils ibd` tools write:

``` r

ibd <- ibd_results(
  per_snp_group  = "ibd_analysis.per_snp_per_group.tsv.gz",
  pairwise_group = "ibd_analysis.per_snp_pairwise_group.tsv.gz",
  selection       = "ibd_selection.per_group.selection_stats.tsv.gz",
  threshold       = "ibd_selection.per_group.threshold.txt",
  reference       = "pf3d7"
)
save_plot("ibd_manhattan.pdf", plot_ibd_sharing_manhattan(ibd), width = 9, height = 4)
```

[`save_plot()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/save_plot.md)
wraps
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)
and, for `.pdf` output, defaults to the cairo PDF device (better font
embedding), falling back to the standard `pdf` device where cairo is
unavailable or unreliable (e.g. Windows). Force a device with
`device =`, or use
[`pdf_device()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pdf_device.md)
directly with `ggsave()`.

It also **sizes the canvas to the drawing**. Plots with a locked panel
shape — the IBD networks and the gene triangles both use `coord_fixed()`
— only fill a canvas of one particular aspect ratio, and on any other
shape the remainder becomes blank margin. Supply one of `width` /
`height` and the other is computed from the built plot’s panel ratio
plus the inches its titles, legends and margins actually need; supply
neither and both are worked out; supply both and they are used verbatim.
`fit = FALSE` disables it. Plots with a free coordinate system (the
Manhattans, tug-of-war, group heatmap) are unaffected — any canvas shape
is legitimate for them, so their attached sizes stand.

Drug-gene triangles read a `genes` track (`name`, `chr`, `start`, `end`)
passed to
[`ibd_results()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_results.md);
a SNP belongs to a gene when its position falls in the gene interval:

``` r

ibd <- ibd_results(pairwise_group = "...pairwise_group.tsv.gz",
                   genes = "drug_resistance_genes.tsv")
plot_pairwise_ibd_for_genes(ibd)
```

## Development

``` r

devtools::load_all()
devtools::test()
devtools::document()
```

## License

GPL-3.
