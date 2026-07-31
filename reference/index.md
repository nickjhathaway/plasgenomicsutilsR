# Package index

## IBD results

Ingest the Python IBD tool outputs into one container.

- [`IbdResults`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
  : IBD post-analysis results
- [`ibd_results()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_results.md)
  : Create an IbdResults object

## Plots

Genome-wide and region-by-region IBD figures.

- [`plot_ibd_manhattan()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_manhattan.md)
  : IBD Manhattan plot
- [`plot_selection_manhattan()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_selection_manhattan.md)
  : IBD selection-statistic Manhattan plot
- [`plot_ibd_tugofwar()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_tugofwar.md)
  : IBD / selection "tug-of-war" mirror plot
- [`plot_ibd_region_heatmap()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_region_heatmap.md)
  : Region-by-region IBD heatmap along the genome
- [`plot_drug_gene_triangles()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_drug_gene_triangles.md)
  : IBD "triangle" panels for genes or specific SNPs

## Population structure

PCA / UMAP and sNMF admixture from genotypes, wrapped in one R6
workspace.

- [`PopStructure`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PopStructure.md)
  : Population-structure workspace (PCA + UMAP + admixture)
- [`run_ld_prune()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_ld_prune.md)
  : LD-prune a VCF and return the genotype matrix
- [`pop_structure()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_structure.md)
  : Compute PCA and UMAP from a genotype matrix
- [`plot_pca()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_pca.md)
  : PCA scatter plot
- [`plot_umap()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_umap.md)
  : UMAP scatter plot
- [`n_pcs_for_variance()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/n_pcs_for_variance.md)
  : Number of PCs explaining a target cumulative variance
- [`run_snmf()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_snmf.md)
  : Run sNMF (LEA) admixture on a genotype matrix
- [`snmf_best_k()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/snmf_best_k.md)
  : Pick the best K from an sNMF fit by cross-entropy
- [`snmf_q()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/snmf_q.md)
  : Best-run Q (ancestry proportion) matrix for a given K
- [`admixture_order()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/admixture_order.md)
  : Sample order for admixture bars
- [`plot_admixture()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_admixture.md)
  : Admixture (STRUCTURE) bar plot from a Q matrix
- [`plot_structure_figure()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_structure_figure.md)
  : Combined UMAP + admixture figure
- [`load_pop_structure()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/load_pop_structure.md)
  : Load a saved PopStructure workspace
- [`example_pop_structure()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/example_pop_structure.md)
  : Public population-structure example datasets

## Population differentiation

Per-SNP allele differentiation between metadata groups (Jost’s D, Nei’s
Gst, Hedrick’s G’st, Hudson’s Fst), a group summary + triangle heatmap,
and marker selection.

- [`pop_diff()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff.md)
  : Population differentiation between metadata groups, per SNP
- [`jost_d()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/jost_d.md)
  : Per-SNP Jost's D between metadata groups
- [`pop_diff_matrix()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff_matrix.md)
  [`jost_d_matrix()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff_matrix.md)
  : Group x group differentiation summary matrix
- [`pop_diff_table()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff_table.md)
  : Group-pair differentiation summary across statistics
- [`plot_diff_heatmap()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_diff_heatmap.md)
  [`plot_jost_d_heatmap()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_diff_heatmap.md)
  : Triangle heatmap of group x group differentiation
- [`top_differentiating_snps()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/top_differentiating_snps.md)
  : The SNPs that most differentiate groups

## Colours

Colour-blind-friendly palettes and a shared level-to-colour map.

- [`meta_colors()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/meta_colors.md)
  : Assign colours to the levels of metadata columns
- [`colorPalette_08`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/color_palettes.md)
  [`colorPalette_12`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/color_palettes.md)
  [`colorPalette_15`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/color_palettes.md)
  : Colour-blind-friendly categorical palettes

## Saving and sizing

- [`save_plot()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/save_plot.md)
  : Save a plot, preferring the cairo PDF device
- [`pdf_device()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pdf_device.md)
  : The preferred PDF graphics device
- [`plot_dims()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_dims.md)
  : Suggested output dimensions for an IBD plot

## Example data

- [`example_ibd_results()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/example_ibd_results.md)
  : Load the bundled example IBD results
- [`EXAMPLE_DRUG_GENES`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/EXAMPLE_DRUG_GENES.md)
  : Pf3D7 drug-resistance gene coordinates

## Reference genome

Species- and assembly-specific facts, namespaced so the tools stay
general.

- [`get_reference()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/get_reference.md)
  : Look up a reference genome's facts by id
- [`available_references()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/available_references.md)
  : List available reference ids
- [`normalise_chr()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/normalise_chr.md)
  : Normalise a chromosome name to a bare number string
- [`DEFAULT_REFERENCE`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/DEFAULT_REFERENCE.md)
  : Default reference id
- [`PF3D7_CORE_CHROM_LENGTHS_BP`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF3D7_CORE_CHROM_LENGTHS_BP.md)
  : Pf3D7 core chromosome lengths (bp)
- [`PF3D7_BP_PER_CM`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF3D7_BP_PER_CM.md)
  : Pf3D7 constant genetic-map rate (bp/cM)
