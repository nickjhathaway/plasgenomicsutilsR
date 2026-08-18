# Package index

## IBD results

Ingest the Python IBD tool outputs into one container.

- [`IbdResults`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/IbdResults.md)
  : IBD post-analysis results
- [`ibd_results()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ibd_results.md)
  : Create an IbdResults object

## Plots

Genome-wide and group-by-group IBD figures.

- [`plot_ibd_sharing_manhattan()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_sharing_manhattan.md)
  : IBD Manhattan plot
- [`plot_selection_manhattan()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_selection_manhattan.md)
  : IBD selection-statistic Manhattan plot
- [`plot_ibd_tugofwar()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_tugofwar.md)
  : IBD / selection "tug-of-war" mirror plot
- [`plot_ibd_locus()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_locus.md)
  : One locus in detail: IBD sharing against a selection scan
- [`plot_ibd_pairwise_group_heatmap()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_pairwise_group_heatmap.md)
  : Group-by-group IBD heatmap along the genome
- [`plot_pairwise_ibd_for_genes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_pairwise_ibd_for_genes.md)
  : IBD "triangle" panels for genes or specific SNPs
- [`plot_ibd_network()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_network.md)
  : IBD network at a gene or locus
- [`plot_ibd_pair_network()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ibd_pair_network.md)
  : Genome-wide IBD relatedness network

## Genome-wide IBD between groups

The per-pair IBD fraction reduced to one row per pair of metadata
groups, over every sample pair spanning them.

- [`pair_fraction_summary()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pair_fraction_summary.md)
  : Summarise genome-wide IBD sharing between metadata groups

## Gene IBD-block overlap

Which pairs share a gene by IBD: the fraction per group pair
(block-based, not SNP-in-gene), and the underlying pair-by-pair
adjacency list with gene coverage.

- [`gene_ibd_overlap()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/gene_ibd_overlap.md)
  : Per-gene IBD-block overlap between groups
- [`gene_ibd_pairs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/gene_ibd_pairs.md)
  : Sample pairs sharing IBD over each gene
- [`add_ibd_clusters()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/add_ibd_clusters.md)
  : Add IBD cluster ids to the stored metadata

## Selection genes

Genes hit by an above-threshold IBD selection signal.

- [`pos_selection_genes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pos_selection_genes.md)
  : Genes under positive selection

## Population structure

PCA / UMAP and sNMF admixture from genotypes, wrapped in one R6
workspace.

- [`PopStructure`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PopStructure.md)
  : Population-structure workspace (PCA + UMAP + admixture)
- [`load_genotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/load_genotypes.md)
  : Load genotypes from a VCF, optionally LD-pruned
- [`run_ld_prune()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_ld_prune.md)
  : Deprecated name for load_genotypes()
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
- [`snmf_cross_entropy()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/snmf_cross_entropy.md)
  : Cross-entropy of every sNMF replicate, summarised per K
- [`plot_snmf_cross_entropy()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_snmf_cross_entropy.md)
  : Cross-entropy elbow plot for choosing K
- [`admixture_order()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/admixture_order.md)
  : Sample order for admixture bars
- [`plot_admixture()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_admixture.md)
  : Admixture (STRUCTURE) bar plot from a Q matrix
- [`plot_admixture_multi_k()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_admixture_multi_k.md)
  : Admixture bar plots across every K, as pages
- [`plot_structure_figure()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_structure_figure.md)
  : Combined UMAP + admixture figure
- [`load_pop_structure()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/load_pop_structure.md)
  : Load a saved PopStructure workspace
- [`example_pop_structure()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/example_pop_structure.md)
  : Public population-structure example datasets

## Population differentiation

Per-SNP allele differentiation between metadata groups (Jost’s D,
Hedrick’s G’st, Hudson’s Fst), a group summary + triangle heatmap, a
genome-wide track, and marker selection.

- [`pop_diff()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff.md)
  : Population differentiation between metadata groups, per SNP
- [`jost_d()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/jost_d.md)
  : Per-SNP Jost's D between metadata groups
- [`pop_diff_snps()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff_snps.md)
  : Per-SNP differentiation in long form
- [`pop_diff_matrix()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff_matrix.md)
  [`jost_d_matrix()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff_matrix.md)
  : Group x group differentiation summary matrix
- [`pop_diff_table()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diff_table.md)
  : Group-pair differentiation summary across statistics
- [`plot_diff_heatmap()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_diff_heatmap.md)
  [`plot_jost_d_heatmap()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_diff_heatmap.md)
  : Triangle heatmap of group x group differentiation
- [`plot_diff_manhattan()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_diff_manhattan.md)
  : Genome-wide differentiation Manhattan plot
- [`top_differentiating_snps()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/top_differentiating_snps.md)
  : The SNPs that most differentiate groups

## Within-population diversity

Nucleotide diversity (per accessible base pair), expected
heterozygosity, Watterson’s theta, Tajima’s D and
haplotype/multilocus-genotype diversity, genome-wide, per gene or in
windows, for each metadata group.

- [`pop_diversity()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pop_diversity.md)
  : Within-population genetic diversity
- [`tajima_d()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/tajima_d.md)
  : Tajima's D from segregating sites and per-site heterozygosity
- [`tajima_d_pvalue()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/tajima_d_pvalue.md)
  : p-value for a Tajima's D
- [`plot_diversity()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_diversity.md)
  : Windowed diversity along the genome

## Linkage disequilibrium

The multilocus index of association as a genome-wide measure of
clonality, and the reader/plot for the r-squared decay curve computed by
`plasgenomicsutils ld_decay`.

- [`ld_index()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ld_index.md)
  : Multilocus linkage disequilibrium: the index of association
- [`read_ld_decay()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/read_ld_decay.md)
  : Read an LD-decay table
- [`plot_ld_decay()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ld_decay.md)
  : Linkage-disequilibrium decay curve

## Selection scans

Recent directional selection from extended haplotype homozygosity (iHS
within a population, Rsb and XP-EHH between two), and long-term
balancing selection from clustered allele frequencies (beta).

- [`parasite_haplotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/parasite_haplotypes.md)
  : Build phased haplotypes for a haplotype-homozygosity scan
- [`subset_haplotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/subset_haplotypes.md)
  : Keep only some of the haplotypes
- [`run_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_ihs.md)
  : Integrated haplotype score (iHS)
- [`run_rsb()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_rsb.md)
  : Cross-population extended haplotype homozygosity (Rsb)
- [`run_xpehh()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/run_xpehh.md)
  : Cross-population extended haplotype homozygosity (XP-EHH)
- [`ihs_genes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ihs_genes.md)
  : Summarise a haplotype scan per gene
- [`plot_ihs()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ihs.md)
  : Manhattan plot of a haplotype-homozygosity scan
- [`plot_ehh()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_ehh.md)
  : EHH decay around one SNP
- [`ehh_candidates()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/ehh_candidates.md)
  : The SNPs an EHH plot had to choose between
- [`plot_region_haplotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_region_haplotypes.md)
  : Genotypes over one region, clustered by sample
- [`beta_score()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/beta_score.md)
  : Beta: balancing selection from clustered allele frequencies
- [`beta_genes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/beta_genes.md)
  : Summarise beta scores per gene
- [`plot_beta()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_beta.md)
  : Manhattan plot of beta scores
- [`selection_peaks()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/selection_peaks.md)
  : Merge a selection scan's significant SNPs into peaks

## Sequencing coverage QC

Read the depth tables written by
`plasgenomicsutils coverage_depth_stats` / `coverage_dropout_regions`,
apply per-sample QC floors, and plot depth per sample, per chromosome,
and the regions almost no sample amplifies.

- [`read_coverage()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/read_coverage.md)
  : Read a coverage table
- [`coverage_qc()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/coverage_qc.md)
  : Per-sample coverage QC verdict
- [`plot_coverage_summary()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_coverage_summary.md)
  : Per-sample coverage overview
- [`plot_coverage_by_chrom()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_coverage_by_chrom.md)
  : Coverage per chromosome, per sample
- [`plot_coverage_dropout()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_coverage_dropout.md)
  : Coverage dropouts along the genome

## Within-host mixtures

Read the per-site table from `plasgenomicsutils wsaf_profile` and draw
each sample’s heterozygous allele fractions, which separate a dominant
clone with minor companions from a mixture of comparable strains – two
things one Fws value cannot tell apart.

- [`plot_wsaf()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_wsaf.md)
  : Per-sample heterozygous allele-fraction distributions

## Colours

Colour-blind-friendly palettes and a shared level-to-colour map.

- [`meta_colors()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/meta_colors.md)
  : Assign colours to the levels of metadata columns
- [`color_palette()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/color_palette.md)
  : A colour-blind-friendly categorical palette

## Saving and sizing

- [`save_plot()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/save_plot.md)
  : Save a plot, preferring the cairo PDF device
- [`pdf_device()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/pdf_device.md)
  : The preferred PDF graphics device
- [`plot_dims()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plot_dims.md)
  : Suggested output dimensions for an IBD plot

## Genomic intervals

Overlap two BED-style interval tables (genes vs. core / paralog regions,
…), label SNPs with what they fall in, and turn amino-acid positions
into genomic intervals. All coordinates in the package follow one
convention.

- [`annotate_snps()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/annotate_snps.md)
  : Which intervals each SNP falls in
- [`aa_intervals()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/aa_intervals.md)
  : Genomic interval of an amino-acid position
- [`snp_aa_positions()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/snp_aa_positions.md)
  : The amino acid a SNP falls in
- [`read_gff_cds()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/read_gff_cds.md)
  : Read the CDS features of a GFF
- [`bed_intersect()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/bed_intersect.md)
  : Intersect two sets of genomic intervals
- [`bed_subtract()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/bed_subtract.md)
  : Subtract one set of genomic intervals from another
- [`write_bed()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/write_bed.md)
  : Write an interval table as a BED file
- [`haplotype_samples()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/haplotype_samples.md)
  : The samples a haplotype set kept
- [`subset_genotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/subset_genotypes.md)
  : Restrict a genotype panel to a set of samples
- [`plasgenomicsutilsR-coordinates`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/plasgenomicsutilsR-coordinates.md)
  : Genomic coordinate conventions

## Example & reference data

- [`example_ibd_results()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/example_ibd_results.md)
  : Load the bundled example IBD results
- [`PF3D7_GENES`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF3D7_GENES.md)
  : Pf3D7 gene coordinates
- [`PF_EXAMPLE_DRUG_GENES`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF_EXAMPLE_DRUG_GENES.md)
  : Pf3D7 drug-resistance / selection example genes
- [`PF3D7_CORE_REGIONS`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF3D7_CORE_REGIONS.md)
  : Pf3D7 core genome regions
- [`PF3D7_PARALOG_GENES`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/PF3D7_PARALOG_GENES.md)
  : Pf3D7 paralogous / hypervariable gene families

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
