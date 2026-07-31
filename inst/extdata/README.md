# Example IBD dataset

A small, **public** fixture for exercising the IBD plots, loaded by
`example_ibd_results()`.

**Provenance.** Derived from publicly available *P. falciparum* whole-genome
samples across five African countries (Tanzania, Kenya, DRC, Ethiopia, Sudan;
205 samples). Identity-by-descent was computed with the `plasgenomicsutils` IBD
tools and downsampled to ~1,357 SNPs for a compact fixture. Values are real (not
synthetic). Regions are the sample countries.

Files (all region-level; no per-sample genotype data):

- `ibd_per_snp_region.tsv.gz` — per-SNP fraction of pairs IBD, per region
  (`chr`, `pos`, `region`, `frac_pairs_ibd`)
- `ibd_pairwise_region.tsv.gz` — per-SNP IBD between region pairs
  (`chr`, `pos`, `region_a`, `region_b`, `frac_pairs_ibd`)
- `ibd_selection_per_region.tsv.gz` (+ `..._threshold.tsv`) — per-region IBD
  selection statistic
- `ibd_selection_global.tsv.gz` (+ `..._threshold.txt`) — genome-wide selection
  statistic
- `sample_regions.tsv` — the public sample → country map

## Example population-structure dataset

`pop_structure_ghana_cambodia.rds`, loaded by `example_pop_structure()`, is a small
**public** genotype matrix for the PCA / UMAP / sNMF-admixture tools.

**Provenance.** 60 publicly available *P. falciparum* samples (30 Ghana, 30 Cambodia)
at 49 biallelic SNPs, extracted from the Pf7 example BCF that also ships with the
companion Python package (`plasgenomicsutils/tests/data/ghana_cambodia.pf7.tiny.bcf`).
Genotypes are real alt-allele dosages (`0`/`1`/`2`, `NA` where uncalled); the two
countries are strongly differentiated, so they separate on PC1 and sNMF picks K = 2.

Contents: a list with `genotype` (a 60 × 49 samples-by-SNP integer matrix, sample ids
as row names) and `meta` (a data frame of `sample`, `country`).

`pop_structure_africa.rds` (`example_pop_structure("africa")`) is a richer, multi-region
**public** dataset for the combined UMAP + admixture figure (`plot_structure_figure()`):
258 published East/Central-African samples (`country`: DRC, Kenya, Tanzania, Uganda;
finer `site`: Kenya_East/West, Tanzania_East/West, three historical Uganda groups, DRC; a
macro `region`). Sample names and metadata come from a published supplemental table of
other-study samples; genotypes are extracted from that study's public combined BCF. To
make the regional structure clear in a compact fixture, the 2,000 SNPs kept are the ones
that **most differentiate the sites** — the top per-pair Jost's D markers
(`top_differentiating_snps()`), selected round-robin across all site pairs from ~20k
unpruned biallelic SNPs (an intentional, documented ascertainment for illustration).
Contents: a list of `genotype` (258 × 2000 integer matrix) and `meta` (`sample`,
`country`, `site`, `region`).
