# Example IBD dataset

A small, **public** fixture for exercising the IBD plots, loaded by
`example_ibd_results()`.

**Provenance.** Derived from publicly available *P. falciparum* whole-genome
samples across five African countries (Tanzania, Kenya, DRC, Ethiopia, Sudan;
205 samples). Identity-by-descent was computed with the `plasgenomicsutils` IBD
tools and downsampled to ~1,357 SNPs for a compact fixture. Values are real (not
synthetic). Regions are the sample countries.

Files (all group-level; no per-sample genotype data). Here the grouping is the sample
country, so the `group` columns hold country names. `pos` is the variant's 0-based
position, matching the package convention (`?plasgenomicsutilsR-coordinates`):

- `ibd_per_snp_group.tsv.gz` — per-SNP fraction of pairs IBD, per group
  (`chr`, `pos`, `group`, `frac_pairs_ibd`)
- `ibd_pairwise_group.tsv.gz` — per-SNP IBD between group pairs
  (`chr`, `pos`, `group_a`, `group_b`, `frac_pairs_ibd`)
- `ibd_selection_per_group.tsv.gz` (+ `..._threshold.tsv`) — per-group IBD
  selection statistic
- `ibd_selection_global.tsv.gz` (+ `..._threshold.txt`) — genome-wide selection
  statistic
- `sample_regions.tsv` — the public sample → country map

## Example population-structure dataset

`pop_structure_ghana_cambodia.rds`, loaded by `example_pop_structure()`, is a small
**public** genotype matrix for the PCA / UMAP / sNMF-admixture tools.

**Provenance.** 60 publicly available *P. falciparum* samples (30 Ghana, 30 Cambodia).
Genotypes are real alt-allele dosages (`0`/`1`/`2`, `NA` where uncalled); the two countries
are strongly differentiated, so they separate on PC1 and sNMF picks K = 2.

**Two SNP panels**, because the two kinds of plot want opposite things:

- `genotype` — 49 biallelic SNPs thinned genome-wide, extracted from the Pf7 example BCF that
  also ships with the companion Python package
  (`plasgenomicsutils/tests/data/ghana_cambodia.pf7.tiny.bcf`). This is what PCA / UMAP /
  admixture read, and the panel those demos were tuned on.
- `full` — the same 49 plus **670** biallelic SNPs (MAF ≥ 2%) covering **pfcrt**, **pfdhps**
  and **pfkelch13** ±30 kb, pulled from the same public Ghana/Cambodia Pf7 callset. 719 SNPs
  in total. Locus-scale plots — `plot_ibd_locus()`, `zoom =`, `plot_region_haplotypes()`,
  `plot_ehh()` — need SNPs *close together* to show anything, and a genome-wide thinning at
  one SNP per few kb leaves 0–2 in a 50 kb window. Registered as the object's `"full"` panel,
  so those plots pick it up on their own while the thinned set stays the primary.

Ids are `chr:pos` with **1-based** VCF positions (SNPRelate's convention, which
`load_genotypes()` follows). Multiallelic sites split by `bcftools norm -m-` can leave two
biallelic records at one position; only the commonest ALT is kept, since a `chr:pos` id can
name just one SNP.

Contents: a list with `genotype` (60 × 49), `full` (60 × 719), `meta` (`sample`, `country`),
and `allele` / `pruned` recording what the dosages count and that the primary panel is
thinned.

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
