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
